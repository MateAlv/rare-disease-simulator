"""Seeded probabilistic patient simulator (v0.4).

Generates synthetic clinical cases from a validated :class:`DiseaseProfile`.
The generative model for one case (``docs/README.md`` documents every knob):

1. **Patient.** Sex from the disease's sex prior (``sex_prior_key``). Onset
   category from the profile's onset distribution (or the configured prior when
   unknown), an onset age uniform in that category's window, and a current age
   = onset + an exponential disease duration.
2. **Truth.** Each profile phenotype gets a per-patient probability
   ``p ~ Beta(k*mean, k*(1-mean))`` around its frequency (ADR-0007, with
   count-based ``n/m`` shrunk toward a Beta prior per ADR-0011; unknown
   frequencies use a configured mean), raised with duration in progressive
   diseases. It is eligible only if its sex restriction matches and its own
   onset (sampled in its category window) is <= the current age; antenatal and
   congenital terms are always eligible. Eligible terms are present with
   probability ``p``.
3. **Reporting** (``reporting.mode: report_model``, ADR-0011). A term budget
   ``k`` is drawn from the report model's histogram. It counts every present
   term a record shows, noise included, so the noise count ``n`` is drawn next
   and the profile terms get ``max(1, k - n)``: the true cardinal terms
   (``cardinal-v1``) are reported first, then true terms are picked without
   replacement with probability proportional to the reporting model's score
   until that many are reported (all of them when it exceeds the true terms).
   ``reporting.mode: observation`` keeps simulator 0.3's model: each present
   term is recorded with the preset's observation rate. Either way,
   recorded terms may be generalized to a parent term, unrecorded true terms
   become ``missing``/``unknown``, and a term whose generalization is already
   recorded becomes ``missing`` (``reason: recorded_as_generalized``), so every
   truly present term appears in the case. Cases with no true term are
   redrawn; after ``max_redraws`` the most probable eligible term is forced in.
4. **Negatives.** 0..k "asked and absent" terms from four sources: the
   disease's own terms the patient lacks (weighted by frequency), the gene's
   other diseases' terms (gene-first mode), hallmark terms of confounder
   diseases the true disease never annotates, and the profile's explicit
   ``NOT`` annotations. A negated term is never equal to, an ancestor of, or a
   descendant of a present term (or of another negated term).
5. **Noise** from an external vocabulary, never related to a negated term.
6. **Covariate missingness** hides sex, age, onset or all negatives.

Every case is fully determined by ``(seed, disease_id, difficulty, index)``
and the configuration, so reruns reproduce identical bytes.

**Gene-first mode** (:func:`simulate_gene_cases`) labels cases with a GNN
gene. Each case draws one of the gene's disease entities by ``sim_weight``,
with the key ``(seed, "gene:" + symbol, difficulty, index)``, and simulates it
from one profile merged from all the entity's profile ids
(``entity_profiles: merged``, ADR-0011) or, as in simulator 0.3, from one of
them drawn uniformly (``uniform``); the entity's inheritance sets the sex
prior.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import TypeVar

from rare_disease_simulator import __version__
from rare_disease_simulator.data_sources.gene_profiles import EntityOption, GeneTarget
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.frequency import beta_shrinkage_mean, pooled_counts
from rare_disease_simulator.profiles.inheritance import SEX_LIMITED_EXPRESSION, derive_sex_bias
from rare_disease_simulator.profiles.merge import merge_entity_profiles
from rare_disease_simulator.profiles.schema import (
    DiseaseProfile,
    OnsetCategory,
    PhenotypeAssociation,
    SexRestriction,
)
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.reporting import Reporting
from rare_disease_simulator.simulation.schema import (
    Age,
    CasePhenotype,
    CaseTarget,
    Difficulty,
    DifficultyPreset,
    FrequencySettings,
    GeneratorMetadata,
    NegativeSource,
    PatientAttributes,
    Sex,
    SexPriorKey,
    SexSettings,
    SimulationConfig,
    SyntheticCase,
)

SIMULATOR_VERSION = "0.4.2"

CARDINAL_ROLES = {"cardinal", "major"}
CONGENITAL_ONSET = "HP:0003577"
MALE_LIMITED = "HP:0001475"
FEMALE_LIMITED = "HP:0034344"
Y_LINKED = "HP:0001450"
X_LINKED_DOMINANT = "HP:0001423"

NEGATIVE_SOURCES: tuple[NegativeSource, ...] = (
    "own_disease",
    "own_gene_other_disease",
    "confounder",
    "not_annotation",
)
NEGATIVE_ORIGINS: dict[NegativeSource, str] = {
    "own_disease": "negative_own_disease",
    "own_gene_other_disease": "negative_own_gene_other_disease",
    "confounder": "negative_confounder",
    "not_annotation": "negative_not_annotation",
}
GENE_ONLY_SOURCES: frozenset[NegativeSource] = frozenset({"own_gene_other_disease"})
EQUIVALENT_PROFILE_REASON = "typical_feature_absent_equivalent_profile"
FORCED_REASON = "forced_min_one"
MERGED_REASON = "recorded_as_generalized"
CARDINAL_REASON = "cardinal"
SPECIALIZED_REASON = "specialized"
RELATED_NOISE_REASON = "related"
VOCABULARY_NOISE_REASON = "nonspecific_finding"
NOT_REPORTED_REASON = "not_reported"

_T = TypeVar("_T")


@dataclass(frozen=True)
class NoiseTerm:
    """A nonspecific phenotype available to inject as clinical noise."""

    hpo_id: str
    label: str
    weight: float = 1.0


class ShrinkageMeanMissing(ValueError):
    """Beta shrinkage was asked for a count-based term without a prior mean."""


def simulation_frequency(
    phenotype: PhenotypeAssociation, settings: FrequencySettings
) -> float | None:
    """The frequency a term is simulated with; None when the profile gives none.

    Count-based terms (``frequency_raw`` ``n/m``) use the configured estimator:
    Beta shrinkage toward ``shrinkage_mean`` with ``shrinkage_strength``
    pseudo-patients (ADR-0011) or the Jeffreys mean the profile stores
    (ADR-0007). Other notations use the profile's estimate.
    """

    counts = pooled_counts(phenotype.frequency_raw)
    if counts is not None and settings.count_estimator == "beta_shrinkage":
        if settings.shrinkage_mean is None:
            raise ShrinkageMeanMissing(
                "frequency.shrinkage_mean is not set: supply the median count-based frequency "
                "of hpoa-recount-v1 (config or calibration key), or set "
                "frequency.count_estimator: jeffreys for simulator 0.3 estimates"
            )
        return round(
            beta_shrinkage_mean(
                counts[0], counts[1], settings.shrinkage_mean, settings.shrinkage_strength
            ),
            4,
        )
    return phenotype.frequency_estimate


def config_hash(config: SimulationConfig) -> str:
    """Return a stable hash of the simulation configuration."""

    payload = json.dumps(config.model_dump(mode="json"), sort_keys=True)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return f"sha256:{digest[:16]}"


def sex_prior_key(profile: DiseaseProfile) -> SexPriorKey:
    """Which ``sex.p_male`` entry applies to a disease.

    Sex-limited expression and Y-linkage fix the sex. Otherwise the profile's
    inheritance-derived ``sex_bias`` decides, then X-linked dominant
    inheritance (affected females outnumber males about 2:1), else 50/50.
    """

    inheritance = {hpo_id for gene in profile.genes for hpo_id in gene.inheritance_hpo_ids}
    bias = profile.sex_bias.value if profile.sex_bias else "unknown"
    return _sex_key(inheritance, bias)


def entity_sex_prior_key(
    inheritance: Sequence[str], profile: DiseaseProfile, ontology: HpoOntology | None
) -> SexPriorKey:
    """Sex prior of a gene-first case: the entity's inheritance, else the profile's.

    The entity's modes of inheritance (ClinGen when curated) replace the
    profile's, but sex-limited expression is a phenotype fact only HPOA
    records, so the profile's male-/female-limited terms are kept.
    """

    if not inheritance:
        return sex_prior_key(profile)
    profile_ids = {hpo_id for gene in profile.genes for hpo_id in gene.inheritance_hpo_ids}
    ids = set(inheritance) | (profile_ids & set(SEX_LIMITED_EXPRESSION))
    bias, _ = derive_sex_bias(sorted(ids), ontology)
    return _sex_key(ids, bias or "unknown")


def _sex_key(inheritance: set[str], bias: str) -> SexPriorKey:
    male_limited = MALE_LIMITED in inheritance or Y_LINKED in inheritance
    female_limited = FEMALE_LIMITED in inheritance
    if male_limited and not female_limited:
        return "male_limited"
    if female_limited and not male_limited:
        return "female_limited"
    if bias == "male":
        return "male_biased"
    if bias == "female":
        return "female_biased"
    if bias == "unknown" and X_LINKED_DOMINANT in inheritance:
        return "x_linked_dominant"
    return "unbiased"


class SexSpecificTerms:
    """Decides whether a term can be observed in a patient of a given sex."""

    def __init__(self, ontology: HpoOntology | None, settings: SexSettings) -> None:
        male = self._subtrees(ontology, settings.male_only_anchors)
        female = self._subtrees(ontology, settings.female_only_anchors)
        self.male_only = frozenset(male - female)
        self.female_only = frozenset(female - male)

    def allowed(self, hpo_id: str, restriction: SexRestriction | None, sex: Sex) -> bool:
        if restriction is not None and restriction != sex:
            return False
        if sex == "male":
            return hpo_id not in self.female_only
        if sex == "female":
            return hpo_id not in self.male_only
        return True

    @staticmethod
    def _subtrees(ontology: HpoOntology | None, anchors: Sequence[str]) -> set[str]:
        if ontology is None:
            return set()
        terms: set[str] = set()
        for anchor in anchors:
            if ontology.is_valid_hpo_id(anchor):
                terms.add(anchor)
                terms |= ontology.get_descendant_set(anchor)
        return terms


@dataclass(frozen=True)
class _Term:
    phenotype: PhenotypeAssociation
    mean: float
    always_eligible: bool
    onset_window: tuple[float, float] | None
    cardinal: bool = False
    report_weight: float = 1.0


@dataclass(frozen=True)
class _Candidate:
    hpo_id: str
    label: str
    weight: float
    sex_restriction: SexRestriction | None = None
    reason: str | None = None
    source_probability: float | None = None


@dataclass
class _DiseaseModel:
    profile: DiseaseProfile
    terms: list[_Term]
    sex_key: SexPriorKey
    disease_term_ids: frozenset[str]
    confounder_terms: list[_Candidate]
    not_terms: list[_Candidate]
    equivalent_terms: list[_Candidate] = field(default_factory=list)
    gene_other_terms: list[_Candidate] = field(default_factory=list)
    gene_first: bool = False
    profile_ids: tuple[str, ...] | None = None
    annotated_closure: frozenset[str] | None = None


@dataclass
class _Patient:
    sex: Sex
    onset_category: OnsetCategory
    onset_age: float
    age: float
    eligible: list[bool]
    probabilities: list[float]


@dataclass
class _Observation:
    positives: list[CasePhenotype] = field(default_factory=list)
    missing: list[CasePhenotype] = field(default_factory=list)
    unknown: list[CasePhenotype] = field(default_factory=list)
    present_ids: set[str] = field(default_factory=set)


def simulate_cases(
    profile: DiseaseProfile,
    config: SimulationConfig,
    *,
    ontology: HpoOntology | None = None,
    noise_vocabulary: Sequence[NoiseTerm] | None = None,
    confounders: ConfounderIndex | None = None,
    gene_label: int | None = None,
    disease_label: int | None = None,
    profile_version: str | None = None,
    source_versions: Mapping[str, str] | None = None,
    sex_terms: SexSpecificTerms | None = None,
    reporting: Reporting | None = None,
) -> list[SyntheticCase]:
    """Simulate all configured cases for a single disease profile.

    Produces ``cases_per_disease_per_difficulty`` cases for each configured
    difficulty, in difficulty-then-index order. ``confounders`` enables the
    confounder negative source; without an ontology, the ancestor/descendant
    checks reduce to identity checks. ``reporting`` is required when
    ``reporting.mode`` is ``report_model``.
    """

    _check_reporting(config, reporting)
    model = _disease_model(profile, config, confounders, reporting, (profile.disease_id,))
    sex_terms = sex_terms or SexSpecificTerms(ontology, config.sex)
    cfg_hash = config_hash(config)
    cases: list[SyntheticCase] = []
    for difficulty in config.difficulties:
        preset = config.presets[difficulty]
        for index in range(config.cases_per_disease_per_difficulty):
            seed = case_seed(config.seed, profile.disease_id, difficulty, index)
            rng = random.Random(seed)
            case, budget = _simulate_one_case(
                model,
                config=config,
                preset=preset,
                rng=rng,
                ontology=ontology,
                noise_vocabulary=noise_vocabulary or (),
                sex_terms=sex_terms,
                reporting=reporting,
            )
            cases.append(
                SyntheticCase(
                    case_id=_case_id(profile.disease_id, difficulty, index),
                    target=CaseTarget(
                        disease_id=profile.disease_id,
                        disease_name=profile.disease_name,
                        gene=primary_gene(profile),
                        gene_label=gene_label,
                        disease_label=disease_label,
                    ),
                    metadata=_metadata(
                        config, cfg_hash, difficulty, seed, model.sex_key,
                        profile_version, source_versions, *budget,
                    ),
                    **case,
                )
            )
    return cases


def simulate_gene_cases(
    target: GeneTarget,
    profiles: Mapping[str, DiseaseProfile],
    config: SimulationConfig,
    *,
    ontology: HpoOntology | None = None,
    noise_vocabulary: Sequence[NoiseTerm] | None = None,
    confounders: ConfounderIndex | None = None,
    profile_version: str | None = None,
    source_versions: Mapping[str, str] | None = None,
    sex_terms: SexSpecificTerms | None = None,
    reporting: Reporting | None = None,
) -> list[SyntheticCase]:
    """Simulate all configured cases for one GNN gene (gene-first mode).

    ``cases_per_disease_per_difficulty`` cases per difficulty. Each case draws
    an entity of ``target`` by ``sim_weight`` and is labelled with the gene's
    v3 symbol and class index. With ``entity_profiles: merged`` the entity is
    simulated from one profile merged from all its profile ids
    (``target.disease_id`` is the entity, ``target.profile_ids`` the ids
    merged); with ``uniform``, from one of its profiles drawn uniformly. Every
    profile id in ``target`` must be in ``profiles``.
    """

    _check_reporting(config, reporting)
    sex_terms = sex_terms or SexSpecificTerms(ontology, config.sex)
    cfg_hash = config_hash(config)
    models: dict[tuple[str, str], _DiseaseModel] = {}
    entity_weights = [entity.sim_weight for entity in target.entities]
    merged = config.entity_profiles == "merged"
    cases: list[SyntheticCase] = []
    for difficulty in config.difficulties:
        preset = config.presets[difficulty]
        for index in range(config.cases_per_disease_per_difficulty):
            seed = case_seed(config.seed, f"gene:{target.symbol}", difficulty, index)
            rng = random.Random(seed)
            entity = target.entities[_weighted_index(entity_weights, rng)]
            profile_id = (
                "+".join(entity.profile_ids)
                if merged
                else entity.profile_ids[rng.randrange(len(entity.profile_ids))]
            )
            model = models.get((entity.entity, profile_id))
            if model is None:
                model = _gene_model(
                    target, entity, profile_id, profiles, config, confounders, ontology, reporting
                )
                models[(entity.entity, profile_id)] = model
            case, budget = _simulate_one_case(
                model,
                config=config,
                preset=preset,
                rng=rng,
                ontology=ontology,
                noise_vocabulary=noise_vocabulary or (),
                sex_terms=sex_terms,
                reporting=reporting,
            )
            profile = model.profile
            cases.append(
                SyntheticCase(
                    case_id=f"synthetic-gene-{target.symbol}-{difficulty}-{index:06d}",
                    target=CaseTarget(
                        disease_id=profile.disease_id,
                        disease_name=profile.disease_name,
                        gene=target.symbol,
                        gene_label=target.index,
                        entity_id=entity.entity,
                        profile_ids=list(model.profile_ids) if model.profile_ids else None,
                    ),
                    metadata=_metadata(
                        config, cfg_hash, difficulty, seed, model.sex_key,
                        profile_version, source_versions, *budget,
                    ),
                    **case,
                )
            )
    return cases


def entity_profile(
    entity: EntityOption,
    profiles: Mapping[str, DiseaseProfile],
    ontology: HpoOntology | None = None,
) -> DiseaseProfile:
    """The merged profile a gene-first case of ``entity`` is simulated from."""

    merged, _ = merge_entity_profiles(
        entity.entity, [profiles[profile_id] for profile_id in entity.profile_ids], ontology
    )
    return merged


def _check_reporting(config: SimulationConfig, reporting: Reporting | None) -> None:
    if config.reporting.mode == "report_model" and reporting is None:
        raise ValueError(
            "reporting.mode is 'report_model' but no report model was given; pass a "
            "report-model-v1 file, or set reporting.mode: observation for simulator 0.3 emission"
        )


def _metadata(
    config: SimulationConfig,
    cfg_hash: str,
    difficulty: Difficulty,
    seed: int,
    sex_key: SexPriorKey,
    profile_version: str | None,
    source_versions: Mapping[str, str] | None,
    report_budget: int | None = None,
    report_budget_profile: int | None = None,
) -> GeneratorMetadata:
    return GeneratorMetadata(
        generator_version=__version__,
        profile_version=profile_version,
        source_versions=dict(source_versions or {}),
        simulator_version=SIMULATOR_VERSION,
        config_hash=cfg_hash,
        seed=config.seed,
        case_seed=seed,
        sex_prior_key=sex_key,
        report_budget=report_budget,
        report_budget_profile=report_budget_profile,
        difficulty=difficulty,
    )


def case_seed(seed: int, key: str, difficulty: Difficulty, index: int) -> int:
    """Per-case RNG seed, stable across runs and machines."""

    digest = hashlib.sha256(f"{seed}|{key}|{difficulty}|{index}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _disease_model(
    profile: DiseaseProfile,
    config: SimulationConfig,
    confounders: ConfounderIndex | None,
    reporting: Reporting | None = None,
    cardinal_ids: Sequence[str] = (),
    confounder_exclude: frozenset[str] | None = None,
) -> _DiseaseModel:
    windows = config.age.onset_years
    cardinal = reporting.cardinal_terms(cardinal_ids) if reporting is not None else frozenset()
    terms: list[_Term] = []
    for phenotype in profile.phenotypes:
        frequency = simulation_frequency(phenotype, config.frequency)
        is_cardinal = phenotype.hpo_id in cardinal
        terms.append(
            _Term(
                phenotype=phenotype,
                mean=frequency if frequency is not None else config.frequency.unknown_frequency,
                always_eligible=phenotype.onset == "antenatal"
                or phenotype.onset_hpo_id == CONGENITAL_ONSET,
                onset_window=(
                    windows.get(phenotype.onset) if phenotype.onset != "unknown" else None
                ),
                cardinal=is_cardinal,
                report_weight=(
                    reporting.scorer.score(phenotype.hpo_id, frequency, is_cardinal)
                    if reporting is not None
                    else 1.0
                ),
            )
        )
    if confounders is None:
        candidates = []
    elif confounder_exclude is None:
        candidates = confounders.candidate_terms(profile.disease_id)
    else:
        candidates = confounders.candidate_terms_for_profile(profile, confounder_exclude)
    confounder_terms = [
        _Candidate(
            hpo_id=term.hpo_id,
            label=term.label,
            weight=term.weight,
            sex_restriction=term.sex_restriction,
            reason=f"confounder:{term.confounder_id}",
        )
        for term in candidates
    ]
    not_terms = [
        _Candidate(
            hpo_id=negative.hpo_id, label=negative.label, weight=1.0, reason="not_annotation"
        )
        for negative in profile.negative_phenotypes
    ]
    return _DiseaseModel(
        profile=profile,
        terms=terms,
        sex_key=sex_prior_key(profile),
        disease_term_ids=frozenset(phenotype.hpo_id for phenotype in profile.phenotypes),
        confounder_terms=confounder_terms,
        not_terms=not_terms,
    )


def _gene_model(
    target: GeneTarget,
    entity: EntityOption,
    profile_id: str,
    profiles: Mapping[str, DiseaseProfile],
    config: SimulationConfig,
    confounders: ConfounderIndex | None,
    ontology: HpoOntology | None,
    reporting: Reporting | None,
) -> _DiseaseModel:
    """A disease model for a gene-first case: entity sex prior and gene-level pools.

    A merged model (``profile_id`` joins the entity's ids with ``+``) never
    takes a confounder from the entity's own profiles, which describe the
    same disease.
    """

    merged = config.entity_profiles == "merged"
    if merged:
        profile = entity_profile(entity, profiles, ontology)
        model = _disease_model(
            profile, config, confounders, reporting, entity.profile_ids,
            confounder_exclude=frozenset(entity.profile_ids),
        )
    else:
        profile = profiles[profile_id]
        model = _disease_model(profile, config, confounders, reporting, entity.profile_ids)
    unknown = config.frequency.unknown_frequency
    own_terms = _entity_terms(entity, profiles, unknown)
    equivalent = [
        _Candidate(
            hpo_id=hpo_id,
            label=label,
            weight=mean,
            sex_restriction=restriction,
            reason=EQUIVALENT_PROFILE_REASON,
            source_probability=mean,
        )
        for hpo_id, (label, mean, restriction) in sorted(own_terms.items())
        if hpo_id not in model.disease_term_ids
    ]

    other: dict[str, tuple[str, float, float, SexRestriction | None, str]] = {}
    for option in target.entities:
        if option.entity == entity.entity:
            continue
        terms = _entity_terms(option, profiles, unknown)
        total = sum(mean for _, mean, _ in terms.values())
        if total <= 0.0:
            continue
        for hpo_id, (label, mean, restriction) in sorted(terms.items()):
            if hpo_id in own_terms or mean <= 0.0:
                continue
            weight = option.sim_weight * mean / total
            previous = other.get(hpo_id)
            if previous is None:
                other[hpo_id] = (label, weight, mean, restriction, option.entity)
                continue
            best = option.entity if weight > previous[1] else previous[4]
            other[hpo_id] = (
                label,
                previous[1] + weight,
                max(previous[2], mean),
                previous[3] or restriction,
                best,
            )
    gene_other = [
        _Candidate(
            hpo_id=hpo_id,
            label=label,
            weight=weight,
            sex_restriction=restriction,
            reason=f"gene_other_disease:{source_entity}",
            source_probability=mean,
        )
        for hpo_id, (label, weight, mean, restriction, source_entity) in sorted(other.items())
    ]
    return replace(
        model,
        sex_key=entity_sex_prior_key(entity.inheritance, profile, ontology),
        disease_term_ids=model.disease_term_ids | frozenset(own_terms),
        equivalent_terms=equivalent,
        gene_other_terms=gene_other,
        gene_first=True,
        profile_ids=entity.profile_ids if merged else None,
    )


def _entity_terms(
    entity: EntityOption, profiles: Mapping[str, DiseaseProfile], unknown_frequency: float
) -> dict[str, tuple[str, float, SexRestriction | None]]:
    """Union of an entity's profile terms: label, highest mean and any sex restriction."""

    terms: dict[str, tuple[str, float, SexRestriction | None]] = {}
    for profile_id in entity.profile_ids:
        for phenotype in profiles[profile_id].phenotypes:
            mean = (
                phenotype.frequency_estimate
                if phenotype.frequency_estimate is not None
                else unknown_frequency
            )
            previous = terms.get(phenotype.hpo_id)
            if previous is None:
                terms[phenotype.hpo_id] = (phenotype.label, mean, phenotype.sex_restriction)
            else:
                terms[phenotype.hpo_id] = (
                    previous[0],
                    max(previous[1], mean),
                    previous[2] or phenotype.sex_restriction,
                )
    return terms


def _simulate_one_case(
    model: _DiseaseModel,
    *,
    config: SimulationConfig,
    preset: DifficultyPreset,
    rng: random.Random,
    ontology: HpoOntology | None,
    noise_vocabulary: Sequence[NoiseTerm],
    sex_terms: SexSpecificTerms,
    reporting: Reporting | None = None,
) -> tuple[dict[str, object], tuple[int | None, int | None]]:
    """One case's fields, and its report budgets ``(k, profile budget)`` (report mode)."""

    patient = _sample_patient(model, config, rng, sex_terms)
    observation: _Observation | None = None
    budget: int | None = None
    profile_budget: int | None = None
    noise_slots: int | None = None
    if config.reporting.mode == "report_model":
        assert reporting is not None
        truth: list[int] = []
        for _ in range(config.max_redraws):
            if not any(patient.eligible):
                patient = _sample_patient(model, config, rng, sex_terms)
                continue
            truth = _sample_truth(patient, rng)
            if truth:
                break
        budget = reporting.budget.draw(rng)
        # The histogram counts every present term a real record shows, noise
        # included, so noise takes its share of the budget before profile terms.
        showable = len(_noise_pool(model, patient, set(), None, noise_vocabulary, sex_terms))
        if config.reporting.noise_count == "budget_share":
            drawn = budget_share_noise_count(
                budget, config.reporting.noise_share, noise_vocabulary, rng
            )
        else:
            drawn = noise_count(preset, noise_vocabulary, rng)
        noise_slots = min(drawn, showable)
        profile_budget = max(1, budget - noise_slots)
        if truth:
            observation = _report(
                model, patient, truth, profile_budget, config.reporting.force_cardinal, preset,
                rng, ontology, sex_terms, config.reporting.specialize_rate,
            )
    else:
        for _ in range(config.max_redraws):
            if not any(patient.eligible):
                patient = _sample_patient(model, config, rng, sex_terms)
                continue
            observation = _observe(model, patient, preset, rng, ontology)
            if observation.positives:
                break
    if observation is None or not observation.positives:
        patient = _presentable_sex(model, patient, sex_terms)
        observation = _forced_observation(model, patient, sex_terms)

    if noise_slots is None:
        negatives = _sample_negatives(
            model, patient, observation.present_ids, config, preset, rng, ontology, sex_terms
        )
        noise = _sample_noise(
            model, patient, observation.present_ids, negatives, preset, noise_vocabulary, rng,
            ontology, sex_terms,
        )
    else:
        # Noise is part of the reported budget, so it is shown before negatives,
        # which then keep clear of it as of any present term.
        noise = _sample_noise(
            model, patient, observation.present_ids, [], preset, noise_vocabulary, rng,
            ontology, sex_terms, count=noise_slots,
            related_share=config.noise.related_share,
            true_ids=_true_profile_ids(observation),
        )
        negatives = _sample_negatives(
            model, patient, observation.present_ids | {term.hpo_id for term in noise},
            config, preset, rng, ontology, sex_terms,
        )

    hide_sex = rng.random() < config.missingness.sex_unknown
    hide_age = rng.random() < config.missingness.age_unknown
    hide_onset = rng.random() < config.missingness.onset_unknown
    hide_negatives = rng.random() < config.missingness.no_negatives

    onset_age = Age(value=round(patient.onset_age, 2), unit="years")
    case = {
        "patient": PatientAttributes(
            sex="unknown" if hide_sex else patient.sex,
            age=None if hide_age else Age(value=round(patient.age, 2), unit="years"),
            age_of_onset=None if hide_onset else onset_age,
            onset_category="unknown" if hide_onset else patient.onset_category,
        ),
        "positive_phenotypes": observation.positives,
        "negative_phenotypes": [] if hide_negatives else negatives,
        "missing_phenotypes": observation.missing,
        "unknown_phenotypes": observation.unknown,
        "noise_phenotypes": noise,
    }
    return case, (budget, profile_budget)


def _sample_patient(
    model: _DiseaseModel,
    config: SimulationConfig,
    rng: random.Random,
    sex_terms: SexSpecificTerms,
) -> _Patient:
    sex: Sex = "male" if rng.random() < config.sex.p_male[model.sex_key] else "female"
    onset_category = _sample_onset_category(model.profile, config, rng)
    low, high = config.age.onset_years[onset_category]
    onset_age = rng.uniform(low, high)
    duration = min(
        rng.expovariate(1.0 / config.age.duration_mean_years), config.age.duration_max_years
    )
    age = max(onset_age, min(onset_age + duration, config.age.max_age_years))
    boost = 0.0
    if model.profile.progression == "progressive" and config.progression.max_boost > 0.0:
        boost = config.progression.max_boost * (
            1.0 - math.exp(-(age - onset_age) / config.progression.timescale_years)
        )

    eligible: list[bool] = []
    probabilities: list[float] = []
    for term in model.terms:
        probability = _patient_probability(term.mean, config.frequency.concentration, rng)
        if boost and not term.always_eligible:
            probability += (1.0 - probability) * boost
        term_onset = rng.uniform(*term.onset_window) if term.onset_window else None
        eligible.append(
            sex_terms.allowed(term.phenotype.hpo_id, term.phenotype.sex_restriction, sex)
            and (term.always_eligible or term_onset is None or term_onset <= age)
        )
        probabilities.append(probability)
    return _Patient(sex, onset_category, onset_age, age, eligible, probabilities)


def _sample_onset_category(
    profile: DiseaseProfile, config: SimulationConfig, rng: random.Random
) -> OnsetCategory:
    windows = config.age.onset_years
    onset = profile.age_of_onset
    if onset is not None:
        weights = [
            (category, weight)
            for category, weight in onset.distribution.items()
            if category in windows and weight > 0.0
        ]
        if weights:
            return _weighted_choice(weights, rng)
        if onset.category in windows:
            return onset.category
    prior = [
        (category, weight)
        for category, weight in config.age.unknown_onset_prior.items()
        if weight > 0.0
    ]
    return _weighted_choice(prior, rng) if prior else "variable"


def _patient_probability(mean: float, concentration: float, rng: random.Random) -> float:
    if mean <= 0.0:
        return 0.0
    if mean >= 1.0:
        return 1.0
    return rng.betavariate(concentration * mean, concentration * (1.0 - mean))


def _observe(
    model: _DiseaseModel,
    patient: _Patient,
    preset: DifficultyPreset,
    rng: random.Random,
    ontology: HpoOntology | None,
) -> _Observation:
    observation = _Observation()
    observed_ids: set[str] = set()
    # A generalized term must not be one the profile restricts to the other sex.
    other_sex_terms = frozenset(
        t.phenotype.hpo_id
        for t in model.terms
        if t.phenotype.sex_restriction is not None and t.phenotype.sex_restriction != patient.sex
    )
    for index, term in enumerate(model.terms):
        if not patient.eligible[index] or rng.random() >= patient.probabilities[index]:
            continue
        phenotype = term.phenotype
        probability = round(patient.probabilities[index], 4)
        observation.present_ids.add(phenotype.hpo_id)
        observe_rate = preset.positive_observation_rate
        if phenotype.diagnostic_role in CARDINAL_ROLES:
            observe_rate = min(1.0, observe_rate + preset.cardinal_observation_boost)
        if rng.random() < observe_rate:
            hpo_id, label = _maybe_generalize(phenotype, preset, ontology, rng, other_sex_terms)
            if hpo_id in observed_ids:
                # The record already shows this term, so the specific one stays
                # unrecorded; keeping it here keeps the case's truth complete.
                if phenotype.hpo_id not in observed_ids:
                    observation.missing.append(
                        _unobserved(
                            phenotype, probability, status="missing", reason=MERGED_REASON
                        )
                    )
                continue
            observed_ids.add(hpo_id)
            observation.present_ids.add(hpo_id)
            generalized = hpo_id != phenotype.hpo_id
            observation.positives.append(
                CasePhenotype(
                    hpo_id=hpo_id,
                    label=label,
                    status="positive",
                    observed=True,
                    source_probability=probability,
                    source_hpo_id=phenotype.hpo_id if generalized else None,
                    simulated_origin="disease_profile",
                    reason="generalized" if generalized else None,
                )
            )
        elif rng.random() < preset.missing_vs_unknown_split:
            observation.missing.append(
                _unobserved(phenotype, probability, status="missing", reason="not_recorded")
            )
        else:
            observation.unknown.append(
                _unobserved(
                    phenotype, probability, status="unknown", reason="not_enough_information"
                )
            )
    return observation


def _sample_truth(patient: _Patient, rng: random.Random) -> list[int]:
    """Indices of the profile terms the patient truly has, in profile order."""

    return [
        index
        for index, eligible in enumerate(patient.eligible)
        if eligible and rng.random() < patient.probabilities[index]
    ]


def _report(
    model: _DiseaseModel,
    patient: _Patient,
    truth: Sequence[int],
    budget: int,
    force_cardinal: bool,
    preset: DifficultyPreset,
    rng: random.Random,
    ontology: HpoOntology | None,
    sex_terms: SexSpecificTerms | None = None,
    specialize_rate: float = 0.0,
) -> _Observation:
    """Record ``budget`` of the true terms: cardinal ones first, the rest by report score.

    Every true cardinal term is reported when the budget is at least 1, even
    beyond the budget; the remaining slots are filled without replacement with
    probability proportional to each term's reporting-model score. Reported
    terms may be generalized to a parent or, when not generalized, shown as a
    descendant with ``specialize_rate``; unreported ones become ``missing``
    (``not_reported``) or ``unknown`` by the preset's split.
    """

    reported: set[int] = set()
    if force_cardinal and budget > 0:
        reported = {index for index in truth if model.terms[index].cardinal}
    pool = [index for index in truth if index not in reported]
    slots = budget - len(reported)
    while slots > 0 and pool:
        chosen = pool.pop(_weighted_index([model.terms[i].report_weight for i in pool], rng))
        reported.add(chosen)
        slots -= 1

    observation = _Observation()
    observed_ids: set[str] = set()
    other_sex_terms = frozenset(
        t.phenotype.hpo_id
        for t in model.terms
        if t.phenotype.sex_restriction is not None and t.phenotype.sex_restriction != patient.sex
    )

    def emit(index: int) -> bool:
        """Record a reported term; False when its own id is already shown, so it fills no slot."""
        term = model.terms[index]
        phenotype = term.phenotype
        probability = round(patient.probabilities[index], 4)
        hpo_id, label = _maybe_generalize(phenotype, preset, ontology, rng, other_sex_terms)
        if hpo_id in observed_ids:
            if phenotype.hpo_id in observed_ids:
                return False
            observation.missing.append(
                _unobserved(phenotype, probability, status="missing", reason=MERGED_REASON)
            )
            return True
        generalized = hpo_id != phenotype.hpo_id
        reason = "generalized" if generalized else (CARDINAL_REASON if term.cardinal else None)
        if (
            not generalized
            and specialize_rate > 0.0
            and ontology is not None
            and sex_terms is not None
            and rng.random() < specialize_rate
        ):
            descendant = _specialize(
                phenotype.hpo_id, ontology, rng, observed_ids | model.disease_term_ids,
                sex_terms, patient.sex,
            )
            if descendant is not None:
                hpo_id = descendant
                label = ontology.get_label(descendant) or descendant
                reason = SPECIALIZED_REASON
        derived = hpo_id != phenotype.hpo_id
        observed_ids.add(hpo_id)
        observation.present_ids.add(hpo_id)
        observation.positives.append(
            CasePhenotype(
                hpo_id=hpo_id,
                label=label,
                status="positive",
                observed=True,
                source_probability=probability,
                source_hpo_id=phenotype.hpo_id if derived else None,
                simulated_origin="disease_profile",
                reason=reason,
            )
        )
        return True

    absorbed = 0
    for index in truth:
        observation.present_ids.add(model.terms[index].phenotype.hpo_id)
        if index in reported and not emit(index):
            absorbed += 1
    # A term whose own id an earlier generalization already shows leaves its
    # slot open; fill it from the remaining true terms.
    while absorbed and pool:
        chosen = pool.pop(_weighted_index([model.terms[i].report_weight for i in pool], rng))
        reported.add(chosen)
        if emit(chosen):
            absorbed -= 1
    unreported = [index for index in truth if index not in reported]
    for index in unreported:
        phenotype = model.terms[index].phenotype
        if phenotype.hpo_id in observed_ids:
            continue
        probability = round(patient.probabilities[index], 4)
        if rng.random() < preset.missing_vs_unknown_split:
            observation.missing.append(
                _unobserved(phenotype, probability, status="missing", reason=NOT_REPORTED_REASON)
            )
        else:
            observation.unknown.append(
                _unobserved(
                    phenotype, probability, status="unknown", reason="not_enough_information"
                )
            )
    return observation


def _specialize(
    hpo_id: str,
    ontology: HpoOntology,
    rng: random.Random,
    blocked: set[str] | frozenset[str],
    sex_terms: SexSpecificTerms,
    sex: Sex,
) -> str | None:
    """A descendant to show instead of ``hpo_id``: a child, then with p 0.5 a grandchild.

    Candidates are phenotypic abnormalities allowed for the patient's sex and
    not in ``blocked`` (terms already shown and the profile's annotated terms,
    so the pick is more specific than any annotation). None when no child
    qualifies; a child without a qualifying child is kept.
    """

    def eligible(ids: Sequence[str]) -> list[str]:
        return [
            term_id
            for term_id in ids
            if ontology.is_phenotypic_abnormality(term_id)
            and term_id not in blocked
            and sex_terms.allowed(term_id, None, sex)
        ]

    children = eligible(ontology.get_direct_children(hpo_id))
    if not children:
        return None
    child = children[rng.randrange(len(children))]
    if rng.random() < 0.5:
        grandchildren = eligible(ontology.get_direct_children(child))
        if grandchildren:
            child = grandchildren[rng.randrange(len(grandchildren))]
    return child


def _true_profile_ids(observation: _Observation) -> list[str]:
    """The profile terms the patient truly has, sorted."""

    ids = {p.source_hpo_id or p.hpo_id for p in observation.positives}
    ids |= {p.hpo_id for p in observation.missing}
    ids |= {p.hpo_id for p in observation.unknown}
    return sorted(ids)


def _presentable_sex(
    model: _DiseaseModel, patient: _Patient, sex_terms: SexSpecificTerms
) -> _Patient:
    """The patient, switched to the other sex when no profile term fits theirs.

    A disease whose every term is restricted to one sex (e.g. only ovarian
    findings) cannot present in the other, so forcing a term would emit a
    sex-inappropriate finding.
    """

    def admits(sex: Sex) -> bool:
        return any(
            sex_terms.allowed(term.phenotype.hpo_id, term.phenotype.sex_restriction, sex)
            for term in model.terms
        )

    other: Sex = "female" if patient.sex == "male" else "male"
    if admits(patient.sex) or not admits(other):
        return patient
    return replace(patient, sex=other)


def _forced_observation(
    model: _DiseaseModel, patient: _Patient, sex_terms: SexSpecificTerms
) -> _Observation:
    """One observed positive when redraws failed: the most probable admissible term."""

    indices = [index for index, eligible in enumerate(patient.eligible) if eligible]
    if not indices:
        indices = [
            index
            for index, term in enumerate(model.terms)
            if sex_terms.allowed(term.phenotype.hpo_id, term.phenotype.sex_restriction, patient.sex)
        ] or list(range(len(model.terms)))
    chosen = max(indices, key=lambda index: (patient.probabilities[index], -index))
    phenotype = model.terms[chosen].phenotype
    observation = _Observation()
    observation.present_ids.add(phenotype.hpo_id)
    observation.positives.append(
        CasePhenotype(
            hpo_id=phenotype.hpo_id,
            label=phenotype.label,
            status="positive",
            observed=True,
            source_probability=round(patient.probabilities[chosen], 4),
            simulated_origin="disease_profile",
            reason=FORCED_REASON,
        )
    )
    return observation


def _maybe_generalize(
    phenotype: PhenotypeAssociation,
    preset: DifficultyPreset,
    ontology: HpoOntology | None,
    rng: random.Random,
    forbidden: frozenset[str] = frozenset(),
) -> tuple[str, str]:
    if ontology is None or preset.ontology_smoothing_rate <= 0.0:
        return phenotype.hpo_id, phenotype.label
    if rng.random() >= preset.ontology_smoothing_rate:
        return phenotype.hpo_id, phenotype.label
    parents = [
        parent
        for parent in ontology.get_direct_parents(phenotype.hpo_id)
        if ontology.is_phenotypic_abnormality(parent) and parent not in forbidden
    ]
    if not parents:
        return phenotype.hpo_id, phenotype.label
    parent_id = rng.choice(parents)
    return parent_id, ontology.get_label(parent_id) or phenotype.label


def _unobserved(
    phenotype: PhenotypeAssociation, probability: float, *, status: str, reason: str
) -> CasePhenotype:
    return CasePhenotype(
        hpo_id=phenotype.hpo_id,
        label=phenotype.label,
        status=status,  # type: ignore[arg-type]
        observed=False if status == "missing" else None,
        source_probability=probability,
        simulated_origin="disease_profile",
        reason=reason,
    )


class _RelatedTerms:
    """Terms plus their ancestors, to test "equal, ancestor or descendant" quickly."""

    def __init__(self, ontology: HpoOntology | None, hpo_ids: set[str] | None = None) -> None:
        self._ontology = ontology
        self.ids: set[str] = set()
        self._ancestors: set[str] = set()
        for hpo_id in hpo_ids or ():
            self.add(hpo_id)

    def add(self, hpo_id: str) -> None:
        self.ids.add(hpo_id)
        if self._ontology is not None:
            self._ancestors |= self._ontology.get_ancestor_set(hpo_id)

    def related(self, hpo_id: str) -> bool:
        if hpo_id in self.ids or hpo_id in self._ancestors:
            return True
        if self._ontology is None:
            return False
        return not self._ontology.get_ancestor_set(hpo_id).isdisjoint(self.ids)


def _sample_negatives(
    model: _DiseaseModel,
    patient: _Patient,
    present_ids: set[str],
    config: SimulationConfig,
    preset: DifficultyPreset,
    rng: random.Random,
    ontology: HpoOntology | None,
    sex_terms: SexSpecificTerms,
) -> list[CasePhenotype]:
    count = negative_count(
        preset.negatives_mean,
        config.negatives.count_dispersion,
        config.negatives.max_per_case,
        rng,
    )
    if count == 0:
        return []

    def admissible(candidates: Sequence[_Candidate]) -> list[_Candidate]:
        return [
            candidate
            for candidate in candidates
            if candidate.weight > 0.0
            and sex_terms.allowed(candidate.hpo_id, candidate.sex_restriction, patient.sex)
        ]

    own = [
        _Candidate(
            hpo_id=term.phenotype.hpo_id,
            label=term.phenotype.label,
            weight=term.mean,
            sex_restriction=term.phenotype.sex_restriction,
            reason="typical_feature_absent",
            source_probability=term.mean,
        )
        for term in model.terms
        if term.phenotype.hpo_id not in present_ids
    ]
    if config.negatives.own_disease_pool == "entity":
        own.extend(term for term in model.equivalent_terms if term.hpo_id not in present_ids)
    pools: dict[NegativeSource, list[_Candidate]] = {
        "own_disease": admissible(own),
        "own_gene_other_disease": admissible(model.gene_other_terms),
        "confounder": admissible(model.confounder_terms),
        "not_annotation": admissible(model.not_terms),
    }
    # A source the mode cannot provide is left out of the draw; a source that
    # exists but runs dry for this case keeps its slots empty (below).
    weighted_sources = [
        (source, config.negatives.source_weights[source])
        for source in NEGATIVE_SOURCES
        if config.negatives.source_weights[source] > 0.0
        and (model.gene_first or source not in GENE_ONLY_SOURCES)
    ]
    if not weighted_sources:
        return []
    # Slots go to sources before any candidate is drawn, and by default a slot
    # its source cannot fill stays empty, so the realized mix follows the weights.
    slots = {source: 0 for source in NEGATIVE_SOURCES}
    for _ in range(count):
        slots[_weighted_choice(weighted_sources, rng)] += 1

    present = _RelatedTerms(ontology, present_ids)
    chosen = _RelatedTerms(ontology)
    negatives: list[CasePhenotype] = []

    def fill(source: NegativeSource) -> bool:
        pool = pools[source]
        while pool:
            candidate = pool.pop(_weighted_index([item.weight for item in pool], rng))
            if present.related(candidate.hpo_id) or chosen.related(candidate.hpo_id):
                continue
            chosen.add(candidate.hpo_id)
            negatives.append(
                CasePhenotype(
                    hpo_id=candidate.hpo_id,
                    label=candidate.label,
                    status="negative",
                    observed=True,
                    source_probability=(
                        round(candidate.source_probability, 4)
                        if candidate.source_probability is not None
                        else None
                    ),
                    simulated_origin=NEGATIVE_ORIGINS[source],
                    reason=candidate.reason,
                )
            )
            return True
        return False

    for source in NEGATIVE_SOURCES:
        for _ in range(slots[source]):
            if not fill(source):
                break
    if config.negatives.unfilled_slots == "redistribute":
        while len(negatives) < count:
            live = [(source, weight) for source, weight in weighted_sources if pools[source]]
            if not live:
                break
            fill(_weighted_choice(live, rng))
    return negatives


def _sample_noise(
    model: _DiseaseModel,
    patient: _Patient,
    present_ids: set[str],
    negatives: Sequence[CasePhenotype],
    preset: DifficultyPreset,
    noise_vocabulary: Sequence[NoiseTerm],
    rng: random.Random,
    ontology: HpoOntology | None,
    sex_terms: SexSpecificTerms,
    *,
    count: int | None = None,
    related_share: float = 0.0,
    true_ids: Sequence[str] = (),
) -> list[CasePhenotype]:
    """Noise terms; ``count`` is drawn here unless the caller drew it already.

    With ``related_share`` > 0 each slot is, with that probability, a term near
    one of ``true_ids`` (see :func:`_related_noise_term`), falling back to the
    vocabulary when none qualifies.
    """

    if count is None:
        count = noise_count(preset, noise_vocabulary, rng)
    if count == 0:
        return []
    negated = _RelatedTerms(ontology, {negative.hpo_id for negative in negatives})
    pool = _noise_pool(model, patient, present_ids, negated, noise_vocabulary, sex_terms)
    related_on = related_share > 0.0 and ontology is not None and bool(true_ids)
    noise: list[CasePhenotype] = []
    while len(noise) < count:
        if related_on and rng.random() < related_share:
            assert ontology is not None
            related = _related_noise_term(
                model, true_ids, ontology, rng,
                present_ids | {item.hpo_id for item in noise}, negated, sex_terms, patient.sex,
            )
            if related is not None:
                hpo_id, seed = related
                pool = [term for term in pool if term.hpo_id != hpo_id]
                noise.append(
                    CasePhenotype(
                        hpo_id=hpo_id,
                        label=ontology.get_label(hpo_id) or hpo_id,
                        status="noise",
                        observed=True,
                        source_hpo_id=seed,
                        simulated_origin="noise",
                        reason=RELATED_NOISE_REASON,
                    )
                )
                continue
        if not pool:
            break
        term = pool.pop(_weighted_index([item.weight for item in pool], rng))
        noise.append(
            CasePhenotype(
                hpo_id=term.hpo_id,
                label=term.label,
                status="noise",
                observed=True,
                simulated_origin="noise",
                reason=VOCABULARY_NOISE_REASON,
            )
        )
    return noise


def _related_noise_term(
    model: _DiseaseModel,
    true_ids: Sequence[str],
    ontology: HpoOntology,
    rng: random.Random,
    shown: set[str],
    negated: _RelatedTerms,
    sex_terms: SexSpecificTerms,
    sex: Sex,
) -> tuple[str, str] | None:
    """A disease-related noise term and the true term it was drawn near.

    From a uniformly chosen true term, go up to a parent (uniform among
    several) or, with probability 0.5, on to one of that parent's parents;
    then pick uniformly among the ancestor's descendants at most 2 levels
    down. The pick is never the ancestor, a term the profile annotates, an
    ancestor or descendant of one (descendants are specialisation, not
    related noise), a shown term, a term related to a negative, or a term
    the patient's sex rules out. None when nothing qualifies.
    """

    if model.annotated_closure is None:
        closure = set(model.disease_term_ids)
        for hpo_id in model.disease_term_ids:
            closure |= ontology.get_ancestor_set(hpo_id)
        model.annotated_closure = frozenset(closure)
    annotated = frozenset(model.disease_term_ids)
    seed = true_ids[rng.randrange(len(true_ids))]
    parents = ontology.get_direct_parents(seed)
    if not parents:
        return None
    anchor = parents[rng.randrange(len(parents))]
    if rng.random() < 0.5:
        grandparents = ontology.get_direct_parents(anchor)
        if grandparents:
            anchor = grandparents[rng.randrange(len(grandparents))]
    children = ontology.get_direct_children(anchor)
    candidates = set(children)
    for child in children:
        candidates.update(ontology.get_direct_children(child))
    eligible = [
        hpo_id
        for hpo_id in sorted(candidates)
        if hpo_id != anchor
        and ontology.is_phenotypic_abnormality(hpo_id)
        and hpo_id not in model.annotated_closure
        and ontology.get_ancestor_set(hpo_id).isdisjoint(annotated)
        and hpo_id not in shown
        and not negated.related(hpo_id)
        and sex_terms.allowed(hpo_id, None, sex)
    ]
    if not eligible:
        return None
    return eligible[rng.randrange(len(eligible))], seed


def _noise_pool(
    model: _DiseaseModel,
    patient: _Patient,
    present_ids: set[str],
    negated: _RelatedTerms | None,
    noise_vocabulary: Sequence[NoiseTerm],
    sex_terms: SexSpecificTerms,
) -> list[NoiseTerm]:
    return [
        term
        for term in noise_vocabulary
        if term.weight > 0.0
        and term.hpo_id not in model.disease_term_ids
        and term.hpo_id not in present_ids
        and (negated is None or not negated.related(term.hpo_id))
        and sex_terms.allowed(term.hpo_id, None, patient.sex)
    ]


def noise_count(
    preset: DifficultyPreset, noise_vocabulary: Sequence[NoiseTerm], rng: random.Random
) -> int:
    """Poisson count of noise terms with the preset's mean; 0 without a vocabulary."""

    if not noise_vocabulary or preset.noise_mean <= 0.0:
        return 0
    return _poisson(preset.noise_mean, rng)


def budget_share_noise_count(
    budget: int, share: float, noise_vocabulary: Sequence[NoiseTerm], rng: random.Random
) -> int:
    """Binomial(budget, share) noise slots, at most budget - 1; 0 without a vocabulary."""

    if not noise_vocabulary or share <= 0.0 or budget <= 1:
        return 0
    return min(sum(1 for _ in range(budget) if rng.random() < share), budget - 1)


def negative_count(
    mean: float, dispersion: float | None, cap: int, rng: random.Random
) -> int:
    """Gamma-Poisson count with the given mean, capped at ``cap``."""

    if mean <= 0.0 or cap <= 0:
        return 0
    rate = rng.gammavariate(dispersion, mean / dispersion) if dispersion else mean
    return min(_poisson(rate, rng), cap)


def _poisson(mean: float, rng: random.Random) -> int:
    if mean <= 0.0:
        return 0
    if mean > 30.0:
        return max(0, round(rng.gauss(mean, math.sqrt(mean))))
    limit = math.exp(-mean)
    count, product = 0, rng.random()
    while product > limit:
        count += 1
        product *= rng.random()
    return count


def _weighted_index(weights: Sequence[float], rng: random.Random) -> int:
    total = sum(weights)
    threshold = rng.random() * total
    cumulative = 0.0
    for index, weight in enumerate(weights):
        cumulative += weight
        if threshold < cumulative:
            return index
    return len(weights) - 1


def _weighted_choice(items: Sequence[tuple[_T, float]], rng: random.Random) -> _T:
    return items[_weighted_index([weight for _, weight in items], rng)][0]


def primary_gene(profile: DiseaseProfile) -> str:
    """Return the gene a case is labelled with: the first causal gene, else the first gene."""

    for gene in profile.genes:
        if gene.association_type == "causal":
            return gene.symbol
    if profile.genes:
        return profile.genes[0].symbol
    return "unknown"


def _case_id(disease_id: str, difficulty: Difficulty, index: int) -> str:
    slug = disease_id.replace(":", "_")
    return f"synthetic-{slug}-{difficulty}-{index:06d}"

