"""Seeded probabilistic patient simulator (v0.2).

Generates synthetic clinical cases from a validated :class:`DiseaseProfile`.
The generative model for one case (``docs/README.md`` documents every knob):

1. **Patient.** Sex from the disease's sex prior (``sex_prior_key``). Onset
   category from the profile's onset distribution (or the configured prior when
   unknown), an onset age uniform in that category's window, and a current age
   = onset + an exponential disease duration.
2. **Truth.** Each profile phenotype gets a per-patient probability
   ``p ~ Beta(k*mean, k*(1-mean))`` around its ADR-0007 frequency estimate
   (unknown frequencies use a configured mean), raised with duration in
   progressive diseases. It is eligible only if its sex restriction matches and
   its own onset (sampled in its category window) is <= the current age;
   antenatal and congenital terms are always eligible. Eligible terms are
   present with probability ``p``.
3. **Observation.** The difficulty preset decides which present terms are
   recorded (the rest become ``missing``/``unknown``) and which are generalized
   to a parent term. Cases with no recorded term are redrawn; after
   ``max_redraws`` the most probable eligible term is forced in.
4. **Negatives.** 0..k "asked and absent" terms from three sources: the
   disease's own terms the patient lacks (weighted by frequency), hallmark
   terms of confounder diseases the true disease never annotates, and the
   profile's explicit ``NOT`` annotations. A negated term is never equal to,
   an ancestor of, or a descendant of a present term (or of another negated
   term).
5. **Noise** from an external vocabulary, never related to a negated term.
6. **Covariate missingness** hides sex, age, onset or all negatives.

Every case is fully determined by ``(seed, disease_id, difficulty, index)``
and the configuration, so reruns reproduce identical bytes.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TypeVar

from rare_disease_simulator import __version__
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import (
    DiseaseProfile,
    OnsetCategory,
    PhenotypeAssociation,
    SexRestriction,
)
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.schema import (
    Age,
    CasePhenotype,
    CaseTarget,
    Difficulty,
    DifficultyPreset,
    GeneratorMetadata,
    NegativeSource,
    PatientAttributes,
    Sex,
    SexPriorKey,
    SexSettings,
    SimulationConfig,
    SyntheticCase,
)

SIMULATOR_VERSION = "0.2.0"

CARDINAL_ROLES = {"cardinal", "major"}
CONGENITAL_ONSET = "HP:0003577"
MALE_LIMITED = "HP:0001475"
FEMALE_LIMITED = "HP:0034344"
Y_LINKED = "HP:0001450"
X_LINKED_DOMINANT = "HP:0001423"

NEGATIVE_SOURCES: tuple[NegativeSource, ...] = ("own_disease", "confounder", "not_annotation")
NEGATIVE_ORIGINS: dict[NegativeSource, str] = {
    "own_disease": "negative_own_disease",
    "confounder": "negative_confounder",
    "not_annotation": "negative_not_annotation",
}
FORCED_REASON = "forced_min_one"

_T = TypeVar("_T")


@dataclass(frozen=True)
class NoiseTerm:
    """A nonspecific phenotype available to inject as clinical noise."""

    hpo_id: str
    label: str
    weight: float = 1.0


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
    male_limited = MALE_LIMITED in inheritance or Y_LINKED in inheritance
    female_limited = FEMALE_LIMITED in inheritance
    if male_limited and not female_limited:
        return "male_limited"
    if female_limited and not male_limited:
        return "female_limited"
    bias = profile.sex_bias.value if profile.sex_bias else "unknown"
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
    profile_ids: frozenset[str]
    confounder_terms: list[_Candidate]
    not_terms: list[_Candidate]


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
) -> list[SyntheticCase]:
    """Simulate all configured cases for a single disease profile.

    Produces ``cases_per_disease_per_difficulty`` cases for each configured
    difficulty, in difficulty-then-index order. ``confounders`` enables the
    confounder negative source; without an ontology, the ancestor/descendant
    checks reduce to identity checks.
    """

    model = _disease_model(profile, config, confounders)
    sex_terms = sex_terms or SexSpecificTerms(ontology, config.sex)
    cfg_hash = config_hash(config)
    cases: list[SyntheticCase] = []
    for difficulty in config.difficulties:
        preset = config.presets[difficulty]
        for index in range(config.cases_per_disease_per_difficulty):
            rng = _case_rng(config.seed, profile.disease_id, difficulty, index)
            case = _simulate_one_case(
                model,
                config=config,
                preset=preset,
                rng=rng,
                ontology=ontology,
                noise_vocabulary=noise_vocabulary or (),
                sex_terms=sex_terms,
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
                    metadata=GeneratorMetadata(
                        generator_version=__version__,
                        profile_version=profile_version,
                        source_versions=dict(source_versions or {}),
                        simulator_version=SIMULATOR_VERSION,
                        config_hash=cfg_hash,
                        seed=config.seed,
                        difficulty=difficulty,
                    ),
                    **case,
                )
            )
    return cases


def _case_rng(seed: int, disease_id: str, difficulty: Difficulty, index: int) -> random.Random:
    """Build a per-case RNG that is stable across runs and machines."""

    key = f"{seed}|{disease_id}|{difficulty}|{index}"
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def _disease_model(
    profile: DiseaseProfile, config: SimulationConfig, confounders: ConfounderIndex | None
) -> _DiseaseModel:
    windows = config.age.onset_years
    terms = [
        _Term(
            phenotype=phenotype,
            mean=(
                phenotype.frequency_estimate
                if phenotype.frequency_estimate is not None
                else config.frequency.unknown_frequency
            ),
            always_eligible=phenotype.onset == "antenatal"
            or phenotype.onset_hpo_id == CONGENITAL_ONSET,
            onset_window=windows.get(phenotype.onset) if phenotype.onset != "unknown" else None,
        )
        for phenotype in profile.phenotypes
    ]
    confounder_terms = (
        [
            _Candidate(
                hpo_id=term.hpo_id,
                label=term.label,
                weight=term.weight,
                sex_restriction=term.sex_restriction,
                reason=f"confounder:{term.confounder_id}",
            )
            for term in confounders.candidate_terms(profile.disease_id)
        ]
        if confounders is not None
        else []
    )
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
        profile_ids=frozenset(phenotype.hpo_id for phenotype in profile.phenotypes),
        confounder_terms=confounder_terms,
        not_terms=not_terms,
    )


def _simulate_one_case(
    model: _DiseaseModel,
    *,
    config: SimulationConfig,
    preset: DifficultyPreset,
    rng: random.Random,
    ontology: HpoOntology | None,
    noise_vocabulary: Sequence[NoiseTerm],
    sex_terms: SexSpecificTerms,
) -> dict[str, object]:
    patient = _sample_patient(model, config, rng, sex_terms)
    observation: _Observation | None = None
    for _ in range(config.max_redraws):
        if not any(patient.eligible):
            patient = _sample_patient(model, config, rng, sex_terms)
            continue
        observation = _observe(model, patient, preset, rng, ontology)
        if observation.positives:
            break
    if observation is None or not observation.positives:
        observation = _forced_observation(model, patient, sex_terms)

    negatives = _sample_negatives(
        model, patient, observation.present_ids, config, preset, rng, ontology, sex_terms
    )
    noise = _sample_noise(
        model, patient, observation.present_ids, negatives, preset, noise_vocabulary, rng,
        ontology, sex_terms,
    )

    hide_sex = rng.random() < config.missingness.sex_unknown
    hide_age = rng.random() < config.missingness.age_unknown
    hide_onset = rng.random() < config.missingness.onset_unknown
    hide_negatives = rng.random() < config.missingness.no_negatives

    onset_age = Age(value=round(patient.onset_age, 2), unit="years")
    return {
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
            hpo_id, label = _maybe_generalize(phenotype, preset, ontology, rng)
            if hpo_id in observed_ids:
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
) -> tuple[str, str]:
    if ontology is None or preset.ontology_smoothing_rate <= 0.0:
        return phenotype.hpo_id, phenotype.label
    if rng.random() >= preset.ontology_smoothing_rate:
        return phenotype.hpo_id, phenotype.label
    parents = [
        parent
        for parent in ontology.get_direct_parents(phenotype.hpo_id)
        if ontology.is_phenotypic_abnormality(parent)
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
    count = min(_poisson(preset.negatives_mean, rng), config.negatives.max_per_case)
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
    pools: dict[NegativeSource, list[_Candidate]] = {
        "own_disease": admissible(own),
        "confounder": admissible(model.confounder_terms),
        "not_annotation": admissible(model.not_terms),
    }
    source_weights = config.negatives.source_weights
    present = _RelatedTerms(ontology, present_ids)
    chosen = _RelatedTerms(ontology)
    negatives: list[CasePhenotype] = []
    while len(negatives) < count:
        sources = [
            (source, source_weights[source])
            for source in NEGATIVE_SOURCES
            if pools[source] and source_weights[source] > 0.0
        ]
        if not sources:
            break
        source = _weighted_choice(sources, rng)
        pool = pools[source]
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
) -> list[CasePhenotype]:
    if not noise_vocabulary or preset.noise_mean <= 0.0:
        return []
    count = _poisson(preset.noise_mean, rng)
    if count == 0:
        return []
    negated = _RelatedTerms(ontology, {negative.hpo_id for negative in negatives})
    pool = [
        term
        for term in noise_vocabulary
        if term.weight > 0.0
        and term.hpo_id not in model.profile_ids
        and term.hpo_id not in present_ids
        and not negated.related(term.hpo_id)
        and sex_terms.allowed(term.hpo_id, None, patient.sex)
    ]
    noise: list[CasePhenotype] = []
    while pool and len(noise) < count:
        term = pool.pop(_weighted_index([item.weight for item in pool], rng))
        noise.append(
            CasePhenotype(
                hpo_id=term.hpo_id,
                label=term.label,
                status="noise",
                observed=True,
                simulated_origin="noise",
                reason="nonspecific_finding",
            )
        )
    return noise


def _poisson(mean: float, rng: random.Random) -> int:
    if mean <= 0.0:
        return 0
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

