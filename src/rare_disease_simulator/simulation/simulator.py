"""Seeded probabilistic patient simulator.

Generates synthetic clinical cases from a validated :class:`DiseaseProfile`.

The generative model for one case is:

1. Sample patient sex (respecting any disease sex bias) and ages.
2. For each profile phenotype, draw a per-patient prevalence from its frequency
   category / probability range and decide via Bernoulli whether the phenotype
   is *truly present* in this patient.
3. Degrade the true picture into an observed picture using the difficulty
   preset: some true positives are observed, the rest split into ``missing`` /
   ``unknown``; explicit profile negatives are asked at the known-negative rate;
   nonspecific noise terms are added; observed specifics are optionally
   generalized to ontology ancestors.

Every case is fully determined by ``(seed, disease_id, difficulty, index)``, so
reruns with the same configuration reproduce identical output.
"""

from __future__ import annotations

import hashlib
import json
import random
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime

from rare_disease_simulator import __version__
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import (
    DiseaseProfile,
    FrequencyCategory,
    OnsetCategory,
    PhenotypeAssociation,
)
from rare_disease_simulator.simulation.difficulty import DifficultyPreset, get_difficulty_preset
from rare_disease_simulator.simulation.schema import (
    Age,
    CasePhenotype,
    CaseTarget,
    Difficulty,
    GeneratorMetadata,
    PatientAttributes,
    Sex,
    SimulationConfig,
    SyntheticCase,
)

SIMULATOR_VERSION = "0.1.0"

FREQUENCY_TO_RANGE: dict[FrequencyCategory, tuple[float, float]] = {
    "obligate": (0.95, 1.0),
    "very_frequent": (0.80, 0.99),
    "frequent": (0.30, 0.79),
    "occasional": (0.05, 0.29),
    "very_rare": (0.01, 0.04),
    "excluded": (0.0, 0.0),
    "unknown": (0.05, 0.50),
}

CARDINAL_ROLES = {"cardinal", "major"}

# Approximate patient age windows (in years) implied by an onset category. The
# simulator samples an onset age within the window and a current age at or after
# onset.
ONSET_TO_AGE_RANGE_YEARS: dict[OnsetCategory, tuple[float, float]] = {
    "antenatal": (0.0, 0.0),
    "neonatal": (0.0, 0.1),
    "infantile": (0.1, 1.0),
    "childhood": (1.0, 10.0),
    "juvenile": (10.0, 16.0),
    "adult": (16.0, 70.0),
    "childhood_or_adolescent": (1.0, 16.0),
    "variable": (0.0, 60.0),
    "unknown": (0.0, 60.0),
}

# A patient can live well beyond onset; current age is onset age plus a sampled
# duration capped here to keep ages plausible.
MAX_YEARS_AFTER_ONSET = 25.0


@dataclass(frozen=True)
class NoiseTerm:
    """A nonspecific phenotype available to inject as clinical noise."""

    hpo_id: str
    label: str


def config_hash(config: SimulationConfig) -> str:
    """Return a stable hash of the simulation configuration."""

    payload = json.dumps(config.model_dump(mode="json"), sort_keys=True)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return f"sha256:{digest[:16]}"


def simulate_cases(
    profile: DiseaseProfile,
    config: SimulationConfig,
    *,
    ontology: HpoOntology | None = None,
    noise_vocabulary: Sequence[NoiseTerm] | None = None,
    gene_label: int | None = None,
    disease_label: int | None = None,
    profile_version: str | None = None,
) -> list[SyntheticCase]:
    """Simulate all configured cases for a single disease profile.

    Produces ``cases_per_disease_per_difficulty`` cases for each configured
    difficulty, in difficulty-then-index order.
    """

    cfg_hash = config_hash(config)
    generated_at = datetime.now(tz=UTC)
    cases: list[SyntheticCase] = []
    for difficulty in config.difficulties:
        preset = get_difficulty_preset(difficulty)
        for index in range(config.cases_per_disease_per_difficulty):
            cases.append(
                _simulate_one_case(
                    profile,
                    config=config,
                    preset=preset,
                    index=index,
                    config_hash_value=cfg_hash,
                    generated_at=generated_at,
                    ontology=ontology,
                    noise_vocabulary=noise_vocabulary,
                    gene_label=gene_label,
                    disease_label=disease_label,
                    profile_version=profile_version,
                )
            )
    return cases


def _case_rng(seed: int, disease_id: str, difficulty: Difficulty, index: int) -> random.Random:
    """Build a per-case RNG that is stable across runs and machines."""

    key = f"{seed}|{disease_id}|{difficulty}|{index}"
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def _simulate_one_case(
    profile: DiseaseProfile,
    *,
    config: SimulationConfig,
    preset: DifficultyPreset,
    index: int,
    config_hash_value: str,
    generated_at: datetime,
    ontology: HpoOntology | None,
    noise_vocabulary: Sequence[NoiseTerm] | None,
    gene_label: int | None,
    disease_label: int | None,
    profile_version: str | None,
) -> SyntheticCase:
    rng = _case_rng(config.seed, profile.disease_id, preset.name, index)

    sex = _sample_sex(profile, rng)
    onset_category = profile.age_of_onset.category if profile.age_of_onset else "unknown"
    onset_age, current_age = _sample_ages(onset_category, rng)

    present_hpo_ids = {
        phenotype.hpo_id
        for phenotype in profile.phenotypes
        if _is_truly_present(phenotype, config, rng)
    }

    positives: list[CasePhenotype] = []
    missing: list[CasePhenotype] = []
    unknown: list[CasePhenotype] = []
    observed_positive_ids: set[str] = set()

    for phenotype in profile.phenotypes:
        if phenotype.hpo_id not in present_hpo_ids:
            continue
        observe_prob = preset.positive_observation_rate
        if phenotype.diagnostic_role in CARDINAL_ROLES:
            observe_prob = min(1.0, observe_prob + preset.cardinal_observation_boost)
        if rng.random() < observe_prob:
            hpo_id, label = _maybe_generalize(phenotype, preset, ontology, rng)
            if hpo_id in observed_positive_ids:
                continue
            observed_positive_ids.add(hpo_id)
            positives.append(
                CasePhenotype(
                    hpo_id=hpo_id,
                    label=label,
                    status="positive",
                    observed=True,
                    source_probability=_phenotype_probability(phenotype),
                    simulated_origin="disease_profile",
                    reason="generalized" if hpo_id != phenotype.hpo_id else None,
                )
            )
        elif rng.random() < preset.missing_vs_unknown_split:
            missing.append(
                _unobserved_phenotype(phenotype, status="missing", reason="not_recorded")
            )
        else:
            unknown.append(
                _unobserved_phenotype(
                    phenotype, status="unknown", reason="not_enough_information"
                )
            )

    negatives = _sample_negatives(profile, preset, observed_positive_ids, rng)
    noise = _sample_noise(
        preset,
        noise_vocabulary,
        excluded_ids=observed_positive_ids
        | {phenotype.hpo_id for phenotype in negatives}
        | {phenotype.hpo_id for phenotype in profile.phenotypes},
        rng=rng,
    )

    return SyntheticCase(
        case_id=_case_id(profile.disease_id, preset.name, index),
        target=CaseTarget(
            disease_id=profile.disease_id,
            disease_name=profile.disease_name,
            gene=primary_gene(profile),
            gene_label=gene_label,
            disease_label=disease_label,
        ),
        patient=PatientAttributes(
            sex=sex,
            age=Age(value=round(current_age, 1), unit="years"),
            age_of_onset=Age(value=round(onset_age, 1), unit="years"),
            onset_category=onset_category,
        ),
        positive_phenotypes=positives,
        negative_phenotypes=negatives,
        missing_phenotypes=missing,
        unknown_phenotypes=unknown,
        noise_phenotypes=noise,
        metadata=GeneratorMetadata(
            generator_version=__version__,
            profile_version=profile_version,
            simulator_version=SIMULATOR_VERSION,
            config_hash=config_hash_value,
            seed=config.seed,
            difficulty=preset.name,
            generated_at=generated_at,
        ),
    )


def _is_truly_present(
    phenotype: PhenotypeAssociation, config: SimulationConfig, rng: random.Random
) -> bool:
    return rng.random() < _phenotype_probability(phenotype)


def _phenotype_probability(phenotype: PhenotypeAssociation) -> float:
    if phenotype.probability_range is not None:
        lower = phenotype.probability_range.lower
        upper = phenotype.probability_range.upper
    else:
        lower, upper = FREQUENCY_TO_RANGE[phenotype.frequency]
    return (lower + upper) / 2.0


def _sample_sex(profile: DiseaseProfile, rng: random.Random) -> Sex:
    bias = profile.sex_bias.value if profile.sex_bias else "unknown"
    if bias == "male":
        return "male" if rng.random() < 0.85 else "female"
    if bias == "female":
        return "female" if rng.random() < 0.85 else "male"
    return "male" if rng.random() < 0.5 else "female"


def _sample_ages(onset_category: OnsetCategory, rng: random.Random) -> tuple[float, float]:
    low, high = ONSET_TO_AGE_RANGE_YEARS.get(onset_category, ONSET_TO_AGE_RANGE_YEARS["unknown"])
    onset_age = rng.uniform(low, high)
    current_age = onset_age + rng.uniform(0.0, MAX_YEARS_AFTER_ONSET)
    return onset_age, current_age


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
    parent_label = ontology.get_label(parent_id) or phenotype.label
    return parent_id, parent_label


def _unobserved_phenotype(
    phenotype: PhenotypeAssociation, *, status: str, reason: str
) -> CasePhenotype:
    return CasePhenotype(
        hpo_id=phenotype.hpo_id,
        label=phenotype.label,
        status=status,  # type: ignore[arg-type]
        observed=False if status == "missing" else None,
        source_probability=_phenotype_probability(phenotype),
        simulated_origin="disease_profile",
        reason=reason,
    )


def _sample_negatives(
    profile: DiseaseProfile,
    preset: DifficultyPreset,
    observed_positive_ids: set[str],
    rng: random.Random,
) -> list[CasePhenotype]:
    negatives: list[CasePhenotype] = []
    for negative in profile.negative_phenotypes:
        if negative.hpo_id in observed_positive_ids:
            continue
        if rng.random() < preset.known_negative_rate:
            negatives.append(
                CasePhenotype(
                    hpo_id=negative.hpo_id,
                    label=negative.label,
                    status="negative",
                    observed=True,
                    simulated_origin="asked_negative",
                    reason="explicitly_absent_in_simulated_visit",
                )
            )
    return negatives


def _sample_noise(
    preset: DifficultyPreset,
    noise_vocabulary: Sequence[NoiseTerm] | None,
    excluded_ids: set[str],
    rng: random.Random,
) -> list[CasePhenotype]:
    if not noise_vocabulary or preset.noise_term_rate <= 0.0:
        return []
    candidates = [term for term in noise_vocabulary if term.hpo_id not in excluded_ids]
    noise: list[CasePhenotype] = []
    for term in candidates:
        if rng.random() < preset.noise_term_rate:
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
