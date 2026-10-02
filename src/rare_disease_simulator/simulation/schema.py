"""Synthetic case and simulation config schemas."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Literal, get_args

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from rare_disease_simulator.profiles.schema import OnsetCategory, Provenance
from rare_disease_simulator.simulation.difficulty import (
    DIFFICULTY_PRESETS,
    Difficulty,
    DifficultyPreset,
)

__all__ = ["Difficulty", "DifficultyPreset"]

Sex = Literal["female", "male", "other", "unknown"]
AgeUnit = Literal["days", "months", "years"]
OnsetAgeSource = Literal["literature", "category"]
PhenotypeObservationStatus = Literal["positive", "negative", "missing", "unknown", "noise"]
SexPriorKey = Literal[
    "male_limited",
    "female_limited",
    "male_biased",
    "female_biased",
    "x_linked_dominant",
    "unbiased",
]
NegativeSource = Literal["own_disease", "own_gene_other_disease", "confounder", "not_annotation"]
OwnDiseasePool = Literal["profile", "entity"]
UnfilledSlots = Literal["drop", "redistribute"]
CountEstimator = Literal["beta_shrinkage", "jeffreys"]
EntityProfiles = Literal["merged", "uniform"]
EmissionMode = Literal["report_model", "independent", "observation"]


class StrictBaseModel(BaseModel):
    """Base model that rejects undeclared fields."""

    model_config = ConfigDict(extra="forbid")


class Age(StrictBaseModel):
    """Patient age or age of onset."""

    value: float
    unit: AgeUnit


class CaseTarget(StrictBaseModel):
    """Training/evaluation target labels for a synthetic case."""

    disease_id: str
    disease_name: str
    gene: str
    gene_label: int | None = None
    disease_label: int | None = None
    entity_id: str | None = Field(
        default=None,
        description="Gene-first mode: the disease entity drawn for the gene; "
        "disease_id is the profile drawn within it (uniform) or the entity (merged).",
    )
    profile_ids: list[str] | None = Field(
        default=None,
        description="Merged gene-first mode: the entity's profile ids the case was "
        "simulated from, merged into one profile.",
    )


class PatientAttributes(StrictBaseModel):
    """Synthetic patient-level attributes."""

    sex: Sex = "unknown"
    age: Age | None = None
    age_of_onset: Age | None = None
    onset_category: OnsetCategory = "unknown"
    region: str = "unknown"


class CasePhenotype(StrictBaseModel):
    """Phenotype observation in a synthetic case."""

    hpo_id: str
    label: str
    status: PhenotypeObservationStatus
    observed: bool | None = None
    source_probability: float | None = Field(default=None, ge=0.0, le=1.0)
    report_probability: float | None = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Independent reporting: the term's report probability q = min(1, "
        "score / frequency), on every true profile term, reported or not.",
    )
    source_hpo_id: str | None = Field(
        default=None,
        description="Profile term this entry was generalized or specialized from; for "
        "related noise, the true term it was drawn near.",
    )
    simulated_origin: str
    reason: str | None = None


class GeneratorMetadata(StrictBaseModel):
    """Reproducibility metadata for generated cases."""

    generator_version: str
    profile_version: str | None = None
    source_versions: dict[str, str] = Field(default_factory=dict)
    prompt_versions: dict[str, str] = Field(default_factory=dict)
    llm_provider: str | None = None
    llm_model: str | None = None
    simulator_version: str
    config_hash: str
    seed: int
    case_seed: int | None = Field(
        default=None, description="Seed of this case's RNG, derived from the run seed and case key."
    )
    sex_prior_key: SexPriorKey | None = Field(
        default=None, description="Sex prior the case was drawn from (sex.p_male key)."
    )
    report_budget: int | None = Field(
        default=None,
        ge=0,
        description="Report-model emission: the term budget k drawn for the case, which "
        "counts profile terms and noise together.",
    )
    report_budget_profile: int | None = Field(
        default=None,
        ge=1,
        description="Report-model emission: the budget left for profile terms, "
        "max(1, k - noise terms drawn).",
    )
    presentation_age_years: tuple[float, float] | None = Field(
        default=None,
        description="Presentation-age range [low, high] the case's age was drawn from "
        "(presentation-ages file); unset when the onset-plus-duration rule set it.",
    )
    onset_age_years: tuple[float, float] | None = Field(
        default=None,
        description="Onset-age range [low, high] the case's onset was drawn from "
        "(onset-ages file); unset when the onset category's window set it.",
    )
    onset_age_source: OnsetAgeSource | None = Field(
        default=None,
        description="Where the onset age came from when an onset-ages file was in use: "
        "'literature' (a row matched) or 'category'; unset without the file.",
    )
    report_count: int | None = Field(
        default=None,
        ge=0,
        description="Budget-normalized independent reporting: true terms whose roll reported "
        "them (1 when forced_min_one).",
    )
    report_scale: float | None = Field(
        default=None,
        ge=0.0,
        description="Budget-normalized independent reporting: the factor c with "
        "sum(min(1, c * q)) = min(k, true terms).",
    )
    difficulty: Difficulty
    generated_at: datetime | None = Field(
        default=None,
        description="Left unset by the simulator so identical runs give identical bytes; "
        "the run summary records the wall-clock time.",
    )
    provenance: list[Provenance] = Field(default_factory=list)


class SyntheticCase(StrictBaseModel):
    """Synthetic patient case generated from a validated disease profile."""

    case_id: str
    target: CaseTarget
    patient: PatientAttributes
    positive_phenotypes: list[CasePhenotype] = Field(default_factory=list)
    negative_phenotypes: list[CasePhenotype] = Field(default_factory=list)
    missing_phenotypes: list[CasePhenotype] = Field(default_factory=list)
    unknown_phenotypes: list[CasePhenotype] = Field(default_factory=list)
    noise_phenotypes: list[CasePhenotype] = Field(default_factory=list)
    metadata: GeneratorMetadata


DEFAULT_P_MALE: dict[SexPriorKey, float] = {
    "male_limited": 1.0,
    "female_limited": 0.0,
    "male_biased": 0.9,
    "female_biased": 0.1,
    "x_linked_dominant": 0.33,
    "unbiased": 0.5,
}


class FrequencySettings(StrictBaseModel):
    """Per-patient phenotype probability around the profile's frequency estimate."""

    concentration: float = Field(
        default=8.0,
        gt=0.0,
        description="Beta concentration k: p ~ Beta(k*mean, k*(1-mean)); larger is tighter.",
    )
    unknown_frequency: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Mean used for terms whose frequency is unknown (ADR-0007: never 'always').",
    )
    count_estimator: CountEstimator = Field(
        default="beta_shrinkage",
        description="Point estimate of a count-based frequency n/m: 'beta_shrinkage' "
        "(ADR-0011: (n + s*mean) / (m + s)) or 'jeffreys' ((n + 0.5) / (m + 1), simulator 0.3).",
    )
    shrinkage_mean: float | None = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Prior mean of the Beta shrinkage: the median count-based frequency of "
        "hpoa-recount-v1, supplied by diagnostic.ar-training. Required when count-based "
        "terms are simulated with 'beta_shrinkage'.",
    )
    shrinkage_strength: float = Field(
        default=2.0,
        gt=0.0,
        description="Prior strength of the Beta shrinkage, in pseudo-patients (ADR-0011: 2).",
    )


class SexSettings(StrictBaseModel):
    """Patient sex prior and sex-specific phenotype anchors."""

    p_male: dict[SexPriorKey, float] = Field(
        default_factory=lambda: dict(DEFAULT_P_MALE),
        description="P(male) per disease sex prior (see simulator.sex_prior_key).",
    )
    male_only_anchors: list[str] = Field(
        default_factory=lambda: ["HP:0010461", "HP:0012874"],
        description="Terms under these (and not under a female anchor) are male-only.",
    )
    female_only_anchors: list[str] = Field(
        default_factory=lambda: ["HP:0010460", "HP:0030012"],
        description="Terms under these (and not under a male anchor) are female-only.",
    )

    @field_validator("p_male")
    @classmethod
    def _complete_probabilities(cls, value: dict[str, float]) -> dict[str, float]:
        merged = {**DEFAULT_P_MALE, **value}
        if any(not 0.0 <= probability <= 1.0 for probability in merged.values()):
            raise ValueError("p_male values must lie in [0, 1]")
        return merged


class AgeSettings(StrictBaseModel):
    """Onset and current-age sampling."""

    onset_years: dict[OnsetCategory, tuple[float, float]] = Field(
        default_factory=lambda: {
            "antenatal": (0.0, 0.0),
            "neonatal": (0.0, 0.0767),
            "infantile": (0.0767, 1.0),
            "childhood": (1.0, 5.0),
            "juvenile": (5.0, 16.0),
            "childhood_or_adolescent": (1.0, 16.0),
            "adult": (16.0, 60.0),
            "variable": (0.0, 60.0),
        },
        description="Uniform onset-age window in years per onset category (HPO definitions).",
    )
    unknown_onset_prior: dict[OnsetCategory, float] = Field(
        default_factory=lambda: {
            "antenatal": 0.11,
            "neonatal": 0.35,
            "infantile": 0.18,
            "childhood": 0.14,
            "juvenile": 0.05,
            "adult": 0.13,
            "variable": 0.04,
        },
        description="Onset category mix for diseases without onset data.",
    )
    duration_mean_years: float = Field(
        default=5.0, gt=0.0, description="Mean of the exponential disease duration."
    )
    duration_mean_by_onset: dict[OnsetCategory, float] = Field(
        default_factory=dict,
        description="Mean exponential duration per onset category; a category listed here "
        "overrides duration_mean_years for patients with that onset.",
    )
    duration_max_years: float = Field(default=40.0, ge=0.0)
    max_age_years: float = Field(default=90.0, gt=0.0)

    @model_validator(mode="after")
    def _check_windows(self) -> AgeSettings:
        for category, (low, high) in self.onset_years.items():
            if not 0.0 <= low <= high:
                raise ValueError(f"onset_years[{category}] must satisfy 0 <= low <= high")
        missing = [
            category
            for category in self.unknown_onset_prior
            if category not in self.onset_years
        ]
        if missing:
            raise ValueError(f"unknown_onset_prior uses categories without a window: {missing}")
        bad = [c for c, mean in self.duration_mean_by_onset.items() if not mean > 0.0]
        if bad:
            raise ValueError(f"duration_mean_by_onset must be positive: {bad}")
        return self


class ProgressionSettings(StrictBaseModel):
    """Raise non-congenital phenotype probability with duration in progressive diseases."""

    max_boost: float = Field(
        default=0.25,
        ge=0.0,
        le=1.0,
        description="p' = p + (1-p) * max_boost * (1 - exp(-duration / timescale)); 0 disables.",
    )
    timescale_years: float = Field(default=10.0, gt=0.0)


class NegativeSettings(StrictBaseModel):
    """Asked-and-absent terms; the per-case count mean lives in the difficulty preset.

    The count is Gamma-Poisson (negative binomial) with the preset mean, capped
    at ``max_per_case``. Own-disease terms dominate the default mix because
    real case reports mostly list the true syndrome's typical signs a patient
    lacks (diagnostic.ar-training EXP-B-001, Addendum 1). The defaults are not
    fitted to any cohort; calibrate them on R1-train only.

    ``own_gene_other_disease`` (the gene's other diseases) only exists in
    gene-first mode; disease-first runs leave it out of the slot draw.
    """

    max_per_case: int = Field(default=15, ge=0)
    count_dispersion: float | None = Field(
        default=2.0,
        gt=0.0,
        description="Gamma shape of the count; smaller is more spread, None is plain Poisson.",
    )
    source_weights: dict[NegativeSource, float] = Field(
        default_factory=lambda: {
            "own_disease": 0.7,
            "own_gene_other_disease": 0.1,
            "confounder": 0.2,
            "not_annotation": 0.1,
        },
        description="Relative odds of each source per negative slot; 0 disables a source.",
    )
    own_disease_pool: OwnDiseasePool = Field(
        default="entity",
        description="Gene-first mode: 'entity' also offers the terms of the entity's other "
        "(equivalent OMIM/ORPHA) profiles as own-disease negatives; 'profile' only the "
        "drawn profile's terms.",
    )
    unfilled_slots: UnfilledSlots = Field(
        default="drop",
        description="A slot whose source has no admissible term left: 'drop' leaves it empty "
        "(the mix follows the weights, the count falls short); 'redistribute' refills it from "
        "the sources that still have terms, by weight (the count follows the draw).",
    )
    confounders_top_n: int = Field(default=10, ge=0)
    min_information_content: float = Field(
        default=2.0,
        ge=0.0,
        description="Only terms with IC >= this (nats) enter the confounder similarity index.",
    )

    @field_validator("source_weights")
    @classmethod
    def _complete_weights(cls, value: dict[str, float]) -> dict[str, float]:
        merged = {source: 0.0 for source in get_args(NegativeSource)}
        merged.update(value)
        if any(weight < 0.0 for weight in merged.values()):
            raise ValueError("source_weights must be non-negative")
        return merged


class ReportingSettings(StrictBaseModel):
    """Which true terms a record shows (ADR-0011 "truth, then reporting")."""

    mode: EmissionMode = Field(
        default="report_model",
        description="'report_model': a term budget and a reporting model pick the reported "
        "terms from the true ones (needs a report-model file); 'independent': each true term "
        "is reported on its own with probability min(1, score / frequency) (simulator 0.5, "
        "needs a report-model file); 'observation': each true term is recorded with the "
        "preset's observation rate (simulator 0.3).",
    )
    force_cardinal: bool = Field(
        default=True,
        description="Report every true cardinal term (cardinal-v1) first, when the budget "
        "is at least 1 (report_model mode; ignored in independent mode).",
    )
    noise_count: Literal["poisson", "budget_share", "proportional"] = Field(
        default="poisson",
        description="How many noise terms a record shows: 'poisson' draws the preset's "
        "noise_mean; 'budget_share' (report_model mode) makes each of the k slots noise with "
        "probability noise_share, at most k - 1; 'proportional' (independent mode) draws "
        "Poisson(m * noise_share / (1 - noise_share)) for m reported profile terms.",
    )
    noise_share: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Share of a real record's present terms that the knowledge base does not "
        "explain; used by noise_count 'budget_share' and 'proportional' (< 1).",
    )
    specialize_rate: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Probability that a reported, not generalized profile term is shown as a "
        "descendant: a child, then with probability 0.5 a grandchild (ADR-0011 amendment 3).",
    )
    q_scale: float = Field(
        default=1.0, gt=0.0,
        description="Independent mode: q = min(1, q_scale * score / frequency); calibrated so "
        "reported explained terms per case match real patients.",
    )
    budget_normalize: bool = Field(
        default=False,
        description="Independent mode: draw a profile-term budget k per case from "
        "profile_budget_histogram and scale every true term's q by one factor c so the "
        "expected reported count is min(k, true terms); terms are still rolled independently.",
    )
    budget_scope: Literal["profile", "record"] = Field(
        default="profile",
        description="What the budget k counts with budget_normalize: 'profile' the reported "
        "profile terms, noise added on top (0.5.2); 'record' every present term, so noise "
        "(noise_count 'budget_share') takes n of the k slots and profile terms aim at k - n.",
    )
    profile_budget_histogram: dict[str, int] | None = Field(
        default=None,
        description="Present terms per real case, {str(k): count}: profile terms with "
        "budget_scope 'profile', all present terms with 'record' (report-model-v2's "
        "term_budget); required by budget_normalize. The 0 bin is ignored and the rest "
        "renormalised.",
    )

    @model_validator(mode="after")
    def _check_noise_count(self) -> ReportingSettings:
        record = self.budget_scope == "record"
        if record and not (self.mode == "independent" and self.budget_normalize):
            raise ValueError("budget_scope 'record' needs independent mode with budget_normalize")
        if record and self.noise_count != "budget_share":
            raise ValueError("budget_scope 'record' takes noise from the budget: use noise_count "
                             "'budget_share'")
        if self.mode == "independent" and self.noise_count == "budget_share" and not record:
            raise ValueError("noise_count 'budget_share' needs a budget: in independent mode use "
                             "budget_scope 'record', or 'proportional' or 'poisson'")
        if self.mode == "report_model" and self.noise_count == "proportional":
            raise ValueError("noise_count 'proportional' is for independent mode")
        if self.noise_count == "proportional" and self.noise_share >= 1.0:
            raise ValueError("noise_share must be below 1 for noise_count 'proportional'")
        if self.budget_normalize:
            if self.mode != "independent":
                raise ValueError("budget_normalize applies only to independent mode")
            if self.profile_budget_histogram is None:
                raise ValueError("budget_normalize needs profile_budget_histogram")
        if self.profile_budget_histogram is not None:
            self.profile_budget()
        return self

    def profile_budget(self) -> dict[int, int]:
        """The profile-budget histogram with integer keys, checked."""

        histogram: dict[int, int] = {}
        for key, count in (self.profile_budget_histogram or {}).items():
            try:
                k = int(key)
            except ValueError as exc:
                raise ValueError(f"profile_budget_histogram key {key!r} is not an int") from exc
            if k < 0 or count < 0:
                raise ValueError("profile_budget_histogram keys and counts must be >= 0")
            histogram[k] = count
        if sum(count for k, count in histogram.items() if k >= 1) <= 0:
            raise ValueError("profile_budget_histogram has no case with at least 1 term")
        return histogram


class NoiseSettings(StrictBaseModel):
    """Where noise terms come from (report-model emission)."""

    related_share: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Probability that a noise slot is a disease-related term (a descendant, "
        "at most 2 levels down, of a true term's parent or grandparent) rather than a "
        "vocabulary term (ADR-0011 amendment 3).",
    )
    related_up_levels: int = Field(
        default=2, ge=1,
        description="Related noise anchor: a true term's parent, then one more level up with "
        "probability 0.5 per level, at most this many levels (2: parent or grandparent; "
        "1: parent only).",
    )
    related_down_levels: int = Field(
        default=2, ge=1,
        description="Related noise candidates: the anchor's descendants at most this many "
        "levels down (1: its children only).",
    )
    related_max_ic: float | None = Field(
        default=None, ge=0.0,
        description="Skip related-noise candidates whose information content (the report "
        "model's information_content feature) is above this; null keeps all.",
    )


class MissingnessSettings(StrictBaseModel):
    """Share of cases whose covariates or negatives are withheld."""

    sex_unknown: float = Field(default=0.10, ge=0.0, le=1.0)
    age_unknown: float = Field(default=0.20, ge=0.0, le=1.0)
    onset_unknown: float = Field(default=0.30, ge=0.0, le=1.0)
    no_negatives: float = Field(default=0.20, ge=0.0, le=1.0)


class SimulationConfig(StrictBaseModel):
    """Simulator configuration; every knob is used by ``simulation/simulator.py``."""

    cases_per_disease_per_difficulty: int = Field(default=100, gt=0)
    difficulties: list[Difficulty] = Field(default_factory=lambda: ["easy", "medium", "hard"])
    seed: int = 42
    presets: dict[Difficulty, DifficultyPreset] = Field(
        default_factory=lambda: dict(DIFFICULTY_PRESETS),
        description="Per-difficulty presets; difficulties not listed keep the built-in preset.",
    )
    frequency: FrequencySettings = Field(default_factory=FrequencySettings)
    sex: SexSettings = Field(default_factory=SexSettings)
    age: AgeSettings = Field(default_factory=AgeSettings)
    progression: ProgressionSettings = Field(default_factory=ProgressionSettings)
    negatives: NegativeSettings = Field(default_factory=NegativeSettings)
    missingness: MissingnessSettings = Field(default_factory=MissingnessSettings)
    entity_profiles: EntityProfiles = Field(
        default="merged",
        description="Gene-first mode: 'merged' simulates an entity from one profile merged "
        "from all its profile ids (ADR-0011); 'uniform' draws one of them per case "
        "(simulator 0.3).",
    )
    reporting: ReportingSettings = Field(default_factory=ReportingSettings)
    noise: NoiseSettings = Field(default_factory=NoiseSettings)
    max_redraws: int = Field(
        default=20,
        ge=0,
        description="Redraws before forcing one observed positive into a case.",
    )

    @field_validator("presets")
    @classmethod
    def _complete_presets(
        cls, value: dict[Difficulty, DifficultyPreset]
    ) -> dict[Difficulty, DifficultyPreset]:
        return {**DIFFICULTY_PRESETS, **value}


def run_config(data: dict[str, Any]) -> SimulationConfig:
    """The config a run summary records, reading a pre-0.4 summary with 0.3 semantics.

    A summary written before simulator 0.4 has no ``entity_profiles``,
    ``reporting`` or ``frequency.count_estimator``; its cases were made with
    uniform profile draws, observation emission and Jeffreys counts, so those
    are filled in instead of the 0.4 defaults.
    """

    config = dict(data)
    config.setdefault("entity_profiles", "uniform")
    config.setdefault("reporting", {"mode": "observation"})
    frequency = dict(config.get("frequency") or {})
    frequency.setdefault("count_estimator", "jeffreys")
    config["frequency"] = frequency
    return SimulationConfig.model_validate(config)
