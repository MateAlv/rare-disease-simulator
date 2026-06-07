"""Synthetic case and simulation config schemas."""

from __future__ import annotations

from datetime import datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from rare_disease_simulator.profiles.schema import OnsetCategory, Provenance

Difficulty = Literal["easy", "medium", "hard"]
Sex = Literal["female", "male", "other", "unknown"]
AgeUnit = Literal["days", "months", "years"]
PhenotypeObservationStatus = Literal["positive", "negative", "missing", "unknown", "noise"]


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
    difficulty: Difficulty
    generated_at: datetime
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


class SimulationConfig(StrictBaseModel):
    """Simulator configuration."""

    cases_per_disease_per_difficulty: int = Field(gt=0)
    difficulties: list[Difficulty]
    seed: int
    positive_observation_rate: float = Field(default=0.75, ge=0.0, le=1.0)
    known_negative_rate: float = Field(default=0.10, ge=0.0, le=1.0)
    missingness_rate: float = Field(default=0.25, ge=0.0, le=1.0)
    ontology_smoothing_rate: float = Field(default=0.15, ge=0.0, le=1.0)
