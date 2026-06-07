"""Disease profile schema definitions."""

from __future__ import annotations

from datetime import date, datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

FrequencyCategory = Literal[
    "obligate",
    "very_frequent",
    "frequent",
    "occasional",
    "very_rare",
    "excluded",
    "unknown",
]
DiagnosticRole = Literal["cardinal", "major", "supportive", "nonspecific", "unknown"]
OnsetCategory = Literal[
    "antenatal",
    "neonatal",
    "infantile",
    "childhood",
    "juvenile",
    "adult",
    "childhood_or_adolescent",
    "variable",
    "unknown",
]
SexBiasValue = Literal["male", "female", "none", "unknown"]
GeneAssociationType = Literal["causal", "susceptibility", "unknown"]
ProgressionPattern = Literal["progressive", "non_progressive", "episodic", "variable", "unknown"]


class StrictBaseModel(BaseModel):
    """Base model that rejects undeclared fields."""

    model_config = ConfigDict(extra="forbid")


class SourceReference(StrictBaseModel):
    """Reference to a source artifact or text section."""

    name: str
    url: str | None = None
    url_or_file: str | None = None
    retrieved_at: date | None = None
    version: str | None = None
    section: str | None = None
    license: str | None = None
    redistribution_allowed: bool = False


class Provenance(StrictBaseModel):
    """Provenance for a source-derived or generated claim."""

    source: SourceReference
    field: str | None = None
    evidence: str | None = None
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    generated_at: datetime | None = None


class MappedDiseaseIds(StrictBaseModel):
    """Cross-source disease identifiers."""

    orpha: str | None = None
    omim: list[str] = Field(default_factory=list)
    mondo: list[str] = Field(default_factory=list)
    medgen: list[str] = Field(default_factory=list)


class DiseaseGene(StrictBaseModel):
    """Gene associated with a disease profile."""

    symbol: str
    association_type: GeneAssociationType = "unknown"
    inheritance: list[str] = Field(default_factory=list)
    provenance: list[Provenance] = Field(default_factory=list)


class ProbabilityRange(StrictBaseModel):
    """Approximate probability range implied by a phenotype frequency category."""

    lower: float = Field(ge=0.0, le=1.0)
    upper: float = Field(ge=0.0, le=1.0)


class PhenotypeAssociation(StrictBaseModel):
    """Positive disease-phenotype association."""

    hpo_id: str
    label: str
    frequency: FrequencyCategory = "unknown"
    probability_range: ProbabilityRange | None = None
    diagnostic_role: DiagnosticRole = "unknown"
    onset: OnsetCategory = "unknown"
    severity: str = "unknown"
    temporality: str = "unknown"
    source: list[str] = Field(default_factory=list)
    provenance: list[Provenance] = Field(default_factory=list)
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)


class NegativePhenotypeAssociation(StrictBaseModel):
    """Explicitly absent disease-phenotype association."""

    hpo_id: str
    label: str
    source: list[str] = Field(default_factory=list)
    provenance: list[Provenance] = Field(default_factory=list)
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)


class AgeOfOnset(StrictBaseModel):
    """Typical age of onset for a disease profile."""

    category: OnsetCategory = "unknown"
    distribution: dict[OnsetCategory, float] = Field(default_factory=dict)
    provenance: list[Provenance] = Field(default_factory=list)
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)


class SexBias(StrictBaseModel):
    """Sex-related disease presentation pattern."""

    value: SexBiasValue = "unknown"
    provenance: list[Provenance] = Field(default_factory=list)
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)


class ClinicalSummary(StrictBaseModel):
    """Short clinical summary with provenance."""

    text: str
    llm_generated: bool = False
    provenance: list[Provenance] = Field(default_factory=list)


class ProfileQuality(StrictBaseModel):
    """Profile quality metadata and warnings."""

    profile_confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    warnings: list[str] = Field(default_factory=list)
    fixture: bool = False


class DiseaseProfile(StrictBaseModel):
    """Validated disease-gene-phenotype profile used by the simulator."""

    disease_id: str
    disease_name: str
    mapped_ids: MappedDiseaseIds = Field(default_factory=MappedDiseaseIds)
    synonyms: list[str] = Field(default_factory=list)
    genes: list[DiseaseGene] = Field(default_factory=list)
    phenotypes: list[PhenotypeAssociation] = Field(default_factory=list)
    negative_phenotypes: list[NegativePhenotypeAssociation] = Field(default_factory=list)
    age_of_onset: AgeOfOnset | None = None
    sex_bias: SexBias | None = None
    progression: ProgressionPattern = "unknown"
    clinical_summary: ClinicalSummary | None = None
    provenance: list[Provenance] = Field(default_factory=list)
    quality: ProfileQuality = Field(default_factory=ProfileQuality)
