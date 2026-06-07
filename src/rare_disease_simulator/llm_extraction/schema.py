"""Schemas for LLM-produced disease profile patches."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from rare_disease_simulator.profiles.schema import OnsetCategory, SourceReference

Confidence = float
HpoStatus = Literal["present", "absent", "uncertain"]


class TextSnippet(BaseModel):
    """Trusted text snippet prepared for LLM extraction."""

    model_config = ConfigDict(extra="forbid")

    snippet_id: str
    disease_id: str
    disease_name: str
    gene: str
    source: SourceReference
    text: str


class PhenotypeExtraction(BaseModel):
    """HPO phenotype claim extracted from source text."""

    model_config = ConfigDict(extra="forbid")

    mention: str
    hpo_id: str
    hpo_label: str
    status: HpoStatus
    frequency: str = "unknown"
    diagnostic_role: str = "unknown"
    onset: OnsetCategory = "unknown"
    evidence_span: str
    confidence: Confidence = Field(ge=0.0, le=1.0)


class EvidenceScalarExtraction(BaseModel):
    """Scalar clinical claim extracted from source text."""

    model_config = ConfigDict(extra="forbid")

    value: str
    evidence_span: str
    confidence: Confidence = Field(ge=0.0, le=1.0)


class DiseaseProfilePatch(BaseModel):
    """Patch proposed by an LLM extractor for later validation and merge."""

    model_config = ConfigDict(extra="forbid")

    disease_id: str
    disease_name: str
    gene: str
    source: SourceReference
    clinical_summary: str | None = None
    phenotypes: list[PhenotypeExtraction] = Field(default_factory=list)
    negative_phenotypes: list[PhenotypeExtraction] = Field(default_factory=list)
    age_of_onset: EvidenceScalarExtraction | None = None
    inheritance: list[EvidenceScalarExtraction] = Field(default_factory=list)
    progression: EvidenceScalarExtraction | None = None
    warnings: list[str] = Field(default_factory=list)
