"""Schemas for LLM-produced disease profile patches."""

from __future__ import annotations

from datetime import date
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

Confidence = float
HpoStatus = Literal["present", "absent", "uncertain"]


class SourceReference(BaseModel):
    """Reference to the source text used for extraction."""

    name: str
    url: str | None = None
    retrieved_at: date | None = None
    section: str | None = None
    license: str | None = None
    redistribution_allowed: bool = False


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
    onset: str = "unknown"
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
