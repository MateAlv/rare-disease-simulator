"""Validation for LLM extraction outputs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from rare_disease_simulator.llm_extraction.schema import (
    DiseaseProfilePatch,
    EvidenceScalarExtraction,
    PhenotypeExtraction,
)
from rare_disease_simulator.validation.hpo import HpoValidator

IssueSeverity = Literal["warning", "review", "error"]
ValidationStatus = Literal["accepted", "needs_review", "rejected"]


@dataclass(frozen=True)
class PatchValidationIssue:
    """Validation issue found in a DiseaseProfilePatch."""

    severity: IssueSeverity
    field: str
    reason: str
    hpo_id: str | None = None
    message: str | None = None


@dataclass(frozen=True)
class PatchValidationPolicy:
    """Thresholds for accepting LLM-extracted claims."""

    accept_confidence: float = 0.75
    review_confidence: float = 0.50
    require_evidence_in_source_text: bool = True


@dataclass(frozen=True)
class PatchValidationResult:
    """Structured validation result for a DiseaseProfilePatch."""

    status: ValidationStatus
    issues: list[PatchValidationIssue]

    @property
    def is_accepted(self) -> bool:
        """Return whether the patch can be merged automatically."""

        return self.status == "accepted"


class DiseaseProfilePatchValidator:
    """Validate LLM-produced profile patches before merge."""

    def __init__(
        self,
        hpo_validator: HpoValidator,
        policy: PatchValidationPolicy | None = None,
    ) -> None:
        self.hpo_validator = hpo_validator
        self.policy = policy or PatchValidationPolicy()

    def validate(
        self,
        patch: DiseaseProfilePatch,
        *,
        source_text: str | None = None,
    ) -> PatchValidationResult:
        """Validate a patch and classify it as accepted, review, or rejected."""

        issues: list[PatchValidationIssue] = []
        for index, phenotype in enumerate(patch.phenotypes):
            issues.extend(
                self._validate_phenotype(
                    phenotype,
                    field=f"phenotypes[{index}]",
                    expected_status="present",
                    source_text=source_text,
                )
            )

        for index, phenotype in enumerate(patch.negative_phenotypes):
            issues.extend(
                self._validate_phenotype(
                    phenotype,
                    field=f"negative_phenotypes[{index}]",
                    expected_status="absent",
                    source_text=source_text,
                )
            )

        if patch.age_of_onset is not None:
            issues.extend(
                self._validate_scalar(
                    patch.age_of_onset,
                    field="age_of_onset",
                    source_text=source_text,
                )
            )
        if patch.progression is not None:
            issues.extend(
                self._validate_scalar(
                    patch.progression,
                    field="progression",
                    source_text=source_text,
                )
            )
        for index, inheritance in enumerate(patch.inheritance):
            issues.extend(
                self._validate_scalar(
                    inheritance,
                    field=f"inheritance[{index}]",
                    source_text=source_text,
                )
            )

        return PatchValidationResult(status=_status_from_issues(issues), issues=issues)

    def _validate_phenotype(
        self,
        phenotype: PhenotypeExtraction,
        *,
        field: str,
        expected_status: Literal["present", "absent"],
        source_text: str | None,
    ) -> list[PatchValidationIssue]:
        issues: list[PatchValidationIssue] = []

        if phenotype.status != expected_status:
            issues.append(
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="wrong_phenotype_status",
                    hpo_id=phenotype.hpo_id,
                    message=f"expected status {expected_status!r}, got {phenotype.status!r}",
                )
            )

        if not self.hpo_validator.is_valid_hpo_id(phenotype.hpo_id):
            issues.append(
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="unknown_hpo_id",
                    hpo_id=phenotype.hpo_id,
                )
            )
        elif not self.hpo_validator.is_phenotypic_abnormality(phenotype.hpo_id):
            issues.append(
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="not_phenotypic_abnormality",
                    hpo_id=phenotype.hpo_id,
                )
            )
        else:
            known_label = self.hpo_validator.get_label(phenotype.hpo_id)
            if known_label and known_label.lower() != phenotype.hpo_label.lower():
                issues.append(
                    PatchValidationIssue(
                        severity="warning",
                        field=field,
                        reason="hpo_label_mismatch",
                        hpo_id=phenotype.hpo_id,
                        message=f"expected label {known_label!r}, got {phenotype.hpo_label!r}",
                    )
                )

        issues.extend(
            self._validate_confidence(
                phenotype.confidence,
                field=field,
                hpo_id=phenotype.hpo_id,
            )
        )
        issues.extend(
            self._validate_evidence_span(
                phenotype.evidence_span,
                field=field,
                source_text=source_text,
                hpo_id=phenotype.hpo_id,
            )
        )
        return issues

    def _validate_scalar(
        self,
        scalar: EvidenceScalarExtraction,
        *,
        field: str,
        source_text: str | None,
    ) -> list[PatchValidationIssue]:
        issues = self._validate_confidence(scalar.confidence, field=field)
        issues.extend(
            self._validate_evidence_span(
                scalar.evidence_span,
                field=field,
                source_text=source_text,
            )
        )
        return issues

    def _validate_confidence(
        self,
        confidence: float,
        *,
        field: str,
        hpo_id: str | None = None,
    ) -> list[PatchValidationIssue]:
        if confidence < self.policy.review_confidence:
            return [
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="confidence_below_reject_threshold",
                    hpo_id=hpo_id,
                    message=f"confidence={confidence}",
                )
            ]
        if confidence < self.policy.accept_confidence:
            return [
                PatchValidationIssue(
                    severity="review",
                    field=field,
                    reason="confidence_requires_review",
                    hpo_id=hpo_id,
                    message=f"confidence={confidence}",
                )
            ]
        return []

    def _validate_evidence_span(
        self,
        evidence_span: str,
        *,
        field: str,
        source_text: str | None,
        hpo_id: str | None = None,
    ) -> list[PatchValidationIssue]:
        if not evidence_span.strip():
            return [
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="missing_evidence_span",
                    hpo_id=hpo_id,
                )
            ]
        if (
            source_text is not None
            and self.policy.require_evidence_in_source_text
            and evidence_span.lower() not in source_text.lower()
        ):
            return [
                PatchValidationIssue(
                    severity="error",
                    field=field,
                    reason="evidence_span_not_found_in_source_text",
                    hpo_id=hpo_id,
                )
            ]
        return []


def _status_from_issues(issues: list[PatchValidationIssue]) -> ValidationStatus:
    if any(issue.severity == "error" for issue in issues):
        return "rejected"
    if any(issue.severity == "review" for issue in issues):
        return "needs_review"
    return "accepted"
