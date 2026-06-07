"""Merge validated LLM patches into disease profiles."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch, PhenotypeExtraction
from rare_disease_simulator.llm_extraction.validator import PatchValidationResult
from rare_disease_simulator.profiles.quality import apply_profile_quality
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    ClinicalSummary,
    DiseaseGene,
    DiseaseProfile,
    FrequencyCategory,
    NegativePhenotypeAssociation,
    OnsetCategory,
    PhenotypeAssociation,
    ProfileQuality,
    ProgressionPattern,
    Provenance,
)

FREQUENCY_VALUES = {
    "obligate",
    "very_frequent",
    "frequent",
    "occasional",
    "very_rare",
    "excluded",
    "unknown",
}
DIAGNOSTIC_ROLE_VALUES = {"cardinal", "major", "supportive", "nonspecific", "unknown"}
ONSET_VALUES = {
    "antenatal",
    "neonatal",
    "infantile",
    "childhood",
    "juvenile",
    "adult",
    "childhood_or_adolescent",
    "variable",
    "unknown",
}
PROGRESSION_VALUES = {"progressive", "non_progressive", "episodic", "variable", "unknown"}


@dataclass(frozen=True)
class MergeResult:
    """Result of merging an accepted patch into a disease profile."""

    profile: DiseaseProfile
    warnings: list[str]


def merge_patch_into_profile(
    profile: DiseaseProfile,
    patch: DiseaseProfilePatch,
    validation: PatchValidationResult,
) -> MergeResult:
    """Merge an accepted LLM patch into a validated disease profile."""

    if not validation.is_accepted:
        raise ValueError(f"cannot merge patch with validation status={validation.status!r}")
    if patch.disease_id != profile.disease_id:
        raise ValueError(
            f"patch disease_id={patch.disease_id!r} does not match "
            f"profile disease_id={profile.disease_id!r}"
        )

    profile_data = profile.model_dump()
    merged = DiseaseProfile.model_validate(profile_data)
    warnings = list(merged.quality.warnings)

    if not _profile_has_gene(merged, patch.gene):
        warnings.append(f"patch gene {patch.gene} not present in profile genes")

    merged.phenotypes = _merge_positive_phenotypes(merged.phenotypes, patch)
    merged.negative_phenotypes = _merge_negative_phenotypes(merged.negative_phenotypes, patch)
    merged.genes = _merge_inheritance(merged.genes, patch)
    merged.age_of_onset = _merge_age_of_onset(merged.age_of_onset, patch, warnings)
    merged.progression = _merge_progression(merged.progression, patch, warnings)
    merged.clinical_summary = _merge_clinical_summary(merged.clinical_summary, patch)
    merged.quality = _merge_quality(merged.quality, patch, warnings)
    merged = apply_profile_quality(merged)

    return MergeResult(profile=merged, warnings=merged.quality.warnings)


def _profile_has_gene(profile: DiseaseProfile, gene: str) -> bool:
    return any(profile_gene.symbol == gene for profile_gene in profile.genes)


def _merge_positive_phenotypes(
    existing: list[PhenotypeAssociation],
    patch: DiseaseProfilePatch,
) -> list[PhenotypeAssociation]:
    by_hpo = {phenotype.hpo_id: phenotype for phenotype in existing}
    ordered_hpo_ids = [phenotype.hpo_id for phenotype in existing]

    for extracted in patch.phenotypes:
        if extracted.hpo_id in by_hpo:
            by_hpo[extracted.hpo_id] = _enrich_positive_phenotype(
                by_hpo[extracted.hpo_id], extracted, patch
            )
        else:
            by_hpo[extracted.hpo_id] = _phenotype_from_extraction(extracted, patch)
            ordered_hpo_ids.append(extracted.hpo_id)

    return [by_hpo[hpo_id] for hpo_id in ordered_hpo_ids]


def _merge_negative_phenotypes(
    existing: list[NegativePhenotypeAssociation],
    patch: DiseaseProfilePatch,
) -> list[NegativePhenotypeAssociation]:
    by_hpo = {phenotype.hpo_id: phenotype for phenotype in existing}
    ordered_hpo_ids = [phenotype.hpo_id for phenotype in existing]

    for extracted in patch.negative_phenotypes:
        provenance = _phenotype_provenance(extracted, field="negative_phenotypes", patch=patch)
        if extracted.hpo_id in by_hpo:
            current = by_hpo[extracted.hpo_id]
            by_hpo[extracted.hpo_id] = current.model_copy(
                update={
                    "source": _append_unique(current.source, patch.source.name),
                    "provenance": [*current.provenance, provenance],
                    "confidence": _max_optional(current.confidence, extracted.confidence),
                }
            )
        else:
            by_hpo[extracted.hpo_id] = NegativePhenotypeAssociation(
                hpo_id=extracted.hpo_id,
                label=extracted.hpo_label,
                source=[patch.source.name],
                provenance=[provenance],
                confidence=extracted.confidence,
            )
            ordered_hpo_ids.append(extracted.hpo_id)

    return [by_hpo[hpo_id] for hpo_id in ordered_hpo_ids]


def _enrich_positive_phenotype(
    current: PhenotypeAssociation,
    extracted: PhenotypeExtraction,
    patch: DiseaseProfilePatch,
) -> PhenotypeAssociation:
    provenance = _phenotype_provenance(extracted, field="phenotypes", patch=patch)
    update = {
        "source": _append_unique(current.source, patch.source.name),
        "provenance": [*current.provenance, provenance],
        "confidence": _max_optional(current.confidence, extracted.confidence),
    }
    if current.diagnostic_role == "unknown":
        update["diagnostic_role"] = _diagnostic_role_or_unknown(extracted.diagnostic_role)
    if current.onset == "unknown":
        update["onset"] = extracted.onset
    if current.frequency == "unknown":
        update["frequency"] = _frequency_or_unknown(extracted.frequency)
    return current.model_copy(update=update)


def _phenotype_from_extraction(
    extracted: PhenotypeExtraction,
    patch: DiseaseProfilePatch,
) -> PhenotypeAssociation:
    return PhenotypeAssociation(
        hpo_id=extracted.hpo_id,
        label=extracted.hpo_label,
        frequency=_frequency_or_unknown(extracted.frequency),
        diagnostic_role=_diagnostic_role_or_unknown(extracted.diagnostic_role),
        onset=extracted.onset,
        source=[patch.source.name],
        provenance=[_phenotype_provenance(extracted, field="phenotypes", patch=patch)],
        confidence=extracted.confidence,
    )


def _merge_inheritance(
    genes: list[DiseaseGene],
    patch: DiseaseProfilePatch,
) -> list[DiseaseGene]:
    inheritance_values = [item.value for item in patch.inheritance]
    if not inheritance_values:
        return genes

    provenance = [
        Provenance(
            source=patch.source,
            field="genes.inheritance",
            evidence=item.evidence_span,
            confidence=item.confidence,
        )
        for item in patch.inheritance
    ]
    merged_genes: list[DiseaseGene] = []
    for gene in genes:
        if gene.symbol != patch.gene:
            merged_genes.append(gene)
            continue
        merged_genes.append(
            gene.model_copy(
                update={
                    "inheritance": _append_many_unique(gene.inheritance, inheritance_values),
                    "provenance": [*gene.provenance, *provenance],
                }
            )
        )
    return merged_genes


def _merge_age_of_onset(
    current: AgeOfOnset | None,
    patch: DiseaseProfilePatch,
    warnings: list[str],
) -> AgeOfOnset | None:
    if patch.age_of_onset is None:
        return current

    category = _onset_or_unknown(patch.age_of_onset.value)
    if category == "unknown" and patch.age_of_onset.value != "unknown":
        warnings.append(f"unknown age_of_onset category: {patch.age_of_onset.value}")

    provenance = Provenance(
        source=patch.source,
        field="age_of_onset",
        evidence=patch.age_of_onset.evidence_span,
        confidence=patch.age_of_onset.confidence,
    )

    if current is None or current.category == "unknown":
        return AgeOfOnset(
            category=category,
            provenance=[provenance],
            confidence=patch.age_of_onset.confidence,
        )
    return current.model_copy(
        update={
            "provenance": [*current.provenance, provenance],
            "confidence": _max_optional(current.confidence, patch.age_of_onset.confidence),
        }
    )


def _merge_progression(
    current: ProgressionPattern,
    patch: DiseaseProfilePatch,
    warnings: list[str],
) -> ProgressionPattern:
    if patch.progression is None:
        return current
    progression = _progression_or_unknown(patch.progression.value)
    if progression == "unknown" and patch.progression.value != "unknown":
        warnings.append(f"unknown progression value: {patch.progression.value}")
    if current == "unknown":
        return progression
    return current


def _merge_clinical_summary(
    current: ClinicalSummary | None,
    patch: DiseaseProfilePatch,
) -> ClinicalSummary | None:
    if patch.clinical_summary is None:
        return current
    provenance = Provenance(
        source=patch.source,
        field="clinical_summary",
        confidence=_mean_optional(_patch_confidences(patch)),
    )
    if current is None:
        return ClinicalSummary(
            text=patch.clinical_summary,
            llm_generated=False,
            provenance=[provenance],
        )
    return current.model_copy(update={"provenance": [*current.provenance, provenance]})


def _merge_quality(
    current: ProfileQuality,
    patch: DiseaseProfilePatch,
    warnings: list[str],
) -> ProfileQuality:
    patch_confidence = _mean_optional(_patch_confidences(patch))
    return current.model_copy(
        update={
            "profile_confidence": _max_optional(current.profile_confidence, patch_confidence),
            "warnings": warnings,
        }
    )


def _phenotype_provenance(
    extracted: PhenotypeExtraction,
    *,
    field: str,
    patch: DiseaseProfilePatch,
) -> Provenance:
    return Provenance(
        source=patch.source,
        field=field,
        evidence=extracted.evidence_span,
        confidence=extracted.confidence,
    )


def _frequency_or_unknown(value: str) -> FrequencyCategory:
    return value if value in FREQUENCY_VALUES else "unknown"  # type: ignore[return-value]


def _diagnostic_role_or_unknown(value: str) -> str:
    return value if value in DIAGNOSTIC_ROLE_VALUES else "unknown"


def _onset_or_unknown(value: str) -> OnsetCategory:
    return value if value in ONSET_VALUES else "unknown"  # type: ignore[return-value]


def _progression_or_unknown(value: str) -> ProgressionPattern:
    return value if value in PROGRESSION_VALUES else "unknown"  # type: ignore[return-value]


def _append_unique(values: list[str], value: str) -> list[str]:
    return values if value in values else [*values, value]


def _append_many_unique(values: list[str], new_values: list[str]) -> list[str]:
    merged = list(values)
    for value in new_values:
        if value not in merged:
            merged.append(value)
    return merged


def _max_optional(current: float | None, new: float | None) -> float | None:
    if current is None:
        return new
    if new is None:
        return current
    return max(current, new)


def _mean_optional(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def _patch_confidences(patch: DiseaseProfilePatch) -> list[float]:
    values = [phenotype.confidence for phenotype in patch.phenotypes]
    values.extend(phenotype.confidence for phenotype in patch.negative_phenotypes)
    values.extend(item.confidence for item in patch.inheritance)
    if patch.age_of_onset is not None:
        values.append(patch.age_of_onset.confidence)
    if patch.progression is not None:
        values.append(patch.progression.confidence)
    return values
