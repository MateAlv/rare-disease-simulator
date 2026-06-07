"""Disease profile quality checks."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_simulator.profiles.schema import DiseaseProfile, ProfileQuality


@dataclass(frozen=True)
class ProfileQualityPolicy:
    """Warning-only quality thresholds for MVP profiles."""

    min_positive_phenotypes: int = 3
    min_profile_confidence: float = 0.75
    warn_missing_negative_phenotypes: bool = True
    warn_missing_age_of_onset: bool = True


def evaluate_profile_quality(
    profile: DiseaseProfile,
    policy: ProfileQualityPolicy | None = None,
) -> ProfileQuality:
    """Return warning-only quality metadata for a disease profile."""

    active_policy = policy or ProfileQualityPolicy()
    warnings = list(profile.quality.warnings)

    if len(profile.phenotypes) < active_policy.min_positive_phenotypes:
        warnings.append("low_phenotype_count")

    if not profile.genes:
        warnings.append("missing_gene")

    if not profile.provenance:
        warnings.append("missing_structured_provenance")

    if not any(phenotype.provenance for phenotype in profile.phenotypes):
        warnings.append("missing_hpo_provenance")

    if active_policy.warn_missing_negative_phenotypes and not profile.negative_phenotypes:
        warnings.append("missing_negative_phenotypes")

    if active_policy.warn_missing_age_of_onset and profile.age_of_onset is None:
        warnings.append("missing_age_of_onset")

    profile_confidence = compute_profile_confidence(profile)
    if profile_confidence is not None and profile_confidence < active_policy.min_profile_confidence:
        warnings.append("low_profile_confidence")

    return ProfileQuality(
        profile_confidence=profile_confidence,
        warnings=_deduplicate_preserving_order(warnings),
        fixture=profile.quality.fixture,
    )


def apply_profile_quality(
    profile: DiseaseProfile,
    policy: ProfileQualityPolicy | None = None,
) -> DiseaseProfile:
    """Return a copy of the profile with refreshed quality metadata."""

    return profile.model_copy(update={"quality": evaluate_profile_quality(profile, policy)})


def compute_profile_confidence(profile: DiseaseProfile) -> float | None:
    """Compute a simple explainable profile confidence from known claim confidences."""

    confidences: list[float] = []

    for phenotype in profile.phenotypes:
        if phenotype.confidence is not None:
            confidences.append(phenotype.confidence)
        confidences.extend(
            provenance.confidence
            for provenance in phenotype.provenance
            if provenance.confidence is not None
        )

    for phenotype in profile.negative_phenotypes:
        if phenotype.confidence is not None:
            confidences.append(phenotype.confidence)
        confidences.extend(
            provenance.confidence
            for provenance in phenotype.provenance
            if provenance.confidence is not None
        )

    for gene in profile.genes:
        confidences.extend(
            provenance.confidence
            for provenance in gene.provenance
            if provenance.confidence is not None
        )

    if profile.age_of_onset is not None:
        if profile.age_of_onset.confidence is not None:
            confidences.append(profile.age_of_onset.confidence)
        confidences.extend(
            provenance.confidence
            for provenance in profile.age_of_onset.provenance
            if provenance.confidence is not None
        )

    if profile.sex_bias is not None:
        if profile.sex_bias.confidence is not None:
            confidences.append(profile.sex_bias.confidence)
        confidences.extend(
            provenance.confidence
            for provenance in profile.sex_bias.provenance
            if provenance.confidence is not None
        )

    if profile.clinical_summary is not None:
        confidences.extend(
            provenance.confidence
            for provenance in profile.clinical_summary.provenance
            if provenance.confidence is not None
        )

    confidences.extend(
        provenance.confidence
        for provenance in profile.provenance
        if provenance.confidence is not None
    )

    if not confidences:
        return None
    return sum(confidences) / len(confidences)


def _deduplicate_preserving_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    deduplicated: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        deduplicated.append(value)
    return deduplicated
