from rare_disease_simulator.profiles.quality import (
    ProfileQualityPolicy,
    apply_profile_quality,
    evaluate_profile_quality,
)
from rare_disease_simulator.profiles.schema import (
    DiseaseGene,
    DiseaseProfile,
    PhenotypeAssociation,
    Provenance,
    SourceReference,
)
from tests.fixtures.readers import read_disease_profile


def test_profile_quality_warns_without_blocking_sparse_profile() -> None:
    profile = DiseaseProfile(
        disease_id="ORPHA:TEST",
        disease_name="Sparse disease",
        phenotypes=[
            PhenotypeAssociation(
                hpo_id="HP:0001250",
                label="Seizure",
                confidence=0.4,
            )
        ],
    )

    quality = evaluate_profile_quality(profile)

    assert quality.profile_confidence == 0.4
    assert quality.warnings == [
        "low_phenotype_count",
        "missing_gene",
        "missing_structured_provenance",
        "missing_hpo_provenance",
        "missing_negative_phenotypes",
        "missing_age_of_onset",
        "low_profile_confidence",
    ]


def test_profile_quality_can_disable_optional_missing_field_warnings() -> None:
    profile = read_disease_profile()

    quality = evaluate_profile_quality(
        profile,
        ProfileQualityPolicy(
            warn_missing_negative_phenotypes=False,
            warn_missing_age_of_onset=False,
        ),
    )

    assert "missing_negative_phenotypes" not in quality.warnings
    assert "missing_age_of_onset" not in quality.warnings


def test_apply_profile_quality_returns_profile_copy_with_quality_metadata() -> None:
    source = SourceReference(name="fixture", redistribution_allowed=True)
    profile = DiseaseProfile(
        disease_id="ORPHA:TEST",
        disease_name="Enough disease",
        genes=[
            DiseaseGene(
                symbol="GENE1",
                provenance=[Provenance(source=source, field="genes", confidence=0.9)],
            )
        ],
        phenotypes=[
            PhenotypeAssociation(
                hpo_id="HP:0001250",
                label="Seizure",
                provenance=[Provenance(source=source, field="phenotypes", confidence=0.9)],
            ),
            PhenotypeAssociation(
                hpo_id="HP:0001251",
                label="Ataxia",
                provenance=[Provenance(source=source, field="phenotypes", confidence=0.9)],
            ),
            PhenotypeAssociation(
                hpo_id="HP:0000511",
                label="Vertical supranuclear gaze palsy",
                provenance=[Provenance(source=source, field="phenotypes", confidence=0.9)],
            ),
        ],
        provenance=[Provenance(source=source, field="profile", confidence=0.9)],
    )

    updated = apply_profile_quality(profile)

    assert updated is not profile
    assert updated.quality.profile_confidence == 0.9
    assert updated.quality.warnings == [
        "missing_negative_phenotypes",
        "missing_age_of_onset",
    ]

