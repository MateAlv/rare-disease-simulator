import pytest

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.llm_extraction.merge_patch import merge_patch_into_profile
from rare_disease_simulator.llm_extraction.validator import DiseaseProfilePatchValidator
from rare_disease_simulator.profiles.schema import (
    DiseaseGene,
    DiseaseProfile,
    MappedDiseaseIds,
    PhenotypeAssociation,
    Provenance,
    SourceReference,
)
from rare_disease_simulator.validation.hpo import HpoValidator
from tests.fixtures.readers import fixture_path, read_profile_patch, read_text_snippets


def _structured_profile() -> DiseaseProfile:
    source = SourceReference(
        name="orphanet_fixture",
        version="fixture-0.1",
        retrieved_at="2026-05-20",
        redistribution_allowed=True,
    )
    provenance = Provenance(source=source, field="profile", confidence=1.0)
    return DiseaseProfile(
        disease_id="ORPHA:646",
        disease_name="Niemann-Pick disease type C1",
        mapped_ids=MappedDiseaseIds(orpha="ORPHA:646", omim=["OMIM:257220"]),
        genes=[
            DiseaseGene(
                symbol="NPC1",
                association_type="causal",
                inheritance=["autosomal recessive"],
                provenance=[Provenance(source=source, field="genes", confidence=1.0)],
            )
        ],
        phenotypes=[
            PhenotypeAssociation(
                hpo_id="HP:0001251",
                label="Ataxia",
                frequency="frequent",
                source=["orphanet_fixture"],
                provenance=[Provenance(source=source, field="phenotypes", confidence=1.0)],
            ),
            PhenotypeAssociation(
                hpo_id="HP:0001433",
                label="Hepatosplenomegaly",
                frequency="occasional",
                source=["orphanet_fixture"],
                provenance=[Provenance(source=source, field="phenotypes", confidence=1.0)],
            ),
            PhenotypeAssociation(
                hpo_id="HP:0000511",
                label="Vertical supranuclear gaze palsy",
                frequency="frequent",
                source=["orphanet_fixture"],
                provenance=[Provenance(source=source, field="phenotypes", confidence=1.0)],
            ),
        ],
        provenance=[provenance],
    )


def _patch_validation():
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")
    validator = DiseaseProfilePatchValidator(HpoValidator(ontology))
    return validator.validate(read_profile_patch(), source_text=read_text_snippets()[0].text)


def test_merge_patch_into_profile_enriches_duplicates_and_adds_new_phenotypes() -> None:
    result = merge_patch_into_profile(
        _structured_profile(),
        read_profile_patch(),
        _patch_validation(),
    )

    profile = result.profile

    assert result.warnings == ["missing_negative_phenotypes"]
    assert [phenotype.hpo_id for phenotype in profile.phenotypes] == [
        "HP:0001251",
        "HP:0001433",
        "HP:0000511",
        "HP:0001250",
    ]

    ataxia = profile.phenotypes[0]
    assert ataxia.frequency == "frequent"
    assert ataxia.diagnostic_role == "major"
    assert ataxia.onset == "childhood"
    assert ataxia.source == ["orphanet_fixture", "synthetic_fixture_text"]
    assert len(ataxia.provenance) == 2

    seizure = profile.phenotypes[-1]
    assert seizure.label == "Seizure"
    assert seizure.source == ["synthetic_fixture_text"]
    assert seizure.confidence == 0.86


def test_merge_patch_into_profile_adds_scalar_clinical_fields() -> None:
    result = merge_patch_into_profile(
        _structured_profile(),
        read_profile_patch(),
        _patch_validation(),
    )
    profile = result.profile

    assert profile.age_of_onset is not None
    assert profile.age_of_onset.category == "childhood_or_adolescent"
    assert profile.progression == "progressive"
    assert profile.clinical_summary is not None
    assert "neurovisceral" in profile.clinical_summary.text
    assert profile.quality.profile_confidence is not None
    assert profile.quality.profile_confidence > 0.85
    assert profile.quality.warnings == ["missing_negative_phenotypes"]


def test_merge_patch_refuses_rejected_validation_result() -> None:
    patch = read_profile_patch()
    phenotype = patch.phenotypes[0].model_copy(update={"hpo_id": "HP:9999999"})
    rejected_patch = patch.model_copy(update={"phenotypes": [phenotype, *patch.phenotypes[1:]]})
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")
    validator = DiseaseProfilePatchValidator(HpoValidator(ontology))
    rejected_validation = validator.validate(
        rejected_patch,
        source_text=read_text_snippets()[0].text,
    )

    with pytest.raises(ValueError, match="cannot merge patch"):
        merge_patch_into_profile(_structured_profile(), rejected_patch, rejected_validation)
