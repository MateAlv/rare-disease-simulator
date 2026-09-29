import pytest
from pydantic import ValidationError

from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.schema_export import core_json_schemas
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase
from tests.fixtures.readers import read_disease_profile, read_profile_patch, read_synthetic_case


def test_disease_profile_round_trip_serialization() -> None:
    profile = read_disease_profile()

    reloaded = DiseaseProfile.model_validate_json(profile.model_dump_json())

    assert reloaded.disease_id == "ORPHA:646"
    assert reloaded.mapped_ids.omim == ["OMIM:257220"]
    assert reloaded.phenotypes[0].provenance[0].source.name == "orphanet_fixture"


def test_profile_patch_round_trip_serialization() -> None:
    patch = read_profile_patch()

    reloaded = DiseaseProfilePatch.model_validate_json(patch.model_dump_json())

    assert reloaded.phenotypes[0].evidence_span
    assert reloaded.phenotypes[0].confidence == 0.92
    assert reloaded.inheritance[0].evidence_span


def test_synthetic_case_round_trip_serialization() -> None:
    synthetic_case = read_synthetic_case()

    reloaded = SyntheticCase.model_validate_json(synthetic_case.model_dump_json())

    assert reloaded.target.gene == "NPC1"
    assert reloaded.metadata.seed == 42
    assert reloaded.positive_phenotypes[0].status == "positive"


def test_simulation_config_schema_accepts_mvp_defaults() -> None:
    config = SimulationConfig(
        cases_per_disease_per_difficulty=100,
        difficulties=["easy", "medium", "hard"],
        seed=42,
    )

    assert set(config.presets) == {"easy", "medium", "hard"}
    assert config.frequency.unknown_frequency == 0.5
    assert config.missingness.sex_unknown == 0.10


def test_core_json_schema_export_contains_public_models() -> None:
    schemas = core_json_schemas()

    assert set(schemas) == {
        "DiseaseProfile",
        "DiseaseProfilePatch",
        "TextSnippet",
        "SyntheticCase",
        "SimulationConfig",
    }
    assert schemas["DiseaseProfile"]["title"] == "DiseaseProfile"
    assert schemas["SyntheticCase"]["title"] == "SyntheticCase"


def test_disease_profile_new_optional_fields_default_to_empty() -> None:
    profile = read_disease_profile()

    phenotype = profile.phenotypes[0]
    assert phenotype.sex_restriction is None
    assert phenotype.onset_hpo_id is None
    assert profile.quality.counters == {}
    assert profile.genes[0].ncbi_gene_id is None


def test_phenotype_sex_restriction_accepts_only_male_or_female() -> None:
    payload = read_disease_profile().model_dump(mode="json")
    payload["phenotypes"][0]["sex_restriction"] = "female"
    assert DiseaseProfile.model_validate(payload).phenotypes[0].sex_restriction == "female"

    payload["phenotypes"][0]["sex_restriction"] = "none"
    with pytest.raises(ValidationError):
        DiseaseProfile.model_validate(payload)
