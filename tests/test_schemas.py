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

    assert config.missingness_rate == 0.25


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
