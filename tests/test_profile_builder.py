from pathlib import Path

from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.profiles.builder import (
    build_profiles_from_fixtures,
    write_profiles_jsonl,
)
from rare_disease_simulator.profiles.schema import DiseaseProfile


def test_build_profiles_from_data_fixtures_merges_validated_patch() -> None:
    result = build_profiles_from_fixtures(Path("data/fixtures"))

    profile = result.profiles[0]

    assert result.patch_validation.status == "accepted"
    assert profile.disease_id == "ORPHA:646"
    assert profile.genes[0].symbol == "NPC1"
    assert [phenotype.hpo_id for phenotype in profile.phenotypes] == [
        "HP:0001251",
        "HP:0001433",
        "HP:0000511",
        "HP:0001250",
    ]
    assert result.warnings == ["missing_negative_phenotypes"]


def test_write_profiles_jsonl_round_trips_built_fixture_profile(tmp_path: Path) -> None:
    result = build_profiles_from_fixtures(Path("data/fixtures"))
    output_path = tmp_path / "profiles.jsonl"

    count = write_profiles_jsonl(output_path, result.profiles)
    reloaded = read_model_jsonl(output_path, DiseaseProfile)

    assert count == 1
    assert reloaded[0].disease_id == "ORPHA:646"
    assert reloaded[0].quality.warnings == ["missing_negative_phenotypes"]

