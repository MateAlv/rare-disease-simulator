import json

from rare_disease_simulator.exports.graphens import write_graphens_json
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.profiles.builder import (
    build_profiles_from_fixtures,
    write_profiles_jsonl,
)
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase
from rare_disease_simulator.simulation.simulator import simulate_cases
from tests.fixtures.readers import FIXTURE_DIR


def test_full_fixture_slice_build_simulate_export(tmp_path) -> None:
    build_result = build_profiles_from_fixtures(FIXTURE_DIR)

    assert build_result.patch_validation.status == "accepted"
    assert len(build_result.profiles) == 1
    profile = build_result.profiles[0]
    assert profile.disease_id == "ORPHA:646"
    # Seizure (HP:0001250) comes only from the LLM patch; merge must fold it in.
    assert any(p.hpo_id == "HP:0001250" for p in profile.phenotypes)

    profiles_path = tmp_path / "profiles.jsonl"
    written = write_profiles_jsonl(profiles_path, build_result.profiles)
    assert written == 1

    config = SimulationConfig(
        cases_per_disease_per_difficulty=10,
        difficulties=["easy"],
        seed=42,
        reporting={"mode": "observation"},
    )
    cases = simulate_cases(profile, config)
    assert len(cases) == 10

    cases_path = tmp_path / "rich_cases.jsonl"
    write_profiles_jsonl(cases_path, cases)
    reloaded = read_model_jsonl(cases_path, SyntheticCase)
    assert len(reloaded) == 10
    assert all(case.target.gene == "NPC1" for case in reloaded)

    graphens_path = tmp_path / "graphens.json"
    write_graphens_json(graphens_path, reloaded)
    export = json.loads(graphens_path.read_text(encoding="utf-8"))
    assert set(export) == {"NPC1"}
    assert len(export["NPC1"]) == 10
