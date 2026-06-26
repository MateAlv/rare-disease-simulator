import json

from rare_disease_simulator.exports.graphens import (
    build_graphens_export,
    write_graphens_json,
)
from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import simulate_cases
from tests.fixtures.readers import read_disease_profile


def _cases():
    profile = read_disease_profile()
    config = SimulationConfig(
        cases_per_disease_per_difficulty=4,
        difficulties=["easy"],
        seed=3,
    )
    return simulate_cases(profile, config)


def test_build_graphens_export_groups_by_gene_with_row_mapping() -> None:
    cases = _cases()

    export, mapping = build_graphens_export(cases)

    assert set(export) == {"NPC1"}
    assert len(export["NPC1"]) == len(cases)
    assert mapping["NPC1"] == [case.case_id for case in cases]
    for row in export["NPC1"]:
        assert all(hpo_id.startswith("HP:") for hpo_id in row)


def test_baseline_export_excludes_negatives() -> None:
    cases = _cases()
    negative_ids = {
        phenotype.hpo_id for case in cases for phenotype in case.negative_phenotypes
    }

    export, _ = build_graphens_export(cases)

    exported_ids = {hpo_id for rows in export.values() for row in rows for hpo_id in row}
    assert exported_ids.isdisjoint(negative_ids)


def test_write_graphens_json_writes_export_and_mapping(tmp_path) -> None:
    cases = _cases()
    out_path = tmp_path / "graphens.json"
    mapping_path = tmp_path / "graphens.mapping.json"

    write_graphens_json(out_path, cases, mapping_path=mapping_path)

    export = json.loads(out_path.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    assert "NPC1" in export
    assert len(mapping["NPC1"]) == len(cases)
