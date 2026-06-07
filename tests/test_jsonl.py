from pathlib import Path

import pytest

from rare_disease_simulator.exports.jsonl import (
    JsonlDecodeError,
    append_jsonl,
    read_jsonl,
    read_model_jsonl,
    write_jsonl,
)
from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch
from rare_disease_simulator.simulation.schema import SyntheticCase
from tests.fixtures.readers import read_profile_patch, read_synthetic_case


def test_write_and_read_model_jsonl_round_trip(tmp_path: Path) -> None:
    patch = read_profile_patch()
    path = tmp_path / "nested" / "profile_patches.jsonl"

    count = write_jsonl(path, [patch])
    reloaded = read_model_jsonl(path, DiseaseProfilePatch)

    assert count == 1
    assert reloaded[0].disease_id == "ORPHA:646"
    assert reloaded[0].phenotypes[0].hpo_id == "HP:0001251"


def test_append_jsonl_keeps_existing_records(tmp_path: Path) -> None:
    case = read_synthetic_case()
    path = tmp_path / "rich_cases.jsonl"

    write_jsonl(path, [case])
    append_jsonl(path, [case.model_copy(update={"case_id": "synthetic-ORPHA_646-000002"})])

    reloaded = read_model_jsonl(path, SyntheticCase)

    assert [record.case_id for record in reloaded] == [
        "synthetic-ORPHA_646-000001",
        "synthetic-ORPHA_646-000002",
    ]


def test_read_jsonl_skips_blank_lines(tmp_path: Path) -> None:
    path = tmp_path / "records.jsonl"
    path.write_text('{"a": 1}\n\n{"b": 2}\n', encoding="utf-8")

    assert read_jsonl(path) == [{"a": 1}, {"b": 2}]


def test_read_jsonl_reports_line_number_for_invalid_json(tmp_path: Path) -> None:
    path = tmp_path / "broken.jsonl"
    path.write_text('{"ok": true}\nnot-json\n', encoding="utf-8")

    with pytest.raises(JsonlDecodeError) as error:
        read_jsonl(path)

    assert "broken.jsonl:2" in str(error.value)
