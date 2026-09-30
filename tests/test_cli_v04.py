import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.cli import app
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.exports.training import iter_dataset
from rare_disease_simulator.simulation.schema import SyntheticCase
from tests.fixtures.readers import fixture_path
from tests.test_cli_gene_first import _text

HPOA = fixture_path("hpoa")
V04 = fixture_path("v04")
CONFIG = ("--config", str(V04 / "config.yaml"))
GENE_OPTIONS = (
    "--gene-profiles", str(HPOA / "gene_profiles.json"),
    "--gnn-genes", str(HPOA / "gnn_genes.json"),
)

NOISE = ("--noise-vocabulary", str(HPOA / "noise_vocabulary.tsv"))


def _invoke(*args: str):
    return CliRunner().invoke(app, [*CONFIG, *args])


def _calibration(directory: Path) -> Path:
    path = directory / "calib.json"
    path.write_text(
        json.dumps({"frequency.shrinkage_mean": 0.4, "_provenance": "fixture, not fitted"}),
        encoding="utf-8",
    )
    return path


@pytest.fixture(scope="module")
def v04_run(tmp_path_factory):
    directory = tmp_path_factory.mktemp("v04")
    profiles = directory / "profiles.jsonl"
    built = _invoke(
        "build-profiles", "--from-hpoa", "--gene-profiles", str(HPOA / "gene_profiles.json"),
        "--output", str(profiles),
    )
    assert built.exit_code == 0, built.output
    cases = directory / "cases.jsonl"
    simulated = _invoke(
        "simulate", "--profiles", str(profiles), "--output", str(cases), *GENE_OPTIONS,
        "--genes", "all", "--cases-per-gene", "20", "--difficulty", "medium",
        "--calibration", str(_calibration(directory)), *NOISE,
    )
    assert simulated.exit_code == 0, simulated.output
    return directory, profiles, cases


def test_build_counts_recounted_rows_and_their_onsets(v04_run) -> None:
    directory, _, _ = v04_run

    summary = json.loads((directory / "profiles.summary.json").read_text("utf-8"))

    assert summary["rows"]["recounted"] == 2
    assert summary["rows"]["recounted_with_onset"] == 2
    assert summary["rows"]["phenotype_rows_used_with_onset"] >= 2


def test_gene_first_v04_run_records_its_artifacts(v04_run) -> None:
    directory, _, cases_path = v04_run

    summary = json.loads((directory / "cases.summary.json").read_text("utf-8"))
    cases = read_model_jsonl(cases_path, SyntheticCase)

    assert summary["simulate"]["simulator_version"] == "0.4.1"
    assert summary["inputs"]["report_model"]["sha256"] == sha256_file(V04 / "report_model.json")
    assert summary["inputs"]["report_model"]["artifact_id"] == "report-model-fixture"
    assert summary["inputs"]["cardinal"]["sha256"] == sha256_file(V04 / "cardinal.tsv")
    assert summary["reporting"]["mode"] == "report_model"
    assert summary["reporting"]["budget_zero_cases_dropped"] == 5
    assert summary["reporting"]["cardinal"]["pairs"] == 2
    assert summary["entity_profiles"]["entities_with_2_profiles"] == 1
    assert summary["config"]["frequency"]["shrinkage_mean"] == 0.4
    merged = [c for c in cases if c.target.entity_id == "OMIM:100008"]
    assert merged and {tuple(c.target.profile_ids) for c in merged} == {
        ("OMIM:100008", "ORPHA:3001")
    }
    assert all(c.metadata.report_budget in {1, 2, 3, 4} for c in cases)
    for case in cases:
        k = case.metadata.report_budget
        assert case.metadata.report_budget_profile == max(1, k - len(case.noise_phenotypes))
    reporting = summary["reporting"]
    total = sum(len(c.positive_phenotypes) + len(c.noise_phenotypes) for c in cases)
    assert reporting["reported_total_mean"] == round(total / len(cases), 4)
    assert reporting["budget_mean_drawn"] == round(
        sum(c.metadata.report_budget for c in cases) / len(cases), 4
    )
    assert reporting["reported_noise_mean"] > 0


def test_v04_simulate_is_byte_identical(v04_run, tmp_path: Path) -> None:
    _, profiles, cases = v04_run
    again = tmp_path / "again.jsonl"

    result = _invoke(
        "simulate", "--profiles", str(profiles), "--output", str(again), *GENE_OPTIONS,
        "--genes", "all", "--cases-per-gene", "20", "--difficulty", "medium",
        "--calibration", str(_calibration(tmp_path)), *NOISE,
    )

    assert result.exit_code == 0, result.output
    assert again.read_bytes() == cases.read_bytes()


def test_v04_export_and_validate_end_to_end(v04_run, tmp_path: Path) -> None:
    _, _, cases = v04_run
    first, second = tmp_path / "a" / "ds.jsonl.gz", tmp_path / "b" / "ds.jsonl.gz"

    for output in (first, second):
        exported = _invoke(
            "export-training", "--cases", str(cases), "--output", str(output), "--allow-dirty"
        )
        assert exported.exit_code == 0, exported.output

    assert first.read_bytes() == second.read_bytes()
    card = json.loads((tmp_path / "a" / "ds.card.json").read_text("utf-8"))
    assert card["dataset_id"].startswith("ds-sim0.4.1-")
    assert card["inputs"]["report_model"]["sha256"] == sha256_file(V04 / "report_model.json")
    assert card["inputs"]["cardinal"]["sha256"] == sha256_file(V04 / "cardinal.tsv")
    assert card["reporting"]["mode"] == "report_model"
    assert card["validate"]["violations"]["count"] == 0
    reporting = card["validate"]["reporting"]
    assert reporting["cardinal"]["reported_terms"] == reporting["cardinal"]["true_terms"]
    for record in iter_dataset(first):
        assert set(record["present"]) <= set(record["true_present"])
        assert not set(record["excluded"]) & set(record["true_present"])

    dataset_check = _invoke(
        "validate", "--dataset", str(first), "--gnn-genes", str(HPOA / "gnn_genes.json"),
        "--hpo-json", str(HPOA / "hp_mini.json"),
    )
    assert dataset_check.exit_code == 0, dataset_check.output
    cases_check = _invoke("validate", "--cases", str(cases))
    assert cases_check.exit_code == 0, cases_check.output
    assert "Reporting (" in cases_check.output
    assert "Invariant violations: 0" in cases_check.output


def test_v04_simulate_refuses_missing_inputs(v04_run, tmp_path: Path) -> None:
    _, profiles, _ = v04_run
    base = ("simulate", "--profiles", str(profiles), "--output", str(tmp_path / "c.jsonl"),
            *GENE_OPTIONS)

    no_mean = _invoke(*base)
    assert no_mean.exit_code != 0
    assert "shrinkage_mean is not set" in _text(no_mean.output)

    bad_model = tmp_path / "bad.json"
    bad_model.write_text('{"format": "report-model"}', encoding="utf-8")
    rejected = _invoke(*base, "--calibration", str(_calibration(tmp_path)),
                       "--report-model", str(bad_model))
    assert rejected.exit_code != 0
    assert "does not match report-model" in _text(rejected.output)

    observation = tmp_path / "observation.json"
    observation.write_text(
        json.dumps({"reporting.mode": "observation", "frequency.shrinkage_mean": 0.4}),
        encoding="utf-8",
    )
    clash = _invoke(*base, "--calibration", str(observation),
                    "--report-model", str(V04 / "report_model.json"))
    assert clash.exit_code != 0
    assert "reporting.mode is observation" in _text(clash.output)
