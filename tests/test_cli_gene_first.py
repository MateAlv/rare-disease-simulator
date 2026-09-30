import json
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from rare_disease_simulator import cli
from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.cli import app
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.exports.training import (
    RECORD_FIELDS,
    card_path_for,
    encode_dataset,
    iter_dataset,
    sort_records,
)
from rare_disease_simulator.simulation.schema import SyntheticCase
from rare_disease_simulator.validation.dataset import validate_dataset
from tests.fixtures.readers import fixture_path

FIXTURES = fixture_path("hpoa")
CONFIG = ("--config", str(FIXTURES / "config.yaml"))
GENE_OPTIONS = (
    "--gene-profiles",
    str(FIXTURES / "gene_profiles.json"),
    "--gnn-genes",
    str(FIXTURES / "gnn_genes.json"),
)


def _text(output: str) -> str:
    """CLI output with Rich's error box and line wrapping removed."""

    return " ".join(re.sub(r"[│╭╮╰╯─]", " ", output).split())


def _invoke(*args: str):
    return CliRunner().invoke(app, [*CONFIG, *args])


def _build(output: Path, *extra: str):
    return _invoke(
        "build-profiles",
        "--from-hpoa",
        "--orphanet-ages",
        str(FIXTURES / "orphanet_ages_mini.xml"),
        "--exclude-pmids",
        str(FIXTURES / "heldout_pmids.txt"),
        "--output",
        str(output),
        *extra,
    )


def _simulate_genes(profiles: Path, output: Path, *extra: str):
    return _invoke(
        "simulate",
        "--profiles",
        str(profiles),
        "--output",
        str(output),
        "--noise-vocabulary",
        str(FIXTURES / "noise_vocabulary.tsv"),
        *GENE_OPTIONS,
        "--cases-per-gene",
        "4",
        *extra,
    )


@pytest.fixture()
def gene_run(tmp_path: Path) -> tuple[Path, Path]:
    profiles = tmp_path / "profiles.jsonl"
    result = _build(profiles, "--gene-profiles", str(FIXTURES / "gene_profiles.json"))
    assert result.exit_code == 0, result.output
    cases = tmp_path / "cases.jsonl"
    result = _simulate_genes(profiles, cases)
    assert result.exit_code == 0, result.output
    return profiles, cases


def test_gene_first_simulate_labels_cases_and_reports_skips(gene_run) -> None:
    _, cases_path = gene_run
    cases = read_model_jsonl(cases_path, SyntheticCase)

    assert {(c.target.gene, c.target.gene_label) for c in cases} == {
        ("GENEA", 1),
        ("GENEB", 2),
        ("OLDC", 3),
    }
    assert len(cases) == 3 * 4 * 3
    assert {c.target.entity_id for c in cases if c.target.gene == "OLDC"} == {"OMIM:100008"}
    assert {c.metadata.sex_prior_key for c in cases if c.target.gene == "GENEB"} == {
        "male_biased"
    }
    summary = json.loads(cases_path.with_name("cases.summary.json").read_text("utf-8"))
    assert summary["simulate"]["mode"] == "gene_first"
    assert summary["genes"]["targets"] == 3
    assert summary["genes"]["skipped"]["unresolved_symbol"] == 2
    assert summary["genes"]["requested"] == "all"
    assert {"gene_profiles", "gnn_genes", "profiles"} <= set(summary["inputs"])
    assert {"held-out reference mask", "gene profiles"} <= set(summary["profile_sources"])


def test_gene_first_rejects_disease_first_options(tmp_path: Path) -> None:
    result = _invoke(
        "simulate",
        "--profiles",
        str(FIXTURES / "gene_profiles.json"),
        "--genes",
        "all",
        "--sample-diseases",
        "2",
    )

    assert result.exit_code != 0
    assert "gene-first mode does not take --sample-diseases" in _text(result.output)


def test_gene_list_and_unknown_symbols(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.jsonl"
    assert _build(profiles, "--gene-profiles", str(FIXTURES / "gene_profiles.json")).exit_code == 0
    genes_file = tmp_path / "genes.txt"
    genes_file.write_text("# wanted\nGENEB\n", encoding="utf-8")

    listed = _simulate_genes(profiles, tmp_path / "a.jsonl", "--genes", "GENEA,OLDC")
    from_file = _simulate_genes(profiles, tmp_path / "b.jsonl", "--genes", str(genes_file))
    unknown = _simulate_genes(profiles, tmp_path / "c.jsonl", "--genes", "NOPE")

    assert listed.exit_code == 0, listed.output
    assert {c.target.gene for c in read_model_jsonl(tmp_path / "a.jsonl", SyntheticCase)} == {
        "GENEA",
        "OLDC",
    }
    assert from_file.exit_code == 0, from_file.output
    assert {c.target.gene for c in read_model_jsonl(tmp_path / "b.jsonl", SyntheticCase)} == {
        "GENEB"
    }
    assert unknown.exit_code != 0 and "not in the GNN vocabulary" in _text(unknown.output)


def test_pinned_gene_profiles_sha_is_enforced(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.jsonl"
    assert _build(profiles).exit_code == 0
    config = tmp_path / "config.yaml"
    config.write_text(
        (FIXTURES / "config.yaml").read_text("utf-8")
        + f"  gene_profiles_path: {FIXTURES / 'gene_profiles.json'}\n"
        + "  gene_profiles_sha256: deadbeef\n"
        + f"  gnn_genes_path: {FIXTURES / 'gnn_genes.json'}\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        app,
        ["--config", str(config), "simulate", "--profiles", str(profiles), "--genes", "all"],
    )

    assert result.exit_code != 0
    assert "not the pinned deadbeef" in _text(result.output)


def test_calibration_file_overrides_knobs_and_is_recorded(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.jsonl"
    assert _build(profiles, "--gene-profiles", str(FIXTURES / "gene_profiles.json")).exit_code == 0
    noise = tmp_path / "noise.tsv"
    noise.write_text("hpo_id\tpatients\tweight\nHP:0000505\t9\t0.5\n", encoding="utf-8")
    calibration = tmp_path / "calib.json"
    calibration.write_text(
        json.dumps(
            {
                "missingness.sex_unknown": 1.0,
                "presets.medium.noise_mean": 3.0,
                "noise.vocabulary_path": "noise.tsv",
            }
        ),
        encoding="utf-8",
    )
    cases_path = tmp_path / "cases.jsonl"

    result = _invoke(
        "simulate", "--profiles", str(profiles), "--output", str(cases_path), *GENE_OPTIONS,
        "--calibration", str(calibration), "--difficulty", "medium",
    )

    assert result.exit_code == 0, result.output
    cases = read_model_jsonl(cases_path, SyntheticCase)
    assert {c.patient.sex for c in cases} == {"unknown"}
    assert {n.hpo_id for c in cases for n in c.noise_phenotypes} <= {"HP:0000505"}
    summary = json.loads((tmp_path / "cases.summary.json").read_text("utf-8"))
    assert summary["calibration"]["sha256"] == sha256_file(calibration)
    assert summary["inputs"]["calibration"]["sha256"] == sha256_file(calibration)
    assert summary["inputs"]["noise_vocabulary"]["path"] == str(noise)
    assert summary["config"]["missingness"]["sex_unknown"] == 1.0

    bad = tmp_path / "bad.json"
    bad.write_text('{"missingness.sex": 0.1}', encoding="utf-8")
    rejected = _invoke(
        "simulate", "--profiles", str(profiles), *GENE_OPTIONS, "--calibration", str(bad)
    )
    assert rejected.exit_code != 0 and "unknown calibration key" in _text(rejected.output)


def test_export_validate_end_to_end_and_byte_identical(gene_run, tmp_path: Path) -> None:
    _, cases = gene_run
    first, second = tmp_path / "x" / "ds.jsonl.gz", tmp_path / "y" / "ds.jsonl.gz"

    result = _invoke("export-training", "--cases", str(cases), "--output", str(first),
                     "--allow-dirty")
    assert result.exit_code == 0, result.output
    assert _invoke("export-training", "--cases", str(cases), "--output", str(second),
                   "--allow-dirty").exit_code == 0

    assert first.read_bytes() == second.read_bytes()
    card = json.loads((tmp_path / "x" / "ds.card.json").read_text("utf-8"))
    assert card["file"]["sha256"] == sha256_file(first)
    assert card["counts"]["cases"] == 36
    assert card["counts"]["per_gene_cases"] == {"GENEA": 12, "GENEB": 12, "OLDC": 12}
    assert card["split"]["fractions"] == {"train": 0.9, "val": 0.05, "test": 0.05}
    assert card["dataset_id"] == "ds-sim0.4.1-hpo2026-02-16-n4-s42"
    assert card["validate"]["violations"]["count"] == 0
    assert card["genes"]["skipped"]["no_gene_record"] == 1
    assert {"gene_profiles", "gnn_genes", "profiles", "profile_sources"} <= set(card["inputs"])
    assert card["simulator"]["allow_dirty"] is True
    assert card["inputs"]["heldout_pmids"]["sha256"] == sha256_file(FIXTURES / "heldout_pmids.txt")

    records = list(iter_dataset(first))
    assert [tuple(record) for record in records] == [RECORD_FIELDS] * 36
    assert records == sort_records(records)

    for record in records:
        assert set(record["present"]) <= set(record["true_present"])
    checked = _invoke("validate", "--dataset", str(first), "--gnn-genes",
                      str(FIXTURES / "gnn_genes.json"), "--hpo-json",
                      str(FIXTURES / "hp_mini.json"))
    assert checked.exit_code == 0, checked.output
    assert "Invariant violations: 0" in checked.output
    cases_check = _invoke("validate", "--cases", str(cases))
    assert cases_check.exit_code == 0, cases_check.output
    assert "and 3 gene(s)" in cases_check.output


def test_export_refuses_a_dirty_tree(gene_run, tmp_path: Path, monkeypatch) -> None:
    _, cases = gene_run
    monkeypatch.setattr(cli, "git_revision", lambda: {"sha": "abc", "dirty": True})

    result = _invoke("export-training", "--cases", str(cases), "--output",
                     str(tmp_path / "ds.jsonl.gz"))

    assert result.exit_code != 0
    assert "--allow-dirty" in _text(result.output)
    assert not (tmp_path / "ds.jsonl.gz").exists()


def test_dataset_validation_catches_tampering(gene_run, tmp_path: Path) -> None:
    _, cases = gene_run
    dataset = tmp_path / "ds.jsonl.gz"
    assert _invoke("export-training", "--cases", str(cases), "--output", str(dataset),
                   "--allow-dirty").exit_code == 0
    card = json.loads(card_path_for(dataset).read_text("utf-8"))
    records = list(iter_dataset(dataset))
    records[0]["excluded"] = sorted(set(records[0]["excluded"]) | {records[0]["present"][0]})
    records[1]["split"] = "val" if records[1]["split"] != "val" else "test"
    dataset.write_bytes(encode_dataset(records))

    report = validate_dataset(dataset, card, ["-", "GENEA", "GENEB", "OLDC"])

    assert report["violations"]["by_type"]["term_present_and_excluded"] == 1
    assert report["violations"]["by_type"]["excluded_true_or_ancestor_of_true"] == 1
    assert report["violations"]["by_type"]["split_not_from_seed"] == 1
    assert report["violations"]["by_type"]["sha256_differs_from_card"] == 1
