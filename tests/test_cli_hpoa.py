import json
from pathlib import Path

from typer.testing import CliRunner

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.cli import app
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.simulation.schema import SyntheticCase
from tests.fixtures.readers import fixture_path

FIXTURES = fixture_path("hpoa")
CONFIG = ("--config", str(FIXTURES / "config.yaml"))


def _build(output: Path, *extra: str):
    return CliRunner().invoke(
        app,
        [
            *CONFIG,
            "build-profiles",
            "--from-hpoa",
            "--hpo-json",
            str(FIXTURES / "hp_mini.json"),
            "--hpoa",
            str(FIXTURES / "phenotype_mini.hpoa"),
            "--genes-to-disease",
            str(FIXTURES / "genes_to_disease_mini.txt"),
            "--orphanet-ages",
            str(FIXTURES / "orphanet_ages_mini.xml"),
            "--output",
            str(output),
            *extra,
        ],
    )


def test_build_profiles_from_hpoa_end_to_end(tmp_path: Path) -> None:
    output = tmp_path / "profiles.jsonl"

    result = _build(output, "--exclude-pmids", str(FIXTURES / "heldout_pmids.txt"))

    assert result.exit_code == 0, result.output
    assert "Wrote 4 profile(s)" in result.output
    assert "Masked 3 HPOA row(s) across 2 disease(s)" in result.output
    profiles = read_model_jsonl(output, DiseaseProfile)
    assert [profile.disease_id for profile in profiles][0] == "OMIM:100001"

    summary = json.loads((tmp_path / "profiles.summary.json").read_text(encoding="utf-8"))
    assert summary["output"]["profiles_written"] == 4
    assert len(summary["output"]["sha256"]) == 64
    assert summary["masking"]["rows_masked"] == 3
    assert set(summary["inputs"]) == {
        "hp_json",
        "phenotype_hpoa",
        "genes_to_disease",
        "orphanet_ages",
        "exclude_pmids",
    }
    assert "sha" in summary["build"]["simulator_git"]


def test_build_profiles_from_hpoa_is_byte_identical_across_runs(tmp_path: Path) -> None:
    first = tmp_path / "a" / "profiles.jsonl"
    second = tmp_path / "b" / "profiles.jsonl"
    mask = ("--exclude-pmids", str(FIXTURES / "heldout_pmids.txt"))

    assert _build(first, *mask).exit_code == 0
    assert _build(second, *mask).exit_code == 0

    assert first.read_bytes() == second.read_bytes()


def test_build_profiles_from_hpoa_warns_without_mask(tmp_path: Path) -> None:
    result = _build(tmp_path / "profiles.jsonl", "--summary", str(tmp_path / "s.json"))

    assert result.exit_code == 0
    assert "no --exclude-pmids given" in result.output
    assert "Wrote 5 profile(s)" in result.output
    assert (tmp_path / "s.json").exists()


def test_build_profiles_from_hpoa_rejects_missing_inputs(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        app,
        [
            *CONFIG,
            "build-profiles",
            "--from-hpoa",
            "--hpo-json",
            str(tmp_path / "missing.json"),
            "--output",
            str(tmp_path / "profiles.jsonl"),
        ],
    )

    assert result.exit_code != 0
    assert "hp.json not found" in result.output


def test_simulate_passes_ontology_noise_vocabulary_and_labels(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.jsonl"
    assert _build(profiles, "--exclude-pmids", str(FIXTURES / "heldout_pmids.txt")).exit_code == 0
    cases_path = tmp_path / "cases.jsonl"

    result = CliRunner().invoke(
        app,
        [
            *CONFIG,
            "simulate",
            "--profiles",
            str(profiles),
            "--output",
            str(cases_path),
            "--hpo-json",
            str(FIXTURES / "hp_mini.json"),
            "--noise-vocabulary",
            str(FIXTURES / "noise_vocabulary.tsv"),
            "--labels",
            str(FIXTURES / "labels.json"),
        ],
    )

    assert result.exit_code == 0, result.output
    cases = read_model_jsonl(cases_path, SyntheticCase)
    alpha = [case for case in cases if case.target.disease_id == "OMIM:100001"]
    assert {case.target.gene_label for case in alpha} == {0}
    assert {case.target.disease_label for case in alpha} == {7}
    gamma = [case for case in cases if case.target.disease_id == "ORPHA:3001"]
    assert {case.target.gene_label for case in gamma} == {None}

    noise_labels = {
        (term.hpo_id, term.label) for case in cases for term in case.noise_phenotypes
    }
    assert ("HP:0001249", "Intellectual disability") in noise_labels
    assert any(
        term.reason == "generalized" for case in cases for term in case.positive_phenotypes
    )


def test_build_profiles_reports_the_omim_orpha_mapping(tmp_path: Path) -> None:
    output = tmp_path / "profiles.jsonl"

    result = _build(
        output,
        "--exclude-pmids",
        str(FIXTURES / "heldout_pmids.txt"),
        "--omim-orpha-map",
        str(FIXTURES / "en_product1_mini.xml"),
    )

    assert result.exit_code == 0, result.output
    summary = json.loads((tmp_path / "profiles.summary.json").read_text(encoding="utf-8"))
    assert summary["inputs"]["omim_orpha_map"]["version"] == "2026-06-23 07:53:50"
    assert summary["omim_orpha_map"]["omim_E_validated"] == 2
    assert summary["age_of_onset"]["gained_via_omim_mapping"] == 1
    assert "Onset coverage: 3/4 disease(s), 1 via OMIM-ORPHA mapping" in result.output


def _simulate(profiles: Path, output: Path, *extra: str):
    return CliRunner().invoke(
        app,
        [
            *CONFIG,
            "simulate",
            "--profiles",
            str(profiles),
            "--output",
            str(output),
            "--hpo-json",
            str(FIXTURES / "hp_mini.json"),
            "--noise-vocabulary",
            str(FIXTURES / "noise_vocabulary.tsv"),
            *extra,
        ],
    )


def test_simulate_samples_diseases_records_provenance_and_is_byte_stable(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.jsonl"
    assert _build(profiles, "--exclude-pmids", str(FIXTURES / "heldout_pmids.txt")).exit_code == 0
    first, second = tmp_path / "a" / "cases.jsonl", tmp_path / "b" / "cases.jsonl"
    options = ("--sample-diseases", "2", "--cases-per-disease", "3", "--difficulty", "medium")

    result = _simulate(profiles, first, *options)
    assert result.exit_code == 0, result.output
    assert _simulate(profiles, second, *options).exit_code == 0

    assert first.read_bytes() == second.read_bytes()
    cases = read_model_jsonl(first, SyntheticCase)
    assert len(cases) == 6
    assert {case.metadata.difficulty for case in cases} == {"medium"}
    versions = cases[0].metadata.source_versions
    assert versions["profiles"] == f"sha256:{sha256_file(profiles)}"
    assert versions["hp_json"].startswith("2026-02-16 sha256:")
    assert "simulator_git" in versions and "HPO disease annotations" in versions

    summary = json.loads((tmp_path / "a" / "cases.summary.json").read_text(encoding="utf-8"))
    assert summary["diseases"]["simulated"] == 2
    assert sum(summary["diseases"]["strata"].values()) == 2
    assert summary["config"]["cases_per_disease_per_difficulty"] == 3
    assert summary["config_hash"] == cases[0].metadata.config_hash
    assert summary["output"]["sha256"] == sha256_file(first)
    assert set(summary["inputs"]) == {"profiles", "hp_json", "noise_vocabulary"}
