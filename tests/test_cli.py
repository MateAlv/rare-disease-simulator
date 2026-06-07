from pathlib import Path

from typer.testing import CliRunner

from rare_disease_simulator.cli import app
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.profiles.schema import DiseaseProfile


def test_cli_help() -> None:
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "fetch-sources" in result.output
    assert "export-graphens" in result.output


def test_cli_validate_loads_mvp_config() -> None:
    result = CliRunner().invoke(app, ["validate"])

    assert result.exit_code == 0
    assert "Config OK: 8 MVP diseases" in result.output
    assert "NPC1" in result.output
    assert "PTPN11" in result.output


def test_cli_fetch_sources_reports_configured_hpo_paths() -> None:
    result = CliRunner().invoke(app, ["fetch-sources"])

    assert result.exit_code == 0
    assert "HPO terms: data/raw/hpo/hpo_terms.tsv" in result.output
    assert "HPO negative annotations" in result.output


def test_cli_build_profiles_from_selected_fixture_dir(tmp_path: Path) -> None:
    output_path = tmp_path / "profiles.jsonl"

    result = CliRunner().invoke(
        app,
        [
            "build-profiles",
            "--fixture-dir",
            "data/fixtures",
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert "Patch validation: accepted" in result.output
    assert "Wrote 1 profile(s)" in result.output

    profiles = read_model_jsonl(output_path, DiseaseProfile)
    assert profiles[0].genes[0].symbol == "NPC1"
