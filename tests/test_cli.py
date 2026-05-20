from typer.testing import CliRunner

from rare_disease_simulator.cli import app


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

