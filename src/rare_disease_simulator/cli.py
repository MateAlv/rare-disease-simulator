"""Command-line entrypoint for rare-disease-simulator."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError

from rare_disease_simulator import __version__
from rare_disease_simulator.config import (
    DEFAULT_CONFIG_PATH,
    AppConfig,
    SimulationSettings,
    load_config,
)
from rare_disease_simulator.exports.graphens import write_graphens_json
from rare_disease_simulator.exports.jsonl import read_model_jsonl, write_jsonl
from rare_disease_simulator.profiles.builder import (
    build_profiles_from_fixtures,
    write_profiles_jsonl,
)
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase
from rare_disease_simulator.simulation.simulator import simulate_cases

app = typer.Typer(
    help="Build disease profiles and simulate rare disease phenotype cases.",
    no_args_is_help=True,
)


def configure_logging(verbose: bool = False) -> None:
    """Configure consistent CLI logging."""

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )


def _load_config_or_exit(config_path: Path) -> AppConfig:
    try:
        return load_config(config_path)
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"config file not found: {config_path}") from exc
    except ValidationError as exc:
        raise typer.BadParameter(f"invalid config file {config_path}: {exc}") from exc


@app.callback()
def callback(
    ctx: typer.Context,
    config: Annotated[
        Path,
        typer.Option(
            "--config",
            "-c",
            help="Path to the YAML configuration file.",
            exists=False,
            dir_okay=False,
            resolve_path=False,
        ),
    ] = DEFAULT_CONFIG_PATH,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Enable debug logging.")] = False,
    version: Annotated[
        bool,
        typer.Option("--version", help="Show package version and exit."),
    ] = False,
) -> None:
    """Load shared CLI configuration."""

    configure_logging(verbose)
    if version:
        typer.echo(__version__)
        raise typer.Exit()
    ctx.obj = {"config_path": config}


def _config_from_context(ctx: typer.Context) -> AppConfig:
    config_path = ctx.obj["config_path"] if ctx.obj else DEFAULT_CONFIG_PATH
    return _load_config_or_exit(config_path)


@app.command("fetch-sources")
def fetch_sources(ctx: typer.Context) -> None:
    """Validate configured source locations for future ingestion."""

    config = _config_from_context(ctx)
    source_paths = [
        ("HPO directory", config.sources.hpo_dir),
        ("HPO terms", config.sources.hpo_terms_path),
        ("HPO phenotype annotations", config.sources.phenotype_annotation_path),
        ("HPO gene-phenotype links", config.sources.genes_to_phenotype_path),
        ("HPO negative annotations", config.sources.negative_phenotype_annotation_path),
        ("Orphadata directory", config.sources.orphadata_dir),
        ("MONDO path", config.sources.mondo_path),
    ]
    for label, path in source_paths:
        status = "found" if path.exists() else "missing"
        typer.echo(f"{label}: {path} [{status}]")


@app.command("extract-profile-patches")
def extract_profile_patches(ctx: typer.Context) -> None:
    """Extract DiseaseProfilePatch artifacts from configured text snippets."""

    config = _config_from_context(ctx)
    typer.echo(
        "LLM extraction is configured "
        f"with provider={config.llm.provider}, model={config.llm.model}."
    )


@app.command("build-profiles")
def build_profiles(
    ctx: typer.Context,
    fixture_dir: Annotated[
        Path | None,
        typer.Option(
            "--fixture-dir",
            help="Run profile building from a local fixture directory.",
            exists=True,
            file_okay=False,
            dir_okay=True,
            resolve_path=False,
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Override the configured profiles JSONL output path.",
            dir_okay=False,
            resolve_path=False,
        ),
    ] = None,
) -> None:
    """Build validated DiseaseProfile artifacts."""

    config = _config_from_context(ctx)
    if fixture_dir is None:
        typer.echo(f"Configured MVP diseases: {len(config.mvp.diseases)}")
        typer.echo(f"Profiles output: {config.exports.profiles_path}")
        typer.echo("No profile source selected. Use --fixture-dir to run the fixture path.")
        return

    result = build_profiles_from_fixtures(fixture_dir)
    output_path = output or config.exports.profiles_path
    written_count = write_profiles_jsonl(output_path, result.profiles)

    typer.echo(f"Patch validation: {result.patch_validation.status}")
    typer.echo(f"Wrote {written_count} profile(s) to {output_path}")
    if result.warnings:
        typer.echo(f"Quality warnings: {', '.join(result.warnings)}")


def _simulation_config(settings: SimulationSettings) -> SimulationConfig:
    return SimulationConfig(
        cases_per_disease_per_difficulty=settings.cases_per_disease_per_difficulty,
        difficulties=settings.difficulties,  # type: ignore[arg-type]
        seed=settings.seed,
    )


@app.command("simulate")
def simulate(
    ctx: typer.Context,
    profiles: Annotated[
        Path | None,
        typer.Option(
            "--profiles",
            help="Override the configured profiles JSONL input path.",
            exists=True,
            dir_okay=False,
            resolve_path=False,
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Override the configured rich cases JSONL output path.",
            dir_okay=False,
            resolve_path=False,
        ),
    ] = None,
) -> None:
    """Simulate synthetic patient cases from validated profiles."""

    config = _config_from_context(ctx)
    profiles_path = profiles or config.exports.profiles_path
    if not profiles_path.exists():
        raise typer.BadParameter(f"profiles file not found: {profiles_path}")

    sim_config = _simulation_config(config.simulation)
    profile_records = read_model_jsonl(profiles_path, DiseaseProfile)
    cases: list[SyntheticCase] = []
    for profile in profile_records:
        cases.extend(simulate_cases(profile, sim_config))

    output_path = output or config.exports.rich_cases_path
    written = write_jsonl(output_path, cases)
    typer.echo(
        f"Simulated {written} case(s) from {len(profile_records)} profile(s) "
        f"({sim_config.cases_per_disease_per_difficulty} per disease/difficulty, "
        f"difficulties={','.join(sim_config.difficulties)}, seed={sim_config.seed})."
    )
    typer.echo(f"Wrote rich cases to {output_path}")


@app.command("export-graphens")
def export_graphens(
    ctx: typer.Context,
    cases: Annotated[
        Path | None,
        typer.Option(
            "--cases",
            help="Override the configured rich cases JSONL input path.",
            exists=True,
            dir_okay=False,
            resolve_path=False,
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Override the configured GraPhens JSON output path.",
            dir_okay=False,
            resolve_path=False,
        ),
    ] = None,
) -> None:
    """Export simulated cases in GraPhens-compatible gene-grouped format."""

    config = _config_from_context(ctx)
    cases_path = cases or config.exports.rich_cases_path
    if not cases_path.exists():
        raise typer.BadParameter(f"rich cases file not found: {cases_path}")

    case_records = read_model_jsonl(cases_path, SyntheticCase)
    output_path = output or config.exports.graphens_path
    mapping_path = output_path.with_name(output_path.stem + ".mapping.json")
    export = write_graphens_json(output_path, case_records, mapping_path=mapping_path)

    total_rows = sum(len(rows) for rows in export.values())
    typer.echo(
        f"Exported {total_rows} case(s) across {len(export)} gene(s) to {output_path}"
    )
    typer.echo(f"Wrote row mapping to {mapping_path}")


@app.command("validate")
def validate(ctx: typer.Context) -> None:
    """Validate the configured MVP setup."""

    config = _config_from_context(ctx)
    genes = ", ".join(disease.gene for disease in config.mvp.diseases)
    typer.echo(f"Config OK: {len(config.mvp.diseases)} MVP diseases ({genes})")


def main() -> None:
    """Run the CLI."""

    app()


if __name__ == "__main__":
    main()
