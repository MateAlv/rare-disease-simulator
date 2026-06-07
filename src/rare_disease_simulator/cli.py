"""Command-line entrypoint for rare-disease-simulator."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError

from rare_disease_simulator import __version__
from rare_disease_simulator.config import DEFAULT_CONFIG_PATH, AppConfig, load_config

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
def build_profiles(ctx: typer.Context) -> None:
    """Build validated DiseaseProfile artifacts."""

    config = _config_from_context(ctx)
    typer.echo(f"Configured MVP diseases: {len(config.mvp.diseases)}")
    typer.echo(f"Profiles output: {config.exports.profiles_path}")


@app.command("simulate")
def simulate(ctx: typer.Context) -> None:
    """Simulate synthetic patient cases from validated profiles."""

    config = _config_from_context(ctx)
    typer.echo(
        "Simulation configured for "
        f"{config.simulation.cases_per_disease_per_difficulty} cases per disease/difficulty "
        f"with seed={config.simulation.seed}."
    )


@app.command("export-graphens")
def export_graphens(ctx: typer.Context) -> None:
    """Export simulated cases in GraPhens-compatible gene-grouped format."""

    config = _config_from_context(ctx)
    typer.echo(f"GraPhens export path: {config.exports.graphens_path}")


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
