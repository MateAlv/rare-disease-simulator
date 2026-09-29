"""Command-line entrypoint for rare-disease-simulator."""

from __future__ import annotations

import json
import logging
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError

from rare_disease_simulator import __version__
from rare_disease_simulator.build_info import git_revision, sha256_file
from rare_disease_simulator.config import (
    DEFAULT_CONFIG_PATH,
    AppConfig,
    MvpDisease,
    load_config,
)
from rare_disease_simulator.data_sources.fetch import DiseaseQuery, fetch_disease_sources
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.data_sources.hpo_annotations import download_hpo_release
from rare_disease_simulator.data_sources.http_client import HttpClient
from rare_disease_simulator.exports.graphens import write_graphens_json
from rare_disease_simulator.exports.jsonl import iter_model_jsonl, read_model_jsonl, write_jsonl
from rare_disease_simulator.profiles.builder import (
    build_profiles_from_fixtures,
    write_profiles_jsonl,
)
from rare_disease_simulator.profiles.hpoa_builder import HpoaBuildInputs, build_profiles_from_hpoa
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.inputs import load_label_maps, load_noise_vocabulary
from rare_disease_simulator.simulation.sampling import sample_diseases as sample_diseases_by_stratum
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase
from rare_disease_simulator.simulation.simulator import (
    SIMULATOR_VERSION,
    SexSpecificTerms,
    config_hash,
    primary_gene,
    simulate_cases,
)
from rare_disease_simulator.validation.cases import format_report, validate_cases

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


def _required_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise typer.BadParameter(f"{label} not found: {path}")
    return path


def _optional_file(path: Path | None, label: str) -> Path | None:
    if path is None:
        return None
    return _required_file(path, label)


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


def _disease_query(disease: MvpDisease) -> DiseaseQuery:
    return DiseaseQuery(
        name=disease.disease,
        gene=disease.gene,
        orpha_id=disease.orpha_id,
        omim=list(disease.omim),
    )


@app.command("download-hpo-release")
def download_hpo_release_command(ctx: typer.Context) -> None:
    """Download the HPO ontology and annotation files into the configured dir."""

    config = _config_from_context(ctx)
    written = download_hpo_release(HttpClient(), config.sources.hpo_dir)
    for name, path in written.items():
        typer.echo(f"{name}: {path}")


@app.command("fetch-disease")
def fetch_disease_command(
    ctx: typer.Context,
    disease: Annotated[
        str | None,
        typer.Option("--disease", help="Substring of a configured disease name to fetch."),
    ] = None,
    all_diseases: Annotated[
        bool,
        typer.Option("--all", help="Fetch every configured MVP disease."),
    ] = False,
    hpoa: Annotated[
        Path | None,
        typer.Option("--hpoa", help="Path to a local phenotype.hpoa for the structured backbone."),
    ] = None,
    output_root: Annotated[
        Path,
        typer.Option("--output-root", help="Directory to write per-disease bundles into."),
    ] = Path("data/raw"),
    retmax: Annotated[
        int,
        typer.Option("--retmax", help="Max PubMed articles to fetch per disease."),
    ] = 5,
    full_text: Annotated[
        bool,
        typer.Option("--full-text/--no-full-text", help="Attempt PMC open-access full text."),
    ] = True,
) -> None:
    """Fetch all available sources for one or more configured diseases."""

    config = _config_from_context(ctx)
    hpoa_path = hpoa or config.sources.phenotype_annotation_path
    resolved_hpoa = hpoa_path if Path(hpoa_path).exists() else None
    if resolved_hpoa is None:
        typer.echo("No phenotype.hpoa found; structured backbone will be skipped.")

    if all_diseases:
        targets = config.mvp.diseases
    elif disease is not None:
        needle = disease.lower()
        targets = [d for d in config.mvp.diseases if needle in d.disease.lower()]
    else:
        raise typer.BadParameter("provide --disease NAME or --all")

    if not targets:
        raise typer.BadParameter(f"no configured disease matches {disease!r}")

    client = HttpClient()
    for target in targets:
        manifest = fetch_disease_sources(
            _disease_query(target),
            client=client,
            output_root=output_root,
            hpoa_path=resolved_hpoa,
            pubmed_retmax=retmax,
            fetch_pmc_full_text=full_text,
        )
        statuses = ", ".join(
            f"{name}={info.get('status', '?')}" for name, info in manifest.sources.items()
        )
        typer.echo(
            f"{target.disease}: {manifest.snippet_count} snippet(s) "
            f"-> {output_root}/{_disease_query(target).slug} [{statuses}]"
        )


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
    from_hpoa: Annotated[
        bool,
        typer.Option(
            "--from-hpoa",
            help="Build one profile per gene-linked disease from the HPO annotation release.",
        ),
    ] = False,
    hpo_json: Annotated[
        Path | None, typer.Option("--hpo-json", help="HPO ontology hp.json.", dir_okay=False)
    ] = None,
    hpoa: Annotated[
        Path | None, typer.Option("--hpoa", help="HPO phenotype.hpoa.", dir_okay=False)
    ] = None,
    genes_to_disease: Annotated[
        Path | None,
        typer.Option("--genes-to-disease", help="HPO genes_to_disease.txt.", dir_okay=False),
    ] = None,
    orphanet_ages: Annotated[
        Path | None,
        typer.Option(
            "--orphanet-ages",
            help="Orphanet en_product9_ages.xml (age-of-onset fallback).",
            dir_okay=False,
        ),
    ] = None,
    omim_orpha_map: Annotated[
        Path | None,
        typer.Option(
            "--omim-orpha-map",
            help="Orphanet en_product1.xml; exact validated OMIM-ORPHA alignments carry onset.",
            dir_okay=False,
        ),
    ] = None,
    exclude_pmids: Annotated[
        Path | None,
        typer.Option(
            "--exclude-pmids",
            help="Held-out references, one per line; HPOA rows citing any are dropped.",
            dir_okay=False,
        ),
    ] = None,
    summary: Annotated[
        Path | None,
        typer.Option(
            "--summary",
            help="Build summary JSON path (default: <output stem>.summary.json).",
            dir_okay=False,
        ),
    ] = None,
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
    if from_hpoa and fixture_dir is not None:
        raise typer.BadParameter("use either --from-hpoa or --fixture-dir, not both")
    if from_hpoa:
        sources = config.sources
        inputs = HpoaBuildInputs(
            hpo_json=_required_file(hpo_json or sources.hpo_json_path, "hp.json"),
            phenotype_hpoa=_required_file(
                hpoa or sources.phenotype_annotation_path, "phenotype.hpoa"
            ),
            genes_to_disease=_required_file(
                genes_to_disease or sources.genes_to_disease_path, "genes_to_disease"
            ),
            orphanet_ages=_optional_file(
                orphanet_ages or sources.orphanet_ages_path, "Orphanet ages"
            ),
            omim_orpha_map=_optional_file(
                omim_orpha_map or sources.omim_orpha_map_path, "OMIM-ORPHA map"
            ),
            exclude_pmids=_optional_file(
                exclude_pmids or sources.exclude_pmids_path, "PMID mask"
            ),
        )
        _build_profiles_from_hpoa(inputs, output or config.exports.profiles_path, summary)
        return
    if fixture_dir is None:
        typer.echo(f"Configured MVP diseases: {len(config.mvp.diseases)}")
        typer.echo(f"Profiles output: {config.exports.profiles_path}")
        typer.echo(
            "No profile source selected. Use --from-hpoa, or --fixture-dir for the fixture path."
        )
        return

    result = build_profiles_from_fixtures(fixture_dir)
    output_path = output or config.exports.profiles_path
    written_count = write_profiles_jsonl(output_path, result.profiles)

    typer.echo(f"Patch validation: {result.patch_validation.status}")
    typer.echo(f"Wrote {written_count} profile(s) to {output_path}")
    if result.warnings:
        typer.echo(f"Quality warnings: {', '.join(result.warnings)}")


def _build_profiles_from_hpoa(
    inputs: HpoaBuildInputs, output_path: Path, summary_path: Path | None
) -> None:
    if inputs.exclude_pmids is None:
        typer.echo(
            "WARNING: no --exclude-pmids given; profiles may include annotations from "
            "evaluation publications.",
            err=True,
        )
    started = time.perf_counter()
    result = build_profiles_from_hpoa(inputs)
    written = write_profiles_jsonl(output_path, result.profiles)
    elapsed = time.perf_counter() - started

    build_summary = {
        "build": {
            "command": "build-profiles --from-hpoa",
            "simulator_version": __version__,
            "simulator_git": git_revision(),
            "generated_at": datetime.now(tz=UTC).isoformat(timespec="seconds"),
            "wall_time_seconds": round(elapsed, 2),
        },
        **result.summary,
        "output": {
            "profiles_path": str(output_path),
            "profiles_written": written,
            "bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
        },
    }
    summary_path = summary_path or output_path.with_name(f"{output_path.stem}.summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(build_summary, indent=2) + "\n", encoding="utf-8")

    masking = result.summary["masking"]
    phenotypes = result.summary["phenotypes"]
    typer.echo(f"Wrote {written} profile(s) to {output_path}")
    typer.echo(
        f"Phenotype annotations: {phenotypes['annotations']}, "
        f"negatives: {result.summary['negatives'].get('annotations', 0)}, "
        f"genes: {result.summary['genes']['unique_symbols']}"
    )
    onset = result.summary["age_of_onset"]
    typer.echo(
        f"Onset coverage: {onset['diseases_with_onset']}/{written} disease(s), "
        f"{onset['gained_via_omim_mapping']} via OMIM-ORPHA mapping"
    )
    if masking["enabled"]:
        typer.echo(
            f"Masked {masking['rows_masked']} HPOA row(s) across "
            f"{masking['diseases_affected']} disease(s); "
            f"{masking['diseases_dropped_zero_positive']} disease(s) left without phenotypes."
        )
    typer.echo(f"Build summary: {summary_path}")


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
    hpo_json: Annotated[
        Path | None,
        typer.Option(
            "--hpo-json",
            help="HPO hp.json enabling ontology generalization (default: configured path).",
            dir_okay=False,
        ),
    ] = None,
    noise_vocabulary: Annotated[
        Path | None,
        typer.Option(
            "--noise-vocabulary",
            help="TSV with hpo_id and optional label columns of nonspecific noise terms.",
            dir_okay=False,
        ),
    ] = None,
    labels: Annotated[
        Path | None,
        typer.Option(
            "--labels",
            help='JSON {"genes": {symbol: int}, "diseases": {disease_id: int}}.',
            dir_okay=False,
        ),
    ] = None,
    disease_ids: Annotated[
        Path | None,
        typer.Option(
            "--disease-ids",
            help="Simulate only these diseases (one id per line); confounders still use all.",
            dir_okay=False,
        ),
    ] = None,
    sample_diseases: Annotated[
        int | None,
        typer.Option(
            "--sample-diseases",
            min=1,
            help="Simulate N diseases spread evenly over (inheritance, onset source) strata.",
        ),
    ] = None,
    cases_per_disease: Annotated[
        int | None,
        typer.Option("--cases-per-disease", min=1, help="Override cases per disease/difficulty."),
    ] = None,
    difficulty: Annotated[
        list[str] | None,
        typer.Option("--difficulty", help="Override the difficulties (repeatable)."),
    ] = None,
    seed: Annotated[int | None, typer.Option("--seed", help="Override the seed.")] = None,
    summary: Annotated[
        Path | None,
        typer.Option(
            "--summary",
            help="Run summary JSON path (default: <output stem>.summary.json).",
            dir_okay=False,
        ),
    ] = None,
) -> None:
    """Simulate synthetic patient cases from validated profiles."""

    config = _config_from_context(ctx)
    profiles_path = profiles or config.exports.profiles_path
    if not profiles_path.exists():
        raise typer.BadParameter(f"profiles file not found: {profiles_path}")
    started = time.perf_counter()

    overrides: dict[str, object] = {}
    if cases_per_disease is not None:
        overrides["cases_per_disease_per_difficulty"] = cases_per_disease
    if difficulty:
        overrides["difficulties"] = difficulty
    if seed is not None:
        overrides["seed"] = seed
    try:
        sim_config = config.simulation.model_validate(
            {**config.simulation.model_dump(), **overrides}
        )
    except ValidationError as exc:
        raise typer.BadParameter(f"invalid simulation override: {exc}") from exc

    hpo_path = hpo_json or config.sources.hpo_json_path
    if hpo_json is not None:
        _required_file(hpo_json, "hp.json")
    ontology = HpoOntology.from_json(hpo_path) if hpo_path.is_file() else None
    if ontology is None:
        typer.echo(
            "No hp.json found; simulating without ontology generalization, sex-specific "
            "anchors or ancestor/descendant checks on negatives."
        )
    noise_path = _optional_file(noise_vocabulary, "noise vocabulary")
    noise_terms = load_noise_vocabulary(noise_path, ontology) if noise_path else None
    label_maps = load_label_maps(_required_file(labels, "labels")) if labels else None

    profile_records = read_model_jsonl(profiles_path, DiseaseProfile)
    targets = profile_records
    strata: dict[str, int] | None = None
    if disease_ids is not None:
        wanted = {
            line.strip()
            for line in _required_file(disease_ids, "disease ids").read_text("utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }
        targets = [profile for profile in targets if profile.disease_id in wanted]
    if sample_diseases is not None:
        targets, strata = sample_diseases_by_stratum(targets, sample_diseases, sim_config.seed)

    confounders = ConfounderIndex(
        profile_records,
        ontology,
        top_n=sim_config.negatives.confounders_top_n,
        min_information_content=sim_config.negatives.min_information_content,
        unknown_frequency=sim_config.frequency.unknown_frequency,
    )
    sex_terms = SexSpecificTerms(ontology, sim_config.sex)
    inputs = {
        "profiles": _input_record(profiles_path),
        "hp_json": _input_record(hpo_path, ontology.version) if ontology else None,
        "noise_vocabulary": _input_record(noise_path) if noise_path else None,
        "labels": _input_record(labels) if labels else None,
        "disease_ids": _input_record(disease_ids) if disease_ids else None,
    }
    git = git_revision()
    source_versions = _source_versions(profile_records, inputs, git)

    cases: list[SyntheticCase] = []
    for profile in targets:
        cases.extend(
            simulate_cases(
                profile,
                sim_config,
                ontology=ontology,
                noise_vocabulary=noise_terms,
                confounders=confounders,
                sex_terms=sex_terms,
                gene_label=label_maps.genes.get(primary_gene(profile)) if label_maps else None,
                disease_label=label_maps.diseases.get(profile.disease_id) if label_maps else None,
                source_versions=source_versions,
            )
        )

    output_path = output or config.exports.rich_cases_path
    written = write_jsonl(output_path, cases)
    run_summary = {
        "simulate": {
            "simulator_version": SIMULATOR_VERSION,
            "package_version": __version__,
            "simulator_git": git,
            "generated_at": datetime.now(tz=UTC).isoformat(timespec="seconds"),
            "wall_time_seconds": round(time.perf_counter() - started, 2),
        },
        "inputs": {key: value for key, value in inputs.items() if value is not None},
        "config_hash": config_hash(sim_config),
        "config": sim_config.model_dump(mode="json"),
        "diseases": {
            "in_profiles": len(profile_records),
            "simulated": len(targets),
            "strata": strata,
        },
        "output": {
            "cases_path": str(output_path),
            "cases_written": written,
            "bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
        },
    }
    summary_path = summary or output_path.with_name(f"{output_path.stem}.summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
    typer.echo(
        f"Simulated {written} case(s) from {len(targets)} of {len(profile_records)} profile(s) "
        f"({sim_config.cases_per_disease_per_difficulty} per disease/difficulty, "
        f"difficulties={','.join(sim_config.difficulties)}, seed={sim_config.seed})."
    )
    typer.echo(f"Wrote rich cases to {output_path}")
    typer.echo(f"Run summary: {summary_path}")


def _input_record(path: Path, version: str | None = None) -> dict[str, object]:
    return {"path": str(path), "version": version, "sha256": sha256_file(path)}


def _source_versions(
    profiles: list[DiseaseProfile],
    inputs: dict[str, dict[str, object] | None],
    git: dict[str, object],
) -> dict[str, str]:
    """Compact per-case provenance: profile sources, run inputs and the simulator SHA."""

    versions: dict[str, str] = {}
    for item in profiles[0].provenance if profiles else []:
        source = item.source
        versions[source.name] = f"{source.version or '-'} sha256:{source.sha256 or '-'}"
    for key, record in inputs.items():
        if record is not None:
            version = f"{record['version']} " if record.get("version") else ""
            versions[key] = f"{version}sha256:{record['sha256']}"
    versions["simulator_git"] = f"{git.get('sha')}{'-dirty' if git.get('dirty') else ''}"
    return dict(sorted(versions.items()))


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
def validate(
    ctx: typer.Context,
    cases: Annotated[
        Path | None,
        typer.Option(
            "--cases",
            help="rich_cases.jsonl to validate; without it only the config is checked.",
            dir_okay=False,
        ),
    ] = None,
    profiles: Annotated[
        Path | None,
        typer.Option(
            "--profiles",
            help="Profiles the cases came from (enables calibration and profile checks).",
            dir_okay=False,
        ),
    ] = None,
    hpo_json: Annotated[
        Path | None,
        typer.Option(
            "--hpo-json",
            help="HPO hp.json for ancestor/descendant checks (default: configured path).",
            dir_okay=False,
        ),
    ] = None,
    noise_vocabulary: Annotated[
        Path | None,
        typer.Option(
            "--noise-vocabulary", help="Noise vocabulary TSV the cases used.", dir_okay=False
        ),
    ] = None,
    report: Annotated[
        Path | None,
        typer.Option(
            "--report",
            help="Report JSON path (default: <cases stem>.validation.json).",
            dir_okay=False,
        ),
    ] = None,
    run_summary: Annotated[
        Path | None,
        typer.Option(
            "--run-summary",
            help="simulate's run summary, whose config is used for the priors "
            "(default: <cases stem>.summary.json when present, else the app config).",
            dir_okay=False,
        ),
    ] = None,
) -> None:
    """Validate the config, or a simulated case set against its profiles and priors.

    Exits with status 1 when any case invariant is violated.
    """

    config = _config_from_context(ctx)
    if cases is None:
        genes = ", ".join(disease.gene for disease in config.mvp.diseases)
        typer.echo(f"Config OK: {len(config.mvp.diseases)} MVP diseases ({genes})")
        return

    cases_path = _required_file(cases, "rich cases")
    hpo_path = hpo_json or config.sources.hpo_json_path
    if hpo_json is not None:
        _required_file(hpo_json, "hp.json")
    ontology = HpoOntology.from_json(hpo_path) if hpo_path.is_file() else None
    if ontology is None:
        typer.echo("No hp.json found; negatives are only checked for identical terms.")
    profile_map = (
        {
            profile.disease_id: profile
            for profile in iter_model_jsonl(_required_file(profiles, "profiles"), DiseaseProfile)
        }
        if profiles is not None
        else None
    )
    noise_ids = (
        {term.hpo_id for term in load_noise_vocabulary(_required_file(noise_vocabulary, "noise"))}
        if noise_vocabulary is not None
        else None
    )

    summary_path = run_summary or cases_path.with_name(f"{cases_path.stem}.summary.json")
    sim_config, config_source = config.simulation, "app config"
    if run_summary is not None or summary_path.is_file():
        run_data = json.loads(_required_file(summary_path, "run summary").read_text("utf-8"))
        sim_config = SimulationConfig.model_validate(run_data["config"])
        config_source = str(summary_path)

    result = validate_cases(
        iter_model_jsonl(cases_path, SyntheticCase),
        profiles=profile_map,
        ontology=ontology,
        config=sim_config,
        noise_vocabulary=noise_ids,
    )
    result = {
        "config_source": config_source,
        "inputs": {
            "cases": _input_record(cases_path),
            "profiles": _input_record(profiles) if profiles else None,
            "hp_json": _input_record(hpo_path, ontology.version) if ontology else None,
            "noise_vocabulary": _input_record(noise_vocabulary) if noise_vocabulary else None,
        },
        **result,
    }
    report_path = report or cases_path.with_name(f"{cases_path.stem}.validation.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    typer.echo(format_report(result))
    typer.echo(f"Validation report: {report_path}")
    if result["violations"]["count"]:
        raise typer.Exit(code=1)


def main() -> None:
    """Run the CLI."""

    app()


if __name__ == "__main__":
    main()
