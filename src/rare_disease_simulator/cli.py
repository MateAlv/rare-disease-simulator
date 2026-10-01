"""Command-line entrypoint for rare-disease-simulator."""

from __future__ import annotations

import json
import logging
import time
from collections import Counter
from collections.abc import Iterable, Iterator
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
from rare_disease_simulator.data_sources.gene_profiles import (
    GenePlan,
    plan_genes,
    read_gene_profiles,
    read_gnn_genes,
)
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.data_sources.hpo_annotations import download_hpo_release
from rare_disease_simulator.data_sources.http_client import HttpClient
from rare_disease_simulator.exports.graphens import write_graphens_json
from rare_disease_simulator.exports.jsonl import iter_model_jsonl, read_model_jsonl, write_jsonl
from rare_disease_simulator.exports.training import (
    EXPORT_FORMAT,
    RECORD_FIELDS,
    SPLIT_RULE,
    ExportError,
    card_path_for,
    dataset_counts,
    encode_dataset,
    parse_split,
    sort_records,
    training_record,
)
from rare_disease_simulator.profiles.builder import (
    build_profiles_from_fixtures,
    write_profiles_jsonl,
)
from rare_disease_simulator.profiles.frequency import pooled_counts
from rare_disease_simulator.profiles.hpoa_builder import HpoaBuildInputs, build_profiles_from_hpoa
from rare_disease_simulator.profiles.merge import MergeStats, merge_entity_profiles
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.simulation.calibration import (
    CalibrationError,
    apply_calibration,
    format_calibration_keys,
    load_calibration,
)
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.inputs import load_label_maps, load_noise_vocabulary
from rare_disease_simulator.simulation.reporting import (
    CardinalIndex,
    PresentationAges,
    Reporting,
    ReportingArtifactError,
    ReportModel,
    load_cardinal,
    load_presentation_ages,
    load_report_model,
)
from rare_disease_simulator.simulation.sampling import sample_diseases as sample_diseases_by_stratum
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase, run_config
from rare_disease_simulator.simulation.simulator import (
    FORCED_REASON,
    MERGED_REASON,
    RELATED_NOISE_REASON,
    SIMULATOR_VERSION,
    SPECIALIZED_REASON,
    SexSpecificTerms,
    config_hash,
    primary_gene,
    simulate_cases,
    simulate_gene_cases,
)
from rare_disease_simulator.validation.cases import format_report, validate_cases
from rare_disease_simulator.validation.dataset import format_dataset_report, validate_dataset

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
    gene_profiles: Annotated[
        Path | None,
        typer.Option(
            "--gene-profiles",
            help="Gene profiles (genes-v1); also build diseases their simulable entities "
            "list that genes_to_disease does not link to a gene.",
            dir_okay=False,
        ),
    ] = None,
    drop_annotations: Annotated[
        Path | None,
        typer.Option(
            "--drop-annotations",
            help="TSV (disease_id, hpo_id) of positive annotations to hold out of the profiles.",
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
            gene_profiles=_optional_file(gene_profiles, "gene profiles"),
            drop_annotations=_optional_file(drop_annotations, "annotation holdout"),
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
    holdout = result.summary["annotation_holdout"]
    if holdout["enabled"]:
        typer.echo(
            f"Held out {holdout['pairs_dropped']} of {holdout['pairs_listed']} listed "
            f"annotation(s) ({holdout['rows_dropped']} row(s)) across "
            f"{holdout['diseases_affected']} disease(s); "
            f"{holdout['diseases_dropped_zero_positive']} disease(s) left without phenotypes."
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
            help="TSV with hpo_id and optional label and weight columns of noise terms.",
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
    genes: Annotated[
        str | None,
        typer.Option(
            "--genes",
            help="Gene-first mode: 'all' GNN genes, a comma-separated list of v3 symbols, "
            "or a file with one symbol per line.",
        ),
    ] = None,
    gene_profiles: Annotated[
        Path | None,
        typer.Option(
            "--gene-profiles",
            help="Gene profiles (genes-v1 genes.json.gz); implies gene-first mode "
            "(default: configured path).",
            dir_okay=False,
        ),
    ] = None,
    gnn_genes: Annotated[
        Path | None,
        typer.Option(
            "--gnn-genes",
            help="GNN gene vocabulary (JSON list; its order is the class index).",
            dir_okay=False,
        ),
    ] = None,
    cases_per_gene: Annotated[
        int | None,
        typer.Option(
            "--cases-per-gene", min=1, help="Gene-first mode: cases per gene and difficulty."
        ),
    ] = None,
    calibration: Annotated[
        Path | None,
        typer.Option(
            "--calibration",
            help="JSON of dotted config paths to override (see calibration-keys).",
            dir_okay=False,
        ),
    ] = None,
    report_model: Annotated[
        Path | None,
        typer.Option(
            "--report-model",
            help="report-model-v1 JSON (term budget and reporting model; default: "
            "sources.report_model_path). Needed when reporting.mode is report_model.",
            dir_okay=False,
        ),
    ] = None,
    cardinal: Annotated[
        Path | None,
        typer.Option(
            "--cardinal",
            help="cardinal-v1 TSV (disease_id, hpo_id, kind, source; default: "
            "sources.cardinal_path). Needed when reporting.force_cardinal is true.",
            dir_okay=False,
        ),
    ] = None,
    presentation_ages: Annotated[
        Path | None,
        typer.Option(
            "--presentation-ages",
            help="TSV (disease_id, age_low_years, age_high_years, source) of age-at-"
            "presentation ranges; a listed disease draws its age from it (default: "
            "sources.presentation_ages_path).",
            dir_okay=False,
        ),
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
    """Simulate synthetic patient cases from validated profiles.

    Disease-first by default; with --genes, --gene-profiles or --cases-per-gene,
    one set of cases per GNN gene, labelled with the gene's class index.
    """

    config = _config_from_context(ctx)
    profiles_path = profiles or config.exports.profiles_path
    if not profiles_path.exists():
        raise typer.BadParameter(f"profiles file not found: {profiles_path}")
    started = time.perf_counter()
    gene_first = genes is not None or gene_profiles is not None or cases_per_gene is not None
    if gene_first:
        clashing = [
            name
            for name, value in (
                ("--labels", labels),
                ("--disease-ids", disease_ids),
                ("--sample-diseases", sample_diseases),
                ("--cases-per-disease", cases_per_disease),
            )
            if value is not None
        ]
        if clashing:
            raise typer.BadParameter(f"gene-first mode does not take {', '.join(clashing)}")

    overrides: dict[str, object] = {}
    per_target = cases_per_gene if gene_first else cases_per_disease
    if per_target is not None:
        overrides["cases_per_disease_per_difficulty"] = per_target
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

    calibration_data = None
    if calibration is not None:
        try:
            calibration_data = load_calibration(_required_file(calibration, "calibration"))
            sim_config = apply_calibration(sim_config, calibration_data.overrides)
        except CalibrationError as exc:
            raise typer.BadParameter(str(exc)) from exc
        if calibration_data.noise_vocabulary is not None and noise_vocabulary is not None:
            raise typer.BadParameter(
                "the calibration file sets noise.vocabulary_path; drop --noise-vocabulary"
            )

    hpo_path = hpo_json or config.sources.hpo_json_path
    if hpo_json is not None:
        _required_file(hpo_json, "hp.json")
    ontology = HpoOntology.from_json(hpo_path) if hpo_path.is_file() else None
    if ontology is None:
        typer.echo(
            "No hp.json found; simulating without ontology generalization, sex-specific "
            "anchors or ancestor/descendant checks on negatives."
        )
    noise_path = _optional_file(noise_vocabulary, "noise vocabulary") or (
        calibration_data.noise_vocabulary if calibration_data else None
    )
    noise_terms = load_noise_vocabulary(noise_path, ontology) if noise_path else None
    label_maps = load_label_maps(_required_file(labels, "labels")) if labels else None

    profile_records = read_model_jsonl(profiles_path, DiseaseProfile)
    _check_count_estimator(sim_config, profile_records)
    reporting, reporting_inputs = _load_reporting(
        config, sim_config, report_model, cardinal, ontology, profile_records
    )
    ages, ages_inputs = _load_presentation_ages(config, presentation_ages)
    reporting_inputs.update(ages_inputs)
    confounders = ConfounderIndex(
        profile_records,
        ontology,
        top_n=sim_config.negatives.confounders_top_n,
        min_information_content=sim_config.negatives.min_information_content,
        unknown_frequency=sim_config.frequency.unknown_frequency,
    )
    sex_terms = SexSpecificTerms(ontology, sim_config.sex)
    inputs: dict[str, dict[str, object] | None] = {
        "profiles": _input_record(profiles_path),
        "hp_json": _input_record(hpo_path, ontology.version) if ontology else None,
        "noise_vocabulary": _input_record(noise_path) if noise_path else None,
        "labels": _input_record(labels) if labels else None,
        "disease_ids": _input_record(disease_ids) if disease_ids else None,
        "calibration": _input_record(calibration) if calibration else None,
        **reporting_inputs,
    }

    plan = None
    profile_map = {profile.disease_id: profile for profile in profile_records}
    targets = profile_records
    strata: dict[str, int] | None = None
    if gene_first:
        plan, gene_inputs = _gene_plan(config, gene_profiles, gnn_genes, genes, profile_records)
        inputs.update(gene_inputs)
    else:
        if disease_ids is not None:
            wanted = {
                line.strip()
                for line in _required_file(disease_ids, "disease ids")
                .read_text("utf-8")
                .splitlines()
                if line.strip() and not line.startswith("#")
            }
            targets = [profile for profile in targets if profile.disease_id in wanted]
        if sample_diseases is not None:
            targets, strata = sample_diseases_by_stratum(targets, sample_diseases, sim_config.seed)

    git = git_revision()
    source_versions = _source_versions(profile_records, inputs, git)
    options = {
        "ontology": ontology,
        "noise_vocabulary": noise_terms,
        "confounders": confounders,
        "sex_terms": sex_terms,
        "source_versions": source_versions,
        "reporting": reporting,
        "presentation_ages": ages,
    }
    if plan is not None:
        cases: Iterable[SyntheticCase] = (
            case
            for target in plan.targets
            for case in simulate_gene_cases(target, profile_map, sim_config, **options)
        )
    else:
        cases = (
            case
            for profile in targets
            for case in simulate_cases(
                profile,
                sim_config,
                gene_label=label_maps.genes.get(primary_gene(profile)) if label_maps else None,
                disease_label=label_maps.diseases.get(profile.disease_id) if label_maps else None,
                **options,
            )
        )

    output_path = output or config.exports.rich_cases_path
    realized = _ReportedTotals(
        sim_config.reporting.mode,
        sim_config.reporting.profile_budget()
        if sim_config.reporting.budget_normalize
        else None,
    )
    written = write_jsonl(output_path, realized.track(cases))
    run_summary: dict[str, object] = {
        "simulate": {
            "mode": "gene_first" if gene_first else "disease_first",
            "simulator_version": SIMULATOR_VERSION,
            "package_version": __version__,
            "simulator_git": git,
            "generated_at": datetime.now(tz=UTC).isoformat(timespec="seconds"),
            "wall_time_seconds": round(time.perf_counter() - started, 2),
        },
        "inputs": {key: value for key, value in inputs.items() if value is not None},
        "profile_sources": _profile_sources(profile_records),
        "calibration": calibration_data.record() if calibration_data else None,
        "config_hash": config_hash(sim_config),
        "config": sim_config.model_dump(mode="json"),
        "diseases": {
            "in_profiles": len(profile_records),
            "simulated": None if gene_first else len(targets),
            "strata": strata,
        },
        "genes": {**plan.summary(), "requested": genes or "all"} if plan is not None else None,
        "entity_profiles": (
            _merge_summary(plan, profile_map, ontology)
            if plan is not None and sim_config.entity_profiles == "merged"
            else None
        ),
        "reporting": (
            {"mode": sim_config.reporting.mode, **reporting.summary(), **realized.summary()}
            if reporting is not None
            else {"mode": sim_config.reporting.mode}
        ),
        "noise": realized.noise_summary() if reporting is not None else None,
        "presentation_ages": (
            {"diseases": len(ages.ranges), "cases": realized.counts["presentation_age_cases"]}
            if ages is not None
            else None
        ),
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
    per = f"{sim_config.cases_per_disease_per_difficulty} per {'gene' if gene_first else 'disease'}"
    run_shape = (
        f"({per}/difficulty, difficulties={','.join(sim_config.difficulties)}, "
        f"seed={sim_config.seed})."
    )
    if plan is not None:
        skipped = sum(len(symbols) for symbols in plan.skipped.values())
        typer.echo(
            f"Simulated {written} case(s) for {len(plan.targets)} gene(s); "
            f"skipped {skipped} of {plan.vocabulary_size - 1} vocabulary gene(s) {run_shape}"
        )
    else:
        typer.echo(
            f"Simulated {written} case(s) from {len(targets)} of {len(profile_records)} "
            f"profile(s) {run_shape}"
        )
    typer.echo(f"Wrote rich cases to {output_path}")
    typer.echo(f"Run summary: {summary_path}")


class _ReportedTotals:
    """Realized report budgets and reported terms of the cases as they are written."""

    def __init__(self, mode: str = "report_model", histogram: dict[int, int] | None = None) -> None:
        self.mode = mode
        positive = {k: n for k, n in (histogram or {}).items() if k >= 1 and n > 0}
        total = sum(positive.values())
        self.histogram = {k: n / total for k, n in positive.items()} if total else None
        self.total_counts: Counter[int] = Counter()
        self.counts: Counter[str] = Counter()
        self.q_sum = 0.0
        self.k_counts: Counter[int] = Counter()
        self.realized_counts: Counter[int] = Counter()
        self.scale_sum = 0.0

    def track(self, cases: Iterable[SyntheticCase]) -> Iterator[SyntheticCase]:
        for case in cases:
            if case.metadata.presentation_age_years is not None:
                self.counts["presentation_age_cases"] += 1
            if case.metadata.report_budget is not None or self.mode == "independent":
                self.counts["cases"] += 1
                self.counts["budget"] += case.metadata.report_budget or 0
                if self.mode == "independent":
                    self._track_independent(case)
                self.counts["profile"] += len(case.positive_phenotypes)
                self.counts["noise"] += len(case.noise_phenotypes)
                self.counts["specialized"] += sum(
                    1 for p in case.positive_phenotypes if p.reason == SPECIALIZED_REASON
                )
                self.counts["related"] += sum(
                    1 for p in case.noise_phenotypes if p.reason == RELATED_NOISE_REASON
                )
            yield case

    def _track_independent(self, case: SyntheticCase) -> None:
        shown = [p for p in case.positive_phenotypes if p.simulated_origin == "disease_profile"]
        merged = [p for p in case.missing_phenotypes if p.reason == MERGED_REASON]
        unreported = [p for p in case.missing_phenotypes if p.reason != MERGED_REASON]
        unreported += case.unknown_phenotypes
        self.counts["true_terms"] += len(shown) + len(merged) + len(unreported)
        self.counts["true_terms_reported"] += len(shown) + len(merged)
        rolled = [
            p.report_probability
            for p in [*shown, *merged, *unreported]
            if p.report_probability is not None
        ]
        self.counts["q_terms"] += len(rolled)
        self.counts["q_capped"] += sum(1 for q in rolled if q >= 1.0)
        self.q_sum += sum(rolled)
        forced = sum(1 for p in case.positive_phenotypes if p.reason == FORCED_REASON)
        self.counts["forced_min_one"] += forced
        if case.metadata.report_budget is not None:
            self.k_counts[case.metadata.report_budget] += 1
            self.realized_counts[len(shown) + len(merged)] += 1
            self.scale_sum += case.metadata.report_scale or 0.0
            if case.metadata.report_budget_profile is not None:
                self.counts["record_cases"] += 1
                self.counts["profile_target"] += case.metadata.report_budget_profile
                self.total_counts[len(shown) + len(merged) + len(case.noise_phenotypes)] += 1

    def summary(self) -> dict[str, object]:
        cases = self.counts["cases"]

        def mean(key: str) -> float | None:
            return round(self.counts[key] / cases, 4) if cases else None

        if self.mode == "independent":
            return {
                "cases": cases,
                "reported_profile_mean": mean("profile"),
                "reported_noise_mean": mean("noise"),
                "reported_total_mean": (
                    round((self.counts["profile"] + self.counts["noise"]) / cases, 4)
                    if cases
                    else None
                ),
                "true_terms_reported_share": self._share("true_terms_reported", "true_terms"),
                "q_mean": (
                    round(self.q_sum / self.counts["q_terms"], 4)
                    if self.counts["q_terms"]
                    else None
                ),
                "q_capped_share": self._share("q_capped", "q_terms"),
                "forced_min_one": self.counts["forced_min_one"],
                "budget_normalized": self._normalized_summary(),
                "specialized": {
                    "terms": self.counts["specialized"],
                    "share_of_profile_terms": self._share("specialized", "profile"),
                },
            }
        return {
            "budget_mean_drawn": mean("budget"),
            "reported_total_mean": (
                round((self.counts["profile"] + self.counts["noise"]) / cases, 4)
                if cases
                else None
            ),
            "reported_profile_mean": mean("profile"),
            "reported_noise_mean": mean("noise"),
            "specialized": {
                "terms": self.counts["specialized"],
                "share_of_profile_terms": self._share("specialized", "profile"),
            },
        }

    def _normalized_summary(self) -> dict[str, object] | None:
        cases = sum(self.k_counts.values())
        if not cases:
            return None
        k_shares = {k: n / cases for k, n in self.k_counts.items()}
        realized = {k: n / cases for k, n in self.realized_counts.items()}
        variation = 0.5 * sum(
            abs(k_shares.get(k, 0.0) - realized.get(k, 0.0)) for k in set(k_shares) | set(realized)
        )
        return {
            "cases": cases,
            "k_mean": round(sum(k * n for k, n in self.k_counts.items()) / cases, 4),
            "realized_profile_mean": round(
                sum(k * n for k, n in self.realized_counts.items()) / cases, 4
            ),
            "total_variation_k_vs_realized": round(variation, 4),
            "scale_mean": round(self.scale_sum / cases, 4),
            "record": self._record_summary(cases),
        }

    def _record_summary(self, cases: int) -> dict[str, object] | None:
        record_cases = self.counts["record_cases"]
        if not record_cases:
            return None
        totals = {k: n / record_cases for k, n in self.total_counts.items()}
        reference = self.histogram or {}
        variation = 0.5 * sum(
            abs(totals.get(k, 0.0) - reference.get(k, 0.0)) for k in set(totals) | set(reference)
        )
        return {
            "k_mean": round(sum(k * n for k, n in self.k_counts.items()) / cases, 4),
            "noise_mean": round(self.counts["noise"] / record_cases, 4),
            "profile_target_mean": round(self.counts["profile_target"] / record_cases, 4),
            "realized_profile_mean": round(
                sum(k * n for k, n in self.realized_counts.items()) / cases, 4
            ),
            "total_shown_mean": round(
                sum(k * n for k, n in self.total_counts.items()) / record_cases, 4
            ),
            "total_variation_total_vs_histogram": round(variation, 4),
        }

    def noise_summary(self) -> dict[str, object]:
        return {
            "terms": self.counts["noise"],
            "related": self.counts["related"],
            "related_share": self._share("related", "noise"),
        }

    def _share(self, part: str, whole: str) -> float | None:
        total = self.counts[whole]
        return round(self.counts[part] / total, 4) if total else None


def _check_count_estimator(config: SimulationConfig, profiles: list[DiseaseProfile]) -> None:
    settings = config.frequency
    if settings.count_estimator != "beta_shrinkage" or settings.shrinkage_mean is not None:
        return
    counted = sum(
        1
        for profile in profiles
        for phenotype in profile.phenotypes
        if pooled_counts(phenotype.frequency_raw) is not None
    )
    if counted:
        raise typer.BadParameter(
            f"{counted} profile term(s) are count-based and frequency.count_estimator is "
            "beta_shrinkage, but frequency.shrinkage_mean is not set: give the median "
            "count-based frequency of hpoa-recount-v1 (config or calibration key), or set "
            "frequency.count_estimator: jeffreys for simulator 0.3 estimates"
        )


def _pinned_record(path: Path, pinned: str | None, from_config: bool, key: str) -> dict:
    record = _input_record(path)
    if from_config and pinned:
        if record["sha256"] != pinned:
            raise typer.BadParameter(
                f"{path} has sha256 {record['sha256']}, not the pinned {pinned} ({key})"
            )
        record["sha256_verified_against"] = key
    return record


def _load_reporting(
    config: AppConfig,
    sim_config: SimulationConfig,
    report_model: Path | None,
    cardinal: Path | None,
    ontology: HpoOntology | None,
    profiles: list[DiseaseProfile],
) -> tuple[Reporting | None, dict[str, dict[str, object] | None]]:
    """Load report-model-v1 and cardinal-v1 when the emission mode needs them."""

    sources = config.sources
    mode = sim_config.reporting.mode
    if mode == "observation":
        if report_model is not None or cardinal is not None:
            raise typer.BadParameter(
                "--report-model/--cardinal given but reporting.mode is observation"
            )
        return None, {}
    model_path = report_model or sources.report_model_path
    if model_path is None:
        raise typer.BadParameter(
            f"reporting.mode is {mode}: pass --report-model (report-model-v1) or set "
            "sources.report_model_path, or set reporting.mode: observation"
        )
    cardinal_path = cardinal or sources.cardinal_path
    if mode == "independent" and sim_config.reporting.force_cardinal:
        typer.echo(
            "Warning: reporting.force_cardinal is ignored in independent mode; cardinal terms "
            "are reported through the report model's cardinal coefficient.",
            err=True,
        )
    if cardinal_path is None and mode == "report_model" and sim_config.reporting.force_cardinal:
        raise typer.BadParameter(
            "reporting.force_cardinal needs --cardinal (cardinal-v1) or sources.cardinal_path; "
            "set reporting.force_cardinal: false to report without cardinal terms"
        )
    model_record = _pinned_record(
        _required_file(model_path, "report model"),
        sources.report_model_sha256,
        report_model is None,
        "sources.report_model_sha256",
    )
    cardinal_record = (
        _pinned_record(
            _required_file(cardinal_path, "cardinal terms"),
            sources.cardinal_sha256,
            cardinal is None,
            "sources.cardinal_sha256",
        )
        if cardinal_path is not None
        else None
    )
    try:
        model = load_report_model(model_path)
        cardinal_index = load_cardinal(cardinal_path, ontology) if cardinal_path else None
        reporting = Reporting.build(
            model,
            ontology=ontology,
            profiles=profiles,
            cardinal=cardinal_index,
            unknown_frequency=sim_config.frequency.unknown_frequency,
        )
    except ReportingArtifactError as exc:
        raise typer.BadParameter(str(exc)) from exc
    if (
        mode == "independent"
        and cardinal_path is None
        and any(feature.kind == "cardinal_flag" for feature in model.features)
    ):
        raise typer.BadParameter(
            "the report model has a cardinal_flag feature: pass --cardinal (cardinal-v1) or "
            "set sources.cardinal_path"
        )
    model_record["artifact_id"] = model.artifact_id
    return reporting, {"report_model": model_record, "cardinal": cardinal_record}


def _load_presentation_ages(
    config: AppConfig, given: Path | None
) -> tuple[PresentationAges | None, dict[str, dict[str, object] | None]]:
    path = given or config.sources.presentation_ages_path
    if path is None:
        return None, {}
    record = _pinned_record(
        _required_file(path, "presentation ages"),
        config.sources.presentation_ages_sha256,
        given is None,
        "sources.presentation_ages_sha256",
    )
    try:
        ages = load_presentation_ages(path)
    except ReportingArtifactError as exc:
        raise typer.BadParameter(str(exc)) from exc
    return ages, {"presentation_ages": record}


def _merge_summary(
    plan: GenePlan, profiles: dict[str, DiseaseProfile], ontology: HpoOntology | None
) -> dict[str, int]:
    """What merging did over the plan's distinct entities (run summary ``entity_profiles``)."""

    total = MergeStats()
    entities = {
        (entity.entity, entity.profile_ids): entity
        for target in plan.targets
        for entity in target.entities
    }
    for (entity_id, profile_ids), _ in sorted(entities.items()):
        _, stats = merge_entity_profiles(
            entity_id, [profiles[pid] for pid in profile_ids], ontology
        )
        total.add(stats)
    return total.summary()


def _gene_plan(
    config: AppConfig,
    gene_profiles: Path | None,
    gnn_genes: Path | None,
    genes: str | None,
    profile_records: list[DiseaseProfile],
) -> tuple[GenePlan, dict[str, dict[str, object] | None]]:
    sources = config.sources
    profiles_path = gene_profiles or sources.gene_profiles_path
    vocabulary_path = gnn_genes or sources.gnn_genes_path
    if profiles_path is None:
        raise typer.BadParameter("gene-first mode needs --gene-profiles")
    if vocabulary_path is None:
        raise typer.BadParameter("gene-first mode needs --gnn-genes")
    _required_file(profiles_path, "gene profiles")
    _required_file(vocabulary_path, "GNN gene vocabulary")
    profiles_record = _input_record(profiles_path)
    if gene_profiles is None and sources.gene_profiles_sha256:
        if profiles_record["sha256"] != sources.gene_profiles_sha256:
            raise typer.BadParameter(
                f"gene profiles {profiles_path} have sha256 {profiles_record['sha256']}, "
                f"not the pinned {sources.gene_profiles_sha256}"
            )
        profiles_record["sha256_verified_against"] = "sources.gene_profiles_sha256"
    try:
        gene_data = read_gene_profiles(profiles_path)
        vocabulary = read_gnn_genes(vocabulary_path)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    wanted = None
    if genes is not None and genes != "all":
        candidate = Path(genes)
        entries = (
            candidate.read_text("utf-8").splitlines() if candidate.is_file() else genes.split(",")
        )
        wanted = {
            entry.strip() for entry in entries if entry.strip() and not entry.startswith("#")
        }
        unknown = sorted(wanted - set(vocabulary))
        if unknown:
            raise typer.BadParameter(f"not in the GNN vocabulary: {', '.join(unknown[:10])}")
    plan = plan_genes(
        gene_data, vocabulary, {profile.disease_id for profile in profile_records}, wanted
    )
    if not plan.targets:
        raise typer.BadParameter("no requested gene has a simulable entity with a profile")
    return plan, {
        "gene_profiles": profiles_record,
        "gnn_genes": _input_record(vocabulary_path),
    }


def _profile_sources(profiles: list[DiseaseProfile]) -> dict[str, dict[str, str | None]]:
    """The sources the profiles were built from (incl. the held-out PMID mask)."""

    sources: dict[str, dict[str, str | None]] = {}
    for profile in profiles:
        for item in profile.provenance:
            source = item.source
            sources.setdefault(
                source.name,
                {"file": source.url_or_file, "version": source.version, "sha256": source.sha256},
            )
    return dict(sorted(sources.items()))


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
    dataset: Annotated[
        Path | None,
        typer.Option(
            "--dataset",
            help="An export-training dataset (.jsonl.gz) to check against its card.",
            dir_okay=False,
        ),
    ] = None,
    card: Annotated[
        Path | None,
        typer.Option(
            "--card", help="Dataset card (default: <dataset>.card.json).", dir_okay=False
        ),
    ] = None,
    gnn_genes: Annotated[
        Path | None,
        typer.Option(
            "--gnn-genes",
            help="GNN gene vocabulary to check gene indices against (dataset mode).",
            dir_okay=False,
        ),
    ] = None,
    report_model: Annotated[
        Path | None,
        typer.Option(
            "--report-model",
            help="report-model-v1 the cases used (default: the run summary's).",
            dir_okay=False,
        ),
    ] = None,
    cardinal: Annotated[
        Path | None,
        typer.Option(
            "--cardinal",
            help="cardinal-v1 the cases used (default: the run summary's).",
            dir_okay=False,
        ),
    ] = None,
) -> None:
    """Validate the config, a simulated case set, or an exported training dataset.

    Exits with status 1 when any case or dataset invariant is violated.
    """

    config = _config_from_context(ctx)
    if dataset is not None:
        if cases is not None:
            raise typer.BadParameter("use either --cases or --dataset, not both")
        hpo_path = hpo_json or config.sources.hpo_json_path
        if hpo_json is not None:
            _required_file(hpo_json, "hp.json")
        ontology = HpoOntology.from_json(hpo_path) if hpo_path.is_file() else None
        if ontology is None:
            typer.echo("No hp.json found; excluded terms are only checked against true terms.")
        _validate_dataset(dataset, card, gnn_genes, report, ontology)
        return
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
    summary_path = run_summary or cases_path.with_name(f"{cases_path.stem}.summary.json")
    sim_config, config_source = config.simulation, "app config"
    run_data: dict | None = None
    if run_summary is not None or summary_path.is_file():
        run_data = json.loads(_required_file(summary_path, "run summary").read_text("utf-8"))
        sim_config = run_config(run_data["config"])
        config_source = str(summary_path)
    if noise_vocabulary is None and run_data is not None:
        recorded = (run_data.get("inputs") or {}).get("noise_vocabulary")
        if recorded and Path(recorded["path"]).is_file():
            noise_vocabulary = Path(recorded["path"])
    noise_ids = (
        {term.hpo_id for term in load_noise_vocabulary(_required_file(noise_vocabulary, "noise"))}
        if noise_vocabulary is not None
        else None
    )
    model, cardinal_index, artifact_inputs = _reporting_artifacts(
        run_data, report_model, cardinal, ontology
    )

    result = validate_cases(
        iter_model_jsonl(cases_path, SyntheticCase),
        profiles=profile_map,
        ontology=ontology,
        config=sim_config,
        noise_vocabulary=noise_ids,
        report_model=model,
        cardinal=cardinal_index,
    )
    result = {
        "config_source": config_source,
        "inputs": {
            "cases": _input_record(cases_path),
            "profiles": _input_record(profiles) if profiles else None,
            "hp_json": _input_record(hpo_path, ontology.version) if ontology else None,
            "noise_vocabulary": _input_record(noise_vocabulary) if noise_vocabulary else None,
            **artifact_inputs,
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


def _reporting_artifacts(
    run_data: dict | None,
    report_model: Path | None,
    cardinal: Path | None,
    ontology: HpoOntology | None,
) -> tuple[ReportModel | None, CardinalIndex | None, dict[str, dict[str, object] | None]]:
    """The report model and cardinal terms to validate against: given, else the run's."""

    recorded = (run_data or {}).get("inputs") or {}

    def pick(given: Path | None, key: str) -> Path | None:
        if given is not None:
            return _required_file(given, key.replace("_", " "))
        record = recorded.get(key)
        return _required_file(Path(record["path"]), key.replace("_", " ")) if record else None

    model_path = pick(report_model, "report_model")
    cardinal_path = pick(cardinal, "cardinal")
    try:
        model = load_report_model(model_path) if model_path else None
        cardinal_index = load_cardinal(cardinal_path, ontology) if cardinal_path else None
    except ReportingArtifactError as exc:
        raise typer.BadParameter(str(exc)) from exc
    return (
        model,
        cardinal_index,
        {
            "report_model": _input_record(model_path) if model_path else None,
            "cardinal": _input_record(cardinal_path) if cardinal_path else None,
        },
    )


def _validate_dataset(
    dataset: Path,
    card: Path | None,
    gnn_genes: Path | None,
    report: Path | None,
    ontology: HpoOntology | None,
) -> None:
    dataset_path = _required_file(dataset, "dataset")
    card_path = _required_file(card or card_path_for(dataset_path), "dataset card")
    card_data = json.loads(card_path.read_text("utf-8"))
    vocabulary = read_gnn_genes(_required_file(gnn_genes, "GNN genes")) if gnn_genes else None
    result = {
        "card": _input_record(card_path),
        "gnn_genes": _input_record(gnn_genes) if gnn_genes else None,
        "hp_json": {"version": ontology.version} if ontology else None,
        **validate_dataset(dataset_path, card_data, vocabulary, ontology),
    }
    report_path = report or card_path.with_name(
        card_path.name.removesuffix(".card.json") + ".dataset-validation.json"
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    typer.echo(format_dataset_report(result))
    typer.echo(f"Validation report: {report_path}")
    if result["violations"]["count"]:
        raise typer.Exit(code=1)


@app.command("export-training")
def export_training(
    ctx: typer.Context,
    cases: Annotated[
        Path, typer.Option("--cases", help="rich_cases.jsonl from simulate.", dir_okay=False)
    ],
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Dataset path (.jsonl.gz).", dir_okay=False),
    ],
    card: Annotated[
        Path | None,
        typer.Option(
            "--card", help="Dataset card path (default: <output>.card.json).", dir_okay=False
        ),
    ] = None,
    run_summary: Annotated[
        Path | None,
        typer.Option(
            "--run-summary",
            help="simulate's run summary (default: <cases stem>.summary.json).",
            dir_okay=False,
        ),
    ] = None,
    profiles: Annotated[
        Path | None,
        typer.Option(
            "--profiles",
            help="Profiles for validation (default: the run summary's).",
            dir_okay=False,
        ),
    ] = None,
    hpo_json: Annotated[
        Path | None,
        typer.Option(
            "--hpo-json", help="hp.json for validation (default: the run summary's).",
            dir_okay=False,
        ),
    ] = None,
    split: Annotated[
        str,
        typer.Option("--split", help="Sim train,val,test fractions, assigned by case seed."),
    ] = "0.90,0.05,0.05",
    allow_dirty: Annotated[
        bool,
        typer.Option(
            "--allow-dirty",
            help="Export even if this tree or the simulate run had uncommitted changes.",
        ),
    ] = False,
) -> None:
    """Write the export-v2 training dataset (gzip JSONL) and its dataset card.

    The cases are validated first; any invariant violation stops the export.
    """

    _config_from_context(ctx)
    started = time.perf_counter()
    cases_path = _required_file(cases, "rich cases")
    summary_path = _required_file(
        run_summary or cases_path.with_name(f"{cases_path.stem}.summary.json"), "run summary"
    )
    run_data = json.loads(summary_path.read_text("utf-8"))
    try:
        fractions = parse_split(split)
    except ExportError as exc:
        raise typer.BadParameter(str(exc)) from exc

    git = git_revision()
    simulate_git = run_data["simulate"].get("simulator_git") or {}
    dirty = bool(git.get("dirty")) or bool(simulate_git.get("dirty"))
    if dirty and not allow_dirty:
        raise typer.BadParameter(
            "the simulator tree (now or at simulate time) had uncommitted changes; "
            "commit them or pass --allow-dirty"
        )

    inputs = run_data.get("inputs", {})
    profiles_path = _required_file(
        profiles or Path(inputs["profiles"]["path"]), "profiles for validation"
    )
    hp_record = inputs.get("hp_json")
    hpo_path = hpo_json or (Path(hp_record["path"]) if hp_record else None)
    ontology = HpoOntology.from_json(_required_file(hpo_path, "hp.json")) if hpo_path else None
    noise_record = inputs.get("noise_vocabulary")
    noise_ids = (
        {term.hpo_id for term in load_noise_vocabulary(Path(noise_record["path"]))}
        if noise_record and Path(noise_record["path"]).is_file()
        else None
    )
    sim_config = run_config(run_data["config"])
    profile_map = {
        profile.disease_id: profile for profile in iter_model_jsonl(profiles_path, DiseaseProfile)
    }
    model, cardinal_index, _ = _reporting_artifacts(run_data, None, None, ontology)
    validation = validate_cases(
        iter_model_jsonl(cases_path, SyntheticCase),
        profiles=profile_map,
        ontology=ontology,
        config=sim_config,
        noise_vocabulary=noise_ids,
        report_model=model,
        cardinal=cardinal_index,
    )
    if validation["violations"]["count"]:
        typer.echo(format_report(validation))
        typer.echo("Not exported: the cases violate invariants.", err=True)
        raise typer.Exit(code=1)

    try:
        records = sort_records(
            training_record(case, fractions)
            for case in iter_model_jsonl(cases_path, SyntheticCase)
        )
    except ExportError as exc:
        raise typer.BadParameter(str(exc)) from exc
    payload = encode_dataset(records)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)

    counts = dataset_counts(records)
    hpo_version = (hp_record or {}).get("version") or "unknown"
    card_data = {
        "dataset_id": (
            f"ds-sim{run_data['simulate']['simulator_version']}-hpo{hpo_version}"
            f"-n{sim_config.cases_per_disease_per_difficulty}-s{sim_config.seed}"
        ),
        "format": EXPORT_FORMAT,
        "record_fields": list(RECORD_FIELDS),
        "file": {
            "path": str(output),
            "bytes": len(payload),
            "sha256": sha256_file(output),
            "records": len(records),
        },
        "split": {"rule": SPLIT_RULE, "fractions": dict(fractions)},
        "counts": counts,
        "genes": run_data.get("genes"),
        "inputs": {
            "cases": _input_record(cases_path),
            "run_summary": _input_record(summary_path),
            **inputs,
            "heldout_pmids": (run_data.get("profile_sources") or {}).get(
                "held-out reference mask"
            ),
            "profile_sources": run_data.get("profile_sources"),
        },
        "calibration": run_data.get("calibration"),
        "entity_profiles": run_data.get("entity_profiles"),
        "reporting": run_data.get("reporting"),
        "simulator": {
            "simulate": run_data["simulate"],
            "export_git": git,
            "allow_dirty": allow_dirty,
            "dirty": dirty,
        },
        "config_hash": run_data["config_hash"],
        "config": run_data["config"],
        "validate": _validation_digest(validation),
        "exported_at": datetime.now(tz=UTC).isoformat(timespec="seconds"),
        "export_wall_time_seconds": round(time.perf_counter() - started, 2),
    }
    card_path = card or card_path_for(output)
    card_path.parent.mkdir(parents=True, exist_ok=True)
    card_path.write_text(json.dumps(card_data, indent=2) + "\n", encoding="utf-8")
    typer.echo(
        f"Exported {len(records)} case(s) for {counts['genes']} gene(s) to {output} "
        f"({len(payload)} bytes); splits "
        + ", ".join(f"{name} {count}" for name, count in counts["by_split"].items())
    )
    typer.echo(f"Dataset card: {card_path}")


def _validation_digest(report: dict) -> dict:
    """The parts of the validate report a dataset card keeps."""

    keys = (
        "cases",
        "diseases",
        "genes",
        "difficulties",
        "config_matches_cases",
        "phenotypes_per_case",
        "positives",
        "negatives",
        "missingness",
        "noise",
        "sex",
        "age",
        "onset",
        "reporting",
        "violations",
    )
    digest = {key: report.get(key) for key in keys}
    calibration = report.get("calibration")
    if calibration is not None:
        digest["calibration"] = {
            view: {
                "pairs": table["pairs"],
                "calibration_error": table["calibration_error"],
                "mean_abs_error_per_pair": table["mean_abs_error_per_pair"],
            }
            for view, table in calibration.items()
        }
    return digest


@app.command("calibration-keys")
def calibration_keys_command() -> None:
    """Print the keys a --calibration file may set, with types and defaults."""

    typer.echo(format_calibration_keys())


def main() -> None:
    """Run the CLI."""

    app()


if __name__ == "__main__":
    main()
