# Implementation Plan

This plan describes the first implementation path for `rare-disease-simulator`.
It is intentionally staged so each milestone produces a usable artifact and a clear validation surface.

## Guiding Principle

```text
LLM = evidence-grounded extraction and normalization
Simulator = seeded, auditable probabilistic patient generation
```

The LLM may propose structured claims, but only validation code can merge them into a `DiseaseProfile`.
The LLM must never generate synthetic patients directly.

## MVP Scope

Initial diseases:

- Niemann-Pick type C
- Usher syndrome type I
- Rett syndrome
- Marfan syndrome
- Duchenne muscular dystrophy
- Phenylketonuria
- Cystic fibrosis
- Tuberous sclerosis complex

Initial sources:

- HPO downloads
- Orphadata API/downloads
- MONDO mappings
- GeneReviews / MedGen / Orphanet text for controlled LLM extraction

Initial outputs:

- `profiles.jsonl`
- `profile_patches.jsonl`
- `rich_cases.jsonl`
- `graphens.json`

Initial training compatibility:

- Rich cases always store both `disease_id` and `gene`.
- First GraPhens export is grouped by gene for compatibility.
- Disease-label export comes later.

## Locked Technical Decisions

These choices are fixed for the first implementation pass unless they prove impractical:

- Python target: `3.11`.
- CLI framework: Typer.
- Schemas: Pydantic v2.
- Config: YAML config files plus CLI overrides.
- Canonical disease ID: `ORPHA:123`.
- HPO MVP inputs: local files in `data/raw`; support `hp.json` or `hp.obo` for ontology plus `phenotype.hpoa` and gene-to-phenotype files where needed.
- Source fetching: manual files in `data/raw` first; automated fetchers later.
- Text snippet storage: JSONL.
- LLM provider: abstract provider interface with dummy provider for tests and a real provider adapter later.
- LLM HPO mapping: LLM may propose HPO candidates; local validation accepts or rejects them. Optional local lexical search can provide candidate context to the LLM.
- LLM confidence policy: auto-merge `>= 0.75`, review `0.50-0.75`, reject `< 0.50`.
- GraPhens coupling: export files only in the MVP. This repo should not call GraPhens directly.

## First Vertical Slice

Before scaling to all sources and all MVP diseases, the first target is a minimal end-to-end run:

```text
one disease
one text snippet
one DiseaseProfilePatch
one validated DiseaseProfile
ten SyntheticCase records
one graphens.json export
```

This vertical slice should work entirely from fixtures and dummy providers before using real downloaded sources or real LLM calls.

## Milestone 0: Repository Foundation

Goal: make the repository installable, testable, and ready for real modules.

Tasks:

- Add `environment.yml` using conda-managed dependencies.
- Add `pyproject.toml` for package metadata, CLI entrypoint, formatting, and test config.
- Fill `configs/mvp.yaml` with real config keys for diseases, source paths, LLM settings, simulation settings, and output paths.
- Add a minimal CLI skeleton with subcommands:
  - `fetch-sources`
  - `extract-profile-patches`
  - `build-profiles`
  - `simulate`
  - `export-graphens`
  - `validate`
- Add structured logging conventions.
- Add basic unit test smoke checks.
- Decide whether CLI config loading uses only YAML or YAML plus CLI overrides.

Deliverables:

- Installable local package.
- `rare-disease-simulator --help` works.
- `pytest` can run.
- Config loading works for `configs/mvp.yaml`.

## Milestone 0.5: Sample Fixtures

Goal: develop the end-to-end pipeline without depending on large downloads or external services.

Tasks:

- Add small fixture files:
  - `tests/fixtures/hpo_terms.tsv`
  - `tests/fixtures/orphadata_mini.json`
  - `tests/fixtures/text_snippets.jsonl`
  - `tests/fixtures/profile_patch.json`
  - `tests/fixtures/disease_profile.json`
- Add fixture readers used only by tests.
- Add a dummy LLM provider that returns a fixed `DiseaseProfilePatch`.
- Add a minimal fixture disease that can exercise the full pipeline.

Deliverables:

- The first vertical slice can run from fixtures only.
- Tests can exercise schemas, validation, merge, simulation, and export without network access.

## Milestone 1: Core Schemas

Goal: define stable typed data contracts before source ingestion.

Modules:

- `profiles/schema.py`
- `llm_extraction/schema.py`
- `simulation/schema.py`

Schemas:

- `Provenance`
- `SourceReference`
- `MappedDiseaseIds`
- `DiseaseGene`
- `PhenotypeAssociation`
- `NegativePhenotypeAssociation`
- `AgeOfOnset`
- `SexBias`
- `DiseaseProfile`
- `DiseaseProfilePatch`
- `SyntheticCase`
- `SimulationConfig`
- `GeneratorMetadata`

Rules:

- Every source-derived field must support provenance.
- Every LLM-derived claim must support evidence span and confidence.
- `negative`, `missing`, and `unknown` must be represented distinctly.
- Cases must store both gene and disease targets.

Deliverables:

- JSON schema export for core models.
- Unit tests for serialization/deserialization.
- Example profile and example case fixtures.

Open decisions:

- Whether enums should be strict or allow `unknown`.
- How to represent probability ranges and calibrated probabilities.

## Milestone 2: HPO Loader And Validation

Goal: load the HPO ontology and validate HPO IDs used by profiles and cases.

Modules:

- `data_sources/hpo.py`
- `validation/hpo.py`

Tasks:

- Download or read configured HPO files.
- Parse ontology terms and labels.
- Parse disease-phenotype annotations when available.
- Parse gene-phenotype links when available.
- Parse negative phenotype annotations.
- Expose:
  - `is_valid_hpo_id`
  - `get_label`
  - `get_ancestors`
  - `get_direct_parents`
  - `is_phenotypic_abnormality`
- Record HPO version in artifacts.

Deliverables:

- HPO cache object.
- CLI/source command can validate configured HPO file paths.
- Unit tests with small fixture files.

Open decisions:

- Whether ontology traversal is needed in this repo or delegated to GraPhens later.

## Milestone 3A: Minimal Structured Profile

Goal: build preliminary disease profiles from the smallest useful structured source subset.

Modules:

- `data_sources/orphadata.py`
- `profiles/builder.py`

Tasks:

- Ingest local Orphadata fixture/raw records for MVP diseases.
- Normalize disease identifiers.
- Extract:
  - ORPHA IDs;
  - disease names;
  - genes;
  - HPO associations;
  - frequency categories;
  - basic provenance.
- Preserve source versions and retrieval dates.

Deliverables:

- `build-profiles --structured-only` produces draft `profiles.jsonl`.
- Profile quality warnings for low coverage, missing genes, and missing HPOs.
- Unit tests for parser fixtures.

Open decisions:

- Exact local Orphadata file format for the MVP.

## Milestone 3B: Extended Structured Profile

Goal: enrich structured profiles after the minimal profile builder is stable.

Modules:

- `data_sources/orphadata.py`
- `data_sources/mondo.py`
- `profiles/builder.py`
- `profiles/quality.py`

Tasks:

- Extract additional structured fields when available:
  - inheritance;
  - age of onset;
  - prevalence;
  - diagnostic criteria;
  - aliases.
- Ingest MONDO mappings for ORPHA/OMIM/MONDO/MedGen.
- Merge HPO and Orphadata phenotype associations.
- Emit conflicts as warnings or review items.
- Preserve extended source versions and retrieval dates.

Deliverables:

- Extended `profiles.jsonl` includes mappings and extra clinical fields when available.
- Conflict and low-coverage reports.

Open decisions:

- How to choose canonical disease profile when Orphadata/MONDO mappings are ambiguous.

## Milestone 4: Text Source Collection

Goal: prepare trusted text snippets for LLM extraction without free-form scraping.

Modules:

- `data_sources/` text helpers if needed.
- `docs/sources.md` source-specific instructions.

Tasks:

- Define accepted text source types:
  - GeneReviews clinical characteristics.
  - GeneReviews diagnosis / suggestive findings.
  - Orphanet clinical description.
  - Orphadata textual fields.
  - MedGen disease summary.
- Define a local input format for text snippets:
  - disease ID;
  - source name;
  - source URL;
  - retrieved date;
  - section name;
  - text content;
  - license/redistribution notes.
- Add examples for 1-2 MVP diseases.
- Ensure restricted source text is not committed unless allowed.

Deliverables:

- Text snippet schema.
- Example snippet fixtures with redistributable or synthetic text.
- Source policy documented in `docs/sources.md`.

Open decisions:

- Whether we manually curate first snippets or build fetchers for public pages.

## Milestone 5: LLM Extraction MVP

Goal: produce validated `DiseaseProfilePatch` objects from trusted text.

Modules:

- `llm_extraction/schema.py`
- `llm_extraction/extractor.py`
- `llm_extraction/validator.py`
- `llm_extraction/prompts/`

Tasks:

- Define strict extraction prompt.
- Define JSON schema for LLM output.
- Implement LLM provider abstraction.
- Implement deterministic extraction settings:
  - low temperature;
  - structured output where available;
  - prompt version;
  - model/provider metadata;
  - extraction schema version;
  - retries with validation errors.
- Add mandatory extraction cache:
  - cache key = hash of `source_text + prompt_version + model + extraction_schema_version`;
  - cache stores raw response, parsed patch, validation result, and metadata;
  - cache is used by default to avoid repeated paid/provider calls.
- Validate:
  - HPO IDs exist;
  - evidence spans are present;
  - confidence is present;
  - unmentioned phenotypes are not treated as negative;
  - no synthetic patients appear in output.
- Write accepted patches to `profile_patches.jsonl`.
- Write rejected/needs-review patches to a review queue.

Review queue outputs:

- `outputs/review/low_confidence_claims.jsonl`
- `outputs/review/conflicts.jsonl`
- `outputs/review/rejected_patches.jsonl`

Deliverables:

- `extract-profile-patches --disease ORPHA:...` works on local snippet input.
- Patch output includes provenance, prompt version, model, confidence, evidence.
- Validation reports for accepted/rejected claims.

Open decisions:

- First LLM provider.

## Milestone 6: Profile Merge And Quality

Goal: merge structured profiles and validated LLM patches into final `DiseaseProfile` artifacts.

Modules:

- `llm_extraction/merge_patch.py`
- `profiles/builder.py`
- `profiles/quality.py`

Tasks:

- Define merge precedence:
  - structured source with clear provenance;
  - LLM extraction with explicit evidence and high confidence;
  - LLM extraction with lower confidence becomes warning/review item.
- Merge phenotype duplicates by HPO ID.
- Merge frequency categories and probability ranges.
- Merge onset, progression, diagnostic role, sex bias, inheritance.
- Emit conflicts as warnings.
- Compute profile-level confidence and source coverage.
- Record all source versions and prompt/model metadata.

Deliverables:

- `build-profiles` produces validated `profiles.jsonl`.
- Quality report summarizing:
  - HPO coverage;
  - phenotype counts;
  - negative phenotype counts;
  - missing source coverage;
  - conflicts;
  - low-confidence claims.

Open decisions:

- How to numerically combine confidence across sources.
- Whether low-confidence fields are retained in profile or stored only as review metadata.

## Milestone 7: Simulation Engine

Goal: generate realistic, seeded synthetic cases from disease profiles.

Modules:

- `simulation/schema.py`
- `simulation/difficulty.py`
- `simulation/simulator.py`

Tasks:

- Implement frequency category to probability range mapping.
- Sample:
  - disease;
  - gene;
  - sex;
  - age;
  - age of onset;
  - true positive phenotypes;
  - observed positive phenotypes after missingness;
  - explicit negatives;
  - missing phenotypes;
  - ontology-smoothed/generalized phenotypes.
- Implement difficulty presets:
  - easy: more cardinal features, high completeness;
  - medium: moderate missingness;
  - hard: fewer positives, more generalized terms, more missingness.
- Preserve seed and config hash per case.
- Prevent impossible states:
  - same HPO cannot be both positive and negative in one case;
  - negative requires explicit source;
  - HPO IDs must validate.

Deliverables:

- `simulate --profiles ...` writes `rich_cases.jsonl`.
- Cases include `disease_id`, `gene`, `gene_label`, `disease_label`.
- Deterministic reruns with same seed/config produce same cases.

Open decisions:

- Exact missingness rates by difficulty.
- Whether ontology smoothing happens before or after missingness.

Simulation v2 candidates:

- Nonspecific clinical noise.
- Differential-diagnosis confounders.
- Asked-negative logic.
- Sex-aware phenotype modulation.

## Milestone 8: Exporters

Goal: make generated cases usable by GraPhens and future tooling.

Modules:

- `exports/jsonl.py`
- `exports/graphens.py`
- `exports/phenopacket.py`

Tasks:

- Rich JSONL writer/reader.
- GraPhens export grouped by gene:

```json
{
  "NPC1": [["HP:...", "HP:..."]]
}
```

- Optional disease-label export grouped by disease.
- Keep negatives and metadata out of GraPhens export for baseline compatibility.
- Preserve a mapping file from export rows back to rich cases.
- Stub Phenopacket export with explicit limitations.

Deliverables:

- `export-graphens --cases outputs/rich_cases.jsonl` writes `outputs/graphens.json`.
- Export report summarizes sample counts per gene/disease.
- GraPhens integration remains file-based only in the MVP.

Open decisions:

- Whether GraPhens export should include noise phenotypes in the positive list.
- How to handle cases with multiple causal genes.

## Milestone 9: Generator Evaluation

Goal: evaluate the generator independently from downstream diagnostic models.

Modules:

- `validation/sanity.py`

Tasks:

- Validate:
  - invalid HPO rate;
  - phenotype count distribution;
  - observed vs expected phenotype frequency;
  - missingness rate;
  - noise rate;
  - negative phenotype rate;
  - disease/gene coverage;
  - difficulty balance;
  - source provenance coverage;
  - warnings and low-confidence extraction rate.
- Produce JSON and human-readable reports.
- Add seed/config reproducibility checks.

Deliverables:

- `validate` command for profiles and cases.
- Evaluation report for MVP generation run.

Open decisions:

- Which thresholds fail the validation command vs only warn.
- How to compare generated cases against Phenopacket Store examples without leakage.

## Milestone 10: Integration With GraPhens

Goal: document and verify file-based integration with the existing GraPhens dataset/training pipeline.

Tasks:

- Export a GraPhens-compatible JSON for MVP cases.
- Run GraPhens JSON -> NPZ conversion externally.
- Verify labels and gene mappings.
- Train or dry-run with a small subset.
- Compare against the current GraPhens simulator output qualitatively:
  - phenotype count;
  - specificity;
  - missingness;
  - noise;
  - disease diversity.

Deliverables:

- Small exported dataset that GraPhens can ingest.
- Integration notes in `docs/design.md` or a dedicated integration doc.

Open decisions:

- Whether to support disease-label exports before or after first GraPhens run.

## Cross-Cutting Requirements

Reproducibility:

- Record HPO version.
- Record Orphadata release/retrieval date.
- Record MONDO version.
- Record prompt version.
- Record LLM provider/model.
- Record simulator version.
- Record config hash.
- Record random seed.
- Record generation timestamp.

Privacy:

- MVP uses public biomedical sources only.
- MVP does not ingest private patient records.
- Generated synthetic cases are not intended to represent real individuals.

Licensing:

- Preserve source attribution.
- Do not include restricted source text in generated artifacts unless redistribution is allowed.
- Avoid automated OMIM ingestion until license constraints are resolved.

Clinical semantics:

- `negative` means explicitly absent.
- `missing` means not observed, not asked, or not recorded.
- `unknown` means not enough information.
- Unmentioned phenotypes are not negative phenotypes.

## Immediate Next Tasks

Suggested next implementation order:

1. Add `environment.yml`, `pyproject.toml`, Typer CLI wiring, and test config.
2. Fill `configs/mvp.yaml` with concrete keys.
3. Implement minimal Pydantic schemas.
4. Add JSONL readers/writers.
5. Add dummy fixtures.
6. Implement minimal HPO validator with fixture data.
7. Implement dummy LLM provider returning a fixed patch.
8. Merge patch into profile.
9. Simulate 10 cases from the profile.
10. Export `graphens.json`.
