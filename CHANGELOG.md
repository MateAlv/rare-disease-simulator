# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.5.1] - 2026-10-01

Candidate 4 (`diagnostic.ar-training` PR #44) reported only 17–19% of true terms in independent mode, drew related noise more specific than real terms, and kept ages far below real ones. Every new knob defaults to 0.5.0 behaviour.

### Added
- `reporting.q_scale`: independent mode uses `q = min(1, q_scale · π / f)`. `validate` and the run summary record the realized `q_mean` next to `q_capped_share`.
- `noise.related_up_levels` and `noise.related_down_levels` (defaults 2 and 2, the 0.5.0 reach), and `noise.related_max_ic`, which skips related-noise candidates above an IC computed with the report model's `information_content` feature and falls back to the vocabulary.
- `--presentation-ages` / `sources.presentation_ages_path` (optional `sources.presentation_ages_sha256` pin): age-at-presentation ranges per disease id. A listed entity, or any of its profile ids, draws its age uniformly in the range, raised to the onset and cut at `max_age_years`. `metadata.presentation_age_years` records the range. The run summary and dataset card record the file's sha256 and the cases that used it, and `validate` checks their ages (`presentation_age_out_of_range`).
- `tests/fixtures/v04/golden_0.5.0_independent.jsonl.gz`: independent-mode cases written by the 0.5.0 code, which the defaults must reproduce.

### Changed
- `SIMULATOR_VERSION` 0.5.1.

## [0.5.0] - 2026-09-30

An independent review of 16 construction-side disease cards found hallmark features rarely reported (10 of 16) and records dominated by unrelated noise (9 of 16), and simulated ages far below real ones. Every new knob defaults to 0.4.2 behaviour.

### Added
- `reporting.mode: independent`: each true profile term is reported on its own with `q = min(1, π/f)` (π the report model's score, f the simulation frequency floored at 1e-3), cardinal terms included (`force_cardinal` is ignored, with a warning). With no successful roll, the highest-q true term is shown as it is (`forced_min_one`). Generalization, specialization and absorbed terms act per term, with no refill. True profile terms carry `CasePhenotype.report_probability`.
- `reporting.noise_count: proportional`: `Poisson(m · ρ/(1 − ρ))` noise terms for `m` reported profile terms, ρ = `reporting.noise_share`, capped at the showable vocabulary. It is for independent mode; `budget_share` stays with `report_model` mode, and the config rejects mismatched pairs.
- `age.duration_mean_by_onset`: a mean exponential duration per onset category, overriding `age.duration_mean_years`.
- `validate` in independent mode: a q calibration table (bins of width 0.1, mean q against the realized report rate), reported terms and noise per case, the share of true terms reported, the capped share and the `forced_min_one` count; the budget checks do not apply.
- Run summary in independent mode: reported profile terms and noise per case, share of true terms reported, `q_capped_share`, `forced_min_one`.

### Changed
- `SIMULATOR_VERSION` 0.5.0.
- The golden-case test compares cases without null fields, as `write_jsonl` writes them, and also checks the new knobs at their defaults.

## [0.4.2] - 2026-09-30

Candidate 2 design of `diagnostic.ar-training` ADR-0011 amendment 3 (evidence: ANALYSIS-002). Both knobs default to 0, and with the defaults a run reproduces 0.4.1 exactly apart from the version and config hash.

### Added
- `reporting.specialize_rate`: a reported, not generalized profile term is shown with this probability as a descendant (a child, then with probability 0.5 a grandchild) that is a phenotypic abnormality, allowed for the patient's sex, not already shown and not annotated to the profile. It fills the same slot, carries `reason: specialized` and `source_hpo_id`, and stays within `true_present`.
- `noise.related_share`: each noise slot is, with this probability, a descendant (at most 2 levels) of a true term's parent or grandparent that the entity neither annotates nor has as an ancestor or descendant of an annotation (descendants count as specialisation, ANALYSIS-002), falling back to the vocabulary when nothing qualifies. Marked `reason: related` with `source_hpo_id` = the true term.
- `validate`: `specialized_not_descendant_of_true_term`, `related_noise_annotated_to_entity`, and the specialised share of positives and related share of noise; related noise is exempt from the vocabulary check.
- Run summary: `reporting.specialized` and a `noise` block (`terms`, `related`, `related_share`).
- `tests/fixtures/v04/golden_0.4.1.jsonl.gz`: cases written by the 0.4.1 code, which the defaults must reproduce.

### Changed
- `SIMULATOR_VERSION` 0.4.2. `CasePhenotype.source_hpo_id` also names the term a specialised positive or related noise came from.

## [0.4.1] - 2026-09-30

### Fixed
- Report mode: a picked true term whose own id an earlier generalization already showed was silently skipped, so its case reported one term fewer than its budget (`reported_count_mismatch`, 622 of 486,800 cases in the first full v0.4 run). Its slot is now refilled from the remaining true terms.

### Added
- `reporting.noise_count: budget_share` with `reporting.noise_share`: each budget slot is noise with that probability, at most `k − 1`, so noise scales with the record instead of being drawn independently of it. The default stays `poisson`.

## [0.4.0] - 2026-09-30

Simulator 0.4.0, the realism rebuild of `diagnostic.ar-training` ADR-0011 (motivated by ANALYSIS-001). Every learned quantity comes from the training repo's artifacts (`hpoa-recount-v1`, `genes-v2`, `report-model-v1`, `cardinal-v1`); nothing is fitted here.

### Added
- Merged entity profiles (`profiles/merge.py`, `entity_profiles: merged`): a gene-first case simulates its entity from one profile merged from all the entity's `profile_ids`. Terms are unioned; counts are pooled and win over percentages, which win over categories (so an Orphanet category never exceeds a disagreeing count-based estimate; capped cases are counted); sex restrictions survive only when all profiles agree; each term keeps its earliest onset; disease onsets are averaged, inheritance is unioned. Cases carry `target.profile_ids`, and `target.disease_id` is the entity. The run summary's `entity_profiles` block counts what merging did.
- Beta shrinkage of count-based frequencies (`frequency.count_estimator: beta_shrinkage`, `frequency.shrinkage_mean`, `frequency.shrinkage_strength` = 2): `(n + s·μ₀)/(m + s)`, with μ₀ the median count-based frequency of `hpoa-recount-v1` supplied by the training repo. There is no built-in μ₀; `simulate` refuses to run without it when a profile term is count-based.
- Truth, then reporting (`reporting.mode: report_model`): the true phenotype is sampled as before; a term budget is drawn from `report-model-v1`'s histogram (its 0 bin dropped), true cardinal terms from `cardinal-v1` are reported first (`reporting.force_cardinal`), and the rest are drawn without replacement in proportion to the reporting model's logistic score. Unreported true terms stay in the case as `missing` (`not_reported`) or `unknown`. Reported cardinal positives carry `reason: cardinal`.
- Noise counts toward the report budget: the histogram's totals include real records' noise, so in `report_model` mode the noise count `n` is drawn right after `k` (capped at the terms the case could show), the profile terms get `max(1, k − n)`, and noise is picked before the negatives, which keep clear of it. Cases record `metadata.report_budget` (`k`) and `metadata.report_budget_profile`. The run summary's `reporting` block adds the realized `budget_mean_drawn`, `reported_total_mean`, `reported_profile_mean` and `reported_noise_mean`. Observation mode is unchanged.
- `simulation/reporting.py`: strict readers for `report-model-v1` (JSON schema `report-model` v1; features `frequency`, `information_content`, `ancestor_count`, `log_reportability`, `cardinal_flag`) and `cardinal-v1` (TSV `disease_id, hpo_id, kind, source`); any mismatch fails at load time.
- `simulate --report-model/--cardinal`, with `sources.report_model_path`, `sources.cardinal_path` and optional `sources.report_model_sha256`/`sources.cardinal_sha256` pins. The run summary and dataset card record both files' sha256; the run summary also records the model id, the effective budget distribution and the cardinal counts (`reporting`).
- `validate` for report-model cases: `budget_outside_histogram`, `profile_budget_mismatch`, `reported_count_mismatch` (reported profile terms ≠ min(true, max(profile budget, true cardinal))), `cardinal_not_reported`, and a `reporting` section comparing drawn budgets and reported totals (profile terms plus noise) with the histogram. `validate --report-model/--cardinal` override the run summary's files.
- Build summary: `rows.phenotype_rows_used(_with_onset)` and `rows.recounted(_with_onset)` for `hpoa-recount-v1` builds (rows citing `R1-TRAIN:`).
- `ConfounderIndex.candidate_terms_for_profile`: confounders of a merged entity, never its own profiles.

### Changed
- `SIMULATOR_VERSION` 0.4.0; the defaults are the v0.4 model. `reporting.mode: observation`, `frequency.count_estimator: jeffreys` and `entity_profiles: uniform` restore simulator 0.3.
- `simulate` needs a report model (and a cardinal file while `force_cardinal` is on) in the default `report_model` mode, and stops before generating when an input is missing.
- `validate` checks merged cases against the profile merged from their `profile_ids`, and calibration expects the frequency the run simulated with (the count estimator), not the stored Jeffreys mean.
- Export v2: with merged profiles `profile_id` is the entity. The record format is unchanged.

## [0.3.0] - 2026-09-29

Simulator 0.3.x as merged on `main` (not tagged).

### Added
- Gene-first simulation (simulator 0.3.0): `simulate --genes all|LIST|FILE --gene-profiles genes.json.gz --gnn-genes genes.json --cases-per-gene N`.
  - One case set per GNN gene, in class-index order.
  - Each case draws an entity by `sim_weight`, then one of its profiles uniformly.
  - The sex prior comes from the entity's inheritance.
  - Cases are labelled with the v3 symbol and class index.
  - Genes that cannot be simulated are skipped and counted by reason in the run summary.
  - The configured `genes-v1` path is pinned by sha256 (`sources.gene_profiles_sha256`).
- `build-profiles --gene-profiles`: builds the diseases a simulable `genes-v1` entity lists even without an HPO gene link (317 on HPO 2026-02-16).
- `build-profiles --drop-annotations holdout.tsv`: an annotation holdout that removes listed (disease, term) positive annotations before profiles are built. The summary records its sha256, the pairs and rows dropped, and the diseases left without phenotypes.
- Negative source `own_gene_other_disease` (gene-first only, default odds 0.1): terms of the gene's other entities.
- `negatives.own_disease_pool` (`entity` by default): own-disease negatives also come from the entity's equivalent profiles.
- `negatives.unfilled_slots` (`drop` by default, or `redistribute`).
- `simulate --calibration calib.json`: overrides knobs by dotted config path.
  - Unknown keys are errors.
  - `noise.vocabulary_path` is accepted, and `_` keys are notes.
  - The file's sha256 is recorded.
  - `calibration-keys` prints the accepted keys, generated from the schema.
- `export-training`: export v2, a deterministic gzip JSONL with `present`, `true_present`, `excluded`, sex, age, onset, difficulty, seed and a seed-hashed sim train/val/test split, plus a dataset card. It validates the cases first and refuses a dirty tree without `--allow-dirty`.
- `validate --dataset`: checks an export against its card, including `present` ⊆ `true_present` and no excluded term being a true term or an ancestor of one.
- Case fields `target.entity_id`, `metadata.case_seed` and `metadata.sex_prior_key`.
- `validate` reports the gene count, checks the per-prior sex share against each case's recorded prior, flags gene-other-disease negatives annotated to the profile, and defaults the noise vocabulary to the run summary's.
- `build-profiles --from-hpoa`. It builds one `DiseaseProfile` per gene-linked disease from HPO `phenotype.hpoa`, `genes_to_disease.txt` and `hp.json`, with an Orphanet age-of-onset fallback. It also writes a build summary JSON with counts, input sha256s and the simulator git SHA.
- `--exclude-pmids`. It drops every HPOA row that cites a held-out publication before profiles are built, and reports the rows dropped, the diseases affected and the diseases left without phenotypes.
- `HpoOntology.from_json` for the official `hp.json`. It covers labels, `is_a` parents, obsolete terms, alt-id / `replaced_by` resolution, and cached ancestor and descendant sets.
- New optional profile fields, all backwards compatible:
  - `PhenotypeAssociation.sex_restriction`, `onset_hpo_id` and `frequency_raw`;
  - `DiseaseGene.ncbi_gene_id` and `inheritance_hpo_ids`;
  - `AgeOfOnset.hpo_ids`;
  - `SourceReference.sha256`;
  - `ProfileQuality.counters`.
- `simulate --hpo-json/--noise-vocabulary/--labels`.
- `PhenotypeAssociation.frequency_estimate`: the ADR-0007 point frequency (`profiles/frequency.py`).
- Build summary: `phenotypes.by_frequency_basis`, `omim_orpha_map` counts by relation and status, and onset coverage with and without the OMIM-ORPHA mapping.
- Simulator v0.2 (`SIMULATOR_VERSION = "0.2.0"`), all knobs in `SimulationConfig`:
  - per-patient phenotype probability ~ Beta around the frequency estimate, with an explicit prior for unknown frequencies;
  - sex from the inheritance-derived prior; sex-restricted terms and the other sex's genital/reproductive subtrees are never emitted;
  - onset from the profile's onset distribution, current age ≥ onset, phenotype onset gating (antenatal/congenital always eligible) and a duration boost for progressive diseases;
  - 0–15 asked-and-absent negatives per case from the disease's own absent terms (default 70%), confounder diseases' hallmark terms from an IC-weighted inverted index (20%) and NOT annotations (10%), never related to a present term;
  - at least one observed positive per case (deterministic redraws, then a forced term);
  - covariate missingness for sex, age, onset and negatives;
  - noise only from an external, optionally weighted vocabulary;
  - `CasePhenotype.source_hpo_id` for generalized positives.
- `simulate --sample-diseases/--disease-ids/--cases-per-disease/--difficulty/--seed/--summary`, and a run summary JSON with the full config, input sha256s and the simulator git SHA.
- `validate --cases [--profiles --hpo-json --noise-vocabulary --run-summary --report]`: counts, negative sources, missingness, noise, sex and onset against the priors, binned frequency calibration and invariant checks; exit status 1 on any violation.

### Changed
- Simulator 0.3.0. Disease-first cases are unchanged apart from the new metadata fields and the fix below.
- The sex bias implied by inheritance is shared by the profile builder and the simulator (`profiles/inheritance.py`).
- HPOA frequencies follow ADR-0007 (`diagnostic.ar-training/docs/decisions/0007-frequency-estimation.md`): pooled counts use the Jeffreys mean `(n+0.5)/(m+1)` with a 95% Jeffreys interval, categories use HPO range midpoints, percentages stay points. `1/1` is no longer "obligate".
- **`0/m` is a low frequency, not a negative phenotype.** Only `NOT` and `Excluded` rows are negatives now.
- The OMIM-ORPHA reader parses Orphanet `en_product1.xml` (2026-06-23): only exact, validated alignments of active disorders transfer onset, and OMIM onset provenance names the alignment.
- `configs/mvp.yaml` builds from the inputs pinned by `diagnostic.ar-training`: HPO 2026-02-16, Orphanet 2026-06-23 ages and alignments, and the `r1-v1` PMID mask.
- Difficulty presets are Pydantic models with `negatives_mean` and `noise_mean` (per-case counts) instead of per-term rates, and can be overridden from the config.
- Cases no longer carry `metadata.generated_at`, so identical runs are byte-identical; `source_versions` now holds input hashes and the simulator SHA. The run summary records the time.

### Removed
- The unused `SimulationConfig` rates (`positive_observation_rate`, `known_negative_rate`, `missingness_rate`, `ontology_smoothing_rate`); the presets hold them.

### Fixed
- Confounder similarities were summed in string-hash order, so near-tied confounders could rank differently from one process to the next and "identical" runs were not byte-identical (3 of 50,500 gene-first cases differed). Term vectors are now built in sorted order.
- When redraws failed and no profile term suited the patient's sex, the forced positive could be a sex-inappropriate term (seen once in 50,500 gene-first cases: a male with ovarian carcinoma). The patient now takes the only sex that can present the disease.
- A recorded term whose generalized parent was already in the case was dropped from the case. It is now kept as `missing` with `reason: recorded_as_generalized`, so the case (and `validate`'s calibration) holds every term the patient has.
- `simulate` now passes the ontology, noise vocabulary and gene/disease labels to the simulator. Before, it dropped them, so ontology generalization and noise injection never ran from the CLI.
