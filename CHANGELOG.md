# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
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
- HPOA frequencies follow ADR-0007 (`diagnostic.ar-training/docs/decisions/0007-frequency-estimation.md`): pooled counts use the Jeffreys mean `(n+0.5)/(m+1)` with a 95% Jeffreys interval, categories use HPO range midpoints, percentages stay points. `1/1` is no longer "obligate".
- **`0/m` is a low frequency, not a negative phenotype.** Only `NOT` and `Excluded` rows are negatives now.
- The OMIM-ORPHA reader parses Orphanet `en_product1.xml` (2026-06-23): only exact, validated alignments of active disorders transfer onset, and OMIM onset provenance names the alignment.
- `configs/mvp.yaml` builds from the inputs pinned by `diagnostic.ar-training`: HPO 2026-02-16, Orphanet 2026-06-23 ages and alignments, and the `r1-v1` PMID mask.
- Difficulty presets are Pydantic models with `negatives_mean` and `noise_mean` (per-case counts) instead of per-term rates, and can be overridden from the config.
- Cases no longer carry `metadata.generated_at`, so identical runs are byte-identical; `source_versions` now holds input hashes and the simulator SHA. The run summary records the time.

### Removed
- The unused `SimulationConfig` rates (`positive_observation_rate`, `known_negative_rate`, `missingness_rate`, `ontology_smoothing_rate`); the presets hold them.

### Fixed
- `simulate` now passes the ontology, noise vocabulary and gene/disease labels to the simulator. Before, it dropped them, so ontology generalization and noise injection never ran from the CLI.
