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

### Fixed
- `simulate` now passes the ontology, noise vocabulary and gene/disease labels to the simulator. Before, it dropped them, so ontology generalization and noise injection never ran from the CLI.
