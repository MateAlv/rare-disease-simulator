# rare-disease-simulator

Synthetic clinical phenotype generation for rare disease diagnosis.

## Status

Early MVP development. The first milestone is an end-to-end pipeline for a small set of diseases: source ingestion, LLM-assisted profile extraction, validated `DiseaseProfile` construction, synthetic case simulation, and GraPhens-compatible export.

## Project Goal

`rare-disease-simulator` will build enriched rare-disease profiles from curated biomedical sources and use them to generate synthetic, auditable clinical cases for phenotype-driven diagnosis.

The project is designed as a higher-fidelity replacement for simple gene-to-HPO synthetic case generation. Instead of sampling clean phenotype lists directly from gene annotations, it will model disease-level clinical profiles with phenotype frequencies, onset, inheritance, sex-related effects, negative phenotypes, missingness, noise, and case difficulty.

The core principle is:

```text
LLM = evidence-grounded extraction and normalization
Simulator = seeded, auditable probabilistic patient generation
```

LLMs are part of the MVP, but they are not allowed to generate patients directly.

The LLM may propose structured claims, but only validation code can merge them into a `DiseaseProfile`.

## Initial Pipeline

```text
HPO / Orphadata / MONDO
        +
trusted public medical text
        |
        v
LLM extractor
        |
        v
DiseaseProfilePatch with evidence
        |
        v
validated DiseaseProfile
        |
        v
probabilistic simulator
        |
        v
synthetic cases
        |
        +--> rich JSONL
        +--> GraPhens-compatible JSON
        +--> future Phenopacket export
```

## MVP Diseases

The MVP will start with 8 varied, well-documented diseases whose target genes are present in the local GraPhens MCRD dataset. The selection is meant to stress different clinical patterns while keeping the first validation loop compatible with MCRD.

| Disease | Gene | Why it is useful for the MVP |
| --- | --- | --- |
| Niemann-Pick disease type C1 | NPC1 | Neurovisceral, progressive, variable onset, strong clinical descriptions. |
| Rett syndrome / atypical Rett syndrome | MECP2 | Neurodevelopmental, sex-related effects, staged/evolving presentation. |
| Duchenne muscular dystrophy | DMD | Neuromuscular, X-linked, clear pediatric progression. |
| Cystic fibrosis | CFTR | Multisystem respiratory/digestive disease with strong documentation. |
| Neuronal ceroid lipofuscinosis 6 | CLN6 | Progressive neurodegenerative lysosomal/storage phenotype. |
| Pompe disease / glycogen storage disease II | GAA | Treatable metabolic-neuromuscular disease with onset and severity variation. |
| Hypophosphatasia | ALPL | Skeletal/metabolic disease with onset variation. |
| Noonan syndrome | PTPN11 | Developmental, craniofacial, cardiac, and growth phenotype. |

Additional candidates for expansion:

- Ataxia-telangiectasia / `ATM`
- PTEN hamartoma tumor syndrome / `PTEN`
- Fibrodysplasia ossificans progressiva / `ACVR1`
- Mucopolysaccharidosis type VI / `ARSB`

## Source Policy

Sources are split into structured sources and textual sources.

### Structured Sources

These are the initial structured reference sources for automatic profile construction.

| Source | Intended use |
| --- | --- |
| [HPO downloads](https://human-phenotype-ontology.github.io/downloads.html) | Ontology, disease-phenotype annotations, gene-phenotype links, negative phenotype annotations. |
| [Orphadata API/downloads](https://www.orphadata.com/_orphadata-api/) | ORPHA IDs, HPO associations, phenotype frequencies, age of onset, inheritance, prevalence, diagnostic criteria when available. |
| [MONDO](https://mondo.monarchinitiative.org/) | Disease ID harmonization across ORPHA, OMIM, MONDO, MedGen, and related sources. |

Structured sources enter the profile builder directly, after schema and ontology validation.

### Textual Sources For LLM Extraction

Trusted text sources are used to enrich disease profiles with information that is often incomplete or absent in tabular data.

Accepted from the beginning:

- GeneReviews clinical characteristics
- GeneReviews diagnosis / suggestive findings
- Orphanet clinical descriptions
- Orphadata textual fields when available
- MedGen disease summaries

Accepted later:

- Open-access review papers
- Open-access clinical guidelines
- Open-access case reports
- Phenopacket Store examples for calibration and evaluation

Avoid for the MVP:

- Wikipedia
- Medical blogs
- Non-curated patient websites
- Free-form web scraping
- OMIM as an automated source until licensing and redistribution constraints are resolved

Useful references:

- [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/) is a clinically oriented expert-authored resource for inherited conditions.
- [MedGen](https://www.ncbi.nlm.nih.gov/medgen/docs/overview/) integrates medical genetics concepts and descriptions from multiple authoritative sources.
- [GA4GH Phenopackets](https://www.ga4gh.org/product/phenopackets/) defines a standard for computable clinical and phenotypic case representation.
- [Phenopacket Store](https://monarch-initiative.github.io/phenopacket-store/collections/) can be used later to inspect and calibrate realistic case structure.

## HPOA Structured Backbone

`build-profiles --from-hpoa` builds one `DiseaseProfile` per OMIM/ORPHA/DECIPHER disease in `phenotype.hpoa` that has at least one gene in `genes_to_disease.txt`. It uses no LLM. Code: `profiles/hpoa_builder.py`.

### Inputs

| Option | File | Required |
| --- | --- | --- |
| `--hpo-json` | HPO `hp.json` (same release as the annotations) | yes |
| `--hpoa` | HPO `phenotype.hpoa` | yes |
| `--genes-to-disease` | HPO `genes_to_disease.txt` | yes |
| `--orphanet-ages` | Orphanet `en_product9_ages.xml` (age-of-onset fallback) | no |
| `--omim-orpha-map` | Orphanet `en_product1.xml` (lets OMIM diseases use the Orphanet fallback through exact, validated alignments) | no |
| `--exclude-pmids` | held-out references, one `PMID:n` (or bare `n`) per line | no, but see below |

Any option not given falls back to `sources.*` in the config. `configs/mvp.yaml` points every input at the copies pinned by `diagnostic.ar-training` (its `data/sources.yaml`), checked out next to this repo:

| Input | Release | sha256 |
| --- | --- | --- |
| `hp.json`, `phenotype.hpoa`, `genes_to_disease.txt` | HPO 2026-02-16 | recorded in the build summary |
| `orphanet_ages_2026_06/en_product9_ages.xml` | Orphanet 2026-06-23 | `38fef7e9e398addb147e57df7bb523e5ef072e03c0208ae81b37a2fd200af1ab` |
| `orphanet_alignments/en_product1.xml` | Orphanet 2026-06-23 | `df8d562a0c6011af36a74eb4000ce81ca7d723e8031010819fb71727c0962bbb` |
| `data/splits/r1-v1.heldout-pmids.txt` | split `r1-v1` | recorded in the build summary |

The 2026-06 ages file supersedes the 2025-12-09 one for new builds. With the default config the build is simply:

```bash
rare-disease-simulator build-profiles --from-hpoa --output outputs/profiles.jsonl
```

The command writes `profiles.jsonl` and `profiles.summary.json` (`--summary` overrides the path). The summary holds:

- counts of diseases, genes, phenotypes (by frequency category and by frequency basis: counts, percent or category), negatives, onset sources, inheritance terms and sex bias;
- onset coverage with and without the OMIM-ORPHA mapping (`age_of_onset.gained_via_omim_mapping`, `fraction_with_onset_without_mapping`);
- the product1 mapping counts (`omim_orpha_map`);
- masking effects;
- each input's path, version and sha256;
- the simulator's git SHA;
- the output's sha256.

Identical inputs produce a byte-identical `profiles.jsonl`.

### Masking rule

With `--exclude-pmids`, the builder drops every HPOA row whose `reference` column (`;`-separated) cites a listed reference. This happens before any other processing and covers all aspects. These publications describe evaluation patients, so their annotations must not reach the simulator.

The summary reports:

- rows dropped;
- diseases affected;
- diseases left with no positive phenotype, which are excluded from the output.

Building without a mask is allowed, but it prints a warning.

### Mapping

- **Terms.** Alternative and obsolete HPO IDs are resolved through `hp.json` (`hasAlternativeId`, `replaced_by`). Rows that still cannot be resolved are dropped and counted.
- **Genes.** Symbol, NCBI id and association type (`MENDELIAN` → `causal`, `POLYGENIC` → `susceptibility`, otherwise `unknown`). Rows with symbol `-` are skipped.
- **Positive phenotypes (aspect `P`).** Frequencies follow the convention below (ADR-0007). An empty frequency stays `unknown`, with no estimate and no range; it never means "always". `frequency_raw` keeps the source value (pooled `n/m`, the percentages, or the category terms).
- **Duplicate rows of the same (disease, term).** Merged per the frequency convention; references are unioned.
- **Negative phenotypes.** Only a `P` row with the `NOT` qualifier or the `Excluded` frequency (HP:0040285) becomes a negative phenotype. A `0/m` count is a low frequency, not a negation. If a term has both positive and negative evidence, it stays positive; the conflict is counted in `quality.counters`.
- **Phenotype onset.** The per-row `onset` becomes the phenotype's `onset` and `onset_hpo_id`. When rows disagree, the earliest onset wins.
- **Sex restriction.** The per-row `sex` becomes `sex_restriction` (`male`/`female`), but only when every row for that term names the same sex.
- **Disease age of onset.** Aspect `C` terms under Onset (HP:0003674) become `age_of_onset`:
  - `distribution` holds the share of the disease's onset terms in each category;
  - `category` is the most supported category, with ties going to the earliest;
  - `hpo_ids` lists the terms;
  - Congenital onset maps to `neonatal`.
- **Orphanet fallback.** Without an HPOA onset, `AverageAgeOfOnset` is used (Infancy → infantile, Adolescent → juvenile, Elderly → adult, All ages → variable). ORPHA diseases match directly. OMIM diseases match only through the mapping below; their onset provenance then names both the ages file and the alignment.
- **Progression.** Aspect `C` pace-of-progression terms set `progression`.
- **Inheritance (aspect `I`).** Modes of inheritance go on every gene of the disease (`inheritance` labels, `inheritance_hpo_ids`).
- **Sex bias.** It is derived only from inheritance, and only where the genetics imply it:
  - male-limited or female-limited expression → that sex;
  - all Mendelian modes X-linked recessive or Y-linked → `male`;
  - all Mendelian modes autosomal → `none`;
  - anything else (X-linked dominant or unspecified, mitochondrial, mixed) is left unset.
- **Other aspects.** Aspects `H` and `M` are not used yet; they are counted in the summary.

### Frequency convention (ADR-0007)

Both this repo and `diagnostic.ar-training` read HPOA frequencies the same way, as decided in `diagnostic.ar-training/docs/decisions/0007-frequency-estimation.md`. Code: `profiles/frequency.py`. Each positive phenotype gets `frequency_estimate` (a point), `probability_range` (a range that always brackets the point) and `frequency` (the HPO category whose range holds the point).

| Notation | Estimate | `probability_range` |
| --- | --- | --- |
| counts `n/m` (pooled over duplicate rows: sum of n, sum of m) | Jeffreys posterior mean `(n + 0.5) / (m + 1)`: `1/1 → 0.75`, `10/10 → 0.95`, `0/3 → 0.125` | 95% equal-tailed Jeffreys interval, Beta(n + ½, m − n + ½) quantiles; the lower bound is 0 when n = 0 and the upper bound is 1 when n = m (Brown, Cai & DasGupta 2001) |
| percentages | the mean of the reported percentages (the cohort size is unknown) | min to max of the reported percentages |
| HPO categories | the mean of the range midpoints: Obligate 1.0, Very frequent 0.895, Frequent 0.545, Occasional 0.17, Very rare 0.025 | the envelope of the categories' HPO ranges |

- **Precedence.** When one (disease, term) has several notations, counts win, then percentages, then categories.
- **Empty.** No frequency means unknown (`frequency_estimate` is null); the simulator uses an explicit prior for it.
- **Negation.** Only `NOT` or `Excluded` negate a term. `0/m` means "not seen in m patients", a low frequency.
- **Category of a point.** ≥ 0.99 obligate, ≥ 0.80 very frequent, ≥ 0.30 frequent, ≥ 0.05 occasional, otherwise very rare. A positive annotation is never `excluded`.
- **Counts above the denominator** (`5/3`) are capped at the denominator; unparseable values are counted in the summary (`rows.unparsed_frequency`) and ignored.

### OMIM ↔ ORPHA mapping

`--omim-orpha-map` reads Orphanet `en_product1.xml` (`data_sources/orphanet_products.py`). Each `Disorder` lists `ExternalReference`s. For `Source = OMIM`, the mapping relation (`E`, `BTNT`, `NTBT`, `ND`) and the validation status are read.

- Only **exact (`E`) and `Validated`** references of **active** disorders become mappings. References of disorders flagged `Inactive`, `Obsolete entity` or `Deprecated entity` are skipped and counted.
- Broader (`BTNT`), narrower (`NTBT`), undefined (`ND`) and not-yet-validated references are counted and ignored. They relate non-equivalent concepts, so no fact moves through them.
- An OMIM id with several exact ORPHA mappings pools their onsets (counted as `omim_ids_with_several_orpha`).
- The mapping is used for two things only: the Orphanet onset fallback of OMIM diseases and `mapped_ids`.

## LLM Extraction Policy

The LLM produces `DiseaseProfilePatch` objects, never synthetic patients.

The extractor should read trusted source text and emit structured claims with evidence:

```json
{
  "disease_id": "ORPHA:123",
  "source": {
    "name": "gene_reviews",
    "url": "...",
    "retrieved_at": "2026-05-20"
  },
  "extractions": {
    "clinical_summary": "...",
    "phenotypes": [
      {
        "mention": "early-onset seizures",
        "hpo_id": "HP:0001250",
        "hpo_label": "Seizure",
        "status": "present",
        "frequency": "frequent",
        "diagnostic_role": "major",
        "onset": "infantile",
        "evidence_span": "patients commonly present with early-onset seizures",
        "confidence": 0.86
      }
    ],
    "negative_phenotypes": [
      {
        "mention": "absence of hepatosplenomegaly",
        "hpo_id": "HP:0001433",
        "hpo_label": "Hepatosplenomegaly",
        "status": "absent",
        "evidence_span": "...",
        "confidence": 0.74
      }
    ],
    "age_of_onset": {
      "category": "infantile",
      "evidence_span": "...",
      "confidence": 0.82
    },
    "inheritance": [
      {
        "value": "autosomal recessive",
        "evidence_span": "...",
        "confidence": 0.9
      }
    ],
    "progression": {
      "value": "progressive",
      "evidence_span": "...",
      "confidence": 0.81
    }
  },
  "warnings": []
}
```

Acceptance rules:

- Every extracted HPO ID must exist in the local HPO version.
- Every extracted claim must have an evidence span unless explicitly marked as inferred.
- Every extracted claim must include confidence.
- Contradictions with structured HPO/Orphadata data become warnings or review items.
- LLM output is merged into profiles only through validation code.
- The same input should produce comparable output by using low temperature and strict schemas.
- Direct LLM-generated patient cases are out of scope.
- Unmentioned phenotypes are not negative phenotypes. A negative phenotype requires explicit textual evidence or a structured negative annotation.

## Disease Profile

`DiseaseProfile` is the consolidated disease-level truth used by the simulator.

It preserves source provenance and separates disease knowledge from simulated cases.

Initial target shape:

```json
{
  "disease_id": "ORPHA:...",
  "disease_name": "...",
  "mapped_ids": {
    "mondo": "MONDO:...",
    "omim": ["OMIM:..."],
    "medgen": ["C..."]
  },
  "genes": [
    {
      "symbol": "NPC1",
      "association_type": "causal",
      "inheritance": ["autosomal recessive"],
      "source": ["orphanet", "hpo"]
    }
  ],
  "phenotypes": [
    {
      "hpo_id": "HP:0001250",
      "label": "Seizure",
      "frequency": "frequent",
      "probability_range": [0.3, 0.79],
      "diagnostic_role": "major",
      "onset": "infantile",
      "severity": "unknown",
      "progression": "unknown",
      "source": ["orphanet", "hpo", "llm_extracted"],
      "confidence": 0.9
    }
  ],
  "negative_phenotypes": [
    {
      "hpo_id": "HP:...",
      "label": "...",
      "source": ["hpo_negative", "llm_extracted"],
      "confidence": 0.8
    }
  ],
  "age_of_onset": {
    "category": "infantile",
    "source": ["orphanet", "llm_extracted"],
    "confidence": 0.8
  },
  "sex_bias": {
    "value": "male|female|none|unknown",
    "confidence": 0.7
  },
  "provenance": [
    {
      "source": "orphanet",
      "url_or_file": "...",
      "retrieved_at": "2026-05-20"
    }
  ],
  "quality": {
    "profile_confidence": 0.86,
    "warnings": []
  }
}
```

## Synthetic Case Design

Every synthetic case stores both disease and gene labels from day one.

Phenotype status is interpreted conservatively:

- `negative`: explicitly absent.
- `missing`: not observed, not asked, or not recorded.
- `unknown`: not enough information.

```json
{
  "case_id": "synthetic-ORPHA_123-000001",
  "target": {
    "disease_id": "ORPHA:123",
    "disease_name": "...",
    "gene": "NPC1",
    "gene_label": 123,
    "disease_label": 456
  },
  "patient": {
    "sex": "female",
    "age": {
      "value": 7,
      "unit": "years"
    },
    "age_of_onset": {
      "value": 1,
      "unit": "years"
    }
  },
  "positive_phenotypes": [],
  "negative_phenotypes": [],
  "missing_phenotypes": [],
  "noise_phenotypes": [],
  "metadata": {
    "generator_version": "0.1.0",
    "profile_version": "...",
    "seed": 1234,
    "difficulty": "easy|medium|hard",
    "completeness": 0.6
  }
}
```

### Label Decision

The rich dataset always keeps both:

- `disease_id`
- `gene`

The first GraPhens-compatible training export will use gene labels because GraPhens currently expects:

```text
phenotype case -> causal gene label
```

The rich dataset keeps `disease_id` because clinical diagnosis is disease-level and because:

- one gene can map to multiple diseases;
- one disease can map to multiple genes;
- two diseases involving the same gene can have different phenotype profiles;
- future reporting and follow-up question modules will likely operate at disease level.

Planned training progression:

1. First training: gene label, for GraPhens compatibility.
2. Second training: disease label export.
3. Later experiment: multitask model with disease logits and gene logits.

## Simulation v0.2

`simulate` turns profiles into synthetic cases. Code: `simulation/simulator.py`, `simulation/confounders.py`, `simulation/sampling.py`; every knob is a field of `SimulationConfig` (`simulation:` in the YAML config, `simulation/schema.py`), and none is fitted to R1. The training repo calibrates them on R1-train only.

### Generative model of one case

1. **Sex.** `P(male)` depends on the disease's sex prior key (`sex_prior_key`):

   | Key | When | Default `sex.p_male` |
   | --- | --- | --- |
   | `male_limited` | male-limited expression (HP:0001475) or Y-linked (HP:0001450) | 1.0 |
   | `female_limited` | female-limited expression (HP:0034344) | 0.0 |
   | `male_biased` | profile `sex_bias = male` (all Mendelian modes X-linked recessive) | 0.9 |
   | `female_biased` | profile `sex_bias = female` | 0.1 |
   | `x_linked_dominant` | no sex bias, X-linked dominant (HP:0001423) among the modes | 0.33 (affected females outnumber males about 2:1) |
   | `unbiased` | anything else | 0.5 |

2. **Onset.** A category is drawn from the profile's `age_of_onset.distribution` (or its `category`); diseases without onset use `age.unknown_onset_prior`, whose default is the onset mix of the 2026-02-16 build (neonatal 0.35, infantile 0.18, childhood 0.14, adult 0.13, antenatal 0.11, juvenile 0.05, variable 0.04). The onset age is uniform in the category's window (`age.onset_years`, HPO definitions):

   | Category | Years |
   | --- | --- |
   | antenatal | 0 |
   | neonatal | 0 – 0.0767 (28 days) |
   | infantile | 0.0767 – 1 |
   | childhood | 1 – 5 |
   | juvenile | 5 – 16 |
   | childhood_or_adolescent | 1 – 16 |
   | adult | 16 – 60 |
   | variable | 0 – 60 |

   Orphanet's "Childhood" (2–11 y) and "Adolescent" (12–18 y) reach these windows through the part-1 category mapping, so their ages are approximate.
3. **Current age** = onset + duration, duration ~ Exponential(mean `age.duration_mean_years` = 5), capped at `age.duration_max_years` = 40, and the age at `age.max_age_years` = 90. The age is never below the onset.
4. **Per-patient probability.** For each profile term, `p ~ Beta(k·μ, k·(1−μ))` with μ the ADR-0007 estimate and k = `frequency.concentration` (default 8). Unknown frequencies use μ = `frequency.unknown_frequency` (default 0.5, close to the 0.505 mean of the known estimates in the 2026-02-16 build). Marginally a term is still present with probability μ; the draw spreads patients around it. μ = 0 or 1 is kept exactly.
5. **Progression.** In diseases with `progression = progressive`, non-congenital terms get `p' = p + (1 − p) · max_boost · (1 − exp(−duration / timescale))` (`progression.max_boost` 0.25, `progression.timescale_years` 10; `max_boost: 0` disables it).
6. **Eligibility.** A term is eligible only if its `sex_restriction` matches the patient, it is not in the other sex's anatomy (terms under `sex.male_only_anchors` HP:0010461, HP:0012874 or `sex.female_only_anchors` HP:0010460, HP:0030012, and not under both), and its own onset, drawn in its category window, is ≤ the current age. Antenatal and congenital (HP:0003577) terms are always eligible; terms without onset are too.
7. **Presence and observation.** Eligible terms are present with probability p. The difficulty preset then records each present term as positive with `positive_observation_rate` (+ `cardinal_observation_boost` for cardinal/major roles), otherwise as `missing` or `unknown` (`missing_vs_unknown_split`), and generalizes recorded terms to a phenotypic parent with `ontology_smoothing_rate` (`source_hpo_id` keeps the original).
8. **At least one observed positive.** A case without a recorded positive is redrawn (a new patient if nothing is eligible), up to `max_redraws` (20); then the most probable admissible term is forced in with `reason: forced_min_one`.
9. **Negatives** (below), then **noise** (below), then **covariate missingness**: sex, age and onset are hidden with `missingness.sex_unknown` 0.10, `age_unknown` 0.20 and `onset_unknown` 0.30, and all negatives with `no_negatives` 0.20. The hidden values still shaped the case.

### Negatives: asked and absent

The count per case is Gamma-Poisson (negative binomial) with mean `negatives_mean` from the difficulty preset and Gamma shape `negatives.count_dispersion` (2; null gives plain Poisson), capped at `negatives.max_per_case` (15). Each slot is assigned to a source with odds `negatives.source_weights` before any term is drawn; a slot its source cannot fill stays empty, so the realized mix follows the weights and a disease with few annotations gets fewer negatives.

| Source | `simulated_origin` | Candidates | Weight within the source | Default share |
| --- | --- | --- | --- | --- |
| (a) own disease | `negative_own_disease` | the disease's terms the patient does not have (including terms not yet eligible at this age) | frequency estimate: a clinician asks about typical features | 0.7 |
| (b) confounders | `negative_confounder` | terms of the `negatives.confounders_top_n` (10) most similar diseases that the true disease does not annotate at any level (not the term, an ancestor or a descendant of an annotated term) | similarity × frequency in the confounder | 0.2 |
| (c) NOT | `negative_not_annotation` | the profile's explicit `NOT`/`Excluded` terms | uniform | 0.1 |

Own-disease terms dominate because real case reports mostly list which typical signs of the syndrome a patient lacks: on R1-dev, 45.1% of excluded terms are annotated to the true gene's diseases against 32.0% for the GNN's top-5 wrong candidates (`diagnostic.ar-training` `docs/experiments/EXP-B-001.md`, Addendum 1). Preset means are 10 (easy), 7 (medium) and 4 (hard), so counts reach the dozen negatives real records carry; with small pools they are lower in practice.

**Confounder index** (`ConfounderIndex`, built once over all profiles). Each disease is its positive terms plus their ancestors; a term's information content is `−ln(share of diseases carrying it)`; only terms with IC ≥ `negatives.min_information_content` (2.0 nats) enter the inverted index. Similarity is the cosine of the IC-weighted vectors. Only the diseases actually simulated are queried, and results are cached.

**Invariant, enforced and checked by `validate`:** a negated term is never equal to, an ancestor of, or a descendant of a term present in the case (recorded, generalized-from, missing, unknown or noise) nor of another negated term; a confounder negative is never annotated to the true disease at any level; sex-inappropriate terms are never negated. Without an ontology the ancestor/descendant parts reduce to identity checks and the anatomy anchors are off.

### Noise

Noise terms come only from `--noise-vocabulary`, a TSV with `hpo_id` and optional `label` and `weight` columns (the training repo builds it from R1-train). The count is Poisson with the preset's `noise_mean` (easy 0, medium 1, hard 2); terms are drawn by weight without replacement, never a profile term, a present term, a term related to a negated one, or the other sex's anatomy.

### Difficulty presets

| Preset field | easy | medium | hard |
| --- | --- | --- | --- |
| `positive_observation_rate` | 0.92 | 0.75 | 0.55 |
| `cardinal_observation_boost` | 0.08 | 0.10 | 0.10 |
| `missing_vs_unknown_split` | 0.8 | 0.6 | 0.5 |
| `negatives_mean` | 10 | 7 | 4 |
| `noise_mean` | 0 | 1 | 2 |
| `ontology_smoothing_rate` | 0.00 | 0.15 | 0.35 |

`simulation.presets` overrides whole presets per difficulty; unlisted difficulties keep the built-in one. A YAML override looks like:

```yaml
simulation:
  cases_per_disease_per_difficulty: 20
  difficulties: [medium]
  seed: 42
  frequency: {concentration: 8.0, unknown_frequency: 0.5}
  negatives: {max_per_case: 15, source_weights: {own_disease: 0.7, confounder: 0.2, not_annotation: 0.1}}
  missingness: {sex_unknown: 0.1, age_unknown: 0.2, onset_unknown: 0.3, no_negatives: 0.2}
```

### Running it

```bash
rare-disease-simulator simulate --profiles outputs/profiles.jsonl \
  --output outputs/rich_cases.jsonl --noise-vocabulary noise.tsv \
  --sample-diseases 300 --cases-per-disease 20 --difficulty medium
```

- `--sample-diseases N` picks N diseases spread evenly over (inheritance class × onset source) strata, seeded by the config seed; `--disease-ids FILE` restricts to listed ids. Confounders always come from all profiles.
- `--cases-per-disease`, `--difficulty` (repeatable) and `--seed` override the config.
- The gene label of a case is unchanged: the first causal gene, else the first gene.

**Provenance and determinism.** Each case's `metadata.source_versions` holds the profile sources (version and sha256), the sha256 of the profiles, `hp.json`, noise vocabulary and labels, and the simulator git SHA (`-dirty` when the tree was dirty); `config_hash` hashes the full `SimulationConfig`. Cases carry no timestamp, so the same inputs, config and seed give byte-identical `rich_cases.jsonl`. The run summary (`<output stem>.summary.json`) holds the full config, its hash, the inputs, the strata quotas, the wall-clock time and the output sha256.

## Validate

`validate --cases rich_cases.jsonl --profiles profiles.jsonl` (code: `validation/cases.py`) writes `<cases stem>.validation.json` and prints a text summary. It takes the priors from the run summary next to the cases (or `--run-summary`), else from the app config. It reports:

- phenotype counts per case (mean, median, min, max, histogram) for positives, negatives, missing, unknown and noise;
- the observed positive rate per difficulty against its preset;
- negatives per case, cases without any, and the share by source;
- missing and unknown rates, and the hidden-sex, hidden-age and hidden-onset shares;
- noise per case and cases with a forced positive;
- sex counts, and the male share per sex prior key against `sex.p_male`;
- age and onset histograms, and the onset-category mix against the profiles' onset distributions;
- **calibration**: per (disease, profile term), the share of the disease's cases where the term is truly present (recorded, generalized from, missing or unknown) against its estimate, in bins [0, 0.05), [0.05, 0.30), [0.30, 0.80), [0.80, 0.99), [0.99, 1]. The error is the case-weighted mean of |observed − expected| per bin. The `ungated_terms` view drops terms that sex, onset or progression act on, so it isolates the frequency model; `all_terms` shows the gating effect;
- **invariant violations**, with examples: no observed positive, duplicate or over-cap negatives, a negative related to a present term, a confounder negative annotated to the true disease, a sex-restricted term for the other sex, age below onset, noise outside the vocabulary, a disease missing from the profiles.

Any violation makes the command exit with status 1. Without `--cases`, `validate` only checks the config, as before. Comparing the same statistics with R1-train is the training repo's job.

## Planned Outputs

Initial exports:

- `profiles.jsonl`: consolidated disease profiles.
- `profile_patches.jsonl`: raw validated LLM extraction patches.
- `rich_cases.jsonl`: full synthetic cases with provenance and metadata.
- `graphens.json`: GraPhens-compatible export grouped by gene.

Future exports:

- Phenopacket JSON.
- Disease-label training export.
- Multitask training export.

## Reproducibility

Generated profiles and cases should record source versions, prompt versions, LLM provider/model, simulator version, configuration hash, random seeds, and generation timestamp.

At minimum, generated artifacts should preserve:

- HPO version.
- Orphadata release or retrieval date.
- MONDO version.
- prompt version.
- LLM provider and model.
- simulator version.
- simulation config hash.
- random seed.
- generation timestamp.

## Privacy

The MVP uses public biomedical sources only. It does not ingest private patient records. Generated synthetic cases are not intended to represent real individuals.

## Planned CLI

```bash
rare-disease-simulator fetch-sources
rare-disease-simulator extract-profile-patches --disease ORPHA:123
rare-disease-simulator build-profiles --from-hpoa --exclude-pmids heldout-pmids.txt
rare-disease-simulator simulate --profiles outputs/profiles.jsonl --sample-diseases 300 --cases-per-disease 20
rare-disease-simulator export-graphens --cases outputs/rich_cases.jsonl
rare-disease-simulator validate --cases outputs/rich_cases.jsonl --profiles outputs/profiles.jsonl
```

## Implementation Plan

The staged implementation plan lives in [plan/implementation_plan.md](plan/implementation_plan.md).

## Proposed Package Layout

```text
rare-disease-simulator/
  configs/
    mvp.yaml
  docs/
    design.md
    plan/
      implementation_plan.md
    schemas.md
    sources.md
  src/rare_disease_simulator/
    data_sources/
    llm_extraction/
      schema.py
      prompts/
      extractor.py
      validator.py
      merge_patch.py
    profiles/
    simulation/
    exports/
    validation/
    cli.py
  tests/
  data/
    raw/
    cache/
    profiles/
  outputs/
```

Generated source data, caches, profiles, and outputs should not be committed unless explicitly curated for a release.

Generated artifacts should preserve source attribution and must not include restricted source text unless redistribution is allowed.

## Non-Goals

- No direct LLM patient generation.
- No free-form web scraping as an MVP source.
- No redistribution of restricted source text without checking license terms.
- No replacement of real curated test sets with synthetic validation only.

## Evaluation Direction

The generator will be evaluated separately from downstream diagnostic models.

Initial checks:

- invalid HPO rate;
- phenotype count distribution;
- observed vs expected phenotype frequency;
- missingness rate;
- noise rate;
- negative phenotype rate;
- disease/gene coverage;
- easy/medium/hard difficulty balance;
- source provenance coverage;
- warnings and low-confidence extraction rate.

Real or curated datasets, including Phenopacket examples, should be reserved for calibration and held-out evaluation. Synthetic train and synthetic validation must be generated with distinct seeds and, when possible, distinct simulator parameters.
