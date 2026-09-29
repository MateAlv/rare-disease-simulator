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
| `--omim-orpha-map` | Orphanet `en_product1` alignments XML (lets OMIM diseases use the Orphanet fallback) | no |
| `--exclude-pmids` | held-out references, one `PMID:n` (or bare `n`) per line | no, but see below |

Any option not given falls back to `sources.*` in the config.

```bash
rare-disease-simulator build-profiles --from-hpoa \
  --hpo-json data/raw/hpo/hp.json \
  --hpoa data/raw/hpo/phenotype.hpoa \
  --genes-to-disease data/raw/hpo/genes_to_disease.txt \
  --orphanet-ages data/raw/orphadata/en_product9_ages.xml \
  --exclude-pmids ../diagnostic.ar-training/data/splits/r1-v1.heldout-pmids.txt \
  --output outputs/profiles.jsonl
```

The command writes `profiles.jsonl` and `profiles.summary.json` (`--summary` overrides the path). The summary holds:

- counts of diseases, genes, phenotypes, negatives, onset sources, inheritance terms and sex bias;
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
- **Positive phenotypes (aspect `P`).** Frequency terms HP:0040280–HP:0040284 map to their HPO-defined ranges (100%, 80–99%, 30–79%, 5–29%, 1–4%). `n/m` and `x%` become point estimates. An empty frequency stays `unknown`, with no range; it never means "always". `frequency_raw` keeps the source value.
- **Duplicate rows of the same (disease, term).** `n/m` counts are pooled (sum n / sum m). Without counts, percentages are averaged. Otherwise, the envelope of the frequency-term ranges is used. References are unioned.
- **Negative phenotypes.** A `P` row becomes a negative phenotype if it has the `NOT` qualifier, the `Excluded` frequency (HP:0040285), or a pooled frequency of `0/m`. If a term has both positive and negative evidence, it stays positive; the conflict is counted in `quality.counters`.
- **Phenotype onset.** The per-row `onset` becomes the phenotype's `onset` and `onset_hpo_id`. When rows disagree, the earliest onset wins.
- **Sex restriction.** The per-row `sex` becomes `sex_restriction` (`male`/`female`), but only when every row for that term names the same sex.
- **Disease age of onset.** Aspect `C` terms under Onset (HP:0003674) become `age_of_onset`:
  - `distribution` holds the share of the disease's onset terms in each category;
  - `category` is the most supported category, with ties going to the earliest;
  - `hpo_ids` lists the terms;
  - Congenital onset maps to `neonatal`.
- **Orphanet fallback.** Without an HPOA onset, `AverageAgeOfOnset` is used (Infancy → infantile, Adolescent → juvenile, Elderly → adult, All ages → variable). ORPHA diseases match directly. OMIM diseases match only through exact, validated alignments in `--omim-orpha-map`.
- **Progression.** Aspect `C` pace-of-progression terms set `progression`.
- **Inheritance (aspect `I`).** Modes of inheritance go on every gene of the disease (`inheritance` labels, `inheritance_hpo_ids`).
- **Sex bias.** It is derived only from inheritance, and only where the genetics imply it:
  - male-limited or female-limited expression → that sex;
  - all Mendelian modes X-linked recessive or Y-linked → `male`;
  - all Mendelian modes autosomal → `none`;
  - anything else (X-linked dominant or unspecified, mitochondrial, mixed) is left unset.
- **Other aspects.** Aspects `H` and `M` are not used yet; they are counted in the summary.

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

## Simulation Decisions

The first simulator will be programmatic and auditable.

It should support:

- phenotype sampling by frequency category;
- cardinal/major phenotype weighting;
- age and onset sampling;
- sex-aware sampling;
- explicit negative phenotypes;
- missingness;
- nonspecific clinical noise;
- ontology smoothing toward more general HPO terms;
- differential-diagnosis confounders;
- difficulty levels: `easy`, `medium`, `hard`;
- deterministic seeds and full provenance.

Frequency categories will be mapped to probability ranges, for example:

```python
FREQUENCY_TO_RANGE = {
    "obligate": (0.95, 1.0),
    "very_frequent": (0.80, 0.99),
    "frequent": (0.30, 0.79),
    "occasional": (0.05, 0.29),
    "very_rare": (0.01, 0.04),
    "excluded": (0.0, 0.0),
    "unknown": (0.05, 0.50),
}
```

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
rare-disease-simulator simulate --profiles data/profiles/profiles.jsonl --cases-per-disease 300
rare-disease-simulator export-graphens --cases outputs/rich_cases.jsonl
rare-disease-simulator validate
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
