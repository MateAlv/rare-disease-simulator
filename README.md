# rare-disease-simulator

Synthetic clinical phenotype generation for rare disease diagnosis.

The full project design and MVP decisions live in [docs/README.md](docs/README.md).

## Quick start: HPOA profiles

Build one disease profile per gene-linked disease from the HPO annotation release, dropping annotations from held-out publications:

```bash
rare-disease-simulator build-profiles --from-hpoa \
  --hpo-json hp.json --hpoa phenotype.hpoa --genes-to-disease genes_to_disease.txt \
  --orphanet-ages en_product9_ages.xml --exclude-pmids heldout-pmids.txt \
  --output outputs/profiles.jsonl
rare-disease-simulator simulate --profiles outputs/profiles.jsonl --hpo-json hp.json
```

Inputs, the masking rule and the mapping rules are described in [docs/README.md](docs/README.md#hpoa-structured-backbone). Changes are listed in [CHANGELOG.md](CHANGELOG.md).
