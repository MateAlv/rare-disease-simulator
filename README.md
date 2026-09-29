# rare-disease-simulator

Synthetic clinical phenotype generation for rare disease diagnosis.

The full project design and MVP decisions live in [docs/README.md](docs/README.md).

## Quick start: HPOA profiles

Build one disease profile per gene-linked disease from the HPO annotation release, dropping annotations from held-out publications:

```bash
rare-disease-simulator build-profiles --from-hpoa \
  --hpo-json hp.json --hpoa phenotype.hpoa --genes-to-disease genes_to_disease.txt \
  --orphanet-ages en_product9_ages.xml --omim-orpha-map en_product1.xml \
  --exclude-pmids heldout-pmids.txt --output outputs/profiles.jsonl
rare-disease-simulator simulate --profiles outputs/profiles.jsonl --hpo-json hp.json \
  --sample-diseases 300 --cases-per-disease 20 --difficulty medium
rare-disease-simulator validate --cases outputs/rich_cases.jsonl --profiles outputs/profiles.jsonl
```

With `configs/mvp.yaml` and `diagnostic.ar-training` checked out next to this repo, every input defaults to the pinned copies, so `build-profiles --from-hpoa` needs no input options.

Inputs, the masking rule, the frequency convention and the mapping rules are described in [docs/README.md](docs/README.md#hpoa-structured-backbone); the simulator model and its knobs in [Simulation v0.2](docs/README.md#simulation-v02), and the report in [Validate](docs/README.md#validate). Changes are listed in [CHANGELOG.md](CHANGELOG.md).
