"""GraPhens-compatible JSON export.

The baseline GraPhens format groups cases by causal gene, where each case is a
flat list of observed (positive) HPO ids:

```json
{"NPC1": [["HP:0001251", "HP:0000511"], ["HP:0001251"]]}
```

Negatives and metadata are intentionally dropped from the baseline export for
compatibility with the current GraPhens ingest. An ``enriched`` export keeps
negatives and case ids so a modified GraPhens can consume the richer signal.
A row-mapping file links every exported row back to its source ``case_id``.
"""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path

from rare_disease_simulator.simulation.schema import SyntheticCase

GraphensExport = dict[str, list[list[str]]]
RowMapping = dict[str, list[str]]


def _observed_positive_ids(case: SyntheticCase, *, include_noise: bool) -> list[str]:
    ids = [phenotype.hpo_id for phenotype in case.positive_phenotypes]
    if include_noise:
        ids.extend(phenotype.hpo_id for phenotype in case.noise_phenotypes)
    return ids


def build_graphens_export(
    cases: Iterable[SyntheticCase],
    *,
    include_noise: bool = True,
) -> tuple[GraphensExport, RowMapping]:
    """Group cases by gene into the baseline GraPhens structure.

    Returns the export mapping ``gene -> [phenotype_id_list, ...]`` and a
    parallel ``gene -> [case_id, ...]`` mapping so each row is traceable.
    """

    export: GraphensExport = defaultdict(list)
    mapping: RowMapping = defaultdict(list)
    for case in cases:
        gene = case.target.gene
        export[gene].append(_observed_positive_ids(case, include_noise=include_noise))
        mapping[gene].append(case.case_id)
    return dict(export), dict(mapping)


def build_enriched_export(cases: Iterable[SyntheticCase]) -> list[dict[str, object]]:
    """Build an enriched per-case export that retains negatives and labels."""

    enriched: list[dict[str, object]] = []
    for case in cases:
        enriched.append(
            {
                "case_id": case.case_id,
                "gene": case.target.gene,
                "disease_id": case.target.disease_id,
                "positive": [phenotype.hpo_id for phenotype in case.positive_phenotypes],
                "negative": [phenotype.hpo_id for phenotype in case.negative_phenotypes],
                "noise": [phenotype.hpo_id for phenotype in case.noise_phenotypes],
            }
        )
    return enriched


def write_graphens_json(
    path: Path | str,
    cases: Sequence[SyntheticCase],
    *,
    include_noise: bool = True,
    mapping_path: Path | str | None = None,
) -> GraphensExport:
    """Write the baseline GraPhens JSON, plus an optional row-mapping file."""

    export, mapping = build_graphens_export(cases, include_noise=include_noise)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(export, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if mapping_path is not None:
        map_path = Path(mapping_path)
        map_path.parent.mkdir(parents=True, exist_ok=True)
        map_path.write_text(
            json.dumps(mapping, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return export
