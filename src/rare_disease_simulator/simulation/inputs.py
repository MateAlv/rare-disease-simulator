"""Optional file inputs for the simulator: noise vocabulary and label maps."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.simulation.simulator import NoiseTerm


@dataclass(frozen=True)
class LabelMaps:
    """Integer class labels for genes and diseases."""

    genes: dict[str, int] = field(default_factory=dict)
    diseases: dict[str, int] = field(default_factory=dict)


def load_noise_vocabulary(
    path: Path | str, ontology: HpoOntology | None = None
) -> list[NoiseTerm]:
    """Read a TSV with an ``hpo_id`` column and optional ``label`` and ``weight`` columns.

    Blank labels are filled from the ontology when one is given. ``weight``
    (default 1) sets how often a term is drawn relative to the others, e.g.
    its frequency among R1-train patients.
    """

    terms: list[NoiseTerm] = []
    with Path(path).open("r", encoding="utf-8", newline="") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            hpo_id = (row.get("hpo_id") or "").strip()
            if not hpo_id:
                continue
            label = (row.get("label") or "").strip()
            if not label and ontology is not None:
                label = ontology.get_label(hpo_id) or ""
            weight = float((row.get("weight") or "").strip() or 1.0)
            if weight < 0.0:
                raise ValueError(f"{path}: negative weight for {hpo_id}")
            terms.append(NoiseTerm(hpo_id=hpo_id, label=label or hpo_id, weight=weight))
    return terms


def load_label_maps(path: Path | str) -> LabelMaps:
    """Read ``{"genes": {symbol: int}, "diseases": {disease_id: int}}`` (both optional)."""

    with Path(path).open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object with 'genes' and/or 'diseases'")
    return LabelMaps(
        genes={str(key): int(value) for key, value in (data.get("genes") or {}).items()},
        diseases={str(key): int(value) for key, value in (data.get("diseases") or {}).items()},
    )
