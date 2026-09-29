"""Reader for HPO's ``genes_to_disease.txt`` gene-disease links."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

MISSING_GENE_SYMBOL = "-"


@dataclass(frozen=True)
class GeneDiseaseLink:
    """One gene-disease association row."""

    disease_id: str
    gene_symbol: str
    ncbi_gene_id: str
    association_type: str
    source: str


@dataclass(frozen=True)
class GenesToDisease:
    """Gene links grouped by disease id, plus rows skipped while reading."""

    links: dict[str, tuple[GeneDiseaseLink, ...]]
    rows_total: int
    rows_missing_symbol: int


def read_genes_to_disease(path: Path | str) -> GenesToDisease:
    """Read ``genes_to_disease.txt``, skipping rows without a gene symbol (``-``).

    Duplicate (disease, gene) rows are collapsed; links are sorted by symbol.
    """

    grouped: dict[str, dict[str, GeneDiseaseLink]] = {}
    rows_total = 0
    rows_missing_symbol = 0
    with Path(path).open("r", encoding="utf-8", newline="") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            rows_total += 1
            symbol = (row.get("gene_symbol") or "").strip()
            disease_id = (row.get("disease_id") or "").strip()
            if not symbol or symbol == MISSING_GENE_SYMBOL or not disease_id:
                rows_missing_symbol += 1
                continue
            link = GeneDiseaseLink(
                disease_id=disease_id,
                gene_symbol=symbol,
                ncbi_gene_id=(row.get("ncbi_gene_id") or "").strip(),
                association_type=(row.get("association_type") or "").strip().upper(),
                source=(row.get("source") or "").strip(),
            )
            grouped.setdefault(disease_id, {}).setdefault(symbol, link)

    links = {
        disease_id: tuple(by_symbol[symbol] for symbol in sorted(by_symbol))
        for disease_id, by_symbol in grouped.items()
    }
    return GenesToDisease(
        links=links, rows_total=rows_total, rows_missing_symbol=rows_missing_symbol
    )
