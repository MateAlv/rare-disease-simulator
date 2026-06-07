"""HPO source loading utilities."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class HpoTerm:
    """Minimal HPO term record used by local validation."""

    hpo_id: str
    label: str
    parents: tuple[str, ...]
    is_phenotypic_abnormality: bool


class HpoOntology:
    """In-memory HPO cache with label and ancestry lookup helpers."""

    def __init__(self, terms: dict[str, HpoTerm], *, version: str | None = None) -> None:
        self.terms = terms
        self.version = version

    @classmethod
    def from_tsv(cls, path: Path | str, *, version: str | None = None) -> HpoOntology:
        """Load a minimal TSV with hpo_id, label, parents, and phenotypic flag columns."""

        tsv_path = Path(path)
        terms: dict[str, HpoTerm] = {}
        with tsv_path.open("r", encoding="utf-8", newline="") as file:
            for row in csv.DictReader(file, delimiter="\t"):
                term = HpoTerm(
                    hpo_id=row["hpo_id"],
                    label=row["label"],
                    parents=_parse_parent_ids(row.get("parents", "")),
                    is_phenotypic_abnormality=_parse_bool(
                        row.get("is_phenotypic_abnormality", "")
                    ),
                )
                terms[term.hpo_id] = term
        return cls(terms, version=version)

    def is_valid_hpo_id(self, hpo_id: str) -> bool:
        """Return whether an HPO ID exists in this ontology cache."""

        return hpo_id in self.terms

    def get_label(self, hpo_id: str) -> str | None:
        """Return a term label, or None when the HPO ID is unknown."""

        term = self.terms.get(hpo_id)
        return term.label if term else None

    def get_direct_parents(self, hpo_id: str) -> tuple[str, ...]:
        """Return direct parent IDs for a term."""

        term = self.terms.get(hpo_id)
        return term.parents if term else ()

    def get_ancestors(self, hpo_id: str) -> tuple[str, ...]:
        """Return known transitive ancestors for a term, preserving traversal order."""

        ancestors: list[str] = []
        seen: set[str] = set()

        def visit(term_id: str) -> None:
            for parent_id in self.get_direct_parents(term_id):
                if parent_id in seen:
                    continue
                seen.add(parent_id)
                ancestors.append(parent_id)
                if parent_id in self.terms:
                    visit(parent_id)

        visit(hpo_id)
        return tuple(ancestors)

    def is_phenotypic_abnormality(self, hpo_id: str) -> bool:
        """Return whether a term is within the phenotypic abnormality branch."""

        term = self.terms.get(hpo_id)
        return bool(term and term.is_phenotypic_abnormality)


def _parse_parent_ids(raw_value: str) -> tuple[str, ...]:
    if not raw_value:
        return ()
    normalized = raw_value.replace("|", ",")
    return tuple(parent.strip() for parent in normalized.split(",") if parent.strip())


def _parse_bool(raw_value: str) -> bool:
    return raw_value.strip().lower() in {"1", "true", "yes", "y"}
