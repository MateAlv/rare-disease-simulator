"""HPO source loading utilities."""

from __future__ import annotations

import csv
import json
import re
from collections import deque
from dataclasses import dataclass, replace
from pathlib import Path

PHENOTYPIC_ABNORMALITY_ID = "HP:0000118"

_OBO_PURL_PREFIX = "http://purl.obolibrary.org/obo/"
_ALT_ID_PRED = "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId"
_REPLACED_BY_PRED = "http://purl.obolibrary.org/obo/IAO_0100001"
_CONSIDER_PRED = "http://www.geneontology.org/formats/oboInOwl#consider"
_RELEASE_VERSION_RE = re.compile(r"/releases/(\d{4}-\d{2}-\d{2})/")


@dataclass(frozen=True)
class HpoTerm:
    """Minimal HPO term record used by local validation."""

    hpo_id: str
    label: str
    parents: tuple[str, ...]
    is_phenotypic_abnormality: bool
    synonyms: tuple[str, ...] = ()
    alt_ids: tuple[str, ...] = ()
    obsolete: bool = False
    replaced_by: str | None = None
    consider: tuple[str, ...] = ()


@dataclass(frozen=True)
class TermResolution:
    """Outcome of mapping a possibly outdated HPO ID onto a current term."""

    original_id: str
    hpo_id: str | None
    status: str


class HpoOntology:
    """In-memory HPO cache with label and ancestry lookup helpers.

    ``terms`` holds only current (non-obsolete) terms. Obsolete terms and
    secondary (alternative) IDs are kept aside so annotations that still cite
    them can be resolved with :meth:`resolve`.
    """

    def __init__(
        self,
        terms: dict[str, HpoTerm],
        *,
        version: str | None = None,
        obsolete_terms: dict[str, HpoTerm] | None = None,
        alt_id_map: dict[str, str] | None = None,
    ) -> None:
        self.terms = terms
        self.version = version
        self.obsolete_terms = obsolete_terms or {}
        self.alt_id_map = alt_id_map or {}
        self._children: dict[str, tuple[str, ...]] | None = None
        self._ancestor_sets: dict[str, frozenset[str]] = {}
        self._descendant_sets: dict[str, frozenset[str]] = {}

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
                    synonyms=_parse_synonyms(row.get("synonyms", "")),
                )
                terms[term.hpo_id] = term
        return cls(terms, version=version)

    @classmethod
    def from_json(cls, path: Path | str, *, version: str | None = None) -> HpoOntology:
        """Load the official obographs ``hp.json`` release.

        Only ``is_a`` edges between HP classes are kept. The release date is read
        from the graph metadata unless ``version`` is given explicitly.
        """

        with Path(path).open("r", encoding="utf-8") as file:
            document = json.load(file)
        graph = document["graphs"][0]

        raw_terms: dict[str, dict[str, object]] = {}
        for node in graph.get("nodes", []):
            hpo_id = _curie(node.get("id", ""))
            if hpo_id is None or node.get("type") != "CLASS":
                continue
            meta = node.get("meta") or {}
            properties = meta.get("basicPropertyValues") or []
            raw_terms[hpo_id] = {
                "label": node.get("lbl") or "",
                "obsolete": bool(meta.get("deprecated", False)),
                "synonyms": tuple(
                    synonym["val"] for synonym in meta.get("synonyms") or [] if synonym.get("val")
                ),
                "alt_ids": tuple(p["val"] for p in properties if p.get("pred") == _ALT_ID_PRED),
                "replaced_by": next(
                    (p["val"] for p in properties if p.get("pred") == _REPLACED_BY_PRED), None
                ),
                "consider": tuple(p["val"] for p in properties if p.get("pred") == _CONSIDER_PRED),
            }

        parents: dict[str, list[str]] = {hpo_id: [] for hpo_id in raw_terms}
        for edge in graph.get("edges", []):
            if edge.get("pred") != "is_a":
                continue
            child = _curie(edge.get("sub", ""))
            parent = _curie(edge.get("obj", ""))
            if child in parents and parent in raw_terms and parent not in parents[child]:
                parents[child].append(parent)

        terms: dict[str, HpoTerm] = {}
        obsolete_terms: dict[str, HpoTerm] = {}
        alt_id_map: dict[str, str] = {}
        for hpo_id, raw in raw_terms.items():
            term = HpoTerm(
                hpo_id=hpo_id,
                label=str(raw["label"]),
                parents=tuple(sorted(parents[hpo_id])),
                is_phenotypic_abnormality=False,
                synonyms=raw["synonyms"],  # type: ignore[arg-type]
                alt_ids=raw["alt_ids"],  # type: ignore[arg-type]
                obsolete=bool(raw["obsolete"]),
                replaced_by=raw["replaced_by"],  # type: ignore[arg-type]
                consider=raw["consider"],  # type: ignore[arg-type]
            )
            if term.obsolete:
                obsolete_terms[hpo_id] = term
            else:
                terms[hpo_id] = term
                for alt_id in term.alt_ids:
                    alt_id_map[alt_id] = hpo_id

        ontology = cls(
            terms,
            version=version or _release_version(graph),
            obsolete_terms=obsolete_terms,
            alt_id_map=alt_id_map,
        )
        phenotypic = ontology.get_descendant_set(PHENOTYPIC_ABNORMALITY_ID)
        for hpo_id in phenotypic:
            terms[hpo_id] = replace(terms[hpo_id], is_phenotypic_abnormality=True)
        return ontology

    def is_valid_hpo_id(self, hpo_id: str) -> bool:
        """Return whether an HPO ID is a current term in this ontology cache."""

        return hpo_id in self.terms

    def resolve(self, hpo_id: str) -> TermResolution:
        """Map an HPO ID onto a current term, following alt IDs and replacements.

        ``status`` is one of ``current``, ``alt_id``, ``replaced`` or
        ``unresolved``; ``hpo_id`` is None only when unresolved.
        """

        if hpo_id in self.terms:
            return TermResolution(hpo_id, hpo_id, "current")
        if hpo_id in self.alt_id_map:
            return TermResolution(hpo_id, self.alt_id_map[hpo_id], "alt_id")

        seen: set[str] = set()
        current = hpo_id
        while current in self.obsolete_terms and current not in seen:
            seen.add(current)
            replacement = self.obsolete_terms[current].replaced_by
            if replacement is None:
                break
            if replacement in self.terms:
                return TermResolution(hpo_id, replacement, "replaced")
            if replacement in self.alt_id_map:
                return TermResolution(hpo_id, self.alt_id_map[replacement], "replaced")
            current = replacement
        return TermResolution(hpo_id, None, "unresolved")

    def get_label(self, hpo_id: str) -> str | None:
        """Return a term label, or None when the HPO ID is unknown."""

        term = self.terms.get(hpo_id)
        return term.label if term else None

    def get_synonyms(self, hpo_id: str) -> tuple[str, ...]:
        """Return known synonyms for a term."""

        term = self.terms.get(hpo_id)
        return term.synonyms if term else ()

    def get_direct_parents(self, hpo_id: str) -> tuple[str, ...]:
        """Return direct parent IDs for a term."""

        term = self.terms.get(hpo_id)
        return term.parents if term else ()

    def get_direct_children(self, hpo_id: str) -> tuple[str, ...]:
        """Return direct child IDs for a term."""

        if self._children is None:
            children: dict[str, list[str]] = {}
            for term in self.terms.values():
                for parent_id in term.parents:
                    children.setdefault(parent_id, []).append(term.hpo_id)
            self._children = {key: tuple(sorted(value)) for key, value in children.items()}
        return self._children.get(hpo_id, ())

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

    def get_ancestor_set(self, hpo_id: str) -> frozenset[str]:
        """Return the (cached) set of transitive ancestors, excluding the term itself."""

        cached = self._ancestor_sets.get(hpo_id)
        if cached is None:
            cached = frozenset(self.get_ancestors(hpo_id))
            self._ancestor_sets[hpo_id] = cached
        return cached

    def get_descendant_set(self, hpo_id: str) -> frozenset[str]:
        """Return the (cached) set of transitive descendants, excluding the term itself."""

        cached = self._descendant_sets.get(hpo_id)
        if cached is None:
            seen: set[str] = set()
            queue = deque(self.get_direct_children(hpo_id))
            while queue:
                child_id = queue.popleft()
                if child_id in seen:
                    continue
                seen.add(child_id)
                queue.extend(self.get_direct_children(child_id))
            cached = frozenset(seen)
            self._descendant_sets[hpo_id] = cached
        return cached

    def is_a(self, hpo_id: str, ancestor_id: str) -> bool:
        """Return whether ``hpo_id`` is ``ancestor_id`` or one of its descendants."""

        return hpo_id == ancestor_id or ancestor_id in self.get_ancestor_set(hpo_id)

    def is_phenotypic_abnormality(self, hpo_id: str) -> bool:
        """Return whether a term is within the phenotypic abnormality branch."""

        term = self.terms.get(hpo_id)
        return bool(term and term.is_phenotypic_abnormality)


def _curie(iri: str) -> str | None:
    if not iri.startswith(_OBO_PURL_PREFIX + "HP_"):
        return None
    return "HP:" + iri[len(_OBO_PURL_PREFIX) + 3 :]


def _release_version(graph: dict[str, object]) -> str | None:
    meta = graph.get("meta") or {}
    version_iri = meta.get("version") if isinstance(meta, dict) else None
    if not isinstance(version_iri, str):
        return None
    match = _RELEASE_VERSION_RE.search(version_iri)
    return match.group(1) if match else version_iri


def _parse_parent_ids(raw_value: str) -> tuple[str, ...]:
    if not raw_value:
        return ()
    normalized = raw_value.replace("|", ",")
    return tuple(parent.strip() for parent in normalized.split(",") if parent.strip())


def _parse_synonyms(raw_value: str) -> tuple[str, ...]:
    if not raw_value:
        return ()
    return tuple(synonym.strip() for synonym in raw_value.split("|") if synonym.strip())


def _parse_bool(raw_value: str) -> bool:
    return raw_value.strip().lower() in {"1", "true", "yes", "y"}
