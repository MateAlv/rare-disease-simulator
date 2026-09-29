"""Readers for gene profiles (``genes-v1``) and the GNN's gene vocabulary.

``genes-v1`` is built by ``diagnostic.ar-training`` (ADR-0008): for each
HGNC-approved gene it lists the disease entities the gene causes, each with
the HPOA ``profile_ids`` that still have phenotypes after the held-out PMID
mask, its inheritance (ClinGen mode of inheritance preferred, else HPOA) and a
``sim_weight`` normalised over the gene's simulable entities. Its
``vocabulary`` maps each GNN v3.0 symbol to the approved symbol (or null).

The GNN vocabulary is a JSON list whose order defines the class index; index
0 is the placeholder ``-``.
"""

from __future__ import annotations

import gzip
import json
from collections import Counter
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

PLACEHOLDER_GENE = "-"

SKIP_REASONS = (
    "unresolved_symbol",
    "no_gene_record",
    "no_simulable_entity",
    "no_profile_available",
)


@dataclass(frozen=True)
class EntityOption:
    """One simulable disease entity of a gene."""

    entity: str
    profile_ids: tuple[str, ...]
    inheritance: tuple[str, ...]
    sim_weight: float


@dataclass(frozen=True)
class GeneTarget:
    """A GNN class to simulate: its label, class index and disease entities."""

    symbol: str
    index: int
    approved_symbol: str
    entities: tuple[EntityOption, ...]


@dataclass
class GenePlan:
    """Genes to simulate in vocabulary order, plus the genes skipped and why."""

    targets: list[GeneTarget]
    vocabulary_size: int
    skipped: dict[str, list[str]] = field(default_factory=dict)
    counters: Counter[str] = field(default_factory=Counter)

    def summary(self) -> dict[str, Any]:
        return {
            "vocabulary_size": self.vocabulary_size,
            "targets": len(self.targets),
            "skipped": {reason: len(self.skipped.get(reason, [])) for reason in SKIP_REASONS},
            "skipped_symbols": {
                reason: sorted(self.skipped[reason])
                for reason in SKIP_REASONS
                if self.skipped.get(reason)
            },
            "entities": dict(sorted(self.counters.items())),
        }


def read_json_maybe_gzip(path: Path | str) -> Any:
    """Read a JSON file, gunzipping it when it starts with the gzip magic bytes."""

    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8"))


def read_gene_profiles(path: Path | str) -> dict[str, Any]:
    """Read a ``genes-v1`` artifact (``{"genes": {...}, "vocabulary": {...}}``)."""

    data = read_json_maybe_gzip(path)
    if not isinstance(data, dict) or not isinstance(data.get("genes"), dict):
        raise ValueError(f"{path}: expected a gene-profiles object with a 'genes' mapping")
    if not isinstance(data.get("vocabulary"), dict):
        raise ValueError(f"{path}: expected a 'vocabulary' mapping of v3 symbols")
    return data


def read_gnn_genes(path: Path | str) -> list[str]:
    """Read the GNN gene vocabulary: a JSON list whose order is the class index."""

    data = read_json_maybe_gzip(path)
    if not isinstance(data, list) or not all(isinstance(symbol, str) for symbol in data):
        raise ValueError(f"{path}: expected a JSON list of gene symbols")
    if len(set(data)) != len(data):
        raise ValueError(f"{path}: duplicate gene symbols in the vocabulary")
    return data


def plan_genes(
    gene_profiles: Mapping[str, Any],
    vocabulary: Sequence[str],
    available_profiles: Collection[str],
    wanted: Collection[str] | None = None,
) -> GenePlan:
    """Pick the GNN genes that can be simulated.

    A gene is simulated when its v3 symbol resolves to an approved symbol with
    a record, and at least one simulable entity with a positive ``sim_weight``
    has a profile in ``available_profiles``. Entity weights are used as given;
    drawing proportionally to them renormalises over the entities kept.
    """

    genes: Mapping[str, Any] = gene_profiles["genes"]
    symbol_map: Mapping[str, str | None] = gene_profiles["vocabulary"]
    plan = GenePlan(targets=[], vocabulary_size=len(vocabulary))
    wanted_set = set(wanted) if wanted is not None else None

    for index, symbol in enumerate(vocabulary):
        if symbol == PLACEHOLDER_GENE:
            continue
        if wanted_set is not None and symbol not in wanted_set:
            continue
        approved = symbol_map.get(symbol)
        if approved is None:
            plan.skipped.setdefault("unresolved_symbol", []).append(symbol)
            continue
        record = genes.get(approved)
        if record is None:
            plan.skipped.setdefault("no_gene_record", []).append(symbol)
            continue
        simulable = [
            entity
            for entity in record.get("entities", [])
            if entity.get("simulable") and (entity.get("sim_weight") or 0.0) > 0.0
        ]
        if not simulable:
            plan.skipped.setdefault("no_simulable_entity", []).append(symbol)
            continue
        options: list[EntityOption] = []
        for entity in simulable:
            listed = entity.get("profile_ids") or []
            present = tuple(sorted(pid for pid in listed if pid in available_profiles))
            plan.counters["profile_ids_listed"] += len(listed)
            plan.counters["profile_ids_missing"] += len(listed) - len(present)
            if not present:
                plan.counters["entities_without_profile"] += 1
                continue
            options.append(
                EntityOption(
                    entity=str(entity["entity"]),
                    profile_ids=present,
                    inheritance=tuple(sorted(entity.get("inheritance") or [])),
                    sim_weight=float(entity["sim_weight"]),
                )
            )
        if not options:
            plan.skipped.setdefault("no_profile_available", []).append(symbol)
            continue
        plan.counters["entities_used"] += len(options)
        plan.targets.append(
            GeneTarget(
                symbol=symbol,
                index=index,
                approved_symbol=approved,
                entities=tuple(sorted(options, key=lambda option: option.entity)),
            )
        )
    return plan


def simulable_profile_links(gene_profiles: Mapping[str, Any]) -> dict[str, list[str]]:
    """Profile id -> approved symbols of the genes whose simulable entities list it."""

    links: dict[str, set[str]] = {}
    for symbol, record in gene_profiles["genes"].items():
        for entity in record.get("entities", []):
            if not entity.get("simulable"):
                continue
            for profile_id in entity.get("profile_ids") or []:
                links.setdefault(profile_id, set()).add(symbol)
    return {profile_id: sorted(symbols) for profile_id, symbols in sorted(links.items())}
