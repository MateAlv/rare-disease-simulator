"""Export v2: compact, feature-carrying training records for the GNN.

One gzip JSONL line per case, in ``(gene_index, case_id)`` order, with keys in
:data:`RECORD_FIELDS` order. The gzip header carries no name and no time, so
the same cases give the same bytes. Each case goes to the sim ``train``,
``val`` or ``test`` split by a hash of its case seed (:func:`split_of`); sim-val
is for early stopping only.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

from rare_disease_simulator.simulation.schema import SyntheticCase

EXPORT_FORMAT = "training-v2"
RECORD_FIELDS: tuple[str, ...] = (
    "case_id",
    "gene",
    "gene_index",
    "entity",
    "profile_id",
    "present",
    "true_present",
    "excluded",
    "sex",
    "age_years",
    "onset_years",
    "difficulty",
    "seed",
    "split",
)
EXPORT_SEXES = ("male", "female", "unknown")
DEFAULT_SPLIT: tuple[tuple[str, float], ...] = (("train", 0.90), ("val", 0.05), ("test", 0.05))
SPLIT_RULE = (
    "u = first 8 bytes of sha256('split|<seed>') / 2**64; "
    "the first split whose cumulative fraction exceeds u"
)


class ExportError(ValueError):
    """A case that cannot be exported."""


def parse_split(text: str) -> tuple[tuple[str, float], ...]:
    """Parse ``"0.9,0.05,0.05"`` into train/val/test fractions summing to 1."""

    try:
        values = [float(part) for part in text.split(",")]
    except ValueError as exc:
        raise ExportError(f"invalid split {text!r}: expected three numbers") from exc
    if len(values) != 3 or any(value < 0.0 for value in values):
        raise ExportError(f"invalid split {text!r}: expected three non-negative fractions")
    if abs(sum(values) - 1.0) > 1e-9:
        raise ExportError(f"invalid split {text!r}: fractions must sum to 1")
    return tuple(zip(("train", "val", "test"), values, strict=True))


def split_of(seed: int, fractions: Sequence[tuple[str, float]] = DEFAULT_SPLIT) -> str:
    """Deterministic split of a case from its seed (see :data:`SPLIT_RULE`)."""

    digest = hashlib.sha256(f"split|{seed}".encode()).digest()
    u = int.from_bytes(digest[:8], "big") / 2**64
    cumulative = 0.0
    for name, fraction in fractions:
        cumulative += fraction
        if u < cumulative:
            return name
    return fractions[-1][0]


def training_record(
    case: SyntheticCase, fractions: Sequence[tuple[str, float]] = DEFAULT_SPLIT
) -> dict[str, Any]:
    """The export-v2 record of one case.

    ``present`` is every observed term the record shows (recorded positives,
    generalized or not, and noise); ``excluded`` the asked-and-absent terms.
    ``missing``/``unknown`` terms are not in the record, as in a real one.
    ``true_present`` is every term the patient has: ``present``, the specific
    terms generalized positives came from, and the ``missing`` and ``unknown``
    terms. It is ground truth for answering simulated questions, not a model
    input.
    """

    target, metadata, patient = case.target, case.metadata, case.patient
    if target.gene_label is None:
        raise ExportError(f"{case.case_id}: no gene_label (simulate in gene-first mode)")
    if metadata.case_seed is None:
        raise ExportError(f"{case.case_id}: no case_seed (simulator older than 0.3.0)")
    if patient.sex not in EXPORT_SEXES:
        raise ExportError(f"{case.case_id}: sex {patient.sex!r} has no export value")
    present = sorted(
        {term.hpo_id for term in case.positive_phenotypes}
        | {term.hpo_id for term in case.noise_phenotypes}
    )
    true_present = sorted(
        set(present)
        | {term.source_hpo_id for term in case.positive_phenotypes if term.source_hpo_id}
        | {term.hpo_id for term in case.missing_phenotypes}
        | {term.hpo_id for term in case.unknown_phenotypes}
    )
    values = {
        "case_id": case.case_id,
        "gene": target.gene,
        "gene_index": target.gene_label,
        "entity": target.entity_id,
        "profile_id": target.disease_id,
        "present": present,
        "true_present": true_present,
        "excluded": sorted({term.hpo_id for term in case.negative_phenotypes}),
        "sex": patient.sex,
        "age_years": _years(patient.age),
        "onset_years": _years(patient.age_of_onset),
        "difficulty": metadata.difficulty,
        "seed": metadata.case_seed,
        "split": split_of(metadata.case_seed, fractions),
    }
    return {name: values[name] for name in RECORD_FIELDS}


def sort_records(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(records, key=lambda record: (record["gene_index"], record["case_id"]))


def encode_dataset(records: Iterable[dict[str, Any]]) -> bytes:
    """Gzip JSONL bytes with a fixed header (no file name, mtime 0)."""

    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, compresslevel=9, mtime=0) as file:
        for record in records:
            line = json.dumps(record, ensure_ascii=False, separators=(",", ":"))
            file.write(line.encode("utf-8") + b"\n")
    return buffer.getvalue()


def iter_dataset(path: Path | str) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as file:
        for line in file:
            if line.strip():
                yield json.loads(line)


def dataset_counts(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Counts the card records and the dataset validator re-derives."""

    per_gene: Counter[str] = Counter(record["gene"] for record in records)
    index_of = {record["gene"]: record["gene_index"] for record in records}
    present = [len(record["present"]) for record in records]
    true_present = [len(record["true_present"]) for record in records]
    excluded = [len(record["excluded"]) for record in records]
    ages = [record["age_years"] for record in records if record["age_years"] is not None]
    onsets = [record["onset_years"] for record in records if record["onset_years"] is not None]
    return {
        "cases": len(records),
        "genes": len(per_gene),
        "by_split": dict(sorted(Counter(record["split"] for record in records).items())),
        "by_difficulty": dict(sorted(Counter(r["difficulty"] for r in records).items())),
        "by_sex": dict(sorted(Counter(record["sex"] for record in records).items())),
        "age_known": len(ages),
        "onset_known": len(onsets),
        "present_per_case_mean": _mean(present),
        "true_present_per_case_mean": _mean(true_present),
        "excluded_per_case_mean": _mean(excluded),
        "cases_with_excluded": sum(1 for count in excluded if count),
        "per_gene_cases": {
            gene: per_gene[gene] for gene in sorted(per_gene, key=lambda g: (index_of[g], g))
        },
    }


def card_path_for(output: Path) -> Path:
    name = output.name
    for suffix in (".jsonl.gz", ".gz", ".jsonl"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return output.with_name(f"{name}.card.json")


def _years(age: Any) -> float | None:
    if age is None:
        return None
    factor = {"years": 1.0, "months": 1.0 / 12.0, "days": 1.0 / 365.25}[age.unit]
    return round(age.value * factor, 4)


def _mean(values: Sequence[float]) -> float | None:
    return round(sum(values) / len(values), 4) if values else None
