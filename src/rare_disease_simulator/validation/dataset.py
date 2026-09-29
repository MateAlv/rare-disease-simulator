"""Checks of an export-v2 training dataset against its dataset card."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.exports.training import (
    EXPORT_SEXES,
    RECORD_FIELDS,
    dataset_counts,
    iter_dataset,
    split_of,
)
from rare_disease_simulator.simulation.schema import Difficulty

DIFFICULTIES: tuple[Difficulty, ...] = ("easy", "medium", "hard")
MAX_EXAMPLES = 20
CARD_COUNT_KEYS = ("cases", "genes", "by_split", "by_difficulty", "by_sex", "per_gene_cases")


def validate_dataset(
    path: Path | str,
    card: Mapping[str, Any],
    vocabulary: Sequence[str] | None = None,
    ontology: HpoOntology | None = None,
) -> dict[str, Any]:
    """Re-derive the card's facts from the file; ``violations.count`` > 0 means failure.

    Without an ontology, "excluded is not a true term or an ancestor of one"
    reduces to "excluded is not a true term".
    """

    path = Path(path)
    violations: Counter[str] = Counter()
    examples: list[dict[str, str]] = []

    def violation(kind: str, case_id: str, detail: str = "") -> None:
        violations[kind] += 1
        if len(examples) < MAX_EXAMPLES:
            examples.append({"type": kind, "case_id": case_id, "detail": detail})

    file_card = card.get("file", {})
    sha256 = sha256_file(path)
    if file_card.get("sha256") != sha256:
        violation("sha256_differs_from_card", "", f"{sha256} != {file_card.get('sha256')}")
    if file_card.get("bytes") != path.stat().st_size:
        violation("bytes_differ_from_card", "", str(path.stat().st_size))

    fractions = [(name, float(value)) for name, value in card["split"]["fractions"].items()]
    records: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    index_of: dict[str, int] = {}
    previous: tuple[int, str] | None = None
    for record in iter_dataset(path):
        case_id = str(record.get("case_id"))
        if tuple(record) != RECORD_FIELDS:
            violation("wrong_fields", case_id, ",".join(record))
            continue
        records.append(record)
        _check_record(record, case_id, fractions, violation, ontology)
        if case_id in seen_ids:
            violation("duplicate_case_id", case_id)
        seen_ids.add(case_id)
        key = (record["gene_index"], case_id)
        if previous is not None and key <= previous:
            violation("not_sorted", case_id)
        previous = key
        known = index_of.setdefault(record["gene"], record["gene_index"])
        if known != record["gene_index"]:
            violation("gene_with_two_indices", case_id, record["gene"])
        if vocabulary is not None and not (
            0 <= record["gene_index"] < len(vocabulary)
            and vocabulary[record["gene_index"]] == record["gene"]
        ):
            violation("gene_index_not_in_vocabulary", case_id, record["gene"])

    counts = dataset_counts(records)
    card_counts = card.get("counts", {})
    for key in CARD_COUNT_KEYS:
        if card_counts.get(key) != counts[key]:
            violation("count_differs_from_card", "", key)
    if len(set(index_of.values())) != len(index_of):
        violation("index_with_two_genes", "", "")

    return {
        "dataset": {"path": str(path), "sha256": sha256, "records": len(records)},
        "counts": {key: value for key, value in counts.items() if key != "per_gene_cases"},
        "violations": {
            "count": sum(violations.values()),
            "by_type": dict(sorted(violations.items())),
            "examples": examples,
        },
    }


def _check_record(
    record: dict[str, Any], case_id: str, fractions, violation, ontology: HpoOntology | None
) -> None:
    present, excluded, truth = record["present"], record["excluded"], record["true_present"]
    for name in ("present", "true_present", "excluded"):
        terms = record[name]
        if not isinstance(terms, list) or terms != sorted(set(terms)):
            violation(f"{name}_not_sorted_unique", case_id)
    if not present:
        violation("no_present_term", case_id)
    overlap = set(present) & set(excluded)
    if overlap:
        violation("term_present_and_excluded", case_id, ",".join(sorted(overlap)))
    unreported = set(present) - set(truth)
    if unreported:
        violation("present_not_in_true_present", case_id, ",".join(sorted(unreported)))
    implied = set(truth)
    if ontology is not None:
        for hpo_id in truth:
            implied |= ontology.get_ancestor_set(hpo_id)
    contradicted = implied & set(excluded)
    if contradicted:
        violation("excluded_true_or_ancestor_of_true", case_id, ",".join(sorted(contradicted)))
    if not isinstance(record["gene_index"], int) or record["gene_index"] <= 0:
        violation("bad_gene_index", case_id, str(record["gene_index"]))
    if record["sex"] not in EXPORT_SEXES:
        violation("bad_sex", case_id, str(record["sex"]))
    if record["difficulty"] not in DIFFICULTIES:
        violation("bad_difficulty", case_id, str(record["difficulty"]))
    age, onset = record["age_years"], record["onset_years"]
    for value in (age, onset):
        if value is not None and (not isinstance(value, int | float) or value < 0):
            violation("bad_age", case_id, str(value))
    if isinstance(age, int | float) and isinstance(onset, int | float) and age < onset:
        violation("age_before_onset", case_id, f"{age}<{onset}")
    if not isinstance(record["seed"], int):
        violation("bad_seed", case_id, str(record["seed"]))
    elif record["split"] != split_of(record["seed"], fractions):
        violation("split_not_from_seed", case_id, str(record["split"]))


def format_dataset_report(report: Mapping[str, Any]) -> str:
    counts = report["counts"]
    lines = [
        f"Dataset: {report['dataset']['records']} record(s), {counts['genes']} gene(s), "
        f"sha256 {report['dataset']['sha256']}",
        "Splits: " + ", ".join(f"{k} {v}" for k, v in counts["by_split"].items()),
        f"Per case: present {counts['present_per_case_mean']}, "
        f"true_present {counts['true_present_per_case_mean']}, "
        f"excluded {counts['excluded_per_case_mean']}",
        f"Invariant violations: {report['violations']['count']}",
    ]
    for kind, count in report["violations"]["by_type"].items():
        lines.append(f"  {kind}: {count}")
    return "\n".join(lines)
