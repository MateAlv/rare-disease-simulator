import gzip
from pathlib import Path

import pytest

from rare_disease_simulator.exports.training import (
    RECORD_FIELDS,
    ExportError,
    card_path_for,
    encode_dataset,
    parse_split,
    split_of,
    training_record,
)
from rare_disease_simulator.simulation.simulator import MERGED_REASON, simulate_cases
from rare_disease_simulator.validation.dataset import validate_dataset
from tests.test_gene_first import _gene_cases, build_ontology
from tests.test_simulation_v02 import _config, _phenotype, _preset, _profile


def test_dataset_validation_checks_the_truth_field(tmp_path: Path) -> None:
    ontology = build_ontology()
    case = next(c for c in _gene_cases(ontology) if c.positive_phenotypes)
    record = training_record(case)
    seizure = {**record, "present": ["HP:0001250"], "true_present": ["HP:0007359"],
               "excluded": ["HP:0000707"]}
    dropped = {**record, "present": ["HP:0001250"], "true_present": ["HP:0001251"],
               "excluded": []}
    card = {"file": {}, "split": {"fractions": {"train": 0.9, "val": 0.05, "test": 0.05}}}

    def violations(records, onto):
        path = tmp_path / "ds.jsonl.gz"
        path.write_bytes(encode_dataset(records))
        return validate_dataset(path, card, None, onto)["violations"]["by_type"]

    with_ontology = violations([seizure], ontology)
    without = violations([seizure], None)
    missing_truth = violations([dropped], ontology)

    assert with_ontology["excluded_true_or_ancestor_of_true"] == 1
    assert "excluded_true_or_ancestor_of_true" not in without
    assert missing_truth["present_not_in_true_present"] == 1


def test_true_present_holds_every_term_the_patient_has() -> None:
    ontology = build_ontology()
    profile = _profile(
        "OMIM:7", [_phenotype("HP:0001250", 1.0), _phenotype("HP:0001251", 1.0),
                   _phenotype("HP:0000505", 0.9)]
    )
    preset = _preset(ontology_smoothing_rate=1.0)
    config = _config(difficulties=["easy"], cases_per_disease_per_difficulty=30,
                     presets={"easy": preset})
    cases = simulate_cases(profile, config, ontology=ontology)
    merged = [m for c in cases for m in c.missing_phenotypes if m.reason == MERGED_REASON]

    assert merged
    for case in cases:
        labelled = case.model_copy(
            update={"target": case.target.model_copy(update={"gene_label": 1}),
                    "metadata": case.metadata.model_copy(update={"case_seed": 1})}
        )
        record = training_record(labelled)
        assert {"HP:0001250", "HP:0001251"} <= set(record["true_present"])
        assert "HP:0000707" in record["present"]
        assert set(record["present"]) <= {"HP:0000707", "HP:0000478"}


def test_training_record_fields_and_split_rule() -> None:
    ontology = build_ontology()
    case = _gene_cases(ontology)[0]

    record = training_record(case)

    assert tuple(record) == RECORD_FIELDS
    assert record["gene"] == "GENE1" and record["gene_index"] == 7
    assert record["present"] == sorted(
        {p.hpo_id for p in case.positive_phenotypes} | {p.hpo_id for p in case.noise_phenotypes}
    )
    assert record["split"] == split_of(case.metadata.case_seed)
    shares = [split_of(seed) for seed in range(20000)]
    assert shares.count("train") / len(shares) == pytest.approx(0.9, abs=0.01)
    assert split_of(5, (("train", 0.0), ("val", 0.0), ("test", 1.0))) == "test"

    target = case.target.model_copy(update={"gene_label": None})
    unlabelled = case.model_copy(update={"target": target})
    with pytest.raises(ExportError, match="gene_label"):
        training_record(unlabelled)


def test_encode_dataset_is_deterministic_gzip() -> None:
    records = [{"case_id": "a", "present": ["HP:1"]}, {"case_id": "b", "present": []}]

    payload = encode_dataset(records)

    assert payload == encode_dataset(records)
    assert payload[4:8] == b"\x00\x00\x00\x00"
    lines = gzip.decompress(payload).decode("utf-8").splitlines()
    assert lines[0] == '{"case_id":"a","present":["HP:1"]}'


def test_parse_split_and_card_path() -> None:
    assert parse_split("0.8,0.1,0.1") == (("train", 0.8), ("val", 0.1), ("test", 0.1))
    for bad in ("0.9,0.1", "0.5,0.5,0.5", "a,b,c", "1.1,-0.05,-0.05"):
        with pytest.raises(ExportError):
            parse_split(bad)
    assert card_path_for(Path("out/ds.jsonl.gz")) == Path("out/ds.card.json")
