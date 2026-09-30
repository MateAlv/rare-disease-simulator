import json
import math
import random
from collections import Counter
from pathlib import Path

import pytest

from rare_disease_simulator.simulation.reporting import (
    BudgetSampler,
    ReportingArtifactError,
    TermScorer,
    load_cardinal,
    load_report_model,
)
from tests.fixtures.readers import fixture_path
from tests.test_simulation_v02 import _phenotype, _profile, build_ontology

V04 = fixture_path("v04")


def _model_data() -> dict:
    return json.loads((V04 / "report_model.json").read_text())


def _write(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "model.json"
    path.write_text(json.dumps(data))
    return path


def test_report_model_fixture_loads() -> None:
    model = load_report_model(V04 / "report_model.json")

    assert model.artifact_id == "report-model-fixture"
    assert [f.kind for f in model.features] == [
        "frequency", "information_content", "ancestor_count", "log_reportability",
        "cardinal_flag",
    ]
    assert model.term_budget.histogram == {0: 5, 1: 10, 2: 20, 3: 15, 4: 5}


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.update({"unexpected": 1}), "unexpected"),
        (lambda d: d.update({"format_version": 2}), "format_version"),
        (lambda d: d["features"].append({"name": "x", "kind": "age", "coefficient": 1}), "kind"),
        (lambda d: d["features"][0].pop("coefficient"), "coefficient"),
        (lambda d: d["features"][3].update({"base": "e"}), "base"),
        (lambda d: d.pop("reportability"), "reportability"),
        (lambda d: d["tables"].clear(), "tables.information_content"),
        (lambda d: d["term_budget"].update({"histogram": {"0": 4}}), "at least 1"),
        (lambda d: d["term_budget"].update({"histogram": {"many": 4}}), "histogram"),
        (lambda d: d["reportability"]["counts"].update({"HP:1": 3}), "counts"),
        (lambda d: d["features"].append(dict(d["features"][0], name="f2")), "unique"),
    ],
)
def test_report_model_mismatches_fail_loudly(tmp_path, mutate, message) -> None:
    data = _model_data()
    mutate(data)

    with pytest.raises(ReportingArtifactError, match=message):
        load_report_model(_write(tmp_path, data))


def test_scorer_computes_features_as_the_model_defines_them() -> None:
    model = load_report_model(V04 / "report_model.json")
    ontology = build_ontology()
    scorer = TermScorer(model, ontology=ontology, unknown_frequency=0.5)

    values = scorer.values("HP:0001250", 0.9, False)
    assert values == {
        "frequency": 0.9,
        "information_content": 2.0,
        "ancestor_count": 2.0,
        "log_reportability": pytest.approx(math.log(50) - math.log(100)),
        "cardinal": 0.0,
    }
    linear = -1.0 + 2.0 * 0.9 + 0.5 * (2.0 - 2.0) / 2.0 - 0.1 * 2 + 0.8 * math.log(0.5)
    assert scorer.score("HP:0001250", 0.9, False) == pytest.approx(1 / (1 + math.exp(-linear)))

    unknown = scorer.values("HP:0000505", None, True)
    assert unknown["frequency"] == 0.3
    assert unknown["information_content"] == 1.0
    assert unknown["log_reportability"] == pytest.approx(math.log(1) - math.log(100))
    assert unknown["cardinal"] == 1.0
    assert scorer.score("HP:0000505", None, True) > scorer.score("HP:0000505", None, False)


def test_scorer_needs_an_ontology_for_ancestor_counts() -> None:
    model = load_report_model(V04 / "report_model.json")

    with pytest.raises(ReportingArtifactError, match="hp.json"):
        TermScorer(model, ontology=None)


def test_information_content_from_profiles_uses_the_declared_base(tmp_path) -> None:
    data = _model_data()
    data["features"][1] = {
        "name": "information_content", "kind": "information_content", "coefficient": 1.0,
        "source": "profiles", "log_base": "2", "missing_value": 9.0,
    }
    model = load_report_model(_write(tmp_path, data))
    profiles = [
        _profile("OMIM:1", [_phenotype("HP:0001250", 0.5)]),
        _profile("OMIM:2", [_phenotype("HP:0001251", 0.5)]),
        _profile("OMIM:3", [_phenotype("HP:0007359", 0.5)]),
        _profile("OMIM:4", [_phenotype("HP:0000518", 0.5)]),
    ]
    scorer = TermScorer(model, ontology=build_ontology(), profiles=profiles)

    # Seizure is carried by OMIM:1 and, through Focal-onset seizure, OMIM:3.
    assert scorer.values("HP:0001250", 0.5, False)["information_content"] == pytest.approx(1.0)
    assert scorer.values("HP:0000707", 0.5, False)["information_content"] == pytest.approx(
        math.log2(4 / 3)
    )
    assert scorer.values("HP:0000964", 0.5, False)["information_content"] == 9.0


def test_budget_sampler_drops_the_zero_bin_and_follows_the_histogram() -> None:
    sampler = BudgetSampler({0: 5, 1: 10, 2: 20, 3: 10})
    rng = random.Random(4)

    draws = Counter(sampler.draw(rng) for _ in range(8000))

    assert sampler.dropped_zero_cases == 5
    assert set(draws) == {1, 2, 3}
    assert draws[2] / 8000 == pytest.approx(0.5, abs=0.02)
    assert sampler.mean() == pytest.approx(2.0)


def test_cardinal_file_is_read_per_disease_and_resolved(tmp_path) -> None:
    index = load_cardinal(V04 / "cardinal.tsv", build_ontology())

    assert index.terms_for(["OMIM:100008", "ORPHA:3001"]) == frozenset({"HP:0001627"})
    assert index.terms_for(["OMIM:100001"]) == frozenset({"HP:0001251"})
    assert index.terms_for(["OMIM:999"]) == frozenset()
    assert index.summary()["rows_by_kind"] == {
        "cardinal_proxy": 1, "diagnostic_criterion": 1, "pathognomonic": 1,
    }


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("disease\thpo_id\tkind\tsource\n", "header"),
        ("disease_id\thpo_id\tkind\tsource\nOMIM:1\tHP:0001250\thallmark\tx\n", "unknown kind"),
        ("disease_id\thpo_id\tkind\tsource\nOMIM:1\tHP:12\tpathognomonic\tx\n", "bad HPO id"),
        ("disease_id\thpo_id\tkind\tsource\nMONDO:1\tHP:0001250\tpathognomonic\tx\n", "disease"),
        ("disease_id\thpo_id\tkind\tsource\nOMIM:1\tHP:0001250\tpathognomonic\n", "columns"),
    ],
)
def test_cardinal_mismatches_fail_loudly(tmp_path, content, message) -> None:
    path = tmp_path / "cardinal.tsv"
    path.write_text(content)

    with pytest.raises(ReportingArtifactError, match=message):
        load_cardinal(path)
