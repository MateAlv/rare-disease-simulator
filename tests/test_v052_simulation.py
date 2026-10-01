import json
from collections import Counter

import pytest
from pydantic import ValidationError

from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import (
    FORCED_REASON,
    MERGED_REASON,
    _disease_model,
    budget_scaled_probabilities,
    simulate_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.test_simulation_v02 import TERMS, _phenotype, _profile, build_ontology
from tests.test_v04_simulation import V04, _reporting
from tests.test_v05_simulation import VOCABULARY, _config

HISTOGRAM = {"0": 7, "1": 10, "2": 20, "3": 30, "4": 25, "5": 15}
RICH_IDS = [i for i in TERMS if i not in {"HP:0000118", "HP:0000078", "HP:0010461",
                                           "HP:0010460", "HP:0000047", "HP:0000028",
                                           "HP:0000008"}]
RICH = _profile("OMIM:7", [_phenotype(i, 1.0) for i in RICH_IDS])


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _normalized(cases=3000, histogram=HISTOGRAM, **reporting):
    return _config(cases=cases, reporting={"budget_normalize": True,
                                           "profile_budget_histogram": histogram, **reporting})


def _run(ontology, config, profile=RICH, vocabulary=()):
    return simulate_cases(profile, config, ontology=ontology,
                          reporting=_reporting(ontology, cardinal=None),
                          noise_vocabulary=list(vocabulary))


def _realized(case) -> int:
    shown = sum(1 for p in case.positive_phenotypes if p.simulated_origin == "disease_profile")
    return shown + sum(1 for p in case.missing_phenotypes if p.reason == MERGED_REASON)


def test_scaled_probabilities_sum_to_the_target() -> None:
    q = [0.05, 0.2, 0.4, 0.9, 3.0]

    for target in (1, 2, 3, 4):
        probs, scale = budget_scaled_probabilities(q, target)
        assert sum(probs) == pytest.approx(target, abs=1e-6)
        assert probs == [min(1.0, scale * x) for x in q]
    probs, _ = budget_scaled_probabilities(q, 5)
    assert probs == [1.0] * 5


def test_zero_q_terms_stay_at_zero_and_all_zero_is_uniform() -> None:
    probs, _ = budget_scaled_probabilities([0.0, 0.3, 0.6], 3)
    assert probs[0] == 0.0 and probs[1:] == [1.0, 1.0]

    uniform, scale = budget_scaled_probabilities([0.0, 0.0, 0.0, 0.0], 2)
    assert uniform == pytest.approx([0.5] * 4)
    assert scale == pytest.approx(0.5)


def test_realized_counts_follow_the_histogram_when_truth_is_rich(ontology) -> None:
    cases = _run(ontology, _normalized(cases=4000))

    expected = {1: 10 / 100, 2: 20 / 100, 3: 30 / 100, 4: 25 / 100, 5: 15 / 100}
    drawn = Counter(c.metadata.report_budget for c in cases)
    realized = Counter(_realized(c) for c in cases)
    for k, share in expected.items():
        assert drawn[k] / len(cases) == pytest.approx(share, abs=0.025)
    mean_k = sum(k * n for k, n in drawn.items()) / len(cases)
    mean_realized = sum(k * n for k, n in realized.items()) / len(cases)
    # Rolls are independent, so counts spread around k but their mean follows it.
    assert mean_realized == pytest.approx(mean_k, abs=0.1)
    assert all(c.metadata.report_count == len(
        [p for p in c.positive_phenotypes if p.simulated_origin == "disease_profile"]
    ) for c in cases)


def test_inclusion_follows_q(ontology) -> None:
    config = _normalized(cases=4000)
    model = _disease_model(RICH, config, None, _reporting(ontology, cardinal=None), ("OMIM:7",))
    q = {t.phenotype.hpo_id: t.report_q_raw for t in model.terms}

    cases = _run(ontology, config)

    included = Counter(p.hpo_id for c in cases for p in c.positive_phenotypes
                       if p.reason != FORCED_REASON)
    ranked = sorted(q, key=q.get)
    low, high = ranked[: len(ranked) // 3], ranked[-(len(ranked) // 3):]
    assert sum(included[t] for t in high) > 2 * sum(included[t] for t in low)
    for case in cases:
        effective = {p.source_hpo_id or p.hpo_id: p.report_probability
                     for p in case.positive_phenotypes + case.missing_phenotypes
                     + case.unknown_phenotypes}
        scale = case.metadata.report_scale
        for term, probability in effective.items():
            assert probability == pytest.approx(min(1.0, scale * q[term]), abs=2e-4)


def test_a_budget_above_the_true_terms_reports_all_of_them(ontology) -> None:
    small = _profile("OMIM:8", [_phenotype("HP:0001250", 1.0), _phenotype("HP:0001251", 1.0)])

    cases = _run(ontology, _normalized(cases=300, histogram={"5": 1}), profile=small)

    for case in cases:
        assert case.metadata.report_budget == 5
        assert {p.report_probability for p in case.positive_phenotypes} == {1.0}
        assert _realized(case) == 2


def test_all_zero_scores_are_reported_uniformly(ontology, tmp_path) -> None:
    data = json.loads((V04 / "report_model.json").read_text())
    data["intercept"] = -800.0
    path = tmp_path / "model.json"
    path.write_text(json.dumps(data))
    reporting = _reporting(ontology, path, cardinal=None)
    config = _normalized(cases=3000, histogram={"2": 1})
    model = _disease_model(RICH, config, None, reporting, ("OMIM:7",))
    assert all(t.report_q_raw == 0.0 for t in model.terms)

    cases = simulate_cases(RICH, config, ontology=ontology, reporting=reporting)

    rate = 2 / len(RICH.phenotypes)
    assert {round(p.report_probability, 4) for c in cases for p in c.positive_phenotypes
            if p.reason != FORCED_REASON} == {round(rate, 4)}
    mean = sum(_realized(c) for c in cases) / len(cases)
    assert mean == pytest.approx(2.0, abs=0.15)


def test_validate_reports_realized_counts_and_effective_calibration(ontology) -> None:
    config = _normalized(cases=3000, noise_count="proportional", noise_share=0.24)
    cases = _run(ontology, config, vocabulary=VOCABULARY)

    report = validate_cases(cases, profiles={"OMIM:7": RICH}, ontology=ontology, config=config,
                            noise_vocabulary={t.hpo_id for t in VOCABULARY})

    assert report["violations"]["count"] == 0, report["violations"]
    normalized = report["reporting"]["budget_normalized"]
    assert normalized["cases"] == len(cases)
    assert normalized["realized_mean"] == pytest.approx(normalized["k_mean"], abs=0.15)
    assert normalized["scale_mean"] > 0
    assert report["reporting"]["q_calibration_error"] < 0.05


@pytest.mark.parametrize(
    "reporting",
    [
        {"mode": "independent", "budget_normalize": True},
        {"mode": "report_model", "budget_normalize": True, "profile_budget_histogram": {"1": 1}},
        {"mode": "independent", "budget_normalize": True, "profile_budget_histogram": {"0": 3}},
        {"mode": "independent", "budget_normalize": True, "profile_budget_histogram": {"x": 3}},
    ],
)
def test_budget_normalize_settings_are_checked(reporting) -> None:
    with pytest.raises(ValidationError):
        SimulationConfig(reporting=reporting)


def test_budget_normalized_cases_are_deterministic(ontology) -> None:
    config = _normalized(cases=200, noise_count="proportional", noise_share=0.24)

    first = [c.model_dump_json() for c in _run(ontology, config, vocabulary=VOCABULARY)]
    again = [c.model_dump_json() for c in _run(ontology, config, vocabulary=VOCABULARY)]

    assert first == again
