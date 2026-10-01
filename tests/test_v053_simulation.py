from math import comb

import pytest
from pydantic import ValidationError

from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import (
    FORCED_REASON,
    MERGED_REASON,
    NoiseTerm,
    simulate_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.test_simulation_v02 import TERMS, _phenotype, _profile, build_ontology
from tests.test_v04_simulation import _reporting
from tests.test_v05_simulation import IDS, _config

TOTALS = {"0": 5, "2": 10, "3": 20, "4": 30, "5": 25, "6": 15}
RICH = _profile("OMIM:7", [_phenotype(i, 1.0) for i in IDS])
SEX_SPECIFIC = {"HP:0010461", "HP:0000047", "HP:0000028", "HP:0010460", "HP:0000008"}
VOCABULARY = [
    NoiseTerm(i, label)
    for i, (label, _) in TERMS.items()
    if i not in IDS and i not in SEX_SPECIFIC and i != "HP:0000118"
]


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _record(cases=3000, histogram=TOTALS, share=0.24, noise=None, **reporting):
    return _config(cases=cases, noise=noise or {}, reporting={
        "budget_normalize": True, "budget_scope": "record", "profile_budget_histogram": histogram,
        "noise_count": "budget_share", "noise_share": share, **reporting,
    })


def _run(ontology, config, profile=RICH, vocabulary=VOCABULARY):
    return simulate_cases(profile, config, ontology=ontology,
                          reporting=_reporting(ontology, cardinal=None),
                          noise_vocabulary=list(vocabulary))


def _profile_shown(case) -> int:
    shown = sum(1 for p in case.positive_phenotypes if p.simulated_origin == "disease_profile")
    return shown + sum(1 for p in case.missing_phenotypes if p.reason == MERGED_REASON)


def _expected_noise(k: int, share: float) -> float:
    return sum(min(n, k - 1) * comb(k, n) * share**n * (1 - share) ** (k - n)
               for n in range(k + 1))


def test_totals_follow_the_histogram_when_truth_is_rich(ontology) -> None:
    cases = _run(ontology, _record(cases=4000))

    weights = {int(k): n for k, n in TOTALS.items() if k != "0"}
    mean_k = sum(k * n for k, n in weights.items()) / sum(weights.values())
    totals = [_profile_shown(c) + len(c.noise_phenotypes) for c in cases]
    assert sum(totals) / len(totals) == pytest.approx(mean_k, abs=0.1)
    for case in cases:
        k = case.metadata.report_budget
        noise = len(case.noise_phenotypes)
        assert noise <= k - 1
        assert case.metadata.report_budget_profile == k - noise


def test_noise_takes_its_share_of_the_slots(ontology) -> None:
    share = 0.24
    cases = _run(ontology, _record(cases=4000, share=share))

    weights = {int(k): n for k, n in TOTALS.items() if k != "0"}
    total_weight = sum(weights.values())
    expected_noise = sum(_expected_noise(k, share) * n for k, n in weights.items()) / total_weight
    mean_k = sum(k * n for k, n in weights.items()) / total_weight
    noise = sum(len(c.noise_phenotypes) for c in cases) / len(cases)
    assert noise == pytest.approx(expected_noise, abs=0.05)
    assert noise / mean_k == pytest.approx(share, abs=0.03)


def test_a_profile_target_above_the_true_terms_reports_all_of_them(ontology) -> None:
    small = _profile("OMIM:8", [_phenotype("HP:0001250", 1.0), _phenotype("HP:0001251", 1.0)])

    cases = _run(ontology, _record(cases=300, histogram={"6": 1}, share=0.1), profile=small)

    for case in cases:
        assert case.metadata.report_budget_profile == 6 - len(case.noise_phenotypes) >= 2
        assert {p.report_probability for p in case.positive_phenotypes} == {1.0}
        assert _profile_shown(case) == 2


def test_a_profile_target_below_the_true_terms_is_the_expected_count(ontology) -> None:
    cases = _run(ontology, _record(cases=3000, histogram={"4": 1}, share=0.5))

    for case in cases:
        target = case.metadata.report_budget_profile
        effective = [p.report_probability for p in
                     case.positive_phenotypes + case.missing_phenotypes + case.unknown_phenotypes
                     if p.simulated_origin == "disease_profile"]
        assert sum(effective) == pytest.approx(target, abs=0.01)
    mean_target = sum(c.metadata.report_budget_profile for c in cases) / len(cases)
    # A forced_min_one term was shown after every roll failed, so it is not a success.
    rolled = [
        0 if any(p.reason == FORCED_REASON for p in c.positive_phenotypes)
        else c.metadata.report_count
        for c in cases
    ]
    assert sum(rolled) / len(cases) == pytest.approx(mean_target, abs=0.08)


def test_validate_compares_totals_with_k(ontology) -> None:
    config = _record(cases=3000, noise={"related_share": 0.4})
    cases = _run(ontology, config)

    report = validate_cases(cases, profiles={"OMIM:7": RICH}, ontology=ontology, config=config,
                            noise_vocabulary={t.hpo_id for t in VOCABULARY})

    assert report["violations"]["count"] == 0, report["violations"]
    record = report["reporting"]["budget_normalized"]["record"]
    assert record["cases"] == len(cases)
    assert record["total_mean"] == pytest.approx(
        report["reporting"]["budget_normalized"]["k_mean"], abs=0.12
    )
    assert record["noise_mean"] > 0 and record["profile_target_mean"] > 0
    assert report["noise"]["related"] > 0


@pytest.mark.parametrize(
    "reporting",
    [
        {"mode": "report_model", "budget_scope": "record"},
        {"mode": "independent", "budget_scope": "record", "noise_count": "budget_share",
         "noise_share": 0.2},
        {"mode": "independent", "budget_scope": "record", "budget_normalize": True,
         "profile_budget_histogram": {"3": 1}, "noise_count": "proportional",
         "noise_share": 0.2},
        {"mode": "independent", "budget_scope": "record", "budget_normalize": True,
         "profile_budget_histogram": {"3": 1}},
        {"mode": "independent", "budget_normalize": True, "profile_budget_histogram": {"3": 1},
         "noise_count": "budget_share", "noise_share": 0.2},
    ],
)
def test_budget_scope_settings_are_checked(reporting) -> None:
    with pytest.raises(ValidationError):
        SimulationConfig(reporting=reporting)


def test_report_model_mode_keeps_budget_share_with_the_default_scope() -> None:
    config = SimulationConfig(
        reporting={"mode": "report_model", "noise_count": "budget_share", "noise_share": 0.3}
    )

    assert config.reporting.budget_scope == "profile"


def test_record_scope_cases_are_deterministic(ontology) -> None:
    config = _record(cases=200)

    first = [c.model_dump_json() for c in _run(ontology, config)]
    again = [c.model_dump_json() for c in _run(ontology, config)]

    assert first == again
