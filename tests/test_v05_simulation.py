import statistics
from collections import Counter

import pytest
from pydantic import ValidationError

from rare_disease_simulator.exports.training import training_record
from rare_disease_simulator.profiles.schema import AgeOfOnset
from rare_disease_simulator.simulation.reporting import CardinalIndex
from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import (
    FORCED_REASON,
    MERGED_REASON,
    REPORT_FREQUENCY_FLOOR,
    NoiseTerm,
    _disease_model,
    simulate_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.test_simulation_v02 import TERMS, _phenotype, _profile, build_ontology
from tests.test_v04_simulation import _config as _v04_config
from tests.test_v04_simulation import _reporting

IDS = ["HP:0001250", "HP:0001251", "HP:0001249", "HP:0000505", "HP:0000518", "HP:0001629",
       "HP:0002650", "HP:0000988"]
FREQUENCIES = [1.0, 0.9, 0.7, 0.5, 0.3, 0.8, 0.6, 0.4]
MIXED = _profile("OMIM:7", [_phenotype(i, f) for i, f in zip(IDS, FREQUENCIES, strict=True)])
VOCABULARY = [NoiseTerm(i, TERMS[i][0]) for i in (
    "HP:0000964", "HP:0001385", "HP:0000008", "HP:0000028", "HP:0000047", "HP:0000707")]


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _config(cases=3000, smoothing=0.0, **overrides):
    easy = {**_v04_config().presets["easy"].model_dump(),
            "ontology_smoothing_rate": smoothing, "negatives_mean": 3.0, "noise_mean": 0.0}
    reporting = {"mode": "independent", **overrides.pop("reporting", {})}
    return _v04_config(
        difficulties=["easy"], cases_per_disease_per_difficulty=cases,
        presets={"easy": easy}, reporting=reporting, **overrides,
    )


def _run(ontology, config, profile=MIXED, cardinal=None, vocabulary=VOCABULARY):
    return simulate_cases(
        profile, config, ontology=ontology, reporting=_reporting(ontology, cardinal=cardinal),
        noise_vocabulary=vocabulary,
    )


def _shown_ids(case) -> set[str]:
    ids = {p.source_hpo_id or p.hpo_id for p in case.positive_phenotypes}
    return ids | {p.hpo_id for p in case.missing_phenotypes if p.reason == MERGED_REASON}


def test_report_probability_is_score_over_frequency_capped_at_one(ontology) -> None:
    reporting = _reporting(ontology, cardinal=None)
    profile = _profile("OMIM:8", [_phenotype("HP:0001250", 0.9), _phenotype("HP:0001251", 0.0),
                                  _phenotype("HP:0001249", None)])
    config = _config()

    model = _disease_model(profile, config, None, reporting, ("OMIM:8",))

    seizure, ataxia, unknown = model.terms
    assert seizure.report_probability == pytest.approx(min(1.0, seizure.report_weight / 0.9))
    assert ataxia.report_probability == pytest.approx(
        min(1.0, ataxia.report_weight / REPORT_FREQUENCY_FLOOR)
    )
    assert unknown.report_probability == pytest.approx(min(1.0, unknown.report_weight / 0.5))


def test_each_true_term_is_reported_with_its_own_probability(ontology) -> None:
    config = _config(cases=6000)
    cases = _run(ontology, config, vocabulary=[])
    model = _disease_model(MIXED, config, None, _reporting(ontology, cardinal=None), ("OMIM:7",))
    q = {t.phenotype.hpo_id: t.report_probability for t in model.terms}

    true = Counter()
    rolled = Counter()
    for case in cases:
        assert case.metadata.report_budget is None
        shown = _shown_ids(case)
        forced = {p.hpo_id for p in case.positive_phenotypes if p.reason == FORCED_REASON}
        for p in case.positive_phenotypes + case.missing_phenotypes + case.unknown_phenotypes:
            if p.simulated_origin != "disease_profile":
                continue
            term = p.source_hpo_id or p.hpo_id
            assert p.report_probability == pytest.approx(q[term], abs=1e-4)
            true[term] += 1
            rolled[term] += term in shown and term not in forced
    for term, count in true.items():
        assert rolled[term] / count == pytest.approx(q[term], abs=0.04), term
    # Marginally a term is reported with min(f, score): P(true) * q.
    seizure = model.terms[0]
    marginal = rolled["HP:0001250"] / len(cases)
    assert marginal == pytest.approx(min(seizure.mean, seizure.report_weight), abs=0.04)


def test_nothing_rolled_forces_the_highest_q_term_as_it_is(ontology) -> None:
    config = _config(cases=800, smoothing=1.0, reporting={"specialize_rate": 1.0})
    cases = _run(ontology, config, vocabulary=[])

    forced = [c for c in cases if any(p.reason == FORCED_REASON for p in c.positive_phenotypes)]
    assert forced
    model = _disease_model(MIXED, config, None, _reporting(ontology, cardinal=None), ("OMIM:7",))
    q = {t.phenotype.hpo_id: t.report_probability for t in model.terms}
    for case in cases:
        assert case.positive_phenotypes
    for case in forced:
        (term,) = case.positive_phenotypes
        assert term.source_hpo_id is None
        truth = {p.hpo_id for p in case.missing_phenotypes + case.unknown_phenotypes}
        assert all(q[term.hpo_id] >= q[other] for other in truth)


def test_cardinal_terms_follow_the_same_rule(ontology) -> None:
    cardinal = CardinalIndex(terms={"OMIM:7": frozenset({"HP:0000518"})})
    config = _config(cases=2000)
    cases = _run(ontology, config, cardinal=cardinal, vocabulary=[])

    true = sum("HP:0000518" in {p.hpo_id for p in c.missing_phenotypes + c.unknown_phenotypes}
               | _shown_ids(c) for c in cases)
    shown = sum("HP:0000518" in _shown_ids(c) for c in cases)
    assert 0 < shown < true


def test_proportional_noise_scales_with_the_reported_terms(ontology) -> None:
    share = 0.24
    config = _config(reporting={"noise_count": "proportional", "noise_share": share})
    cases = _run(ontology, config)

    reported = sum(len(c.positive_phenotypes) for c in cases)
    noise = sum(len(c.noise_phenotypes) for c in cases)
    assert noise / reported == pytest.approx(share / (1 - share), abs=0.03)
    assert noise / (noise + reported) == pytest.approx(share, abs=0.02)

    small = _run(ontology, config, vocabulary=VOCABULARY[:1])
    assert max(len(c.noise_phenotypes) for c in small) == 1


@pytest.mark.parametrize(
    "reporting",
    [
        {"mode": "independent", "noise_count": "budget_share", "noise_share": 0.3},
        {"mode": "report_model", "noise_count": "proportional", "noise_share": 0.3},
        {"mode": "independent", "noise_count": "proportional", "noise_share": 1.0},
    ],
)
def test_noise_count_modes_are_checked_against_the_reporting_mode(reporting) -> None:
    with pytest.raises(ValidationError):
        SimulationConfig(reporting=reporting)


def test_duration_mean_can_depend_on_onset(ontology) -> None:
    adult = _profile("OMIM:9", [_phenotype("HP:0001250", 1.0)],
                     age_of_onset=AgeOfOnset(category="adult", distribution={"adult": 1.0}))
    infantile = adult.model_copy(update={
        "disease_id": "OMIM:10",
        "age_of_onset": AgeOfOnset(category="infantile", distribution={"infantile": 1.0}),
    })
    config = _config(cases=3000, age={"duration_mean_by_onset": {"adult": 20.0},
                                      "duration_max_years": 200.0, "max_age_years": 500.0})

    def mean_duration(profile):
        cases = _run(ontology, config, profile=profile, vocabulary=[])
        return statistics.fmean(c.patient.age.value - c.patient.age_of_onset.value
                                for c in cases)

    assert mean_duration(adult) == pytest.approx(20.0, rel=0.08)
    assert mean_duration(infantile) == pytest.approx(5.0, rel=0.08)


def test_validate_calibrates_q_and_skips_budget_checks(ontology) -> None:
    config = _config(
        cases=3000, smoothing=0.3,
        reporting={"noise_count": "proportional", "noise_share": 0.24, "specialize_rate": 0.3},
        noise={"related_share": 0.3},
    )
    cases = _run(ontology, config)

    report = validate_cases(
        cases, profiles={"OMIM:7": MIXED}, ontology=ontology, config=config,
        noise_vocabulary={t.hpo_id for t in VOCABULARY},
    )

    assert report["violations"]["count"] == 0, report["violations"]
    reporting = report["reporting"]
    assert reporting["mode"] == "independent"
    assert reporting["q_calibration"]
    for row in reporting["q_calibration"]:
        if row["terms"] >= 300:
            assert row["reported_share"] == pytest.approx(row["mean_q"], abs=0.06)
    assert reporting["q_calibration_error"] < 0.05
    for case in cases:
        record = training_record(case.model_copy(
            update={"target": case.target.model_copy(update={"gene_label": 1})}
        ))
        assert set(record["present"]) <= set(record["true_present"])


def test_independent_cases_are_deterministic(ontology) -> None:
    config = _config(
        cases=200, smoothing=0.3,
        reporting={"noise_count": "proportional", "noise_share": 0.24, "specialize_rate": 0.3},
        noise={"related_share": 0.3}, age={"duration_mean_by_onset": {"infantile": 9.8}},
    )

    first = [c.model_dump_json() for c in _run(ontology, config)]
    again = [c.model_dump_json() for c in _run(ontology, config)]

    assert first == again
