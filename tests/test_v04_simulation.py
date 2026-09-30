import json
from collections import Counter
from pathlib import Path

import pytest

from rare_disease_simulator.data_sources.gene_profiles import (
    EntityOption,
    GeneTarget,
    plan_genes,
    read_gene_profiles,
    read_gnn_genes,
)
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.reporting import (
    CardinalIndex,
    Reporting,
    load_report_model,
)
from rare_disease_simulator.simulation.schema import FrequencySettings
from rare_disease_simulator.simulation.simulator import (
    CARDINAL_REASON,
    MERGED_REASON,
    ShrinkageMeanMissing,
    simulate_cases,
    simulate_gene_cases,
    simulation_frequency,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.fixtures.readers import fixture_path
from tests.test_gene_first import ALL_PROFILES, PROFILE_MAP, TARGET
from tests.test_simulation_v02 import (
    TRUE,
    _phenotype,
    _profile,
    build_ontology,
)
from tests.test_simulation_v02 import _config as _base_config

V04 = fixture_path("v04")
CARDINAL = CardinalIndex(
    terms={"OMIM:1": frozenset({"HP:0001251"}), "ORPHA:1": frozenset({"HP:0000518"})}
)


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _config(**overrides):
    return _base_config(**{"reporting": {"mode": "report_model"}, **overrides})


def _reporting(ontology, model_path: Path = V04 / "report_model.json", cardinal=CARDINAL):
    return Reporting.build(
        load_report_model(model_path), ontology=ontology, profiles=ALL_PROFILES,
        cardinal=cardinal,
    )


def _phenotype_raw(hpo_id, raw, estimate):
    return _phenotype(hpo_id, estimate, frequency_raw=raw)


def _true_terms(case) -> set[str]:
    ids = {p.source_hpo_id or p.hpo_id for p in case.positive_phenotypes}
    ids |= {p.hpo_id for p in case.missing_phenotypes}
    ids |= {p.hpo_id for p in case.unknown_phenotypes}
    return ids


def _reported(case) -> int:
    return len(case.positive_phenotypes) + sum(
        1 for p in case.missing_phenotypes if p.reason == MERGED_REASON
    )


def _gene_cases(ontology, config=None, reporting=None, target=TARGET):
    index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)
    return simulate_gene_cases(
        target, PROFILE_MAP, config or _config(), ontology=ontology, confounders=index,
        reporting=reporting or _reporting(ontology),
    )


def test_simulation_frequency_applies_the_configured_count_estimator() -> None:
    counted = _phenotype_raw("HP:0001250", "1/1", 0.75)
    category = _phenotype_raw("HP:0001251", "HP:0040281", 0.895)
    unknown = _phenotype("HP:0001249", None)
    shrink = FrequencySettings(shrinkage_mean=0.4)
    jeffreys = FrequencySettings(count_estimator="jeffreys")

    assert simulation_frequency(counted, shrink) == pytest.approx(0.6)
    assert simulation_frequency(counted, jeffreys) == 0.75
    assert simulation_frequency(category, shrink) == 0.895
    assert simulation_frequency(unknown, shrink) is None
    with pytest.raises(ShrinkageMeanMissing, match="shrinkage_mean"):
        simulation_frequency(counted, FrequencySettings())
    assert simulation_frequency(category, FrequencySettings()) == 0.895


def test_report_mode_needs_a_report_model(ontology) -> None:
    with pytest.raises(ValueError, match="report model"):
        simulate_cases(TRUE, _config(), ontology=ontology)


def test_reported_terms_follow_budget_and_cardinal_rule(ontology) -> None:
    reporting = _reporting(ontology)
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=150), reporting)

    support = set(reporting.budget.budgets)
    for case in cases:
        budget = case.metadata.report_budget
        assert budget in support
        truth = _true_terms(case)
        cardinal = truth & CARDINAL.terms_for(case.target.profile_ids or [])
        assert _reported(case) == min(len(truth), max(budget, len(cardinal)))
        unreported = {p.hpo_id for p in case.missing_phenotypes if p.reason != MERGED_REASON}
        unreported |= {p.hpo_id for p in case.unknown_phenotypes}
        assert not cardinal & unreported
        assert all(p.reason in (None, "generalized", CARDINAL_REASON, "forced_min_one")
                   for p in case.positive_phenotypes)
    assert any(p.reason == CARDINAL_REASON for c in cases for p in c.positive_phenotypes)


def test_reported_count_follows_the_budget_when_truth_is_rich(ontology, tmp_path) -> None:
    ids = ["HP:0001250", "HP:0001251", "HP:0001249", "HP:0000505", "HP:0000518",
           "HP:0001629", "HP:0002650", "HP:0000988"]
    rich = _profile("OMIM:7", [_phenotype(hpo_id, 1.0) for hpo_id in ids])
    config = _config(
        difficulties=["easy"], cases_per_disease_per_difficulty=3000,
        presets={"easy": {**_config().presets["easy"].model_dump(),
                          "ontology_smoothing_rate": 0.0, "noise_mean": 0.0}},
    )

    cases = simulate_cases(rich, config, ontology=ontology, reporting=_reporting(ontology))

    counts = Counter(len(case.positive_phenotypes) for case in cases)
    assert set(counts) == {1, 2, 3, 4}
    for k, expected in {1: 10 / 50, 2: 20 / 50, 3: 15 / 50, 4: 5 / 50}.items():
        assert counts[k] / len(cases) == pytest.approx(expected, abs=0.025)
    assert all(len(_true_terms(case)) == len(ids) for case in cases)


def test_report_score_decides_which_true_terms_are_reported(ontology, tmp_path) -> None:
    data = json.loads((V04 / "report_model.json").read_text())
    data["term_budget"]["histogram"] = {"1": 1}
    data["reportability"]["counts"] = {"HP:0001250": 98, "HP:0001251": 0}
    path = tmp_path / "model.json"
    path.write_text(json.dumps(data))
    reporting = _reporting(ontology, path, cardinal=None)
    pair = _profile("OMIM:8", [_phenotype("HP:0001250", 1.0), _phenotype("HP:0001251", 1.0)])
    config = _config(
        difficulties=["easy"], cases_per_disease_per_difficulty=2000,
        presets={"easy": {**_config().presets["easy"].model_dump(),
                          "ontology_smoothing_rate": 0.0}},
    )

    cases = simulate_cases(pair, config, ontology=ontology, reporting=reporting)

    chosen = Counter(case.positive_phenotypes[0].hpo_id for case in cases)
    seizure = reporting.scorer.score("HP:0001250", 1.0, False)
    ataxia = reporting.scorer.score("HP:0001251", 1.0, False)
    assert chosen["HP:0001250"] / len(cases) == pytest.approx(
        seizure / (seizure + ataxia), abs=0.03
    )
    assert all(len(case.positive_phenotypes) == 1 for case in cases)


def test_truth_stays_age_gated_in_report_mode(ontology) -> None:
    cases = simulate_cases(
        TRUE, _config(cases_per_disease_per_difficulty=200), ontology=ontology,
        reporting=_reporting(ontology),
    )

    for case in cases:
        truth = _true_terms(case)
        if case.patient.age.value < 16.0:
            assert "HP:0000505" not in truth
        if case.patient.sex == "female":
            assert "HP:0000047" not in truth


def test_merged_entity_cases_use_every_profile_of_the_entity(ontology) -> None:
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=200))

    merged = [case for case in cases if case.target.entity_id == "OMIM:1"]
    assert {case.target.disease_id for case in merged} == {"OMIM:1"}
    assert {tuple(case.target.profile_ids) for case in merged} == {("OMIM:1", "ORPHA:1")}
    truth = set().union(*(_true_terms(case) for case in merged))
    assert "HP:0000518" in truth  # only ORPHA:1 annotates it
    assert "HP:0001251" in truth  # only OMIM:1 annotates it
    reasons = {n.reason for case in merged for n in case.negative_phenotypes}
    assert "typical_feature_absent_equivalent_profile" not in reasons
    assert "confounder:ORPHA:1" not in reasons
    single = [case for case in cases if case.target.entity_id == "OMIM:4"]
    assert {tuple(case.target.profile_ids) for case in single} == {("OMIM:4",)}


def test_merged_cardinal_comes_from_any_profile_id(ontology) -> None:
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=200))

    for case in cases:
        if case.target.entity_id != "OMIM:1":
            continue
        truth = _true_terms(case)
        shown = {p.source_hpo_id or p.hpo_id for p in case.positive_phenotypes}
        shown |= {p.hpo_id for p in case.missing_phenotypes if p.reason == MERGED_REASON}
        assert truth & {"HP:0000518", "HP:0001251"} <= shown


def test_v04_cases_are_byte_identical_across_runs(ontology) -> None:
    config = _config(cases_per_disease_per_difficulty=30)

    first = [case.model_dump_json() for case in _gene_cases(ontology, config)]
    again = [case.model_dump_json() for case in _gene_cases(ontology, config)]

    assert first == again


def test_uniform_entity_profiles_keep_the_simulator_03_draw(ontology) -> None:
    cases = _gene_cases(ontology, _config(entity_profiles="uniform"))

    assert {case.target.disease_id for case in cases if case.target.entity_id == "OMIM:1"} == {
        "OMIM:1",
        "ORPHA:1",
    }
    assert all(case.target.profile_ids is None for case in cases)


def test_single_entity_target_matches_disease_first_truth(ontology) -> None:
    target = GeneTarget("G", 1, "G", (EntityOption("OMIM:1", ("OMIM:1",), (), 1.0),))

    cases = _gene_cases(ontology, target=target)

    assert {case.target.profile_ids[0] for case in cases} == {"OMIM:1"}
    assert all(case.metadata.report_budget is not None for case in cases)


def _validate(ontology, cases):
    return validate_cases(
        cases,
        profiles=PROFILE_MAP,
        ontology=ontology,
        config=_config(cases_per_disease_per_difficulty=120),
        report_model=load_report_model(V04 / "report_model.json"),
        cardinal=CARDINAL,
    )


def test_validate_accepts_v04_cases_and_reports_the_reporting_rule(ontology) -> None:
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=120))

    report = _validate(ontology, cases)

    assert report["violations"]["count"] == 0, report["violations"]
    reporting = report["reporting"]
    assert reporting["cases"] == len(cases)
    assert reporting["budget_expected"] == {"1": 0.2, "2": 0.4, "3": 0.3, "4": 0.1}
    assert reporting["total_variation_drawn_vs_expected"] < 0.05
    assert reporting["cardinal"]["reported_terms"] == reporting["cardinal"]["true_terms"] > 0
    assert report["diseases"] == 2
    assert report["calibration"]["all_terms"]["pairs"] > 0


def test_validate_flags_broken_reporting(ontology) -> None:
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=120))
    case = next(
        c for c in cases
        if any(p.reason == CARDINAL_REASON for p in c.positive_phenotypes)
        and len(c.positive_phenotypes) > 1
    )
    cardinal = next(p for p in case.positive_phenotypes if p.reason == CARDINAL_REASON)
    hidden = case.model_copy(
        update={
            "positive_phenotypes": [p for p in case.positive_phenotypes if p is not cardinal],
            "missing_phenotypes": [
                *case.missing_phenotypes,
                cardinal.model_copy(update={"status": "missing", "reason": "not_reported"}),
            ],
        }
    )
    outside = case.model_copy(
        update={"metadata": case.metadata.model_copy(update={"report_budget": 9})}
    )

    report = _validate(ontology, [hidden, outside])

    assert set(report["violations"]["by_type"]) >= {
        "cardinal_not_reported",
        "reported_count_mismatch",
        "budget_outside_histogram",
    }


def test_genes_v2_fields_are_read_like_genes_v1(tmp_path) -> None:
    data = read_gene_profiles(fixture_path("hpoa") / "gene_profiles.json")
    for record in data["genes"].values():
        for entity in record["entities"]:
            entity["r1_train_patients"] = 3
    somatic = data["genes"]["GENEA"]["entities"][1]
    somatic.update({"sim_weight": 0.0, "weight_zero_reason": "somatic_only"})
    path = tmp_path / "genes-v2.json"
    path.write_text(json.dumps(data))

    plan = plan_genes(
        read_gene_profiles(path),
        read_gnn_genes(fixture_path("hpoa") / "gnn_genes.json"),
        {"OMIM:100001", "OMIM:100002", "OMIM:100006", "OMIM:100008", "ORPHA:3001"},
    )

    genea = next(target for target in plan.targets if target.symbol == "GENEA")
    assert [entity.entity for entity in genea.entities] == ["OMIM:100001"]
