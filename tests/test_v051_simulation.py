import json
import statistics

import pytest

from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.reporting import (
    PresentationAges,
    ReportingArtifactError,
    load_presentation_ages,
)
from rare_disease_simulator.simulation.simulator import (
    RELATED_NOISE_REASON,
    _disease_model,
    simulate_cases,
    simulate_gene_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.fixtures.readers import fixture_path
from tests.test_gene_first import ALL_PROFILES, PROFILE_MAP, TARGET
from tests.test_simulation_v02 import UNRELATED, build_ontology
from tests.test_v04_simulation import V04, _reporting
from tests.test_v05_simulation import MIXED, VOCABULARY, _config

AGES = fixture_path("v05") / "presentation_ages.tsv"


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _run(ontology, config, profile=MIXED, vocabulary=VOCABULARY, reporting=None, **options):
    return simulate_cases(
        profile, config, ontology=ontology,
        reporting=reporting or _reporting(ontology, cardinal=None),
        noise_vocabulary=vocabulary, **options,
    )


def _profile_terms_per_case(cases) -> float:
    return statistics.fmean(len(c.positive_phenotypes) for c in cases)


def test_q_scale_multiplies_the_report_probability(ontology) -> None:
    reporting = _reporting(ontology, cardinal=None)
    base = _disease_model(MIXED, _config(), None, reporting, ("OMIM:7",))
    scaled = _disease_model(MIXED, _config(reporting={"q_scale": 2.0}), None, reporting,
                            ("OMIM:7",))

    for one, two in zip(base.terms, scaled.terms, strict=True):
        assert two.report_probability == pytest.approx(
            min(1.0, 2.0 * one.report_weight / max(one.mean, 1e-3))
        )
    single = _run(ontology, _config(cases=2000), vocabulary=[])
    double = _run(ontology, _config(cases=2000, reporting={"q_scale": 2.0}), vocabulary=[])
    assert _profile_terms_per_case(double) > _profile_terms_per_case(single) + 0.5
    report = validate_cases(double, profiles={"OMIM:7": MIXED}, ontology=ontology,
                            config=_config(cases=2000, reporting={"q_scale": 2.0}))
    assert report["violations"]["count"] == 0
    assert report["reporting"]["q_mean"] > 0
    assert report["reporting"]["q_capped_share"] > 0


def _related(cases):
    return [n for c in cases for n in c.noise_phenotypes if n.reason == RELATED_NOISE_REASON]


def _noise_config(**noise):
    return _config(
        cases=3000,
        reporting={"noise_count": "proportional", "noise_share": 0.4, "q_scale": 2.0},
        noise={"related_share": 1.0, **noise},
    )


def test_one_level_up_and_down_draws_siblings_of_a_true_term(ontology) -> None:
    near = _related(_run(ontology, _noise_config(related_up_levels=1, related_down_levels=1)))
    wide = _related(_run(ontology, _noise_config()))

    def sibling(noise) -> bool:
        parents = set(ontology.get_direct_parents(noise.source_hpo_id))
        return bool(parents & set(ontology.get_direct_parents(noise.hpo_id)))

    assert near and all(sibling(n) for n in near)
    assert wide and not all(sibling(n) for n in wide)


def test_up_levels_bound_the_anchor_and_down_levels_the_depth(ontology) -> None:
    up_one = _related(_run(ontology, _noise_config(related_up_levels=1)))

    for noise in up_one:
        parents = ontology.get_direct_parents(noise.source_hpo_id)
        assert any(
            noise.hpo_id in ontology.get_direct_children(p)
            or any(noise.hpo_id in ontology.get_direct_children(c)
                   for c in ontology.get_direct_children(p))
            for p in parents
        )


def test_related_max_ic_skips_informative_candidates(ontology, tmp_path) -> None:
    # The fixture model's IC table gives Seizure 2, Ataxia 3 and heart morphology 4,
    # and every other term its missing value 1.
    capped = _related(_run(ontology, _noise_config(related_max_ic=1.5)))
    assert capped and not {n.hpo_id for n in capped} & {"HP:0001250", "HP:0001251",
                                                       "HP:0001627"}

    nothing = _run(ontology, _noise_config(related_max_ic=0.5))
    assert not _related(nothing)
    assert any(c.noise_phenotypes for c in nothing)

    data = json.loads((V04 / "report_model.json").read_text())
    data["features"] = [f for f in data["features"] if f["kind"] != "information_content"]
    path = tmp_path / "model.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="information_content"):
        _run(ontology, _noise_config(related_max_ic=1.5),
             reporting=_reporting(ontology, path, cardinal=None))


def test_presentation_ages_set_the_age_of_listed_entities(ontology) -> None:
    ages = load_presentation_ages(AGES)
    config = _config(cases=400)
    index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)

    cases = simulate_gene_cases(
        TARGET, PROFILE_MAP, config, ontology=ontology, confounders=index,
        reporting=_reporting(ontology), presentation_ages=ages,
    )

    listed = [c for c in cases if c.target.entity_id == "OMIM:1"]  # ORPHA:1 is a profile id
    other = [c for c in cases if c.target.entity_id == "OMIM:4"]
    assert listed and other
    for case in listed:
        assert case.metadata.presentation_age_years == (20.0, 40.0)
        age, onset = case.patient.age.value, case.patient.age_of_onset.value
        assert max(20.0, onset) - 0.01 <= age <= max(40.0, onset) + 0.01
    assert statistics.fmean(c.patient.age.value for c in listed) > 25
    assert all(c.metadata.presentation_age_years is None for c in other)
    report = validate_cases(cases, profiles=PROFILE_MAP, ontology=ontology, config=config,
                            report_model=_reporting(ontology).model)
    assert report["violations"]["count"] == 0, report["violations"]
    assert report["age"]["presentation_age_cases"] == len(listed)

    tampered = listed[0].model_copy(update={"patient": listed[0].patient.model_copy(
        update={"age": listed[0].patient.age.model_copy(update={"value": 70.0})}
    )})
    flagged = validate_cases([tampered], profiles=PROFILE_MAP, ontology=ontology, config=config)
    assert flagged["violations"]["by_type"]["presentation_age_out_of_range"] == 1


def test_presentation_age_is_raised_to_onset_and_cut_at_max_age(ontology) -> None:
    profile = UNRELATED.model_copy(update={"age_of_onset": None})
    ages = PresentationAges(ranges={"OMIM:3": (0.0, 0.5)})
    config = _config(cases=400, age={"unknown_onset_prior": {"adult": 1.0}})

    cases = _run(ontology, config, profile=profile, vocabulary=[], presentation_ages=ages)

    for case in cases:
        assert case.patient.age.value == case.patient.age_of_onset.value
    old = PresentationAges(ranges={"OMIM:3": (80.0, 200.0)})
    capped = _run(ontology, _config(cases=200, age={"max_age_years": 85.0}), profile=profile,
                  vocabulary=[], presentation_ages=old)
    assert max(c.patient.age.value for c in capped) <= 85.0


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("disease\tage_low_years\tage_high_years\tsource\n", "header"),
        ("disease_id\tage_low_years\tage_high_years\tsource\nOMIM:1\t5\t2\tx\n", "low <= high"),
        ("disease_id\tage_low_years\tage_high_years\tsource\nOMIM:1\ta\t2\tx\n", "numbers"),
        ("disease_id\tage_low_years\tage_high_years\tsource\nOMIM:1\t1\t2\tx\nOMIM:1\t1\t3\tx\n",
         "duplicate"),
    ],
)
def test_presentation_age_file_mismatches_fail_loudly(tmp_path, content, message) -> None:
    path = tmp_path / "ages.tsv"
    path.write_text(content)

    with pytest.raises(ReportingArtifactError, match=message):
        load_presentation_ages(path)


def test_v051_knobs_are_deterministic(ontology) -> None:
    config = _noise_config(related_up_levels=1, related_down_levels=1, related_max_ic=1.5)
    ages = PresentationAges(ranges={"OMIM:7": (10.0, 30.0)})

    first = [c.model_dump_json() for c in _run(ontology, config, presentation_ages=ages)]
    again = [c.model_dump_json() for c in _run(ontology, config, presentation_ages=ages)]

    assert first == again
