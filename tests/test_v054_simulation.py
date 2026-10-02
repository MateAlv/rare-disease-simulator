import json
import random

import pytest

from rare_disease_simulator.profiles.schema import AgeOfOnset
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.reporting import (
    OnsetAges,
    PresentationAges,
    ReportingArtifactError,
    load_onset_ages,
    load_presentation_ages,
)
from rare_disease_simulator.simulation.schema import AgeSettings
from rare_disease_simulator.simulation.simulator import (
    SexSpecificTerms,
    _category_for_onset_age,
    _disease_model,
    _sample_patient,
    simulate_cases,
    simulate_gene_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.fixtures.readers import fixture_path
from tests.test_gene_first import ALL_PROFILES, PROFILE_MAP, TARGET
from tests.test_simulation_v02 import _phenotype, _profile, build_ontology
from tests.test_v04_simulation import _reporting
from tests.test_v05_simulation import _config
from tests.v04_golden import (
    GOLDEN_053,
    PRESENTATION_AGES,
    SCENARIOS_053,
    golden_file_lines,
    golden_lines,
)

ONSET_AGES = fixture_path("v05") / "onset_ages.tsv"
ADULT = _profile(
    "OMIM:9",
    [_phenotype("HP:0001250", 1.0), _phenotype("HP:0000505", 1.0, onset="adult")],
    age_of_onset=AgeOfOnset(category="adult", distribution={"adult": 1.0}),
)
HEADER = "disease_id\tage_low_years\tage_high_years\tsource\n"


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _run(ontology, config, profile=ADULT, **options):
    return simulate_cases(
        profile, config, ontology=ontology, reporting=_reporting(ontology, cardinal=None),
        noise_vocabulary=[], **options,
    )


def _patients(config, ages, n=300):
    model = _disease_model(ADULT, config, None, None, ("OMIM:9",))
    if ages is not None:
        model.onset_age_range = ages.range_for(("OMIM:9",))
    sex_terms = SexSpecificTerms(None, config.sex)
    rng = random.Random(7)
    return model, [_sample_patient(model, config, rng, sex_terms) for _ in range(n)]


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("disease\tage_low_years\tage_high_years\tsource\n", "header"),
        (HEADER + "OMIM:1\t5\t2\tx\n", "low <= high"),
        (HEADER + "OMIM:1\ta\t2\tx\n", "numbers"),
        (HEADER + "bad id\t1\t2\tx\n", "bad disease id"),
        (HEADER + "OMIM:1\t1\t2\t\n", "empty source"),
        (HEADER + "OMIM:1\t1\t2\tx\nOMIM:1\t1\t3\tx\n", "duplicate"),
    ],
)
def test_onset_age_file_mismatches_fail_loudly(tmp_path, content, message) -> None:
    path = tmp_path / "onset.tsv"
    path.write_text(content)

    with pytest.raises(ReportingArtifactError, match=message):
        load_onset_ages(path)


def test_onset_ages_lookup_takes_the_entity_before_its_profile_ids() -> None:
    ages = load_onset_ages(ONSET_AGES)

    assert ages.range_for(("OMIM:1", "ORPHA:1")) == (2.0, 4.0)
    assert ages.range_for(("OMIM:1", "ORPHA:1", "OMIM:3")) == (2.0, 4.0)
    both = OnsetAges({"OMIM:1": (9.0, 9.5), "ORPHA:1": (2.0, 4.0)})
    assert both.range_for(("OMIM:1", "ORPHA:1")) == (9.0, 9.5)
    assert ages.range_for(("OMIM:77",)) is None


def test_onset_is_drawn_within_the_literature_range_not_the_category_window(ontology) -> None:
    ages = OnsetAges({"OMIM:9": (2.0, 4.0)})
    config = _config(cases=300, missingness={"onset_unknown": 0.0, "age_unknown": 0.0})

    cases = _run(ontology, config, onset_ages=ages)

    onsets = [c.patient.age_of_onset.value for c in cases]
    assert all(2.0 - 0.006 <= onset <= 4.0 + 0.006 for onset in onsets)
    assert max(onsets) - min(onsets) > 1.0
    assert all(c.metadata.onset_age_source == "literature" for c in cases)
    assert all(c.metadata.onset_age_years == (2.0, 4.0) for c in cases)
    assert all(c.patient.onset_category == "childhood" for c in cases)
    assert all(c.patient.age.value >= c.patient.age_of_onset.value for c in cases)


def test_category_is_the_narrowest_window_holding_the_onset() -> None:
    windows = AgeSettings().onset_years

    assert _category_for_onset_age(0.0, windows) == "antenatal"
    assert _category_for_onset_age(0.05, windows) == "neonatal"
    assert _category_for_onset_age(0.5, windows) == "infantile"
    assert _category_for_onset_age(1.0, windows) == "infantile"
    assert _category_for_onset_age(3.0, windows) == "childhood"
    assert _category_for_onset_age(10.0, windows) == "juvenile"
    assert _category_for_onset_age(30.0, windows) == "adult"
    assert _category_for_onset_age(75.0, windows) == "adult"


def test_gating_and_duration_follow_the_mapped_category() -> None:
    config = _config(age={"duration_max_years": 1.0, "duration_mean_by_onset": {"childhood": 0.01}})
    model, early = _patients(config, OnsetAges({"OMIM:9": (2.0, 4.0)}))
    adult_term = next(i for i, term in enumerate(model.terms) if term.onset_window)

    _, category = _patients(config, None)

    assert all(not p.eligible[adult_term] for p in early)
    assert any(p.eligible[adult_term] for p in category)
    assert all(p.onset_category == "childhood" for p in early)
    assert all(p.age - p.onset_age < 0.2 for p in early)


def test_presentation_age_still_sets_the_age(ontology) -> None:
    onsets = OnsetAges({"OMIM:9": (0.1, 0.5)})
    late = PresentationAges({"OMIM:9": (30.0, 40.0)})
    config = _config(cases=200, missingness={"onset_unknown": 0.0, "age_unknown": 0.0})

    cases = _run(ontology, config, onset_ages=onsets, presentation_ages=late)

    for case in cases:
        assert 0.1 - 0.006 <= case.patient.age_of_onset.value <= 0.5 + 0.006
        assert 30.0 - 0.006 <= case.patient.age.value <= 40.0 + 0.006
        assert case.metadata.onset_age_source == "literature"

    raised = _run(
        ontology, config, onset_ages=OnsetAges({"OMIM:9": (10.0, 12.0)}),
        presentation_ages=PresentationAges({"OMIM:9": (0.0, 0.2)}),
    )
    assert all(c.patient.age.value == c.patient.age_of_onset.value for c in raised)


def test_validate_checks_onset_ages(ontology) -> None:
    ages = OnsetAges({"OMIM:9": (2.0, 4.0)})
    config = _config(cases=100, missingness={"onset_unknown": 0.0})
    cases = _run(ontology, config, onset_ages=ages)

    report = validate_cases(
        cases, profiles={"OMIM:9": ADULT}, ontology=ontology, config=config,
        report_model=_reporting(ontology, cardinal=None).model,
    )

    assert report["violations"]["count"] == 0, report["violations"]
    assert report["age"]["onset_age_cases"] == 100
    onset = cases[0].patient.age_of_onset
    tampered = cases[0].model_copy(update={"patient": cases[0].patient.model_copy(
        update={"age_of_onset": onset.model_copy(update={"value": 30.0})}
    )})
    flagged = validate_cases(
        [tampered], profiles={"OMIM:9": ADULT}, ontology=ontology, config=config
    )
    assert flagged["violations"]["by_type"]["onset_age_out_of_range"] == 1


def test_onset_ages_are_deterministic(ontology) -> None:
    ages = OnsetAges({"OMIM:9": (2.0, 4.0)})
    config = _config(cases=50)

    first = [c.model_dump_json() for c in _run(ontology, config, onset_ages=ages)]
    again = [c.model_dump_json() for c in _run(ontology, config, onset_ages=ages)]

    assert first == again


def _without_source(lines: list[str]) -> list[str]:
    out = []
    for line in lines:
        data = json.loads(line)
        data["metadata"].pop("onset_age_source", None)
        out.append(json.dumps(data, sort_keys=True))
    return out


def test_defaults_reproduce_simulator_053_cases() -> None:
    ages = load_presentation_ages(PRESENTATION_AGES)

    assert golden_lines(SCENARIOS_053, ages) == golden_file_lines(GOLDEN_053)


def test_onset_file_matching_nothing_reproduces_simulator_053_cases() -> None:
    ages = load_presentation_ages(PRESENTATION_AGES)
    unmatched = OnsetAges({"OMIM:99999": (1.0, 2.0)})

    lines = golden_lines(SCENARIOS_053, ages, unmatched)

    assert all('"onset_age_source": "category"' in line for line in lines)
    assert _without_source(lines) == _without_source(golden_file_lines(GOLDEN_053))


def test_gene_first_cases_use_the_entity_then_profile_ids(ontology) -> None:
    config = _config(cases=200)
    index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)
    ages = load_onset_ages(ONSET_AGES)

    cases = simulate_gene_cases(
        TARGET, PROFILE_MAP, config, ontology=ontology, confounders=index,
        reporting=_reporting(ontology), onset_ages=ages,
    )

    listed = [c for c in cases if c.target.entity_id == "OMIM:1"]
    other = [c for c in cases if c.target.entity_id != "OMIM:1"]
    assert listed and other
    assert all(c.metadata.onset_age_years == (2.0, 4.0) for c in listed)
    assert all(c.metadata.onset_age_source == "category" for c in other)
