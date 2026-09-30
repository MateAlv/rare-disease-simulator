from collections import Counter

import pytest

from rare_disease_simulator.data_sources.hpo import HpoOntology, HpoTerm
from rare_disease_simulator.exports.training import training_record
from rare_disease_simulator.simulation.simulator import (
    RELATED_NOISE_REASON,
    SPECIALIZED_REASON,
    NoiseTerm,
    simulate_cases,
)
from rare_disease_simulator.validation.cases import validate_cases
from tests.test_simulation_v02 import TERMS, _phenotype, _profile, build_ontology
from tests.test_v04_simulation import _config as _v04_config
from tests.test_v04_simulation import _reporting

RICH = _profile(
    "OMIM:7",
    [
        _phenotype("HP:0001250", 1.0),  # Seizure: one child, Focal-onset seizure
        _phenotype("HP:0000078", 1.0),  # genital system: male and female subtrees
        _phenotype("HP:0001627", 1.0),  # heart morphology: its only child is annotated
        _phenotype("HP:0001629", 1.0),  # VSD
        _phenotype("HP:0000988", 1.0),  # Skin rash: no child
        _phenotype("HP:0000924", 1.0),  # skeletal system: Scoliosis, Hip dysplasia
    ],
)
MALE_SUBTREE = {"HP:0010461", "HP:0000047", "HP:0000028"}
FEMALE_SUBTREE = {"HP:0010460", "HP:0000008"}
VOCABULARY = [NoiseTerm("HP:0000964", "Eczema"), NoiseTerm("HP:0000505", "Visual impairment")]


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _config(smoothing=0.0, **overrides):
    easy = {**_v04_config().presets["easy"].model_dump(),
            "ontology_smoothing_rate": smoothing, "negatives_mean": 3.0}
    return _v04_config(
        difficulties=["easy"], cases_per_disease_per_difficulty=1500,
        presets={"easy": easy}, **overrides,
    )


def _run(ontology, config, profile=RICH, vocabulary=VOCABULARY):
    return simulate_cases(
        profile, config, ontology=ontology, reporting=_reporting(ontology, cardinal=None),
        noise_vocabulary=vocabulary,
    )


def _validate(ontology, cases, config, profile=RICH, vocabulary=VOCABULARY):
    return validate_cases(
        cases, profiles={profile.disease_id: profile}, ontology=ontology, config=config,
        noise_vocabulary={term.hpo_id for term in vocabulary},
    )


def test_specialized_terms_are_eligible_descendants_that_fill_the_slot(ontology) -> None:
    config = _config(reporting={"mode": "report_model", "specialize_rate": 1.0})

    cases = _run(ontology, config)

    specialized = Counter()
    for case in cases:
        assert len(case.positive_phenotypes) == case.metadata.report_budget_profile
        for p in case.positive_phenotypes:
            if p.reason != SPECIALIZED_REASON:
                # No eligible descendant: the leaf, and the heart term whose child is annotated.
                assert p.hpo_id in {"HP:0000988", "HP:0001627", "HP:0001629"}
                continue
            specialized[(p.source_hpo_id, p.hpo_id)] += 1
            assert ontology.is_a(p.hpo_id, p.source_hpo_id) and p.hpo_id != p.source_hpo_id
            assert p.hpo_id not in {q.hpo_id for q in RICH.phenotypes}
            shown = [q.hpo_id for q in case.positive_phenotypes + case.noise_phenotypes]
            assert len(shown) == len(set(shown))
            if case.patient.sex == "female":
                assert p.hpo_id not in MALE_SUBTREE
            if case.patient.sex == "male":
                assert p.hpo_id not in FEMALE_SUBTREE
        record = training_record(case.model_copy(
            update={"target": case.target.model_copy(update={"gene_label": 1})}
        ))
        assert set(record["present"]) <= set(record["true_present"])
        for p in case.positive_phenotypes:
            if p.reason == SPECIALIZED_REASON:
                assert p.source_hpo_id in record["true_present"]
    seizure = {child for source, child in specialized if source == "HP:0001250"}
    assert seizure == {"HP:0007359"}
    genital = Counter({c: n for (s, c), n in specialized.items() if s == "HP:0000078"})
    grandchildren = sum(n for c, n in genital.items() if c in {"HP:0000047", "HP:0000028",
                                                              "HP:0000008"})
    # A child is taken, then a grandchild with probability 0.5 (the female child
    # has one, the male child two).
    assert grandchildren / sum(genital.values()) == pytest.approx(0.5, abs=0.06)
    report = _validate(ontology, cases, config)
    assert report["violations"]["count"] == 0, report["violations"]
    assert report["positives"]["specialized"] == sum(specialized.values())


def test_generalized_terms_are_never_specialized(ontology) -> None:
    config = _config(
        smoothing=0.6, reporting={"mode": "report_model", "specialize_rate": 1.0}
    )

    cases = _run(ontology, config)

    reasons = Counter(p.reason for case in cases for p in case.positive_phenotypes)
    assert reasons["generalized"] and reasons[SPECIALIZED_REASON]
    for case in cases:
        for p in case.positive_phenotypes:
            if p.reason == "generalized":
                assert ontology.is_a(p.source_hpo_id, p.hpo_id)
    assert _validate(ontology, cases, config)["violations"]["count"] == 0


def test_specialize_rate_sets_the_share(ontology) -> None:
    config = _config(reporting={"mode": "report_model", "specialize_rate": 0.3})

    cases = _run(ontology, config)

    eligible = specialized = 0
    for case in cases:
        for p in case.positive_phenotypes:
            source = p.source_hpo_id or p.hpo_id
            if source in {"HP:0001250", "HP:0000078", "HP:0000924"}:
                eligible += 1
                specialized += p.reason == SPECIALIZED_REASON
    assert specialized / eligible == pytest.approx(0.3, abs=0.04)


def test_related_noise_is_near_the_disease_and_never_annotated(ontology) -> None:
    config = _config(
        reporting={"mode": "report_model", "noise_count": "budget_share", "noise_share": 0.5},
        noise={"related_share": 1.0},
    )

    cases = _run(ontology, config)

    annotated = {p.hpo_id for p in RICH.phenotypes}
    closure = annotated.union(*(ontology.get_ancestor_set(i) for i in annotated))
    related = [n for case in cases for n in case.noise_phenotypes
               if n.reason == RELATED_NOISE_REASON]
    assert related
    for case in cases:
        k = case.metadata.report_budget
        assert len(case.positive_phenotypes) + len(case.noise_phenotypes) == k
        shown = {p.hpo_id for p in case.positive_phenotypes}
        for noise in case.noise_phenotypes:
            assert noise.simulated_origin == "noise"
            if noise.reason != RELATED_NOISE_REASON:
                continue
            assert noise.hpo_id not in closure and noise.hpo_id not in shown
            assert ontology.get_ancestor_set(noise.hpo_id).isdisjoint(annotated)
            seed = noise.source_hpo_id
            assert seed in annotated
            anchors = set(ontology.get_direct_parents(seed))
            anchors |= {g for a in anchors for g in ontology.get_direct_parents(a)}
            assert any(
                noise.hpo_id in ontology.get_direct_children(a)
                or any(noise.hpo_id in ontology.get_direct_children(c)
                       for c in ontology.get_direct_children(a))
                for a in anchors
            )
            if case.patient.sex == "female":
                assert noise.hpo_id not in MALE_SUBTREE
            if case.patient.sex == "male":
                assert noise.hpo_id not in FEMALE_SUBTREE
        negatives = {n.hpo_id for n in case.negative_phenotypes}
        assert not negatives & {n.hpo_id for n in case.noise_phenotypes}
    report = _validate(ontology, cases, config)
    assert report["violations"]["count"] == 0, report["violations"]
    assert report["noise"]["related"] == len(related)


def test_related_share_sets_the_share_of_noise(ontology) -> None:
    config = _config(
        reporting={"mode": "report_model", "noise_count": "budget_share", "noise_share": 0.5},
        noise={"related_share": 0.3},
    )

    cases = _run(ontology, config, vocabulary=[NoiseTerm(i, TERMS[i][0]) for i in (
        "HP:0000964", "HP:0000505", "HP:0000518", "HP:0001385", "HP:0002650")])

    noise = [n for case in cases for n in case.noise_phenotypes]
    related = sum(n.reason == RELATED_NOISE_REASON for n in noise)
    # Some related draws find nothing eligible and fall back to the vocabulary.
    assert 0.15 < related / len(noise) <= 0.34


def test_related_noise_falls_back_to_the_vocabulary(ontology) -> None:
    tiny = HpoOntology(
        {
            "HP:0000118": HpoTerm("HP:0000118", "Phenotypic abnormality", (), True),
            "HP:0000001": HpoTerm("HP:0000001", "A", ("HP:0000118",), True),
            "HP:0000002": HpoTerm("HP:0000002", "B", ("HP:0000001",), True),
            "HP:0000003": HpoTerm("HP:0000003", "C", ("HP:0000001",), True),
            # Outside Phenotypic abnormality, so only the vocabulary can offer it.
            "HP:0000004": HpoTerm("HP:0000004", "D", (), False),
        },
        version="tiny",
    )
    profile = _profile("OMIM:9", [])
    profile = profile.model_copy(update={"phenotypes": [
        _phenotype("HP:0001250", 1.0).model_copy(update={"hpo_id": i, "label": i})
        for i in ("HP:0000002", "HP:0000003")
    ]})
    vocabulary = [NoiseTerm("HP:0000004", "D")]
    config = _config(
        reporting={"mode": "report_model", "noise_count": "budget_share", "noise_share": 0.5},
        noise={"related_share": 1.0},
    )
    cases = simulate_cases(
        profile, config, ontology=tiny,
        reporting=_reporting(build_ontology(), cardinal=None),
        noise_vocabulary=vocabulary,
    )

    noise = [n for case in cases for n in case.noise_phenotypes]
    assert noise
    assert {(n.hpo_id, n.reason) for n in noise} == {("HP:0000004", "nonspecific_finding")}


def test_v042_knobs_are_deterministic(ontology) -> None:
    config = _config(
        smoothing=0.3,
        reporting={"mode": "report_model", "specialize_rate": 0.5,
                   "noise_count": "budget_share", "noise_share": 0.5},
        noise={"related_share": 0.5},
    )

    first = [case.model_dump_json() for case in _run(ontology, config)]
    again = [case.model_dump_json() for case in _run(ontology, config)]

    assert first == again


def test_validate_flags_bad_specializations_and_annotated_related_noise(ontology) -> None:
    config = _config(
        reporting={"mode": "report_model", "specialize_rate": 1.0,
                   "noise_count": "budget_share", "noise_share": 0.5},
        noise={"related_share": 1.0},
    )
    cases = _run(ontology, config)
    case = next(
        c for c in cases
        if any(p.reason == SPECIALIZED_REASON for p in c.positive_phenotypes)
        and any(n.reason == RELATED_NOISE_REASON for n in c.noise_phenotypes)
    )
    positives = [
        p.model_copy(update={"source_hpo_id": "HP:0000988"}) if p.reason == SPECIALIZED_REASON
        else p
        for p in case.positive_phenotypes
    ]
    noise = [
        n.model_copy(update={"hpo_id": "HP:0000707"}) if n.reason == RELATED_NOISE_REASON else n
        for n in case.noise_phenotypes
    ]
    tampered = case.model_copy(
        update={"positive_phenotypes": positives, "noise_phenotypes": noise}
    )

    report = _validate(ontology, [tampered], config)

    assert set(report["violations"]["by_type"]) >= {
        "specialized_not_descendant_of_true_term",
        "related_noise_annotated_to_entity",
    }
