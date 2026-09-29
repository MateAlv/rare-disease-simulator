import gzip
import json
import random
from collections import Counter
from pathlib import Path

import pytest

from rare_disease_simulator.data_sources.gene_profiles import (
    EntityOption,
    GeneTarget,
    plan_genes,
    read_gene_profiles,
    read_gnn_genes,
    simulable_profile_links,
)
from rare_disease_simulator.profiles.inheritance import derive_sex_bias
from rare_disease_simulator.profiles.schema import DiseaseGene
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.simulator import (
    EQUIVALENT_PROFILE_REASON,
    NoiseTerm,
    SexSpecificTerms,
    _disease_model,
    _Patient,
    _sample_noise,
    entity_sex_prior_key,
    simulate_cases,
    simulate_gene_cases,
)
from tests.fixtures.readers import fixture_path
from tests.test_simulation_v02 import (
    CONFOUNDER,
    TRUE,
    UNRELATED,
    _config,
    _phenotype,
    _present_ids,
    _profile,
    _related,
    build_ontology,
)

FIXTURES = fixture_path("hpoa")

EQUIVALENT = _profile(
    "ORPHA:1",
    [
        _phenotype("HP:0001250", 0.9),
        _phenotype("HP:0000518", 0.9),
        _phenotype("HP:0001385", 0.7),
    ],
)
OTHER = _profile("OMIM:4", [_phenotype("HP:0000964", 0.8), _phenotype("HP:0000028", 0.6)])
ALL_PROFILES = [TRUE, CONFOUNDER, UNRELATED, EQUIVALENT, OTHER]
PROFILE_MAP = {profile.disease_id: profile for profile in ALL_PROFILES}
TARGET = GeneTarget(
    symbol="GENE1",
    index=7,
    approved_symbol="GENE1",
    entities=(
        EntityOption("OMIM:1", ("OMIM:1", "ORPHA:1"), ("HP:0000006",), 0.8),
        EntityOption("OMIM:4", ("OMIM:4",), (), 0.2),
    ),
)


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _gene_cases(ontology, config=None, target=TARGET):
    config = config or _config()
    index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)
    return simulate_gene_cases(
        target, PROFILE_MAP, config, ontology=ontology, confounders=index
    )


def _negatives_config(**negatives):
    return _config(
        difficulties=["easy"],
        cases_per_disease_per_difficulty=200,
        negatives={"max_per_case": 15, **negatives},
    )


def test_read_gene_profiles_accepts_gzip_and_plain(tmp_path: Path) -> None:
    plain = FIXTURES / "gene_profiles.json"
    packed = tmp_path / "genes.json.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))

    assert read_gene_profiles(packed) == read_gene_profiles(plain)


def test_read_gnn_genes_rejects_duplicates_and_non_lists(tmp_path: Path) -> None:
    assert read_gnn_genes(FIXTURES / "gnn_genes.json")[:2] == ["-", "GENEA"]
    duplicated = tmp_path / "dup.json"
    duplicated.write_text('["-", "A", "A"]')
    with pytest.raises(ValueError, match="duplicate"):
        read_gnn_genes(duplicated)
    mapping = tmp_path / "map.json"
    mapping.write_text('{"A": 1}')
    with pytest.raises(ValueError, match="list"):
        read_gnn_genes(mapping)


def test_plan_genes_keeps_vocabulary_order_and_counts_skips() -> None:
    data = read_gene_profiles(FIXTURES / "gene_profiles.json")
    vocabulary = read_gnn_genes(FIXTURES / "gnn_genes.json")
    available = {"OMIM:100001", "OMIM:100002", "OMIM:100008", "ORPHA:3001"}

    plan = plan_genes(data, vocabulary, available)

    assert [(t.symbol, t.index, t.approved_symbol) for t in plan.targets] == [
        ("GENEA", 1, "GENEA"),
        ("GENEB", 2, "GENEB"),
        ("OLDC", 3, "GENEC"),
    ]
    assert [e.entity for e in plan.targets[0].entities] == ["OMIM:100001"]
    assert plan.targets[2].entities[0].profile_ids == ("OMIM:100008", "ORPHA:3001")
    summary = plan.summary()
    assert summary["skipped"] == {
        "unresolved_symbol": 2,
        "no_gene_record": 1,
        "no_simulable_entity": 1,
        "no_profile_available": 1,
    }
    assert summary["skipped_symbols"]["unresolved_symbol"] == ["GENEU", "GENEX"]
    assert summary["entities"]["entities_without_profile"] == 2

    wanted = plan_genes(data, vocabulary, available | {"OMIM:100006"}, wanted={"GENEA"})
    assert [t.symbol for t in wanted.targets] == ["GENEA"]
    assert [e.entity for e in wanted.targets[0].entities] == ["OMIM:100001", "OMIM:100006"]


def test_simulable_profile_links_list_approved_symbols() -> None:
    links = simulable_profile_links(read_gene_profiles(FIXTURES / "gene_profiles.json"))

    assert links["OMIM:100006"] == ["GENEA"]
    assert links["ORPHA:3001"] == ["GENEC"]
    assert "OMIM:100009" not in links and "OMIM:100010" not in links


def test_derive_sex_bias_without_ontology_uses_identity() -> None:
    assert derive_sex_bias(["HP:0001419"], None) == ("male", ["HP:0001419"])
    assert derive_sex_bias(["HP:0000006", "HP:0000007"], None)[0] == "none"
    assert derive_sex_bias(["HP:0001423"], None) == (None, [])
    assert derive_sex_bias(["HP:0000006", "HP:0001475"], None) == ("male", ["HP:0001475"])


def test_entity_inheritance_sets_the_sex_prior() -> None:
    base = _profile("OMIM:9", [_phenotype("HP:0001250", 0.5)])
    male_limited = base.model_copy(
        update={
            "genes": [DiseaseGene(symbol="G", inheritance_hpo_ids=["HP:0001423", "HP:0001475"])]
        }
    )
    x_dominant = base.model_copy(
        update={"genes": [DiseaseGene(symbol="G", inheritance_hpo_ids=["HP:0001423"])]}
    )

    assert entity_sex_prior_key(["HP:0001419"], x_dominant, None) == "male_biased"
    assert entity_sex_prior_key(["HP:0000006"], x_dominant, None) == "unbiased"
    assert entity_sex_prior_key([], x_dominant, None) == "x_linked_dominant"
    # Sex-limited expression only comes from HPOA, so it survives an entity override.
    assert entity_sex_prior_key(["HP:0000006"], male_limited, None) == "male_limited"


def test_gene_cases_are_labelled_with_the_gene_and_byte_stable(ontology) -> None:
    cases = _gene_cases(ontology)

    assert len(cases) == 180
    assert {case.target.gene for case in cases} == {"GENE1"}
    assert {case.target.gene_label for case in cases} == {7}
    assert cases[0].case_id == "synthetic-gene-GENE1-easy-000000"
    assert all(case.metadata.case_seed is not None for case in cases)
    assert {case.metadata.sex_prior_key for case in cases if case.target.entity_id == "OMIM:1"} == {
        "unbiased"
    }
    again = _gene_cases(ontology)
    assert [c.model_dump_json() for c in cases] == [c.model_dump_json() for c in again]


def test_gene_cases_draw_entities_by_weight_and_profiles_uniformly(ontology) -> None:
    cases = _gene_cases(ontology, _config(cases_per_disease_per_difficulty=400))

    entities = Counter(case.target.entity_id for case in cases)
    assert entities["OMIM:1"] / len(cases) == pytest.approx(0.8, abs=0.03)
    profiles = Counter(
        case.target.disease_id for case in cases if case.target.entity_id == "OMIM:1"
    )
    assert set(profiles) == {"OMIM:1", "ORPHA:1"}
    assert profiles["OMIM:1"] / sum(profiles.values()) == pytest.approx(0.5, abs=0.04)
    assert {
        case.target.disease_id for case in cases if case.target.entity_id == "OMIM:4"
    } == {"OMIM:4"}


def test_gene_other_disease_negatives_come_from_the_genes_other_entities(ontology) -> None:
    config = _negatives_config(
        source_weights={"own_gene_other_disease": 1.0}, unfilled_slots="drop"
    )
    cases = [c for c in _gene_cases(ontology, config) if c.target.entity_id == "OMIM:1"]

    negatives = [n for case in cases for n in case.negative_phenotypes]
    assert negatives
    assert {n.simulated_origin for n in negatives} == {"negative_own_gene_other_disease"}
    assert {n.hpo_id for n in negatives} <= {"HP:0000964", "HP:0000028"}
    assert {n.reason for n in negatives} == {"gene_other_disease:OMIM:4"}
    for case in cases:
        present = _present_ids(case)
        for negative in case.negative_phenotypes:
            assert not any(_related(ontology, negative.hpo_id, p) for p in present)
        if case.patient.sex == "female":
            assert "HP:0000028" not in {n.hpo_id for n in case.negative_phenotypes}


def test_gene_other_disease_source_is_absent_in_disease_first_mode(ontology) -> None:
    config = _negatives_config(source_weights={"own_gene_other_disease": 1.0})
    index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)

    cases = simulate_cases(TRUE, config, ontology=ontology, confounders=index)

    assert not any(case.negative_phenotypes for case in cases)


def test_entity_pool_offers_the_equivalent_profiles_terms(ontology) -> None:
    only_own = {"own_disease": 1.0}
    entity = _gene_cases(ontology, _negatives_config(source_weights=only_own))
    profile = _gene_cases(
        ontology, _negatives_config(source_weights=only_own, own_disease_pool="profile")
    )

    def reasons(cases):
        return Counter(
            (n.hpo_id, n.reason)
            for case in cases
            if case.target.disease_id == "OMIM:1"
            for n in case.negative_phenotypes
        )

    equivalent = {
        hpo_id for hpo_id, reason in reasons(entity) if reason == EQUIVALENT_PROFILE_REASON
    }
    assert equivalent and equivalent <= {"HP:0000518", "HP:0001385"}
    assert all(reason != EQUIVALENT_PROFILE_REASON for _, reason in reasons(profile))


def test_redistribute_refills_slots_a_source_cannot_fill(ontology) -> None:
    easy = {**_config().presets["easy"].model_dump(), "negatives_mean": 6.0}

    def run(unfilled_slots: str):
        config = _config(
            difficulties=["easy"],
            cases_per_disease_per_difficulty=200,
            presets={"easy": easy},
            negatives={
                "source_weights": {"not_annotation": 1.0, "confounder": 1.0},
                "count_dispersion": None,
                "unfilled_slots": unfilled_slots,
            },
        )
        return _gene_cases(ontology, config)

    dropped, refilled = run("drop"), run("redistribute")

    def mean(cases):
        return sum(len(case.negative_phenotypes) for case in cases) / len(cases)

    assert mean(refilled) > mean(dropped) + 1.0
    origins = Counter(n.simulated_origin for case in refilled for n in case.negative_phenotypes)
    assert origins["negative_confounder"] > origins["negative_not_annotation"]
    for case in refilled:
        ids = [n.hpo_id for n in case.negative_phenotypes]
        assert len(ids) == len(set(ids))


def test_noise_is_drawn_by_weight(ontology) -> None:
    rng = random.Random(3)
    vocabulary = [
        NoiseTerm("HP:0000988", "Skin rash", weight=20.0),
        NoiseTerm("HP:0000505", "Visual impairment", weight=1.0),
        NoiseTerm("HP:0000964", "Eczema", weight=0.0),
    ]
    model = _disease_model(CONFOUNDER, _config(), None)
    patient = _Patient("female", "infantile", 0.5, 1.0, [], [])
    sex_terms = SexSpecificTerms(ontology, _config().sex)
    preset = _config().presets["hard"].model_copy(update={"noise_mean": 1.0})
    first: Counter[str] = Counter()
    for _ in range(400):
        noise = _sample_noise(
            model, patient, set(), [], preset, vocabulary, rng, ontology, sex_terms
        )
        if noise:
            first[noise[0].hpo_id] += 1

    assert first["HP:0000964"] == 0
    assert first["HP:0000988"] > 5 * first["HP:0000505"]


def test_gene_cases_pass_invariants_with_every_source(ontology) -> None:
    config = _negatives_config(
        source_weights={
            "own_disease": 1.0,
            "own_gene_other_disease": 1.0,
            "confounder": 1.0,
            "not_annotation": 1.0,
        },
        unfilled_slots="redistribute",
    )
    for case in _gene_cases(ontology, config):
        present = _present_ids(case)
        ids = [n.hpo_id for n in case.negative_phenotypes]
        assert len(ids) <= 15 and len(ids) == len(set(ids))
        for negative in ids:
            assert not any(_related(ontology, negative, p) for p in present)
        profile_terms = {p.hpo_id for p in PROFILE_MAP[case.target.disease_id].phenotypes}
        for negative in case.negative_phenotypes:
            if negative.simulated_origin == "negative_own_gene_other_disease":
                assert negative.hpo_id not in profile_terms


def test_json_dump_of_gene_case_has_entity_and_seed(ontology) -> None:
    record = json.loads(_gene_cases(ontology)[0].model_dump_json(exclude_none=True))

    assert record["target"]["entity_id"] in {"OMIM:1", "OMIM:4"}
    assert isinstance(record["metadata"]["case_seed"], int)
