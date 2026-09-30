import pytest

from rare_disease_simulator.profiles.frequency import beta_shrinkage_mean, pooled_counts
from rare_disease_simulator.profiles.merge import merge_entity_profiles
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiseaseGene,
    NegativePhenotypeAssociation,
)
from tests.test_simulation_v02 import _phenotype, _profile, build_ontology


def _counted(hpo_id, raw, estimate, **extra):
    return _phenotype(hpo_id, estimate, frequency_raw=raw, **extra)


def test_pooled_counts_reads_only_count_notation() -> None:
    assert pooled_counts("3/7") == (3, 7)
    assert pooled_counts("45%") is None
    assert pooled_counts("HP:0040281") is None
    assert pooled_counts(None) is None
    assert pooled_counts("1/0") is None


def test_beta_shrinkage_moves_small_counts_toward_the_prior() -> None:
    assert beta_shrinkage_mean(1, 1, 0.4, 2.0) == pytest.approx((1 + 0.8) / 3)
    assert beta_shrinkage_mean(0, 3, 0.4, 2.0) == pytest.approx(0.8 / 5)
    assert beta_shrinkage_mean(90, 100, 0.4, 2.0) == pytest.approx(90.8 / 102)


def _entity_profiles():
    omim = _profile(
        "OMIM:10",
        [
            _counted("HP:0001250", "2/10", 0.2273, source=["PMID:1"]),
            _counted("HP:0001251", "45%", 0.45),
            _counted("HP:0000047", "3/4", 0.7, sex_restriction="male"),
            _counted("HP:0000505", "1/2", 0.5, onset="adult", onset_hpo_id="HP:0003581"),
            _counted("HP:0002650", "4/4", 0.9, source=["PMID:7"]),
        ],
        negative_phenotypes=[
            NegativePhenotypeAssociation(hpo_id="HP:0001385", label="Hip dysplasia"),
            NegativePhenotypeAssociation(hpo_id="HP:0000518", label="Cataract"),
        ],
        age_of_onset=AgeOfOnset(category="infantile", distribution={"infantile": 1.0}),
    ).model_copy(
        update={"genes": [DiseaseGene(symbol="G1", inheritance_hpo_ids=["HP:0000006"])]}
    )
    orpha = _profile(
        "ORPHA:20",
        [
            _counted("HP:0001250", "HP:0040281", 0.895),
            _counted("HP:0001251", "HP:0040283", 0.17),
            _counted("HP:0000047", "HP:0040282", 0.545),
            _counted("HP:0000505", "HP:0040282", 0.545, onset="childhood",
                     onset_hpo_id="HP:0011463"),
            _counted("HP:0001627", "HP:0040280", 1.0),
            _counted("HP:0000518", "HP:0040283", 0.17),
            _counted("HP:0002650", "4/4", 0.9, source=["PMID:7"]),
        ],
        negative_phenotypes=[
            NegativePhenotypeAssociation(hpo_id="HP:0001385", label="Hip dysplasia")
        ],
        age_of_onset=AgeOfOnset(category="childhood", distribution={"childhood": 1.0}),
    ).model_copy(
        update={"genes": [DiseaseGene(symbol="G1", inheritance_hpo_ids=["HP:0000007"])]}
    )
    return omim, orpha


def test_merge_takes_the_union_and_applies_frequency_precedence() -> None:
    omim, orpha = _entity_profiles()

    merged, stats = merge_entity_profiles("OMIM:10", [omim, orpha], build_ontology())
    terms = {p.hpo_id: p for p in merged.phenotypes}

    assert merged.disease_id == "OMIM:10"
    assert merged.disease_name == omim.disease_name
    assert set(terms) == {
        "HP:0001250", "HP:0001251", "HP:0000047", "HP:0000505", "HP:0001627",
        "HP:0000518", "HP:0002650",
    }
    # Counts win over the Orphanet bin, which is capped (0.895 vs 0.23).
    assert terms["HP:0001250"].frequency_raw == "2/10"
    assert terms["HP:0001250"].frequency_estimate == pytest.approx(2.5 / 11, abs=1e-4)
    # Percentages win over categories.
    assert terms["HP:0001251"].frequency_raw == "45%"
    # A category only one source has keeps its value.
    assert terms["HP:0001627"].frequency_estimate == 1.0
    # The same count row repeated by both profiles is counted once.
    assert terms["HP:0002650"].frequency_raw == "4/4"
    assert stats.counters["count_rows_deduplicated"] == 1
    assert stats.counters["category_capped"] == 1
    assert stats.counters["terms_added_over_largest_profile"] == 0


def test_merge_is_conservative_on_sex_onset_and_negatives() -> None:
    omim, orpha = _entity_profiles()

    merged, stats = merge_entity_profiles("OMIM:10", [omim, orpha], None)
    terms = {p.hpo_id: p for p in merged.phenotypes}

    # One profile restricts the term to males and the other does not: no restriction.
    assert terms["HP:0000047"].sex_restriction is None
    assert stats.counters["sex_restriction_dropped"] == 1
    # The earliest onset any profile gives the term.
    assert (terms["HP:0000505"].onset, terms["HP:0000505"].onset_hpo_id) == (
        "childhood",
        "HP:0011463",
    )
    # A NOT term another profile annotates as present is dropped.
    assert [n.hpo_id for n in merged.negative_phenotypes] == ["HP:0001385"]
    assert merged.age_of_onset.distribution == {"infantile": 0.5, "childhood": 0.5}
    assert merged.age_of_onset.category == "infantile"
    assert merged.genes[0].inheritance_hpo_ids == ["HP:0000006", "HP:0000007"]
    assert merged.sex_bias is not None and merged.sex_bias.value == "none"
    assert merged.mapped_ids.orpha == "ORPHA:20" and merged.mapped_ids.omim == ["OMIM:10"]


def test_merging_one_profile_keeps_it() -> None:
    omim, _ = _entity_profiles()

    merged, _ = merge_entity_profiles("OMIM:10", [omim], build_ontology())

    assert merged.phenotypes == omim.phenotypes
    assert merged.age_of_onset == omim.age_of_onset


def test_merge_is_deterministic_and_order_insensitive_for_terms() -> None:
    omim, orpha = _entity_profiles()

    first, _ = merge_entity_profiles("OMIM:10", [omim, orpha], build_ontology())
    again, _ = merge_entity_profiles("OMIM:10", [omim, orpha], build_ontology())

    assert first.model_dump_json() == again.model_dump_json()
    assert [p.hpo_id for p in first.phenotypes] == sorted(p.hpo_id for p in first.phenotypes)
