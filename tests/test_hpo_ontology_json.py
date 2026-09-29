from rare_disease_simulator.data_sources.hpo import HpoOntology
from tests.fixtures.readers import fixture_path


def _ontology() -> HpoOntology:
    return HpoOntology.from_json(fixture_path("hpoa/hp_mini.json"))


def test_from_json_reads_version_labels_and_parents() -> None:
    ontology = _ontology()

    assert ontology.version == "2026-02-16"
    assert ontology.get_label("HP:0001250") == "Seizure"
    assert ontology.get_synonyms("HP:0001250") == ("Epileptic seizure",)
    assert ontology.get_direct_parents("HP:0001250") == ("HP:0000707",)
    assert ontology.get_ancestors("HP:0001250") == ("HP:0000707", "HP:0000118", "HP:0000001")


def test_from_json_flags_phenotypic_abnormality_branch() -> None:
    ontology = _ontology()

    assert ontology.is_phenotypic_abnormality("HP:0001250")
    assert not ontology.is_phenotypic_abnormality("HP:0000118")
    assert not ontology.is_phenotypic_abnormality("HP:0000007")
    assert not ontology.is_phenotypic_abnormality("HP:0003593")


def test_from_json_keeps_obsolete_terms_out_of_current_terms() -> None:
    ontology = _ontology()

    assert not ontology.is_valid_hpo_id("HP:0000999")
    assert "HP:0000999" in ontology.obsolete_terms
    assert ontology.obsolete_terms["HP:0000999"].replaced_by == "HP:0001251"
    assert all(not term.obsolete for term in ontology.terms.values())


def test_resolve_follows_alt_ids_and_replacements() -> None:
    ontology = _ontology()

    assert ontology.resolve("HP:0001250").status == "current"
    alt = ontology.resolve("HP:0001275")
    assert (alt.hpo_id, alt.status) == ("HP:0001250", "alt_id")
    replaced = ontology.resolve("HP:0000999")
    assert (replaced.hpo_id, replaced.status) == ("HP:0001251", "replaced")
    assert ontology.resolve("HP:0000489").hpo_id is None
    assert ontology.resolve("HP:1234567").status == "unresolved"


def test_ancestor_and_descendant_sets() -> None:
    ontology = _ontology()

    assert ontology.get_descendant_set("HP:0001417") == frozenset({"HP:0001419", "HP:0001423"})
    assert "HP:0003674" in ontology.get_ancestor_set("HP:0003593")
    assert ontology.is_a("HP:0003593", "HP:0410280")
    assert ontology.is_a("HP:0003593", "HP:0003593")
    assert not ontology.is_a("HP:0003581", "HP:0410280")
    assert ontology.get_direct_children("HP:0410280") == (
        "HP:0003593",
        "HP:0003621",
        "HP:0011463",
    )


def test_from_json_explicit_version_overrides_metadata() -> None:
    ontology = HpoOntology.from_json(fixture_path("hpoa/hp_mini.json"), version="custom")

    assert ontology.version == "custom"
