from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.validation.hpo import HpoValidator
from tests.fixtures.readers import fixture_path


def test_hpo_ontology_loads_fixture_terms() -> None:
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")

    assert ontology.version == "fixture-0.1"
    assert ontology.is_valid_hpo_id("HP:0001251")
    assert ontology.get_label("HP:0001251") == "Ataxia"
    assert ontology.get_label("HP:9999999") is None


def test_hpo_ontology_returns_direct_parents_and_ancestors() -> None:
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"))

    assert ontology.get_direct_parents("HP:0002066") == ("HP:0001251",)
    assert ontology.get_ancestors("HP:0002066") == (
        "HP:0001251",
        "HP:0001311",
        "HP:0000707",
        "HP:0000118",
    )


def test_hpo_validator_checks_phenotypic_abnormality_status() -> None:
    validator = HpoValidator(HpoOntology.from_tsv(fixture_path("hpo_terms.tsv")))

    assert validator.is_phenotypic_abnormality("HP:0001250")
    assert not validator.is_phenotypic_abnormality("HP:0000007")
    assert not validator.is_phenotypic_abnormality("HP:9999999")


def test_hpo_validator_reports_unknown_and_non_phenotypic_ids() -> None:
    validator = HpoValidator(HpoOntology.from_tsv(fixture_path("hpo_terms.tsv")))

    issues = validator.validate_phenotypic_hpo_ids(
        ["HP:0001250", "HP:0000007", "HP:9999999"]
    )

    assert [(issue.hpo_id, issue.reason) for issue in issues] == [
        ("HP:0000007", "not_phenotypic_abnormality"),
        ("HP:9999999", "unknown_hpo_id"),
    ]

