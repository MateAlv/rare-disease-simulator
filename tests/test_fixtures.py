from rare_disease_simulator.llm_extraction.extractor import DummyLlmProvider
from tests.fixtures.readers import (
    fixture_path,
    read_disease_profile,
    read_hpo_terms,
    read_orphadata_mini,
    read_profile_patch,
    read_synthetic_case,
    read_text_snippets,
)


def test_fixture_readers_load_minimal_vertical_slice_inputs() -> None:
    hpo_terms = read_hpo_terms()
    orphadata = read_orphadata_mini()
    snippets = read_text_snippets()
    patch = read_profile_patch()
    profile = read_disease_profile()
    synthetic_case = read_synthetic_case()

    assert {term["hpo_id"] for term in hpo_terms} >= {"HP:0001251", "HP:0000511"}
    assert orphadata["diseases"][0]["orpha_id"] == "ORPHA:646"
    assert snippets[0].disease_id == "ORPHA:646"
    assert patch.gene == "NPC1"
    assert profile.genes[0].symbol == "NPC1"
    assert synthetic_case.target.gene == "NPC1"


def test_dummy_llm_provider_returns_fixed_patch_for_matching_snippet() -> None:
    snippet = read_text_snippets()[0]
    provider = DummyLlmProvider.from_patch_file(fixture_path("profile_patch.json"))

    patch = provider.extract(snippet)

    assert patch.disease_id == snippet.disease_id
    assert patch.gene == snippet.gene
    assert [phenotype.hpo_id for phenotype in patch.phenotypes] == [
        "HP:0001251",
        "HP:0000511",
        "HP:0001250",
    ]
    assert patch.inheritance[0].value == "autosomal recessive"


def test_synthetic_case_fixture_represents_negative_missing_and_unknown_distinctly() -> None:
    synthetic_case = read_synthetic_case()

    assert synthetic_case.negative_phenotypes[0].status == "negative"
    assert synthetic_case.missing_phenotypes[0].status == "missing"
    assert synthetic_case.unknown_phenotypes[0].status == "unknown"
    assert synthetic_case.target.disease_id == "ORPHA:646"
