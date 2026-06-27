from rare_disease_simulator.data_sources.hpo_annotations import _frequency_category, parse_hpoa
from tests.fixtures.readers import fixture_path


def _annotations():
    return parse_hpoa(fixture_path("phenotype_mini.hpoa"), {"OMIM:257220"})


def test_parse_hpoa_filters_to_requested_disease() -> None:
    result = _annotations()

    assert set(result) == {"OMIM:257220"}
    annotations = result["OMIM:257220"]
    assert annotations.disease_name == "Niemann-Pick disease type C1"


def test_parse_hpoa_splits_positive_negative_inheritance_course() -> None:
    annotations = _annotations()["OMIM:257220"]

    positive_ids = {p.hpo_id for p in annotations.phenotypes}
    assert positive_ids == {"HP:0001251", "HP:0001433", "HP:0000511"}
    assert [p.hpo_id for p in annotations.negative_phenotypes] == ["HP:0001250"]
    assert annotations.inheritance_terms == ["HP:0000007"]
    assert annotations.clinical_course_terms == ["HP:0003593"]


def test_parse_hpoa_maps_frequency_and_onset() -> None:
    annotations = _annotations()["OMIM:257220"]
    by_id = {p.hpo_id: p for p in annotations.phenotypes}

    assert by_id["HP:0001251"].frequency_category == "frequent"  # HP:0040282
    assert by_id["HP:0001433"].frequency_category == "frequent"  # 12/20 = 0.6
    assert by_id["HP:0000511"].frequency_category == "occasional"  # 20%
    assert by_id["HP:0000511"].onset_category == "infantile"  # HP:0003593


def test_frequency_category_handles_terms_ratios_and_percentages() -> None:
    assert _frequency_category("HP:0040281") == "very_frequent"
    assert _frequency_category("9/10") == "very_frequent"
    assert _frequency_category("2%") == "very_rare"
    assert _frequency_category(None) is None
    assert _frequency_category("not-a-frequency") is None
