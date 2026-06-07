from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.llm_extraction.validator import DiseaseProfilePatchValidator
from rare_disease_simulator.validation.hpo import HpoValidator
from tests.fixtures.readers import fixture_path, read_profile_patch, read_text_snippets


def _validator() -> DiseaseProfilePatchValidator:
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")
    return DiseaseProfilePatchValidator(HpoValidator(ontology))


def test_profile_patch_validator_accepts_fixture_patch() -> None:
    patch = read_profile_patch()
    source_text = read_text_snippets()[0].text

    result = _validator().validate(patch, source_text=source_text)

    assert result.status == "accepted"
    assert result.issues == []


def test_profile_patch_validator_marks_low_confidence_claim_for_review() -> None:
    patch = read_profile_patch()
    source_text = read_text_snippets()[0].text
    phenotype = patch.phenotypes[0].model_copy(update={"confidence": 0.60})
    patch = patch.model_copy(update={"phenotypes": [phenotype, *patch.phenotypes[1:]]})

    result = _validator().validate(patch, source_text=source_text)

    assert result.status == "needs_review"
    assert [(issue.severity, issue.reason) for issue in result.issues] == [
        ("review", "confidence_requires_review")
    ]


def test_profile_patch_validator_rejects_unknown_hpo_id() -> None:
    patch = read_profile_patch()
    source_text = read_text_snippets()[0].text
    phenotype = patch.phenotypes[0].model_copy(update={"hpo_id": "HP:9999999"})
    patch = patch.model_copy(update={"phenotypes": [phenotype, *patch.phenotypes[1:]]})

    result = _validator().validate(patch, source_text=source_text)

    assert result.status == "rejected"
    assert ("error", "unknown_hpo_id") in [
        (issue.severity, issue.reason) for issue in result.issues
    ]


def test_profile_patch_validator_rejects_evidence_not_found_in_source_text() -> None:
    patch = read_profile_patch()
    source_text = read_text_snippets()[0].text
    phenotype = patch.phenotypes[0].model_copy(update={"evidence_span": "not in source"})
    patch = patch.model_copy(update={"phenotypes": [phenotype, *patch.phenotypes[1:]]})

    result = _validator().validate(patch, source_text=source_text)

    assert result.status == "rejected"
    assert ("error", "evidence_span_not_found_in_source_text") in [
        (issue.severity, issue.reason) for issue in result.issues
    ]

