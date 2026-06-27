from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.data_sources.hpo_match import HpoMatcher, normalize_mention
from tests.fixtures.readers import fixture_path


def _matcher() -> HpoMatcher:
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")
    return HpoMatcher(ontology)


def test_normalize_mention_strips_punctuation_and_case() -> None:
    assert normalize_mention("  Vertical-Gaze   Palsy! ") == "vertical gaze palsy"


def test_exact_label_match() -> None:
    match = _matcher().match("ataxia")

    assert match is not None
    assert match.hpo_id == "HP:0001251"
    assert match.method == "exact"
    assert match.score == 1.0


def test_synonym_match_resolves_to_canonical_id() -> None:
    match = _matcher().match("vertical gaze palsy")

    assert match is not None
    assert match.hpo_id == "HP:0000511"
    assert match.label == "Vertical supranuclear gaze palsy"


def test_plural_synonym_match() -> None:
    match = _matcher().match("Seizures")

    assert match is not None
    assert match.hpo_id == "HP:0001250"


def test_fuzzy_match_above_threshold() -> None:
    match = _matcher().match("hepatosplenomegally")

    assert match is not None
    assert match.hpo_id == "HP:0001433"
    assert match.method == "fuzzy"


def test_unmatchable_mention_returns_none() -> None:
    assert _matcher().match("polydactyly of the left hand") is None
