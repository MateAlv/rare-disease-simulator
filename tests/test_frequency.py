import pytest

from rare_disease_simulator.profiles.frequency import (
    category_for_probability,
    estimate_frequency,
    jeffreys_interval,
    parse_frequencies,
)


def _estimate(*values: str | None):
    return estimate_frequency(parse_frequencies(values))


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        (("1/1",), 0.75),
        (("10/10",), 0.9545),
        (("0/3",), 0.125),
        (("1/4", "3/4"), 0.5),
    ],
)
def test_adr_0007_count_vectors_use_the_jeffreys_mean(values, expected) -> None:
    estimate = _estimate(*values)

    assert estimate.basis == "counts"
    assert estimate.estimate == expected
    assert estimate.lower <= estimate.estimate <= estimate.upper


def test_adr_0007_precedence_counts_then_percent_then_category() -> None:
    assert _estimate("HP:0040281", "30%", "1/1").basis == "counts"
    assert _estimate("HP:0040281", "30%", "1/1").estimate == 0.75

    percent = _estimate("HP:0040281", "30%", "50%")
    assert (percent.basis, percent.estimate, percent.lower, percent.upper) == (
        "percent",
        0.4,
        0.3,
        0.5,
    )

    category = _estimate("HP:0040281", "HP:0040283")
    assert (category.basis, category.estimate) == ("category", 0.5325)


@pytest.mark.parametrize(
    ("term", "midpoint"),
    [
        ("HP:0040280", 1.0),
        ("HP:0040281", 0.895),
        ("HP:0040282", 0.545),
        ("HP:0040283", 0.17),
        ("HP:0040284", 0.025),
    ],
)
def test_categories_use_hpo_range_midpoints(term, midpoint) -> None:
    estimate = _estimate(term)

    assert estimate.estimate == midpoint
    assert estimate.lower <= midpoint <= estimate.upper


def test_empty_frequency_is_unknown_never_always() -> None:
    assert _estimate(None, "") is None


def test_unparsed_values_are_counted_and_ignored() -> None:
    parsed = parse_frequencies(["often", "2/0", "3/4"])

    assert parsed.unparsed == 2
    assert parsed.ratios == ((3, 4),)


def test_counts_above_the_denominator_are_capped() -> None:
    assert parse_frequencies(["5/3"]).ratios == ((3, 3),)


def test_jeffreys_interval_matches_reference_values() -> None:
    assert jeffreys_interval(0, 10) == (0.0, 0.2172)
    assert jeffreys_interval(5, 10) == (0.2235, 0.7765)
    assert jeffreys_interval(10, 10)[1] == 1.0


def test_positive_probabilities_are_never_excluded() -> None:
    assert category_for_probability(0.0) == "very_rare"
    assert category_for_probability(0.75) == "frequent"
    assert category_for_probability(0.99) == "obligate"
