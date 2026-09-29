"""Phenotype frequency estimation from HPOA notations (ADR-0007).

The convention is shared with ``diagnostic.ar-training``
(``docs/decisions/0007-frequency-estimation.md``):

- counts ``n/m`` pool across duplicate rows and use the Jeffreys posterior mean
  ``(n + 0.5) / (m + 1)``; the range is the 95% Jeffreys interval;
- percentages are points; the range spans the reported percentages;
- HPO categories use the midpoint of their HPO-defined range; the range is the
  envelope of the categories' ranges;
- precedence: counts, then percentages, then categories;
- no frequency means unknown (``None``), never "always".
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache

HPO_FREQUENCY_RANGES: dict[str, tuple[float, float]] = {
    "HP:0040280": (1.0, 1.0),
    "HP:0040281": (0.80, 0.99),
    "HP:0040282": (0.30, 0.79),
    "HP:0040283": (0.05, 0.29),
    "HP:0040284": (0.01, 0.04),
}

HPO_FREQUENCY_MIDPOINTS: dict[str, float] = {
    "HP:0040280": 1.0,
    "HP:0040281": 0.895,
    "HP:0040282": 0.545,
    "HP:0040283": 0.17,
    "HP:0040284": 0.025,
}

HPO_FREQUENCY_CATEGORIES: dict[str, str] = {
    "HP:0040280": "obligate",
    "HP:0040281": "very_frequent",
    "HP:0040282": "frequent",
    "HP:0040283": "occasional",
    "HP:0040284": "very_rare",
}

JEFFREYS_INTERVAL_LEVEL = 0.95

_RATIO_RE = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
_PERCENT_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*%\s*$")
_ROUND_DIGITS = 4


@dataclass(frozen=True)
class FrequencyEstimate:
    """A phenotype frequency read from one or more HPOA rows."""

    category: str
    estimate: float
    lower: float
    upper: float
    raw: str
    basis: str


@dataclass(frozen=True)
class ParsedFrequencies:
    """HPOA frequency values of one (disease, term), split by notation."""

    ratios: tuple[tuple[int, int], ...] = ()
    percents: tuple[float, ...] = ()
    categories: tuple[str, ...] = ()
    unparsed: int = 0


def parse_frequencies(values: Iterable[str | None]) -> ParsedFrequencies:
    """Split raw HPOA frequency values by notation; empty values are skipped."""

    ratios: list[tuple[int, int]] = []
    percents: list[float] = []
    categories: list[str] = []
    unparsed = 0
    for value in values:
        if not value:
            continue
        if value in HPO_FREQUENCY_MIDPOINTS:
            categories.append(value)
            continue
        ratio = _RATIO_RE.match(value)
        if ratio and int(ratio.group(2)) > 0:
            numerator, denominator = int(ratio.group(1)), int(ratio.group(2))
            ratios.append((min(numerator, denominator), denominator))
            continue
        percent = _PERCENT_RE.match(value)
        if percent:
            percents.append(min(float(percent.group(1)) / 100.0, 1.0))
            continue
        unparsed += 1
    return ParsedFrequencies(tuple(ratios), tuple(percents), tuple(categories), unparsed)


def estimate_frequency(parsed: ParsedFrequencies) -> FrequencyEstimate | None:
    """Apply the ADR-0007 precedence; None means the frequency is unknown."""

    if parsed.ratios:
        numerator = sum(n for n, _ in parsed.ratios)
        denominator = sum(m for _, m in parsed.ratios)
        estimate = jeffreys_mean(numerator, denominator)
        lower, upper = jeffreys_interval(numerator, denominator)
        return FrequencyEstimate(
            category=category_for_probability(estimate),
            estimate=_round(estimate),
            lower=lower,
            upper=upper,
            raw=f"{numerator}/{denominator}",
            basis="counts",
        )
    if parsed.percents:
        estimate = sum(parsed.percents) / len(parsed.percents)
        return FrequencyEstimate(
            category=category_for_probability(estimate),
            estimate=_round(estimate),
            lower=_round(min(parsed.percents)),
            upper=_round(max(parsed.percents)),
            raw=";".join(f"{_round(p * 100):g}%" for p in sorted(set(parsed.percents))),
            basis="percent",
        )
    if parsed.categories:
        distinct = sorted(set(parsed.categories))
        estimate = sum(HPO_FREQUENCY_MIDPOINTS[term] for term in parsed.categories) / len(
            parsed.categories
        )
        category = (
            HPO_FREQUENCY_CATEGORIES[distinct[0]]
            if len(distinct) == 1
            else category_for_probability(estimate)
        )
        return FrequencyEstimate(
            category=category,
            estimate=_round(estimate),
            lower=min(HPO_FREQUENCY_RANGES[term][0] for term in distinct),
            upper=max(HPO_FREQUENCY_RANGES[term][1] for term in distinct),
            raw=";".join(distinct),
            basis="category",
        )
    return None


def category_for_probability(probability: float) -> str:
    """HPO frequency category whose range contains a positive phenotype's probability.

    A positive annotation is never ``excluded``: probabilities below the
    Very rare range still map to ``very_rare``.
    """

    if probability >= 0.99:
        return "obligate"
    if probability >= 0.80:
        return "very_frequent"
    if probability >= 0.30:
        return "frequent"
    if probability >= 0.05:
        return "occasional"
    return "very_rare"


def jeffreys_mean(successes: int, trials: int) -> float:
    """Posterior mean of a Binomial proportion under the Jeffreys Beta(1/2, 1/2) prior."""

    return (successes + 0.5) / (trials + 1.0)


@cache
def jeffreys_interval(
    successes: int, trials: int, level: float = JEFFREYS_INTERVAL_LEVEL
) -> tuple[float, float]:
    """Equal-tailed Jeffreys interval, with the usual 0/1 bounds at the extremes.

    Following Brown, Cai & DasGupta (2001), the lower bound is 0 when
    ``successes == 0`` and the upper bound is 1 when ``successes == trials``.
    """

    alpha = (1.0 - level) / 2.0
    a = successes + 0.5
    b = trials - successes + 0.5
    lower = 0.0 if successes == 0 else _beta_quantile(alpha, a, b)
    upper = 1.0 if successes == trials else _beta_quantile(1.0 - alpha, a, b)
    return _round(lower), _round(upper)


def _beta_quantile(probability: float, a: float, b: float) -> float:
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = (low + high) / 2.0
        if _regularized_incomplete_beta(middle, a, b) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _regularized_incomplete_beta(x: float, a: float, b: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    front = math.exp(log_front)
    # The continued fraction converges fast only below the mean; use symmetry above it.
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(x, a, b) / a
    return 1.0 - front * _beta_continued_fraction(1.0 - x, b, a) / b


def _beta_continued_fraction(x: float, a: float, b: float) -> float:
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    result = d
    for m in range(1, 300):
        m2 = 2 * m
        numerator = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + numerator * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + numerator / c
        c = c if abs(c) > tiny else tiny
        result *= d * c
        numerator = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + numerator * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + numerator / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return result


def _round(value: float) -> float:
    return round(value, _ROUND_DIGITS)
