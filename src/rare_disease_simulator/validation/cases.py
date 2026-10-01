"""Validation report over simulated cases (``rich_cases.jsonl``).

It reports the shape of the dataset (phenotype counts, negative sources,
missingness, noise, sex and age distributions) against the priors that
generated it, calibrates observed phenotype frequencies against the profiles'
estimates, and lists invariant violations. Any violation makes the report
fail.

"Truly present" in a case means recorded as positive (or generalized from),
``missing`` or ``unknown``. Calibration compares, per (disease, profile
term), the share of that disease's cases where the term is truly present
with its frequency estimate. The ``ungated`` view keeps only terms that no
mechanism other than frequency acts on (no sex restriction or sex-specific
anchor, no onset gating, not in a progressive disease), so it should match
the estimate up to sampling noise and the at-least-one-positive redraw.

Merged gene-first cases (``target.profile_ids``) are checked against the
profile merged from those ids, as the simulator built it. Report-model cases
(``metadata.report_budget``) are also checked for the reporting rule: the
reported profile-term count is ``min(true terms, max(profile budget, true
cardinal terms))``, with the profile budget ``max(1, k - noise drawn)``, every
true cardinal term is reported when the budget is at least 1, and ``k`` lies
in the report model's histogram; the report compares the drawn budgets and
the reported counts, profile terms plus noise, with that histogram.
Independent-mode cases (a report probability ``q`` on every true profile term)
have no budget; the report bins the true terms by ``q`` and compares each
bin's mean ``q`` with the share of its terms that were reported.
"""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.merge import merge_entity_profiles
from rare_disease_simulator.profiles.schema import DiseaseProfile, PhenotypeAssociation
from rare_disease_simulator.simulation.reporting import BudgetSampler, CardinalIndex, ReportModel
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase
from rare_disease_simulator.simulation.simulator import (
    CONGENITAL_ONSET,
    FORCED_REASON,
    MERGED_REASON,
    RELATED_NOISE_REASON,
    SPECIALIZED_REASON,
    SexSpecificTerms,
    ShrinkageMeanMissing,
    config_hash,
    sex_prior_key,
    simulation_frequency,
)

BUCKETS = (
    "positive_phenotypes",
    "negative_phenotypes",
    "missing_phenotypes",
    "unknown_phenotypes",
    "noise_phenotypes",
)
CALIBRATION_BINS: tuple[tuple[float, float], ...] = (
    (0.0, 0.05),
    (0.05, 0.30),
    (0.30, 0.80),
    (0.80, 0.99),
    (0.99, 1.0001),
)
AGE_BINS_YEARS: tuple[tuple[float, float, str], ...] = (
    (0.0, 0.0767, "<28d"),
    (0.0767, 1.0, "28d-1y"),
    (1.0, 5.0, "1-5y"),
    (5.0, 16.0, "5-16y"),
    (16.0, 40.0, "16-40y"),
    (40.0, 60.0, "40-60y"),
    (60.0, float("inf"), ">=60y"),
)
MAX_EXAMPLES = 20


@dataclass
class _Accumulator:
    cases: int = 0
    diseases: Counter[str] = field(default_factory=Counter)
    genes: Counter[str] = field(default_factory=Counter)
    difficulties: Counter[str] = field(default_factory=Counter)
    config_hashes: Counter[str] = field(default_factory=Counter)
    bucket_counts: dict[str, list[int]] = field(default_factory=lambda: defaultdict(list))
    negative_origins: Counter[str] = field(default_factory=Counter)
    positives_generalized: int = 0
    positives_specialized: int = 0
    positives_profile: int = 0
    noise_related: int = 0
    forced_cases: int = 0
    observed_by_difficulty: dict[str, Counter[str]] = field(
        default_factory=lambda: defaultdict(Counter)
    )
    sex: Counter[str] = field(default_factory=Counter)
    sex_by_prior: dict[str, Counter[str]] = field(default_factory=lambda: defaultdict(Counter))
    ages: list[float] = field(default_factory=list)
    onsets: list[float] = field(default_factory=list)
    onset_categories: Counter[str] = field(default_factory=Counter)
    expected_onset: Counter[str] = field(default_factory=Counter)
    present_counts: dict[str, Counter[str]] = field(default_factory=lambda: defaultdict(Counter))
    case_profiles: dict[str, DiseaseProfile] = field(default_factory=dict)
    budgets: Counter[int] = field(default_factory=Counter)
    reported: Counter[int] = field(default_factory=Counter)
    reported_total: Counter[int] = field(default_factory=Counter)
    q_rolls: list[tuple[float, bool]] = field(default_factory=list)
    presentation_age_cases: int = 0
    independent: Counter[str] = field(default_factory=Counter)
    cardinal: Counter[str] = field(default_factory=Counter)
    violations: Counter[str] = field(default_factory=Counter)
    examples: list[dict[str, str]] = field(default_factory=list)

    def violation(self, kind: str, case_id: str, detail: str) -> None:
        self.violations[kind] += 1
        if len(self.examples) < MAX_EXAMPLES:
            self.examples.append({"type": kind, "case_id": case_id, "detail": detail})


def validate_cases(
    cases: Iterable[SyntheticCase],
    *,
    profiles: Mapping[str, DiseaseProfile] | None = None,
    ontology: HpoOntology | None = None,
    config: SimulationConfig | None = None,
    noise_vocabulary: set[str] | None = None,
    report_model: ReportModel | None = None,
    cardinal: CardinalIndex | None = None,
) -> dict[str, Any]:
    """Build the validation report; ``report["violations"]["count"]`` > 0 means failure."""

    acc = _Accumulator()
    sex_terms = SexSpecificTerms(ontology, config.sex) if config is not None else None
    closures: dict[str, frozenset[str]] = {}
    merged: dict[tuple[str, ...], DiseaseProfile | None] = {}
    budget = BudgetSampler(report_model.term_budget.histogram) if report_model else None
    force_cardinal = config.reporting.force_cardinal if config is not None else True
    for case in cases:
        profile = _case_profile(case, profiles, ontology, merged)
        if profiles is not None and profile is None:
            acc.violation("unknown_disease", case.case_id, _profile_key(case))
        if profile is not None:
            acc.case_profiles.setdefault(_profile_key(case), profile)
        _count_case(acc, case, profile, config)
        _check_case(acc, case, profile, ontology, config, sex_terms, noise_vocabulary, closures)
        if case.metadata.report_budget is not None:
            _check_reporting(acc, case, budget, cardinal, force_cardinal)
        elif config is not None and config.reporting.mode == "independent":
            _collect_report_rolls(acc, case)

    return {
        "cases": acc.cases,
        "diseases": len(acc.diseases),
        "genes": len(acc.genes),
        "difficulties": dict(sorted(acc.difficulties.items())),
        "config_hashes": dict(sorted(acc.config_hashes.items())),
        "config_matches_cases": (
            set(acc.config_hashes) == {config_hash(config)} if config is not None else None
        ),
        "phenotypes_per_case": {
            bucket.removesuffix("_phenotypes"): _distribution(acc.bucket_counts[bucket])
            for bucket in BUCKETS
        },
        "positives": {
            "generalized": acc.positives_generalized,
            "specialized": acc.positives_specialized,
            "specialized_share": _share(acc.positives_specialized, acc.positives_profile),
            "cases_with_forced_positive": acc.forced_cases,
            "observation_rate_by_difficulty": _observation_rates(acc, config),
        },
        "negatives": _negative_report(acc),
        "missingness": _missingness_report(acc),
        "noise": {
            "per_case_mean": _mean(acc.bucket_counts["noise_phenotypes"]),
            "cases_with_noise_share": _share(
                sum(1 for count in acc.bucket_counts["noise_phenotypes"] if count), acc.cases
            ),
            "related": acc.noise_related,
            "related_share": _share(
                acc.noise_related, sum(acc.bucket_counts["noise_phenotypes"])
            ),
        },
        "sex": _sex_report(acc, config),
        "age": {
            **_age_report(acc.ages, acc.cases),
            "presentation_age_cases": acc.presentation_age_cases,
        },
        "onset": {
            **_age_report(acc.onsets, acc.cases),
            "by_category": _category_report(acc.onset_categories, acc.expected_onset),
        },
        "reporting": (
            _independent_report(acc)
            if acc.independent["cases"]
            else _reporting_report(acc, budget, cardinal)
        ),
        "calibration": (
            _calibration(acc, ontology, config) if profiles is not None else None
        ),
        "violations": {
            "count": sum(acc.violations.values()),
            "by_type": dict(sorted(acc.violations.items())),
            "examples": acc.examples,
        },
    }


def _profile_key(case: SyntheticCase) -> str:
    if case.target.profile_ids:
        return "+".join(case.target.profile_ids)
    return case.target.disease_id


def _case_profile(
    case: SyntheticCase,
    profiles: Mapping[str, DiseaseProfile] | None,
    ontology: HpoOntology | None,
    merged: dict[tuple[str, ...], DiseaseProfile | None],
) -> DiseaseProfile | None:
    if profiles is None:
        return None
    ids = case.target.profile_ids
    if not ids:
        return profiles.get(case.target.disease_id)
    key = tuple(ids)
    if key not in merged:
        members = [profiles.get(profile_id) for profile_id in ids]
        merged[key] = (
            merge_entity_profiles(
                case.target.entity_id or case.target.disease_id,
                [member for member in members if member is not None],
                ontology,
            )[0]
            if all(member is not None for member in members)
            else None
        )
    return merged[key]


def _check_reporting(
    acc: _Accumulator,
    case: SyntheticCase,
    budget: BudgetSampler | None,
    cardinal: CardinalIndex | None,
    force_cardinal: bool,
) -> None:
    """The reporting rule of report-model cases (see the module docstring)."""

    k = case.metadata.report_budget
    assert k is not None
    acc.budgets[k] += 1
    if budget is not None and k not in budget.probabilities:
        acc.violation("budget_outside_histogram", case.case_id, str(k))
    profile_budget = case.metadata.report_budget_profile
    if profile_budget is None:
        profile_budget = k
    elif not 1 <= profile_budget <= max(1, k - len(case.noise_phenotypes)):
        # max(1, k - noise drawn); fewer noise terms are shown when the pool runs dry.
        acc.violation("profile_budget_mismatch", case.case_id, f"{profile_budget} for k={k}")
    profile_positives = [p for p in case.positive_phenotypes if p.simulated_origin != "noise"]
    reported_ids = {p.source_hpo_id or p.hpo_id for p in profile_positives}
    reported_ids |= {p.hpo_id for p in case.missing_phenotypes if p.reason == MERGED_REASON}
    unreported = {p.hpo_id for p in case.missing_phenotypes if p.reason != MERGED_REASON}
    unreported |= {p.hpo_id for p in case.unknown_phenotypes}
    reported = len(profile_positives) + sum(
        1 for p in case.missing_phenotypes if p.reason == MERGED_REASON
    )
    acc.reported[reported] += 1
    acc.reported_total[reported + len(case.noise_phenotypes)] += 1
    true_count = reported + len(unreported - reported_ids)
    if cardinal is None:
        if reported < min(profile_budget, true_count):
            acc.violation(
                "reported_below_budget", case.case_id, f"{reported}<{profile_budget}"
            )
        return
    ids = case.target.profile_ids or [case.target.disease_id]
    true_cardinal = (reported_ids | unreported) & cardinal.terms_for(ids)
    acc.cardinal["cases"] += 1
    acc.cardinal["true_terms"] += len(true_cardinal)
    acc.cardinal["reported_terms"] += len(true_cardinal & reported_ids)
    if true_cardinal:
        acc.cardinal["cases_with_true_cardinal"] += 1
    forced = len(true_cardinal) if force_cardinal and profile_budget > 0 else 0
    expected = min(true_count, max(profile_budget, forced))
    if reported != expected:
        acc.violation("reported_count_mismatch", case.case_id, f"{reported}!={expected}")
    if forced:
        missed = sorted(true_cardinal & (unreported - reported_ids))
        if missed:
            acc.violation("cardinal_not_reported", case.case_id, ",".join(missed))


def _collect_report_rolls(acc: _Accumulator, case: SyntheticCase) -> None:
    """Each true profile term's report probability and whether its roll reported it.

    A ``forced_min_one`` term was shown after every roll of its case failed,
    so it counts as a failed roll.
    """

    acc.independent["cases"] += 1
    rolls: list[tuple[float | None, bool]] = []
    for p in case.positive_phenotypes:
        if p.simulated_origin == "disease_profile":
            rolls.append((p.report_probability, p.reason != FORCED_REASON))
            acc.independent["reported"] += 1
            acc.independent["forced_min_one"] += p.reason == FORCED_REASON
    for p in case.missing_phenotypes:
        merged = p.reason == MERGED_REASON
        rolls.append((p.report_probability, merged))
        acc.independent["reported"] += merged
    rolls.extend((p.report_probability, False) for p in case.unknown_phenotypes)
    acc.independent["true_terms"] += len(rolls)
    acc.independent["noise"] += len(case.noise_phenotypes)
    acc.q_rolls.extend((q, reported) for q, reported in rolls if q is not None)


def _independent_report(acc: _Accumulator) -> dict[str, Any]:
    cases = acc.independent["cases"]
    bins = []
    for decile in range(10):
        low, high = decile / 10, (decile + 1) / 10
        members = [
            (q, reported)
            for q, reported in acc.q_rolls
            if low <= q < high or (decile == 9 and q == 1.0)
        ]
        if not members:
            continue
        bins.append(
            {
                "q_range": [low, high],
                "terms": len(members),
                "mean_q": round(sum(q for q, _ in members) / len(members), 4),
                "reported_share": round(sum(r for _, r in members) / len(members), 4),
            }
        )
    rolls = len(acc.q_rolls)
    error = sum(abs(b["mean_q"] - b["reported_share"]) * b["terms"] for b in bins)
    return {
        "mode": "independent",
        "cases": cases,
        "q_calibration": bins,
        "q_calibration_error": round(error / rolls, 4) if rolls else None,
        "reported_per_case_mean": round(acc.independent["reported"] / cases, 3),
        "noise_per_case_mean": round(acc.independent["noise"] / cases, 3),
        "true_terms_reported_share": _share(
            acc.independent["reported"], acc.independent["true_terms"]
        ),
        "q_mean": round(sum(q for q, _ in acc.q_rolls) / rolls, 4) if rolls else None,
        "q_capped_share": _share(sum(1 for q, _ in acc.q_rolls if q >= 1.0), rolls),
        "forced_min_one": acc.independent["forced_min_one"],
    }


def _reporting_report(
    acc: _Accumulator, budget: BudgetSampler | None, cardinal: CardinalIndex | None
) -> dict[str, Any] | None:
    cases = sum(acc.budgets.values())
    if not cases:
        return None
    drawn = {k: count / cases for k, count in acc.budgets.items()}
    reported = {k: count / cases for k, count in acc.reported.items()}
    total = {k: count / cases for k, count in acc.reported_total.items()}
    expected = budget.probabilities if budget is not None else None
    return {
        "cases": cases,
        "budget_expected": _share_table(expected) if expected is not None else None,
        "budget_drawn": _share_table(drawn),
        "reported_per_case": _share_table(reported),
        "reported_total_per_case": _share_table(total),
        "reported_mean": round(sum(k * c for k, c in acc.reported.items()) / cases, 3),
        "reported_total_mean": round(
            sum(k * c for k, c in acc.reported_total.items()) / cases, 3
        ),
        "budget_mean": round(sum(k * c for k, c in acc.budgets.items()) / cases, 3),
        "total_variation_drawn_vs_expected": (
            _total_variation(drawn, expected) if expected is not None else None
        ),
        "total_variation_reported_vs_expected": (
            _total_variation(reported, expected) if expected is not None else None
        ),
        "total_variation_reported_total_vs_expected": (
            _total_variation(total, expected) if expected is not None else None
        ),
        "cardinal": (
            {
                "cases_with_true_cardinal": acc.cardinal["cases_with_true_cardinal"],
                "true_terms": acc.cardinal["true_terms"],
                "reported_terms": acc.cardinal["reported_terms"],
                "reported_share": _share(
                    acc.cardinal["reported_terms"], acc.cardinal["true_terms"]
                ),
            }
            if cardinal is not None
            else None
        ),
    }


def _share_table(shares: Mapping[int, float]) -> dict[str, float]:
    return {str(k): round(shares[k], 4) for k in sorted(shares)}


def _total_variation(a: Mapping[int, float], b: Mapping[int, float]) -> float:
    return round(0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b)), 4)


def _count_case(
    acc: _Accumulator,
    case: SyntheticCase,
    profile: DiseaseProfile | None,
    config: SimulationConfig | None,
) -> None:
    acc.cases += 1
    acc.diseases[_profile_key(case)] += 1
    acc.genes[case.target.gene] += 1
    acc.difficulties[case.metadata.difficulty] += 1
    acc.config_hashes[case.metadata.config_hash] += 1
    for bucket in BUCKETS:
        acc.bucket_counts[bucket].append(len(getattr(case, bucket)))
    acc.negative_origins.update(p.simulated_origin for p in case.negative_phenotypes)
    acc.positives_generalized += sum(
        1 for p in case.positive_phenotypes
        if p.source_hpo_id and p.reason != SPECIALIZED_REASON
    )
    acc.positives_specialized += sum(
        1 for p in case.positive_phenotypes if p.reason == SPECIALIZED_REASON
    )
    acc.positives_profile += len(case.positive_phenotypes)
    acc.noise_related += sum(
        1 for p in case.noise_phenotypes if p.reason == RELATED_NOISE_REASON
    )
    if any(p.reason == FORCED_REASON for p in case.positive_phenotypes):
        acc.forced_cases += 1
    observed = acc.observed_by_difficulty[case.metadata.difficulty]
    observed["positive"] += len(case.positive_phenotypes)
    observed["unobserved"] += len(case.missing_phenotypes) + len(case.unknown_phenotypes)

    patient = case.patient
    acc.sex[patient.sex] += 1
    prior_key = case.metadata.sex_prior_key or (
        sex_prior_key(profile) if profile is not None else None
    )
    if prior_key is not None:
        acc.sex_by_prior[prior_key][patient.sex] += 1
    if patient.age is not None:
        acc.ages.append(patient.age.value)
    if patient.age_of_onset is not None:
        acc.onsets.append(patient.age_of_onset.value)
    if patient.onset_category != "unknown":
        acc.onset_categories[patient.onset_category] += 1
        if profile is not None and config is not None:
            for category, weight in _onset_prior(profile, config).items():
                acc.expected_onset[category] += weight

    for hpo_id in _truly_present(case):
        acc.present_counts[_profile_key(case)][hpo_id] += 1


def _check_case(
    acc: _Accumulator,
    case: SyntheticCase,
    profile: DiseaseProfile | None,
    ontology: HpoOntology | None,
    config: SimulationConfig | None,
    sex_terms: SexSpecificTerms | None,
    noise_vocabulary: set[str] | None,
    closures: dict[str, frozenset[str]],
) -> None:
    case_id = case.case_id
    if not case.positive_phenotypes:
        acc.violation("no_observed_positive", case_id, "")

    negatives = [p.hpo_id for p in case.negative_phenotypes]
    if len(negatives) != len(set(negatives)):
        acc.violation("duplicate_negative", case_id, ",".join(negatives))
    if config is not None and len(negatives) > config.negatives.max_per_case:
        acc.violation("negatives_above_cap", case_id, str(len(negatives)))

    present = _truly_present(case) | {p.hpo_id for p in case.noise_phenotypes}
    for negative in case.negative_phenotypes:
        related = next((p for p in sorted(present) if _related(negative.hpo_id, p, ontology)), None)
        if related is not None:
            acc.violation("negative_related_to_present", case_id, f"{negative.hpo_id}~{related}")
        if negative.simulated_origin == "negative_confounder" and profile is not None:
            closure = closures.get(profile.disease_id)
            if closure is None:
                closure = _annotation_closure(profile, ontology)
                closures[profile.disease_id] = closure
            if negative.hpo_id in closure:
                acc.violation("confounder_negative_annotated_to_disease", case_id, negative.hpo_id)
        if (
            negative.simulated_origin == "negative_own_gene_other_disease"
            and profile is not None
            and any(p.hpo_id == negative.hpo_id for p in profile.phenotypes)
        ):
            acc.violation("gene_other_negative_annotated_to_profile", case_id, negative.hpo_id)

    patient = case.patient
    age_range = case.metadata.presentation_age_years
    if age_range is not None:
        acc.presentation_age_cases += 1
        if patient.age is not None:
            _check_presentation_age(acc, case, age_range, config)
    if (
        patient.age is not None
        and patient.age_of_onset is not None
        and patient.age.value < patient.age_of_onset.value
    ):
        acc.violation(
            "age_before_onset", case_id, f"{patient.age.value}<{patient.age_of_onset.value}"
        )

    if patient.sex in {"male", "female"} and profile is not None:
        restrictions = {p.hpo_id: p.sex_restriction for p in profile.phenotypes}
        for hpo_id in sorted(present | set(negatives)):
            restriction = restrictions.get(hpo_id)
            allowed = (
                sex_terms.allowed(hpo_id, restriction, patient.sex)
                if sex_terms is not None
                else restriction in (None, patient.sex)
            )
            if not allowed:
                acc.violation("sex_restricted_term_for_other_sex", case_id, hpo_id)

    if noise_vocabulary is not None:
        for noise in case.noise_phenotypes:
            if noise.reason != RELATED_NOISE_REASON and noise.hpo_id not in noise_vocabulary:
                acc.violation("noise_not_in_vocabulary", case_id, noise.hpo_id)

    if profile is not None:
        profile_ids = {p.hpo_id for p in profile.phenotypes}
        for positive in case.positive_phenotypes:
            if positive.reason != SPECIALIZED_REASON:
                continue
            source = positive.source_hpo_id
            descends = source is not None and (
                ontology.is_a(positive.hpo_id, source) if ontology is not None else True
            )
            if source not in profile_ids or positive.hpo_id == source or not descends:
                acc.violation(
                    "specialized_not_descendant_of_true_term", case_id,
                    f"{positive.hpo_id}<-{source}",
                )
        related = [n for n in case.noise_phenotypes if n.reason == RELATED_NOISE_REASON]
        if related:
            key = "annotated:" + profile.disease_id
            annotated = closures.get(key)
            if annotated is None:
                annotated = frozenset(
                    profile_ids.union(
                        *(ontology.get_ancestor_set(i) for i in profile_ids)
                    )
                    if ontology is not None
                    else profile_ids
                )
                closures[key] = annotated
            for noise in related:
                below = (
                    ontology is not None
                    and not ontology.get_ancestor_set(noise.hpo_id).isdisjoint(profile_ids)
                )
                if noise.hpo_id in annotated or below:
                    acc.violation("related_noise_annotated_to_entity", case_id, noise.hpo_id)


AGE_ROUNDING = 0.006


def _check_presentation_age(
    acc: _Accumulator,
    case: SyntheticCase,
    age_range: tuple[float, float],
    config: SimulationConfig | None,
) -> None:
    """The age is uniform in [low, high], raised to the onset and cut at max_age_years."""

    assert case.patient.age is not None
    age = case.patient.age.value
    low, high = age_range
    max_age = config.age.max_age_years if config is not None else float("inf")
    lowest, highest = min(low, max_age), min(high, max_age)
    onset = case.patient.age_of_onset
    if onset is not None:
        lowest, highest = max(onset.value, lowest), max(onset.value, highest)
        too_high = age > highest + AGE_ROUNDING
    else:
        too_high = age > max(highest, max_age) + AGE_ROUNDING
    if age < lowest - AGE_ROUNDING or too_high:
        acc.violation(
            "presentation_age_out_of_range", case.case_id, f"{age} not in [{low}, {high}]"
        )


def _truly_present(case: SyntheticCase) -> set[str]:
    ids = {p.source_hpo_id or p.hpo_id for p in case.positive_phenotypes}
    ids |= {p.hpo_id for p in case.positive_phenotypes}
    ids |= {p.hpo_id for p in case.missing_phenotypes}
    ids |= {p.hpo_id for p in case.unknown_phenotypes}
    return ids


def _related(a: str, b: str, ontology: HpoOntology | None) -> bool:
    if a == b:
        return True
    if ontology is None:
        return False
    return ontology.is_a(a, b) or ontology.is_a(b, a)


def _annotation_closure(profile: DiseaseProfile, ontology: HpoOntology | None) -> frozenset[str]:
    closure: set[str] = set()
    for phenotype in profile.phenotypes:
        closure.add(phenotype.hpo_id)
        if ontology is not None:
            closure |= ontology.get_ancestor_set(phenotype.hpo_id)
            closure |= ontology.get_descendant_set(phenotype.hpo_id)
    return frozenset(closure)


def _onset_prior(profile: DiseaseProfile, config: SimulationConfig) -> dict[str, float]:
    windows = config.age.onset_years
    onset = profile.age_of_onset
    if onset is not None:
        weights = {c: w for c, w in onset.distribution.items() if c in windows and w > 0.0}
        if weights:
            return _normalize(weights)
        if onset.category in windows:
            return {onset.category: 1.0}
    return _normalize({c: w for c, w in config.age.unknown_onset_prior.items() if w > 0.0})


def _normalize(weights: Mapping[str, float]) -> dict[str, float]:
    total = sum(weights.values())
    return {key: value / total for key, value in weights.items()} if total else {}


def _calibration(
    acc: _Accumulator,
    ontology: HpoOntology | None,
    config: SimulationConfig | None,
) -> dict[str, Any]:
    unknown = config.frequency.unknown_frequency if config is not None else 0.5
    sex_terms = SexSpecificTerms(ontology, config.sex) if config is not None else None
    anchored = (sex_terms.male_only | sex_terms.female_only) if sex_terms is not None else set()
    all_pairs: list[tuple[float, float, int]] = []
    ungated: list[tuple[float, float, int]] = []
    for disease_id, cases in acc.diseases.items():
        profile = acc.case_profiles.get(disease_id)
        if profile is None:
            continue
        present = acc.present_counts[disease_id]
        for phenotype in profile.phenotypes:
            frequency = _simulated_frequency(phenotype, config)
            expected = frequency if frequency is not None else unknown
            pair = (expected, present[phenotype.hpo_id] / cases, cases)
            all_pairs.append(pair)
            gated = (
                phenotype.sex_restriction is not None
                or phenotype.hpo_id in anchored
                or (
                    phenotype.onset not in ("unknown", "antenatal")
                    and phenotype.onset_hpo_id != CONGENITAL_ONSET
                )
                or profile.progression == "progressive"
            )
            if not gated:
                ungated.append(pair)
    return {
        "all_terms": _calibration_table(all_pairs),
        "ungated_terms": _calibration_table(ungated),
    }


def _simulated_frequency(
    phenotype: PhenotypeAssociation, config: SimulationConfig | None
) -> float | None:
    """The frequency the run simulated the term with; the stored estimate without a prior mean."""

    if config is None:
        return phenotype.frequency_estimate
    try:
        return simulation_frequency(phenotype, config.frequency)
    except ShrinkageMeanMissing:
        return phenotype.frequency_estimate


def _calibration_table(pairs: Sequence[tuple[float, float, int]]) -> dict[str, Any]:
    bins = []
    total_weight = sum(weight for _, _, weight in pairs)
    weighted_error = 0.0
    for low, high in CALIBRATION_BINS:
        members = [pair for pair in pairs if low <= pair[0] < high]
        weight = sum(w for _, _, w in members)
        if not weight:
            continue
        expected = sum(e * w for e, _, w in members) / weight
        observed = sum(o * w for _, o, w in members) / weight
        weighted_error += weight * abs(observed - expected)
        bins.append(
            {
                "range": [low, min(high, 1.0)],
                "pairs": len(members),
                "mean_expected": round(expected, 4),
                "mean_observed": round(observed, 4),
            }
        )
    return {
        "pairs": len(pairs),
        "bins": bins,
        "calibration_error": round(weighted_error / total_weight, 4) if total_weight else None,
        "mean_abs_error_per_pair": (
            round(sum(abs(o - e) * w for e, o, w in pairs) / total_weight, 4)
            if total_weight
            else None
        ),
    }


def _observation_rates(
    acc: _Accumulator, config: SimulationConfig | None
) -> dict[str, dict[str, float | None]]:
    rates = {}
    for difficulty, counts in sorted(acc.observed_by_difficulty.items()):
        total = counts["positive"] + counts["unobserved"]
        rates[difficulty] = {
            "observed": _share(counts["positive"], total),
            "preset": (
                config.presets[difficulty].positive_observation_rate  # type: ignore[index]
                if config is not None
                else None
            ),
        }
    return rates


def _negative_report(acc: _Accumulator) -> dict[str, Any]:
    counts = acc.bucket_counts["negative_phenotypes"]
    total = sum(acc.negative_origins.values())
    return {
        "per_case_mean": _mean(counts),
        "cases_with_none_share": _share(sum(1 for count in counts if count == 0), acc.cases),
        "by_source": dict(sorted(acc.negative_origins.items())),
        "by_source_share": {
            origin: _share(count, total) for origin, count in sorted(acc.negative_origins.items())
        },
    }


def _missingness_report(acc: _Accumulator) -> dict[str, float | None]:
    positives = sum(acc.bucket_counts["positive_phenotypes"])
    missing = sum(acc.bucket_counts["missing_phenotypes"])
    unknown = sum(acc.bucket_counts["unknown_phenotypes"])
    present = positives + missing + unknown
    return {
        "missing_rate": _share(missing, present),
        "unknown_rate": _share(unknown, present),
        "sex_unknown_share": _share(acc.sex["unknown"], acc.cases),
        "age_unknown_share": _share(acc.cases - len(acc.ages), acc.cases),
        "onset_unknown_share": _share(acc.cases - len(acc.onsets), acc.cases),
    }


def _sex_report(acc: _Accumulator, config: SimulationConfig | None) -> dict[str, Any]:
    by_prior = {}
    for key, counts in sorted(acc.sex_by_prior.items()):
        known = counts["male"] + counts["female"]
        by_prior[key] = {
            "cases_with_known_sex": known,
            "male_share": _share(counts["male"], known),
            "expected_male_share": config.sex.p_male[key] if config is not None else None,  # type: ignore[index]
        }
    return {"counts": dict(sorted(acc.sex.items())), "by_prior": by_prior}


def _age_report(values: Sequence[float], cases: int) -> dict[str, Any]:
    histogram = {label: 0 for _, _, label in AGE_BINS_YEARS}
    for value in values:
        for low, high, label in AGE_BINS_YEARS:
            if low <= value < high:
                histogram[label] += 1
                break
    return {
        "known": len(values),
        "known_share": _share(len(values), cases),
        "mean_years": _mean(values),
        "median_years": round(statistics.median(values), 2) if values else None,
        "histogram": histogram,
    }


def _category_report(observed: Counter[str], expected: Counter[str]) -> dict[str, Any]:
    observed_total = sum(observed.values())
    expected_total = sum(expected.values())
    return {
        category: {
            "observed_share": _share(observed[category], observed_total),
            "expected_share": (
                round(expected[category] / expected_total, 4) if expected_total else None
            ),
        }
        for category in sorted(set(observed) | set(expected))
    }


def _distribution(counts: Sequence[int]) -> dict[str, Any]:
    histogram = Counter(min(count, 15) for count in counts)
    return {
        "mean": _mean(counts),
        "median": statistics.median(counts) if counts else None,
        "min": min(counts) if counts else None,
        "max": max(counts) if counts else None,
        "histogram": {
            ("15+" if value == 15 else str(value)): histogram[value] for value in sorted(histogram)
        },
    }


def _mean(values: Sequence[float]) -> float | None:
    return round(statistics.fmean(values), 3) if values else None


def _share(part: int, total: int) -> float | None:
    return round(part / total, 4) if total else None


def format_report(report: Mapping[str, Any]) -> str:
    """Readable text rendering of :func:`validate_cases` output."""

    lines = [
        f"Cases: {report['cases']} across {report['diseases']} disease(s) "
        f"and {report['genes']} gene(s)"
    ]
    if report["config_matches_cases"] is False:
        lines.append("WARNING: the config hash differs from the cases' config hash")
    per_case = report["phenotypes_per_case"]
    lines.append(
        "Per case (mean): "
        + ", ".join(f"{name} {stats['mean']}" for name, stats in per_case.items())
    )
    negatives = report["negatives"]
    lines.append(
        f"Negatives: {negatives['per_case_mean']} per case, "
        f"{negatives['cases_with_none_share']} of cases without any; by source "
        + ", ".join(f"{k} {v}" for k, v in negatives["by_source_share"].items())
    )
    missingness = report["missingness"]
    lines.append(
        "Missingness: "
        + ", ".join(f"{key} {value}" for key, value in missingness.items())
    )
    lines.append(
        f"Noise: {report['noise']['per_case_mean']} per case "
        f"(related share {report['noise']['related_share']}); "
        f"specialized share of positives {report['positives']['specialized_share']}; "
        f"forced positives in {report['positives']['cases_with_forced_positive']} case(s)"
    )
    sex = report["sex"]
    lines.append("Sex: " + ", ".join(f"{k} {v}" for k, v in sex["counts"].items()))
    for key, stats in sex["by_prior"].items():
        lines.append(
            f"  {key}: male share {stats['male_share']} "
            f"(expected {stats['expected_male_share']}, n={stats['cases_with_known_sex']})"
        )
    for name in ("age", "onset"):
        stats = report[name]
        lines.append(
            f"{name.capitalize()}: known {stats['known_share']}, mean {stats['mean_years']} y, "
            f"median {stats['median_years']} y; "
            + ", ".join(f"{k} {v}" for k, v in stats["histogram"].items())
        )
    for category, stats in report["onset"]["by_category"].items():
        lines.append(
            f"  onset {category}: observed {stats['observed_share']}, "
            f"expected {stats['expected_share']}"
        )
    calibration = report["calibration"]
    if calibration is not None:
        for view in ("all_terms", "ungated_terms"):
            table = calibration[view]
            lines.append(
                f"Calibration ({view}, {table['pairs']} pairs): "
                f"error {table['calibration_error']}, "
                f"per-pair MAE {table['mean_abs_error_per_pair']}"
            )
            for row in table["bins"]:
                lines.append(
                    f"  [{row['range'][0]:.2f}, {row['range'][1]:.2f}): "
                    f"expected {row['mean_expected']}, observed {row['mean_observed']} "
                    f"(n={row['pairs']})"
                )
    reporting = report.get("reporting")
    if reporting is not None and reporting.get("mode") == "independent":
        lines.append(
            f"Independent reporting ({reporting['cases']} cases): "
            f"{reporting['reported_per_case_mean']} profile terms and "
            f"{reporting['noise_per_case_mean']} noise per case; "
            f"{reporting['true_terms_reported_share']} of true terms reported; "
            f"q capped {reporting['q_capped_share']}; "
            f"forced_min_one {reporting['forced_min_one']}; "
            f"q calibration error {reporting['q_calibration_error']}"
        )
        for row in reporting["q_calibration"]:
            lines.append(
                f"  q [{row['q_range'][0]:.1f}, {row['q_range'][1]:.1f}): mean q "
                f"{row['mean_q']}, reported {row['reported_share']} (n={row['terms']})"
            )
    elif reporting is not None:
        lines.append(
            f"Reporting ({reporting['cases']} report-model cases): budget mean "
            f"{reporting['budget_mean']}, reported mean {reporting['reported_total_mean']} "
            f"({reporting['reported_mean']} profile terms + noise); TV vs histogram: "
            f"drawn {reporting['total_variation_drawn_vs_expected']}, reported "
            f"{reporting['total_variation_reported_total_vs_expected']}"
        )
        cardinal = reporting["cardinal"]
        if cardinal is not None:
            lines.append(
                f"  cardinal: {cardinal['reported_terms']} of {cardinal['true_terms']} true "
                f"cardinal terms reported, in {cardinal['cases_with_true_cardinal']} case(s)"
            )
    violations = report["violations"]
    lines.append(f"Invariant violations: {violations['count']}")
    for kind, count in violations["by_type"].items():
        lines.append(f"  {kind}: {count}")
    return "\n".join(lines)
