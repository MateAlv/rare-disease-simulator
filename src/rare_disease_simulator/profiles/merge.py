"""One merged profile per disease entity (ADR-0011, simulator v0.4).

A ``genes-v1``/``genes-v2`` entity lists the exactly-equivalent OMIM and ORPHA
ids of one disease. Simulating it from one of them at random halves every
feature only one source records, so the entity is simulated from a single
profile merged from all of them:

- **terms:** the union;
- **frequency:** the ADR-0007 notations of every profile are pooled and the
  usual precedence applies, counts (summed ``n/m``) over percentages over
  categories. A count row that two profiles repeat with the same references is
  counted once. Under this precedence an Orphanet category never overrides a
  count-based estimate, so the ADR's cap (a category is capped at the
  count-based estimate when they differ by more than 0.3) always holds; the
  merge counts how often it bound (``category_capped``);
- **sex restriction:** kept only when every profile with the term names the
  same sex, the rule the builder applies to HPOA rows;
- **phenotype onset:** the earliest onset any profile gives the term, the
  builder's rule for rows, so no source's early presentation is gated away;
- **disease onset:** the mean of the profiles' onset distributions;
- **inheritance:** the union of the profiles' modes, with the sex bias derived
  again from it; **progression** is the shared value, ``variable`` when the
  profiles disagree;
- **negatives:** the union of ``NOT`` terms that no profile annotates as present.

Frequencies stay ADR-0007 estimates here; the simulator applies its count
estimator (Beta shrinkage or Jeffreys) to the pooled ``n/m``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.frequency import (
    HPO_FREQUENCY_MIDPOINTS,
    ParsedFrequencies,
    estimate_frequency,
    parse_frequencies,
    pooled_counts,
)
from rare_disease_simulator.profiles.inheritance import derive_sex_bias
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiagnosticRole,
    DiseaseGene,
    DiseaseProfile,
    MappedDiseaseIds,
    NegativePhenotypeAssociation,
    OnsetCategory,
    PhenotypeAssociation,
    ProbabilityRange,
    ProfileQuality,
    ProgressionPattern,
    Provenance,
    SexBias,
)

CATEGORY_CAP_THRESHOLD = 0.3

ONSET_ORDER: tuple[OnsetCategory, ...] = (
    "antenatal",
    "neonatal",
    "infantile",
    "childhood",
    "childhood_or_adolescent",
    "juvenile",
    "adult",
    "variable",
)
ROLE_ORDER: tuple[DiagnosticRole, ...] = ("cardinal", "major", "supportive", "nonspecific")


@dataclass
class MergeStats:
    """What merging an entity's profiles did, summed over entities."""

    counters: Counter[str] = field(default_factory=Counter)

    def add(self, other: MergeStats) -> None:
        self.counters.update(other.counters)

    def summary(self) -> dict[str, int]:
        return dict(sorted(self.counters.items()))


def merge_entity_profiles(
    entity_id: str,
    profiles: Sequence[DiseaseProfile],
    ontology: HpoOntology | None = None,
) -> tuple[DiseaseProfile, MergeStats]:
    """Merge the profiles of one entity into a profile whose ``disease_id`` is the entity."""

    if not profiles:
        raise ValueError(f"entity {entity_id} has no profile to merge")
    stats = MergeStats()
    stats.counters["entities"] += 1
    stats.counters[f"entities_with_{len(profiles)}_profiles"] += 1
    if len(profiles) == 1:
        stats.counters["terms"] += len(profiles[0].phenotypes)
        for phenotype in profiles[0].phenotypes:
            stats.counters[f"basis_{_basis(phenotype)}"] += 1
        return profiles[0].model_copy(update={"disease_id": entity_id}), stats

    by_term: dict[str, list[PhenotypeAssociation]] = {}
    for profile in profiles:
        for phenotype in profile.phenotypes:
            by_term.setdefault(phenotype.hpo_id, []).append(phenotype)
    phenotypes = [
        _merge_phenotype(by_term[hpo_id], stats) for hpo_id in sorted(by_term)
    ]
    stats.counters["terms"] += len(phenotypes)
    stats.counters["terms_in_several_profiles"] += sum(
        1 for items in by_term.values() if len(items) > 1
    )
    if len(profiles) > 1:
        largest = max(len(profile.phenotypes) for profile in profiles)
        stats.counters["terms_added_over_largest_profile"] += len(phenotypes) - largest

    positive_ids = set(by_term)
    negatives: dict[str, NegativePhenotypeAssociation] = {}
    for profile in profiles:
        for negative in profile.negative_phenotypes:
            if negative.hpo_id in positive_ids:
                stats.counters["negatives_dropped_positive_elsewhere"] += 1
                continue
            previous = negatives.get(negative.hpo_id)
            sources = sorted(set(negative.source) | set(previous.source if previous else ()))
            negatives[negative.hpo_id] = negative.model_copy(update={"source": sources})

    genes = _merge_genes(profiles)
    inheritance = sorted({hpo_id for gene in genes for hpo_id in gene.inheritance_hpo_ids})
    names = {profile.disease_id: profile.disease_name for profile in profiles}
    return (
        DiseaseProfile(
            disease_id=entity_id,
            disease_name=names.get(entity_id, profiles[0].disease_name),
            mapped_ids=_merge_mapped_ids(profiles),
            synonyms=sorted({name for name in names.values()} - {names.get(entity_id)}),
            genes=genes,
            phenotypes=phenotypes,
            negative_phenotypes=[negatives[hpo_id] for hpo_id in sorted(negatives)],
            age_of_onset=_merge_onset(profiles),
            sex_bias=_merge_sex_bias(profiles, inheritance, ontology),
            progression=_merge_progression(profiles),
            provenance=_merge_provenance(profiles),
            quality=ProfileQuality(
                warnings=sorted({w for profile in profiles for w in profile.quality.warnings}),
                counters={
                    "merged_profiles": len(profiles),
                    **{
                        key: value
                        for key, value in sorted(stats.counters.items())
                        if key.startswith(("category_capped", "sex_restriction_dropped"))
                    },
                },
            ),
        ),
        stats,
    )


def _merge_phenotype(
    items: Sequence[PhenotypeAssociation], stats: MergeStats
) -> PhenotypeAssociation:
    first = items[0]
    if len(items) == 1:
        stats.counters[f"basis_{_basis(first)}"] += 1
        return first

    ratios: dict[tuple[tuple[int, int], tuple[str, ...]], None] = {}
    percents: list[float] = []
    categories: list[str] = []
    for item in items:
        counts = pooled_counts(item.frequency_raw)
        if counts is not None:
            key = (counts, tuple(sorted(item.source)))
            if key in ratios:
                stats.counters["count_rows_deduplicated"] += 1
            ratios[key] = None
            continue
        parsed = parse_frequencies((item.frequency_raw or "").split(";"))
        if parsed.percents:
            percents.extend(parsed.percents)
        elif parsed.categories:
            categories.extend(parsed.categories)
        elif item.frequency_estimate is not None:
            # A profile without HPOA notation (fixture or curated) contributes its point.
            percents.append(item.frequency_estimate)
    pooled = ParsedFrequencies(
        ratios=tuple(counts for counts, _ in ratios),
        percents=tuple(percents),
        categories=tuple(categories),
    )
    estimate = estimate_frequency(pooled)
    if pooled.ratios and categories and estimate is not None:
        category_mean = sum(HPO_FREQUENCY_MIDPOINTS[c] for c in categories) / len(categories)
        if category_mean - estimate.estimate > CATEGORY_CAP_THRESHOLD:
            stats.counters["category_capped"] += 1

    restrictions = {item.sex_restriction for item in items}
    restriction = next(iter(restrictions)) if len(restrictions) == 1 else None
    if len(restrictions) > 1 and restrictions - {None}:
        stats.counters["sex_restriction_dropped"] += 1

    onsets = [
        (ONSET_ORDER.index(item.onset), item.onset_hpo_id or "", item.onset)
        for item in items
        if item.onset in ONSET_ORDER
    ]
    onset: OnsetCategory = "unknown"
    onset_hpo_id = None
    if onsets:
        _, onset_id, onset = min(onsets)
        onset_hpo_id = onset_id or None
        if len({category for _, _, category in onsets}) > 1:
            stats.counters["onset_conflicts_earliest_kept"] += 1

    roles = [item.diagnostic_role for item in items if item.diagnostic_role in ROLE_ORDER]
    merged = PhenotypeAssociation(
        hpo_id=first.hpo_id,
        label=first.label,
        frequency=estimate.category if estimate else "unknown",  # type: ignore[arg-type]
        frequency_raw=estimate.raw if estimate else None,
        frequency_estimate=estimate.estimate if estimate else None,
        probability_range=(
            ProbabilityRange(lower=estimate.lower, upper=estimate.upper) if estimate else None
        ),
        diagnostic_role=min(roles, key=ROLE_ORDER.index) if roles else "unknown",
        onset=onset,
        onset_hpo_id=onset_hpo_id,
        sex_restriction=restriction,
        source=sorted({reference for item in items for reference in item.source}),
    )
    stats.counters[f"basis_{_basis(merged)}"] += 1
    return merged


def _basis(phenotype: PhenotypeAssociation) -> str:
    if pooled_counts(phenotype.frequency_raw) is not None:
        return "counts"
    parsed = parse_frequencies((phenotype.frequency_raw or "").split(";"))
    if parsed.percents:
        return "percent"
    if parsed.categories:
        return "category"
    return "estimate_only" if phenotype.frequency_estimate is not None else "unknown"


def _merge_genes(profiles: Sequence[DiseaseProfile]) -> list[DiseaseGene]:
    genes: dict[str, DiseaseGene] = {}
    for profile in profiles:
        for gene in profile.genes:
            previous = genes.get(gene.symbol)
            if previous is None:
                genes[gene.symbol] = gene
                continue
            genes[gene.symbol] = previous.model_copy(
                update={
                    "ncbi_gene_id": previous.ncbi_gene_id or gene.ncbi_gene_id,
                    "association_type": (
                        "causal"
                        if "causal" in (previous.association_type, gene.association_type)
                        else previous.association_type
                    ),
                    "inheritance": list(dict.fromkeys([*previous.inheritance, *gene.inheritance])),
                    "inheritance_hpo_ids": sorted(
                        set(previous.inheritance_hpo_ids) | set(gene.inheritance_hpo_ids)
                    ),
                }
            )
    return list(genes.values())


def _merge_mapped_ids(profiles: Sequence[DiseaseProfile]) -> MappedDiseaseIds:
    ids = [profile.disease_id for profile in profiles]
    orpha = next((i for i in ids if i.startswith("ORPHA:")), None)
    omim = sorted({i for i in ids if i.startswith("OMIM:")})
    for profile in profiles:
        orpha = orpha or profile.mapped_ids.orpha
        omim = sorted(set(omim) | set(profile.mapped_ids.omim))
    return MappedDiseaseIds(orpha=orpha, omim=omim)


def _merge_onset(profiles: Sequence[DiseaseProfile]) -> AgeOfOnset | None:
    known = [profile.age_of_onset for profile in profiles if profile.age_of_onset is not None]
    if not known:
        return None
    if len(known) == 1:
        return known[0]
    totals: Counter[OnsetCategory] = Counter()
    for onset in known:
        distribution = onset.distribution or (
            {onset.category: 1.0} if onset.category != "unknown" else {}
        )
        weight = sum(distribution.values())
        for category, share in distribution.items():
            totals[category] += share / weight / len(known) if weight else 0.0
    distribution = {
        category: round(totals[category], 4) for category in ONSET_ORDER if totals[category] > 0
    }
    category = (
        min(distribution, key=lambda c: (-distribution[c], ONSET_ORDER.index(c)))
        if distribution
        else "unknown"
    )
    return AgeOfOnset(
        category=category,
        distribution=distribution,
        hpo_ids=sorted({hpo_id for onset in known for hpo_id in onset.hpo_ids}),
        provenance=[item for onset in known for item in onset.provenance],
    )


def _merge_sex_bias(
    profiles: Sequence[DiseaseProfile], inheritance: Sequence[str], ontology: HpoOntology | None
) -> SexBias | None:
    value, basis = derive_sex_bias(inheritance, ontology)
    if value is None:
        return None
    return SexBias(
        value=value,
        provenance=[
            Provenance(
                source=profiles[0].provenance[0].source,
                field="sex_bias",
                evidence="derived from the merged profiles' inheritance: " + "; ".join(basis),
            )
        ]
        if profiles[0].provenance
        else [],
    )


def _merge_progression(profiles: Sequence[DiseaseProfile]) -> ProgressionPattern:
    known = {profile.progression for profile in profiles} - {"unknown"}
    if not known:
        return "unknown"
    if len(known) > 1:
        return "variable"
    return next(iter(known))


def _merge_provenance(profiles: Sequence[DiseaseProfile]) -> list[Provenance]:
    seen: dict[tuple[str, str | None, str | None], Provenance] = {}
    for profile in profiles:
        for item in profile.provenance:
            key = (item.source.name, item.source.sha256, item.field)
            seen.setdefault(key, item)
    return list(seen.values())
