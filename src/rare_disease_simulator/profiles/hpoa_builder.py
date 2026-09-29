"""Build :class:`DiseaseProfile` records from the HPO annotation release.

This is the structured, LLM-free backbone: one profile per OMIM/ORPHA/DECIPHER
disease in ``phenotype.hpoa`` that has at least one gene in
``genes_to_disease.txt``. Mapping rules (see ``docs/README.md``):

- aspect ``P`` rows become phenotypes, their frequency read per ADR-0007
  (``profiles/frequency.py``); only ``NOT`` rows and the ``Excluded`` frequency
  term become negative phenotypes (a ``0/m`` count is a low frequency);
- aspect ``I`` rows become gene inheritance and, where sound, a sex bias;
- aspect ``C`` onset rows become the disease age of onset, with Orphanet average
  age of onset as fallback; ``C`` pace-of-progression rows set ``progression``;
- rows citing a held-out reference are dropped before anything else.

Output is deterministic: diseases, genes and terms are emitted in sorted order
and no timestamps are written into profiles.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.data_sources.hpo_annotations import (
    HpoaRow,
    iter_hpoa_rows,
    read_hpoa_header,
    read_reference_mask,
)
from rare_disease_simulator.data_sources.hpo_genes import GeneDiseaseLink, read_genes_to_disease
from rare_disease_simulator.data_sources.orphanet_products import (
    OmimOrphaMap,
    OrphanetOnsets,
    read_omim_orpha_map,
    read_orphanet_onsets,
)
from rare_disease_simulator.profiles.frequency import (
    FrequencyEstimate,
    estimate_frequency,
    parse_frequencies,
)
from rare_disease_simulator.profiles.inheritance import derive_sex_bias
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiseaseGene,
    DiseaseProfile,
    GeneAssociationType,
    MappedDiseaseIds,
    NegativePhenotypeAssociation,
    OnsetCategory,
    PhenotypeAssociation,
    ProbabilityRange,
    ProfileQuality,
    ProgressionPattern,
    Provenance,
    SexBias,
    SourceReference,
)

ONSET_ROOT = "HP:0003674"
MODE_OF_INHERITANCE_ROOT = "HP:0000005"
EXCLUDED_FREQUENCY = "HP:0040285"

# Checked in order; the first anchor the term is (a descendant of) wins, so the
# specific pediatric anchors precede their parent "Pediatric onset".
ONSET_ANCHORS: tuple[tuple[str, OnsetCategory], ...] = (
    ("HP:0030674", "antenatal"),
    # Congenital (present at birth) has no own category; neonatal is the closest age window.
    ("HP:0003577", "neonatal"),
    ("HP:0003623", "neonatal"),
    ("HP:0003593", "infantile"),
    ("HP:0011463", "childhood"),
    ("HP:0003621", "juvenile"),
    ("HP:0410280", "childhood_or_adolescent"),
    ("HP:0003581", "adult"),
    ("HP:4000040", "adult"),
    ("HP:6000314", "adult"),
    ("HP:6000315", "adult"),
)

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

ORPHANET_ONSET_CATEGORIES: dict[str, OnsetCategory] = {
    "Antenatal": "antenatal",
    "Neonatal": "neonatal",
    "Infancy": "infantile",
    "Childhood": "childhood",
    "Adolescent": "juvenile",
    "Adult": "adult",
    "Elderly": "adult",
    "All ages": "variable",
}

PROGRESSION_ANCHORS: tuple[tuple[str, ProgressionPattern], ...] = (
    ("HP:0003676", "progressive"),
    ("HP:0003680", "non_progressive"),
    ("HP:0025303", "episodic"),
)

ASSOCIATION_TYPES: dict[str, GeneAssociationType] = {
    "MENDELIAN": "causal",
    "POLYGENIC": "susceptibility",
}


SEX_RESTRICTIONS = {"MALE": "male", "FEMALE": "female"}

_ID_NUMBER_RE = re.compile(r"(\d+)$")
_ROUND_DIGITS = 4


@dataclass(frozen=True)
class HpoaBuildInputs:
    """Paths to every file the HPOA build reads."""

    hpo_json: Path
    phenotype_hpoa: Path
    genes_to_disease: Path
    orphanet_ages: Path | None = None
    omim_orpha_map: Path | None = None
    exclude_pmids: Path | None = None


@dataclass
class HpoaBuildResult:
    """Built profiles plus a JSON-serializable build summary."""

    profiles: list[DiseaseProfile]
    summary: dict[str, Any]


@dataclass
class _SourceContext:
    hpoa: SourceReference
    ontology: SourceReference
    genes: SourceReference
    orphanet: SourceReference | None
    alignments: SourceReference | None
    mask: SourceReference | None

    def base_provenance(self) -> list[Provenance]:
        sources = [self.hpoa, self.genes, self.ontology]
        if self.mask is not None:
            sources.append(self.mask)
        return [Provenance(source=source) for source in sources]


@dataclass
class _Stats:
    rows: Counter[str] = field(default_factory=Counter)
    rows_by_aspect: Counter[str] = field(default_factory=Counter)
    diseases: Counter[str] = field(default_factory=Counter)
    diseases_by_prefix: Counter[str] = field(default_factory=Counter)
    phenotypes: Counter[str] = field(default_factory=Counter)
    frequency_categories: Counter[str] = field(default_factory=Counter)
    frequency_basis: Counter[str] = field(default_factory=Counter)
    negatives: Counter[str] = field(default_factory=Counter)
    onset_sources: Counter[str] = field(default_factory=Counter)
    onset_categories: Counter[str] = field(default_factory=Counter)
    inheritance_terms: Counter[str] = field(default_factory=Counter)
    sex_bias: Counter[str] = field(default_factory=Counter)
    progression: Counter[str] = field(default_factory=Counter)
    masked_diseases: set[str] = field(default_factory=set)
    genes: set[str] = field(default_factory=set)
    gene_links: int = 0
    multi_gene_diseases: int = 0
    unique_terms: set[str] = field(default_factory=set)


def build_profiles_from_hpoa(
    inputs: HpoaBuildInputs, *, ontology: HpoOntology | None = None
) -> HpoaBuildResult:
    """Build one profile per gene-linked HPOA disease and summarize the build."""

    ontology = ontology or HpoOntology.from_json(inputs.hpo_json)
    hpoa_header = read_hpoa_header(inputs.phenotype_hpoa)
    genes = read_genes_to_disease(inputs.genes_to_disease)
    mask = read_reference_mask(inputs.exclude_pmids) if inputs.exclude_pmids else frozenset()
    orphanet = read_orphanet_onsets(inputs.orphanet_ages) if inputs.orphanet_ages else None
    omim_orpha = read_omim_orpha_map(inputs.omim_orpha_map) if inputs.omim_orpha_map else None

    input_records = _input_records(inputs, ontology, hpoa_header, orphanet, omim_orpha)
    sources = _source_context(input_records)
    stats = _Stats()

    rows_by_disease: dict[str, list[HpoaRow]] = defaultdict(list)
    names: dict[str, str] = {}
    positives_before_mask: set[str] = set()
    for row in iter_hpoa_rows(inputs.phenotype_hpoa):
        stats.rows["total"] += 1
        stats.rows_by_aspect[row.aspect] += 1
        names.setdefault(row.disease_id, row.disease_name)
        if row.aspect == "P" and not row.negated:
            positives_before_mask.add(row.disease_id)
        if mask and any(reference in mask for reference in row.references):
            stats.rows["masked"] += 1
            if row.disease_id in genes.links:
                stats.rows["masked_gene_linked"] += 1
            stats.masked_diseases.add(row.disease_id)
            continue
        if row.disease_id not in genes.links:
            stats.rows["disease_without_gene"] += 1
            continue
        rows_by_disease[row.disease_id].append(row)

    stats.diseases["in_hpoa"] = len(names)
    stats.diseases["with_gene"] = sum(1 for disease_id in names if disease_id in genes.links)
    stats.diseases["gene_linked_not_in_hpoa"] = sum(
        1 for disease_id in genes.links if disease_id not in names
    )

    profiles: list[DiseaseProfile] = []
    gene_linked = [disease_id for disease_id in names if disease_id in genes.links]
    for disease_id in sorted(gene_linked, key=_disease_sort_key):
        profile = _build_profile(
            disease_id,
            names[disease_id],
            rows_by_disease.get(disease_id, []),
            genes.links[disease_id],
            ontology=ontology,
            orphanet=orphanet,
            omim_orpha=omim_orpha,
            sources=sources,
            stats=stats,
        )
        if profile is None:
            stats.diseases["dropped_zero_positive"] += 1
            if disease_id in positives_before_mask and disease_id in stats.masked_diseases:
                stats.diseases["dropped_zero_positive_due_to_mask"] += 1
            continue
        profiles.append(profile)
        _record_profile_stats(profile, stats)

    summary = _summary(
        stats,
        profiles,
        input_records=input_records,
        omim_orpha=omim_orpha,
        mask_size=len(mask),
        genes_rows_total=genes.rows_total,
        genes_rows_missing_symbol=genes.rows_missing_symbol,
    )
    return HpoaBuildResult(profiles=profiles, summary=summary)


def _build_profile(
    disease_id: str,
    disease_name: str,
    rows: Sequence[HpoaRow],
    gene_links: Sequence[GeneDiseaseLink],
    *,
    ontology: HpoOntology,
    orphanet: OrphanetOnsets | None,
    omim_orpha: OmimOrphaMap | None,
    sources: _SourceContext,
    stats: _Stats,
) -> DiseaseProfile | None:
    counters: Counter[str] = Counter()
    positive_rows: dict[str, list[HpoaRow]] = defaultdict(list)
    negative_rows: dict[str, list[HpoaRow]] = defaultdict(list)
    negative_reasons: dict[str, set[str]] = defaultdict(set)
    inheritance_ids: set[str] = set()
    onset_ids: set[str] = set()
    progression_ids: set[str] = set()

    for row in rows:
        resolution = ontology.resolve(row.hpo_id)
        if resolution.hpo_id is None:
            counters["rows_unresolved_term"] += 1
            stats.rows["unresolved_term"] += 1
            continue
        if resolution.status != "current":
            counters[f"rows_term_{resolution.status}"] += 1
            stats.rows[f"term_{resolution.status}"] += 1
        hpo_id = resolution.hpo_id

        if row.aspect == "P":
            stats.rows["processed"] += 1
            if row.negated:
                negative_rows[hpo_id].append(row)
                negative_reasons[hpo_id].add("not_qualifier")
            elif row.frequency == EXCLUDED_FREQUENCY:
                negative_rows[hpo_id].append(row)
                negative_reasons[hpo_id].add("excluded_frequency")
            else:
                positive_rows[hpo_id].append(row)
        elif row.aspect == "I":
            stats.rows["processed"] += 1
            if ontology.is_a(hpo_id, MODE_OF_INHERITANCE_ROOT):
                inheritance_ids.add(hpo_id)
            else:
                stats.rows["inheritance_not_mode_of_inheritance"] += 1
        elif row.aspect == "C":
            stats.rows["processed"] += 1
            if ontology.is_a(hpo_id, ONSET_ROOT):
                onset_ids.add(hpo_id)
            elif _progression_for(hpo_id, ontology) is not None:
                progression_ids.add(hpo_id)
            else:
                stats.rows["clinical_course_other"] += 1
        else:
            stats.rows[f"skipped_aspect_{row.aspect or 'blank'}"] += 1

    phenotypes = [
        _merge_positive(hpo_id, positive_rows[hpo_id], ontology, stats)
        for hpo_id in sorted(positive_rows)
    ]

    if not phenotypes:
        return None

    positive_ids = {phenotype.hpo_id for phenotype in phenotypes}
    negatives: list[NegativePhenotypeAssociation] = []
    for hpo_id in sorted(negative_rows):
        if hpo_id in positive_ids:
            counters["positive_negative_conflicts"] += 1
            stats.negatives["dropped_conflict_with_positive"] += 1
            continue
        for reason in negative_reasons[hpo_id]:
            stats.negatives[f"from_{reason}"] += 1
        negatives.append(
            NegativePhenotypeAssociation(
                hpo_id=hpo_id,
                label=ontology.get_label(hpo_id) or hpo_id,
                source=_references(negative_rows[hpo_id]),
            )
        )

    counters["phenotypes_unknown_frequency"] = sum(
        1 for phenotype in phenotypes if phenotype.frequency == "unknown"
    )
    inheritance = sorted(inheritance_ids)
    age_of_onset, onset_source = _age_of_onset(
        disease_id, onset_ids, ontology, orphanet, omim_orpha, sources
    )
    stats.onset_sources[onset_source] += 1

    warnings: list[str] = []
    if age_of_onset is None:
        warnings.append("missing_age_of_onset")
    if not inheritance:
        warnings.append("missing_inheritance")

    return DiseaseProfile(
        disease_id=disease_id,
        disease_name=disease_name,
        mapped_ids=_mapped_ids(disease_id, omim_orpha),
        genes=[_disease_gene(link, inheritance, ontology) for link in gene_links],
        phenotypes=phenotypes,
        negative_phenotypes=negatives,
        age_of_onset=age_of_onset,
        sex_bias=_sex_bias(inheritance, ontology, sources.hpoa),
        progression=_progression(progression_ids, ontology),
        provenance=sources.base_provenance(),
        quality=ProfileQuality(
            warnings=warnings,
            counters={key: value for key, value in sorted(counters.items()) if value},
        ),
    )


def _merge_positive(
    hpo_id: str,
    rows: Sequence[HpoaRow],
    ontology: HpoOntology,
    stats: _Stats,
) -> PhenotypeAssociation:
    """Merge all positive rows of one (disease, term)."""

    frequency = _merge_frequency(rows, stats)
    onset_category, onset_hpo_id = _merge_phenotype_onset(rows, ontology, stats)
    return PhenotypeAssociation(
        hpo_id=hpo_id,
        label=ontology.get_label(hpo_id) or hpo_id,
        frequency=frequency.category if frequency else "unknown",  # type: ignore[arg-type]
        frequency_raw=frequency.raw if frequency else None,
        frequency_estimate=frequency.estimate if frequency else None,
        probability_range=(
            ProbabilityRange(lower=frequency.lower, upper=frequency.upper) if frequency else None
        ),
        onset=onset_category,
        onset_hpo_id=onset_hpo_id,
        sex_restriction=_merge_sex(rows, stats),  # type: ignore[arg-type]
        source=_references(rows),
    )


def _merge_frequency(rows: Sequence[HpoaRow], stats: _Stats) -> FrequencyEstimate | None:
    """Combine the frequency of duplicate rows per ADR-0007 (see ``profiles/frequency.py``)."""

    parsed = parse_frequencies(row.frequency for row in rows)
    if parsed.unparsed:
        stats.rows["unparsed_frequency"] += parsed.unparsed
    estimate = estimate_frequency(parsed)
    if estimate is not None:
        stats.frequency_basis[estimate.basis] += 1
    return estimate


def _merge_phenotype_onset(
    rows: Sequence[HpoaRow], ontology: HpoOntology, stats: _Stats
) -> tuple[OnsetCategory, str | None]:
    """Pick the earliest onset any source reports for the phenotype."""

    candidates: list[tuple[int, str, OnsetCategory]] = []
    for row in rows:
        if row.onset is None:
            continue
        resolved = ontology.resolve(row.onset).hpo_id
        category = _onset_category(resolved, ontology) if resolved else None
        if resolved is None or category is None:
            stats.rows["invalid_phenotype_onset"] += 1
            continue
        candidates.append((ONSET_ORDER.index(category), resolved, category))
    if not candidates:
        return "unknown", None
    _, hpo_id, category = min(candidates)
    return category, hpo_id


def _merge_sex(rows: Sequence[HpoaRow], stats: _Stats) -> str | None:
    """A restriction holds only when every row reports the same single sex."""

    values = {SEX_RESTRICTIONS.get(row.sex or "") for row in rows}
    if len(values) == 1:
        (value,) = values
        return value
    if values - {None}:
        stats.phenotypes["sex_restriction_conflicts"] += 1
    return None


def _onset_category(hpo_id: str, ontology: HpoOntology) -> OnsetCategory | None:
    for anchor, category in ONSET_ANCHORS:
        if ontology.is_a(hpo_id, anchor):
            return category
    return None


def _age_of_onset(
    disease_id: str,
    onset_ids: Iterable[str],
    ontology: HpoOntology,
    orphanet: OrphanetOnsets | None,
    omim_orpha: OmimOrphaMap | None,
    sources: _SourceContext,
) -> tuple[AgeOfOnset | None, str]:
    categorized = sorted(
        (hpo_id, category)
        for hpo_id in onset_ids
        if (category := _onset_category(hpo_id, ontology)) is not None
    )
    if categorized:
        return (
            AgeOfOnset(
                category=_modal_category([category for _, category in categorized]),
                distribution=_distribution([category for _, category in categorized]),
                hpo_ids=[hpo_id for hpo_id, _ in categorized],
                provenance=[Provenance(source=sources.hpoa, field="age_of_onset")],
            ),
            "hpoa",
        )

    if orphanet is None or sources.orphanet is None:
        return None, "none"
    via = "orphanet_direct"
    orpha_ids: tuple[str, ...] = (disease_id,) if disease_id.startswith("ORPHA:") else ()
    if not orpha_ids and omim_orpha is not None:
        orpha_ids = omim_orpha.omim_to_orpha.get(disease_id, ())
        via = "orphanet_via_omim"
    labels = [
        (orpha_id, label)
        for orpha_id in orpha_ids
        for label in orphanet.onsets.get(orpha_id, ())
        if label in ORPHANET_ONSET_CATEGORIES
    ]
    categories = [ORPHANET_ONSET_CATEGORIES[label] for _, label in labels]
    if not categories:
        return None, "none"
    evidence = "Orphanet AverageAgeOfOnset: " + "; ".join(
        f"{label} ({orpha_id})" for orpha_id, label in labels
    )
    provenance = [Provenance(source=sources.orphanet, field="age_of_onset", evidence=evidence)]
    if via == "orphanet_via_omim" and sources.alignments is not None:
        provenance.append(
            Provenance(
                source=sources.alignments,
                field="age_of_onset",
                evidence="exact validated OMIM-ORPHA alignment: "
                + ", ".join(sorted({orpha_id for orpha_id, _ in labels})),
            )
        )
    return (
        AgeOfOnset(
            category=_modal_category(categories),
            distribution=_distribution(categories),
            provenance=provenance,
        ),
        via,
    )


def _modal_category(categories: Sequence[OnsetCategory]) -> OnsetCategory:
    """Most supported category; ties go to the earliest onset."""

    counts = Counter(categories)
    return min(counts, key=lambda category: (-counts[category], ONSET_ORDER.index(category)))


def _distribution(categories: Sequence[OnsetCategory]) -> dict[OnsetCategory, float]:
    counts = Counter(categories)
    total = sum(counts.values())
    return {
        category: _round(counts[category] / total)
        for category in ONSET_ORDER
        if category in counts
    }


def _progression_for(hpo_id: str, ontology: HpoOntology) -> ProgressionPattern | None:
    for anchor, pattern in PROGRESSION_ANCHORS:
        if ontology.is_a(hpo_id, anchor):
            return pattern
    return None


def _progression(progression_ids: Iterable[str], ontology: HpoOntology) -> ProgressionPattern:
    patterns = {_progression_for(hpo_id, ontology) for hpo_id in progression_ids} - {None}
    if not patterns:
        return "unknown"
    if len(patterns) > 1:
        return "variable"
    (pattern,) = patterns
    return pattern  # type: ignore[return-value]


def _sex_bias(
    inheritance: Sequence[str], ontology: HpoOntology, source: SourceReference
) -> SexBias | None:
    """Derive a sex bias from inheritance only where the genetics imply it."""

    value, basis = derive_sex_bias(inheritance, ontology)
    if value is None:
        return None
    evidence = "derived from inheritance: " + "; ".join(
        f"{ontology.get_label(hpo_id)} ({hpo_id})" for hpo_id in basis
    )
    return SexBias(
        value=value,
        provenance=[Provenance(source=source, field="sex_bias", evidence=evidence)],
    )


def _disease_gene(
    link: GeneDiseaseLink, inheritance: Sequence[str], ontology: HpoOntology
) -> DiseaseGene:
    return DiseaseGene(
        symbol=link.gene_symbol,
        ncbi_gene_id=link.ncbi_gene_id or None,
        association_type=ASSOCIATION_TYPES.get(link.association_type, "unknown"),
        inheritance=[ontology.get_label(hpo_id) or hpo_id for hpo_id in inheritance],
        inheritance_hpo_ids=list(inheritance),
    )


def _mapped_ids(disease_id: str, omim_orpha: OmimOrphaMap | None) -> MappedDiseaseIds:
    if disease_id.startswith("OMIM:"):
        orpha_ids = omim_orpha.omim_to_orpha.get(disease_id, ()) if omim_orpha else ()
        return MappedDiseaseIds(
            omim=[disease_id], orpha=orpha_ids[0] if len(orpha_ids) == 1 else None
        )
    if disease_id.startswith("ORPHA:"):
        omim_ids = omim_orpha.orpha_to_omim.get(disease_id, ()) if omim_orpha else ()
        return MappedDiseaseIds(orpha=disease_id, omim=list(omim_ids))
    return MappedDiseaseIds()


def _references(rows: Iterable[HpoaRow]) -> list[str]:
    return sorted({reference for row in rows for reference in row.references})


def _round(value: float) -> float:
    return round(value, _ROUND_DIGITS)


def _disease_sort_key(disease_id: str) -> tuple[str, int, str]:
    prefix = disease_id.split(":", 1)[0]
    match = _ID_NUMBER_RE.search(disease_id)
    return prefix, int(match.group(1)) if match else -1, disease_id


def _input_records(
    inputs: HpoaBuildInputs,
    ontology: HpoOntology,
    hpoa_header: dict[str, str],
    orphanet: OrphanetOnsets | None,
    omim_orpha: OmimOrphaMap | None,
) -> dict[str, dict[str, Any]]:
    entries: list[tuple[str, Path | None, str | None]] = [
        ("hp_json", inputs.hpo_json, ontology.version),
        ("phenotype_hpoa", inputs.phenotype_hpoa, hpoa_header.get("version")),
        ("genes_to_disease", inputs.genes_to_disease, None),
        ("orphanet_ages", inputs.orphanet_ages, orphanet.version if orphanet else None),
        ("omim_orpha_map", inputs.omim_orpha_map, omim_orpha.version if omim_orpha else None),
        ("exclude_pmids", inputs.exclude_pmids, None),
    ]
    return {
        name: {
            "file": path.name,
            "path": str(path),
            "sha256": sha256_file(path),
            "version": version,
        }
        for name, path, version in entries
        if path is not None
    }


def _source_context(records: dict[str, dict[str, Any]]) -> _SourceContext:
    def reference(key: str, name: str, license_name: str | None) -> SourceReference | None:
        record = records.get(key)
        if record is None:
            return None
        return SourceReference(
            name=name,
            url_or_file=record["file"],
            version=record["version"],
            license=license_name,
            sha256=record["sha256"],
        )

    hpoa = reference("phenotype_hpoa", "HPO disease annotations", "HPO license")
    ontology = reference("hp_json", "Human Phenotype Ontology", "HPO license")
    genes = reference("genes_to_disease", "HPO genes_to_disease", "HPO license")
    assert hpoa is not None and ontology is not None and genes is not None
    return _SourceContext(
        hpoa=hpoa,
        ontology=ontology,
        genes=genes,
        orphanet=reference("orphanet_ages", "Orphanet average age of onset", "CC-BY-4.0"),
        alignments=reference("omim_orpha_map", "Orphanet OMIM-ORPHA alignments", "CC-BY-4.0"),
        mask=reference("exclude_pmids", "held-out reference mask", None),
    )


def _record_profile_stats(profile: DiseaseProfile, stats: _Stats) -> None:
    stats.diseases["built"] += 1
    stats.diseases_by_prefix[profile.disease_id.split(":", 1)[0]] += 1
    stats.gene_links += len(profile.genes)
    stats.genes.update(gene.symbol for gene in profile.genes)
    if len(profile.genes) > 1:
        stats.multi_gene_diseases += 1
    for phenotype in profile.phenotypes:
        stats.phenotypes["annotations"] += 1
        stats.unique_terms.add(phenotype.hpo_id)
        stats.frequency_categories[phenotype.frequency] += 1
        if phenotype.onset != "unknown":
            stats.phenotypes["with_onset"] += 1
        if phenotype.sex_restriction is not None:
            stats.phenotypes[f"sex_restricted_{phenotype.sex_restriction}"] += 1
    stats.negatives["annotations"] += len(profile.negative_phenotypes)
    if profile.negative_phenotypes:
        stats.diseases["with_negatives"] += 1
    if profile.age_of_onset is not None:
        stats.onset_categories[profile.age_of_onset.category] += 1
    inheritance = profile.genes[0].inheritance if profile.genes else []
    if inheritance:
        stats.diseases["with_inheritance"] += 1
    stats.inheritance_terms.update(inheritance)
    stats.sex_bias[profile.sex_bias.value if profile.sex_bias else "unset"] += 1
    stats.progression[profile.progression] += 1


def _summary(
    stats: _Stats,
    profiles: Sequence[DiseaseProfile],
    *,
    input_records: dict[str, dict[str, Any]],
    omim_orpha: OmimOrphaMap | None,
    mask_size: int,
    genes_rows_total: int,
    genes_rows_missing_symbol: int,
) -> dict[str, Any]:
    built = len(profiles)
    rows_total = stats.rows["total"]
    with_onset = built - stats.onset_sources["none"]
    via_mapping = stats.onset_sources["orphanet_via_omim"]
    return {
        "inputs": input_records,
        "masking": {
            "enabled": "exclude_pmids" in input_records,
            "references_in_mask": mask_size,
            "rows_masked": stats.rows["masked"],
            "rows_masked_fraction": _round(stats.rows["masked"] / rows_total) if rows_total else 0,
            "diseases_affected": len(stats.masked_diseases),
            "diseases_dropped_zero_positive": stats.diseases["dropped_zero_positive_due_to_mask"],
        },
        "rows": {
            **dict(sorted(stats.rows.items())),
            "by_aspect": dict(sorted(stats.rows_by_aspect.items())),
        },
        "diseases": {
            **dict(sorted(stats.diseases.items())),
            "by_prefix": dict(sorted(stats.diseases_by_prefix.items())),
        },
        "genes": {
            "unique_symbols": len(stats.genes),
            "disease_gene_links": stats.gene_links,
            "diseases_with_multiple_genes": stats.multi_gene_diseases,
            "genes_to_disease_rows": genes_rows_total,
            "genes_to_disease_rows_without_symbol": genes_rows_missing_symbol,
        },
        "phenotypes": {
            **dict(sorted(stats.phenotypes.items())),
            "unique_terms": len(stats.unique_terms),
            "by_frequency": dict(sorted(stats.frequency_categories.items())),
            "by_frequency_basis": dict(sorted(stats.frequency_basis.items())),
        },
        "negatives": dict(sorted(stats.negatives.items())),
        "age_of_onset": {
            "diseases_with_onset": with_onset,
            "fraction_with_onset": _round(with_onset / built) if built else 0,
            "gained_via_omim_mapping": via_mapping,
            "fraction_with_onset_without_mapping": (
                _round((with_onset - via_mapping) / built) if built else 0
            ),
            "by_source": dict(sorted(stats.onset_sources.items())),
            "by_category": dict(sorted(stats.onset_categories.items())),
        },
        "inheritance": {
            "diseases_with_inheritance": stats.diseases["with_inheritance"],
            "by_term": dict(
                sorted(stats.inheritance_terms.items(), key=lambda item: (-item[1], item[0]))
            ),
        },
        "sex_bias": dict(sorted(stats.sex_bias.items())),
        "progression": dict(sorted(stats.progression.items())),
        "omim_orpha_map": omim_orpha.stats if omim_orpha else None,
    }
