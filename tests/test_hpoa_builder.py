from dataclasses import replace
from pathlib import Path

import pytest

from rare_disease_simulator.data_sources.hpo_annotations import read_reference_mask
from rare_disease_simulator.data_sources.hpo_genes import read_genes_to_disease
from rare_disease_simulator.data_sources.orphanet_products import (
    read_omim_orpha_map,
    read_orphanet_onsets,
)
from rare_disease_simulator.profiles.hpoa_builder import (
    HpoaBuildInputs,
    HpoaBuildResult,
    build_profiles_from_hpoa,
)
from rare_disease_simulator.profiles.schema import DiseaseProfile
from tests.fixtures.readers import fixture_path

FIXTURES = fixture_path("hpoa")


def _inputs(**overrides: Path | None) -> HpoaBuildInputs:
    inputs = HpoaBuildInputs(
        hpo_json=FIXTURES / "hp_mini.json",
        phenotype_hpoa=FIXTURES / "phenotype_mini.hpoa",
        genes_to_disease=FIXTURES / "genes_to_disease_mini.txt",
        orphanet_ages=FIXTURES / "orphanet_ages_mini.xml",
        omim_orpha_map=FIXTURES / "omim_orpha_mini.xml",
        exclude_pmids=FIXTURES / "heldout_pmids.txt",
    )
    return replace(inputs, **overrides)


@pytest.fixture(scope="module")
def result() -> HpoaBuildResult:
    return build_profiles_from_hpoa(_inputs())


def _profile(result: HpoaBuildResult, disease_id: str) -> DiseaseProfile:
    return next(profile for profile in result.profiles if profile.disease_id == disease_id)


def _phenotype(profile: DiseaseProfile, hpo_id: str):
    return next(phenotype for phenotype in profile.phenotypes if phenotype.hpo_id == hpo_id)


def test_builds_one_profile_per_gene_linked_disease_in_sorted_order(result) -> None:
    assert [profile.disease_id for profile in result.profiles] == [
        "OMIM:100001",
        "OMIM:100002",
        "OMIM:100008",
        "ORPHA:3001",
    ]


def test_genes_keep_symbol_ncbi_id_and_association_type(result) -> None:
    alpha = _profile(result, "OMIM:100001")
    assert [(g.symbol, g.ncbi_gene_id, g.association_type) for g in alpha.genes] == [
        ("GENEA", "NCBIGene:1", "causal")
    ]
    gamma = _profile(result, "ORPHA:3001")
    assert [gene.symbol for gene in gamma.genes] == ["GENEC", "GENED"]
    assert {gene.association_type for gene in gamma.genes} == {"unknown"}


def test_frequency_terms_map_to_hpo_range_midpoints(result) -> None:
    alpha = _profile(result, "OMIM:100001")
    seizure = _phenotype(alpha, "HP:0001250")
    assert seizure.frequency == "frequent"
    assert seizure.frequency_estimate == 0.545
    assert (seizure.probability_range.lower, seizure.probability_range.upper) == (0.30, 0.79)

    gamma = _profile(result, "ORPHA:3001")
    heart = _phenotype(gamma, "HP:0001627")
    assert heart.frequency == "obligate"
    assert heart.frequency_estimate == 1.0
    assert (heart.probability_range.lower, heart.probability_range.upper) == (1.0, 1.0)


def test_ratio_and_percent_frequencies_become_point_estimates(result) -> None:
    alpha = _profile(result, "OMIM:100001")
    vision = _phenotype(alpha, "HP:0000505")
    assert vision.frequency_raw == "20%"
    assert vision.frequency == "occasional"
    assert vision.frequency_estimate == 0.2
    assert vision.probability_range.lower == vision.probability_range.upper == 0.2


def test_duplicate_rows_pool_ratio_counts_across_references(result) -> None:
    ataxia = _phenotype(_profile(result, "OMIM:100001"), "HP:0001251")

    assert ataxia.frequency_raw == "4/8"
    assert ataxia.frequency_estimate == 0.5
    lower, upper = ataxia.probability_range.lower, ataxia.probability_range.upper
    assert (lower, upper) == (0.199, 0.801)
    assert ataxia.source == ["OMIM:100001", "PMID:2", "PMID:3"]


def test_small_counts_are_shrunk_so_n_of_n_is_not_obligate(result) -> None:
    hypospadias = _phenotype(_profile(result, "OMIM:100002"), "HP:0000047")

    assert hypospadias.frequency_raw == "2/3"
    assert hypospadias.frequency_estimate == 0.625
    assert hypospadias.frequency == "frequent"
    assert result.summary["phenotypes"]["by_frequency_basis"] == {
        "category": 4,
        "counts": 6,
        "percent": 1,
    }


def test_duplicate_frequency_terms_use_the_envelope_of_their_ranges(result) -> None:
    ataxia = _phenotype(_profile(result, "ORPHA:3001"), "HP:0001251")

    assert ataxia.frequency_raw == "HP:0040281;HP:0040283"
    assert ataxia.frequency_estimate == 0.5325
    assert (ataxia.probability_range.lower, ataxia.probability_range.upper) == (0.05, 0.99)


def test_missing_frequency_stays_unknown_not_always(result) -> None:
    disability = _phenotype(_profile(result, "OMIM:100001"), "HP:0001249")

    assert disability.frequency == "unknown"
    assert disability.frequency_estimate is None
    assert disability.probability_range is None
    assert disability.frequency_raw is None


def test_negatives_come_only_from_not_and_excluded_frequency(result) -> None:
    alpha = _profile(result, "OMIM:100001")

    assert [negative.hpo_id for negative in alpha.negative_phenotypes] == [
        "HP:0000047",
        "HP:0001627",
    ]
    positive_ids = {phenotype.hpo_id for phenotype in alpha.phenotypes}
    assert positive_ids.isdisjoint({"HP:0000047", "HP:0001627"})
    assert result.summary["negatives"] == {
        "annotations": 2,
        "dropped_conflict_with_positive": 1,
        "from_excluded_frequency": 1,
        "from_not_qualifier": 1,
    }


def test_zero_count_is_a_low_frequency_not_a_negative(result) -> None:
    ptosis = _phenotype(_profile(result, "OMIM:100001"), "HP:0000140")

    assert ptosis.frequency_raw == "0/5"
    assert ptosis.frequency_estimate == 0.0833
    assert ptosis.frequency == "occasional"
    assert ptosis.probability_range.lower == 0.0


def test_positive_evidence_wins_over_a_conflicting_not_row(result) -> None:
    theta = _profile(result, "OMIM:100008")

    assert "HP:0001249" in {phenotype.hpo_id for phenotype in theta.phenotypes}
    assert theta.negative_phenotypes == []
    assert theta.quality.counters["positive_negative_conflicts"] == 1


def test_per_row_onset_maps_to_phenotype_onset_earliest_wins(result) -> None:
    beta = _profile(result, "OMIM:100002")

    hypospadias = _phenotype(beta, "HP:0000047")
    assert (hypospadias.onset, hypospadias.onset_hpo_id) == ("neonatal", "HP:0003577")
    seizure = _phenotype(beta, "HP:0001250")
    assert (seizure.onset, seizure.onset_hpo_id) == ("infantile", "HP:0003593")
    assert _phenotype(_profile(result, "OMIM:100001"), "HP:0000505").onset == "infantile"


def test_sex_column_becomes_phenotype_sex_restriction(result) -> None:
    beta = _profile(result, "OMIM:100002")

    assert _phenotype(beta, "HP:0000047").sex_restriction == "male"
    assert _phenotype(beta, "HP:0000140").sex_restriction == "female"
    assert _phenotype(beta, "HP:0001250").sex_restriction is None
    assert result.summary["phenotypes"]["sex_restriction_conflicts"] == 1


def test_clinical_course_onset_sets_disease_age_of_onset(result) -> None:
    onset = _profile(result, "OMIM:100001").age_of_onset

    assert onset is not None
    assert onset.category == "infantile"
    assert onset.distribution == {"infantile": 0.5, "childhood": 0.5}
    assert onset.hpo_ids == ["HP:0003593", "HP:0011463"]
    assert onset.provenance[0].source.name == "HPO disease annotations"


def test_orphanet_onset_fallback_direct_and_via_omim_mapping(result) -> None:
    gamma = _profile(result, "ORPHA:3001").age_of_onset
    assert gamma is not None
    assert gamma.category == "neonatal"
    assert gamma.distribution == {"neonatal": 0.5, "infantile": 0.5}
    assert gamma.provenance[0].source.name == "Orphanet average age of onset"

    beta = _profile(result, "OMIM:100002").age_of_onset
    assert beta is not None
    assert (beta.category, beta.distribution) == ("adult", {"adult": 1.0})
    assert "ORPHA:2002" in (beta.provenance[0].evidence or "")
    assert result.summary["age_of_onset"]["by_source"] == {
        "hpoa": 1,
        "none": 1,
        "orphanet_direct": 1,
        "orphanet_via_omim": 1,
    }


def test_without_omim_orpha_map_only_orpha_diseases_get_orphanet_onset() -> None:
    result = build_profiles_from_hpoa(_inputs(omim_orpha_map=None))

    assert _profile(result, "OMIM:100002").age_of_onset is None
    assert _profile(result, "ORPHA:3001").age_of_onset is not None


def test_progression_from_clinical_course(result) -> None:
    assert _profile(result, "OMIM:100001").progression == "progressive"
    assert _profile(result, "OMIM:100002").progression == "unknown"


def test_inheritance_terms_are_attached_to_every_gene_with_labels(result) -> None:
    gamma = _profile(result, "ORPHA:3001")

    for gene in gamma.genes:
        assert gene.inheritance_hpo_ids == ["HP:0001423"]
        assert gene.inheritance == ["X-linked dominant inheritance"]


def test_sex_bias_is_derived_from_inheritance_only_where_sound(result) -> None:
    assert _profile(result, "OMIM:100001").sex_bias.value == "none"
    assert _profile(result, "OMIM:100002").sex_bias.value == "male"
    theta = _profile(result, "OMIM:100008").sex_bias
    assert theta.value == "male"
    assert "Male-limited expression" in theta.provenance[0].evidence
    assert _profile(result, "ORPHA:3001").sex_bias is None


def test_alt_and_obsolete_ids_are_resolved_and_unresolvable_rows_dropped(result) -> None:
    alpha = _profile(result, "OMIM:100001")
    ids = {phenotype.hpo_id for phenotype in alpha.phenotypes}

    assert "HP:0000505" in ids
    assert "HP:0000504" not in ids
    assert "HP:0000999" not in ids
    assert "HP:0000489" not in ids
    assert alpha.quality.counters["rows_term_alt_id"] == 1
    assert alpha.quality.counters["rows_term_replaced"] == 1
    assert alpha.quality.counters["rows_unresolved_term"] == 1
    assert result.summary["rows"]["unresolved_term"] == 1


def test_pmid_mask_drops_rows_and_reports_affected_diseases(result) -> None:
    masking = result.summary["masking"]

    assert masking["rows_masked"] == 3
    assert masking["diseases_affected"] == 2
    assert masking["diseases_dropped_zero_positive"] == 1
    assert "HP:0000478" not in {p.hpo_id for p in _profile(result, "OMIM:100001").phenotypes}
    assert "OMIM:100004" not in {profile.disease_id for profile in result.profiles}
    assert len(result.summary["inputs"]["exclude_pmids"]["sha256"]) == 64


def test_without_mask_the_masked_rows_are_kept() -> None:
    result = build_profiles_from_hpoa(_inputs(exclude_pmids=None))

    assert "OMIM:100004" in {profile.disease_id for profile in result.profiles}
    alpha = _profile(result, "OMIM:100001")
    assert "HP:0000478" in {phenotype.hpo_id for phenotype in alpha.phenotypes}
    assert result.summary["masking"]["enabled"] is False


def test_diseases_without_positive_phenotypes_are_excluded_and_counted(result) -> None:
    ids = {profile.disease_id for profile in result.profiles}

    assert "OMIM:100005" not in ids
    assert result.summary["diseases"]["dropped_zero_positive"] == 2


def test_diseases_without_usable_gene_are_skipped(result) -> None:
    ids = {profile.disease_id for profile in result.profiles}

    assert ids.isdisjoint({"OMIM:100006", "OMIM:100007", "DECIPHER:1"})
    assert result.summary["genes"]["genes_to_disease_rows_without_symbol"] == 1
    assert result.summary["diseases"]["gene_linked_not_in_hpoa"] == 1


def test_profiles_record_source_versions_and_checksums(result) -> None:
    sources = {item.source.name: item.source for item in result.profiles[0].provenance}

    assert set(sources) == {
        "HPO disease annotations",
        "HPO genes_to_disease",
        "Human Phenotype Ontology",
        "held-out reference mask",
    }
    assert sources["HPO disease annotations"].version == "2026-02-16"
    assert sources["Human Phenotype Ontology"].version == "2026-02-16"
    assert all(source.sha256 and len(source.sha256) == 64 for source in sources.values())


def test_other_aspects_are_skipped_and_counted(result) -> None:
    rows = result.summary["rows"]

    assert rows["skipped_aspect_H"] == 1
    assert rows["skipped_aspect_M"] == 1


def test_reference_mask_accepts_prefixed_and_bare_pmids() -> None:
    assert read_reference_mask(FIXTURES / "heldout_pmids.txt") == {"PMID:900", "PMID:901"}


def test_genes_to_disease_reader_skips_dash_and_deduplicates() -> None:
    genes = read_genes_to_disease(FIXTURES / "genes_to_disease_mini.txt")

    assert "OMIM:100007" not in genes.links
    assert [link.gene_symbol for link in genes.links["ORPHA:3001"]] == ["GENEC", "GENED"]


def test_orphanet_readers_keep_exact_validated_mappings_and_real_onsets() -> None:
    onsets = read_orphanet_onsets(FIXTURES / "orphanet_ages_mini.xml")
    assert onsets.onsets == {
        "ORPHA:2002": ("Adult", "Elderly"),
        "ORPHA:3001": ("Infancy", "Neonatal"),
    }

    mapping = read_omim_orpha_map(FIXTURES / "omim_orpha_mini.xml")
    assert mapping.omim_to_orpha == {"OMIM:100002": ("ORPHA:2002",)}
    assert mapping.orpha_to_omim == {"ORPHA:2002": ("OMIM:100002",)}
