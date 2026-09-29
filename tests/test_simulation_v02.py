from pathlib import Path

import pytest

from rare_disease_simulator.data_sources.hpo import HpoOntology, HpoTerm
from rare_disease_simulator.exports.jsonl import write_jsonl
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiseaseGene,
    DiseaseProfile,
    NegativePhenotypeAssociation,
    PhenotypeAssociation,
    SexBias,
)
from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.difficulty import DifficultyPreset
from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import (
    FORCED_REASON,
    MERGED_REASON,
    NoiseTerm,
    sex_prior_key,
    simulate_cases,
)

TERMS = {
    "HP:0000118": ("Phenotypic abnormality", ()),
    "HP:0000707": ("Abnormality of the nervous system", ("HP:0000118",)),
    "HP:0001250": ("Seizure", ("HP:0000707",)),
    "HP:0007359": ("Focal-onset seizure", ("HP:0001250",)),
    "HP:0001251": ("Ataxia", ("HP:0000707",)),
    "HP:0001249": ("Intellectual disability", ("HP:0000707",)),
    "HP:0000478": ("Abnormality of the eye", ("HP:0000118",)),
    "HP:0000505": ("Visual impairment", ("HP:0000478",)),
    "HP:0000518": ("Cataract", ("HP:0000478",)),
    "HP:0000078": ("Abnormality of the genital system", ("HP:0000118",)),
    "HP:0010461": ("Abnormality of the male genitalia", ("HP:0000078",)),
    "HP:0000047": ("Hypospadias", ("HP:0010461",)),
    "HP:0000028": ("Cryptorchidism", ("HP:0010461",)),
    "HP:0010460": ("Abnormality of the female genitalia", ("HP:0000078",)),
    "HP:0000008": ("Abnormal morphology of female internal genitalia", ("HP:0010460",)),
    "HP:0001627": ("Abnormal heart morphology", ("HP:0000118",)),
    "HP:0001629": ("Ventricular septal defect", ("HP:0001627",)),
    "HP:0000924": ("Abnormality of the skeletal system", ("HP:0000118",)),
    "HP:0002650": ("Scoliosis", ("HP:0000924",)),
    "HP:0001385": ("Hip dysplasia", ("HP:0000924",)),
    "HP:0001574": ("Abnormality of the integument", ("HP:0000118",)),
    "HP:0000988": ("Skin rash", ("HP:0001574",)),
    "HP:0000964": ("Eczema", ("HP:0001574",)),
}


def build_ontology() -> HpoOntology:
    return HpoOntology(
        {
            hpo_id: HpoTerm(
                hpo_id=hpo_id, label=label, parents=parents, is_phenotypic_abnormality=True
            )
            for hpo_id, (label, parents) in TERMS.items()
        },
        version="test",
    )


@pytest.fixture(scope="module")
def ontology() -> HpoOntology:
    return build_ontology()


def _phenotype(hpo_id: str, estimate: float | None, **extra) -> PhenotypeAssociation:
    return PhenotypeAssociation(
        hpo_id=hpo_id, label=TERMS[hpo_id][0], frequency_estimate=estimate, **extra
    )


def _profile(disease_id: str, phenotypes, **extra) -> DiseaseProfile:
    return DiseaseProfile(
        disease_id=disease_id,
        disease_name=f"{disease_id} disease",
        genes=[DiseaseGene(symbol=f"G{disease_id[-1]}", association_type="causal")],
        phenotypes=phenotypes,
        **extra,
    )


TRUE = _profile(
    "OMIM:1",
    [
        _phenotype("HP:0001250", 0.9),
        _phenotype("HP:0001251", 0.5),
        _phenotype("HP:0001249", None),
        _phenotype("HP:0000047", 0.8, sex_restriction="male"),
        _phenotype("HP:0000505", 0.6, onset="adult", onset_hpo_id="HP:0003581"),
        _phenotype("HP:0001627", 1.0, onset="neonatal", onset_hpo_id="HP:0003577"),
    ],
    negative_phenotypes=[NegativePhenotypeAssociation(hpo_id="HP:0002650", label="Scoliosis")],
    age_of_onset=AgeOfOnset(
        category="infantile", distribution={"infantile": 0.5, "childhood": 0.5}
    ),
    sex_bias=SexBias(value="none"),
    progression="progressive",
)
CONFOUNDER = _profile(
    "OMIM:2",
    [
        _phenotype("HP:0001250", 0.9),
        _phenotype("HP:0001251", 0.9),
        _phenotype("HP:0007359", 0.9),
        _phenotype("HP:0001629", 0.9),
        _phenotype("HP:0000518", 0.9),
        _phenotype("HP:0001385", 0.7),
    ],
)
UNRELATED = _profile("OMIM:3", [_phenotype("HP:0000988", 0.5), _phenotype("HP:0002650", 0.5)])
PROFILES = [TRUE, CONFOUNDER, UNRELATED]
NOISE = [
    NoiseTerm("HP:0000988", "Skin rash"),
    NoiseTerm("HP:0000964", "Eczema", weight=0.0),
    NoiseTerm("HP:0001385", "Hip dysplasia"),
]


def _config(**overrides) -> SimulationConfig:
    base = {
        "cases_per_disease_per_difficulty": 60,
        "difficulties": ["easy", "medium", "hard"],
        "seed": 11,
        "missingness": {"sex_unknown": 0.0, "age_unknown": 0.0, "onset_unknown": 0.0,
                        "no_negatives": 0.0},
    }
    base.update(overrides)
    return SimulationConfig.model_validate(base)


def _preset(**overrides) -> DifficultyPreset:
    base = {
        "positive_observation_rate": 1.0,
        "cardinal_observation_boost": 0.0,
        "missing_vs_unknown_split": 0.5,
        "negatives_mean": 0.0,
        "noise_mean": 0.0,
        "ontology_smoothing_rate": 0.0,
    }
    base.update(overrides)
    return DifficultyPreset(**base)


def _simulate(ontology, config=None, profile=TRUE, **kwargs):
    config = config or _config()
    index = ConfounderIndex(
        PROFILES, ontology, top_n=config.negatives.confounders_top_n, min_information_content=0.0
    )
    return simulate_cases(
        profile, config, ontology=ontology, noise_vocabulary=NOISE, confounders=index, **kwargs
    )


def _present_ids(case) -> set[str]:
    ids = set()
    for phenotype in case.positive_phenotypes:
        ids.add(phenotype.hpo_id)
        if phenotype.source_hpo_id:
            ids.add(phenotype.source_hpo_id)
    for bucket in (case.missing_phenotypes, case.unknown_phenotypes, case.noise_phenotypes):
        ids.update(phenotype.hpo_id for phenotype in bucket)
    return ids


def _related(ontology, a: str, b: str) -> bool:
    return a == b or ontology.is_a(a, b) or ontology.is_a(b, a)


def test_same_seed_and_config_give_identical_bytes(ontology, tmp_path: Path) -> None:
    first, second = tmp_path / "a.jsonl", tmp_path / "b.jsonl"

    write_jsonl(first, _simulate(ontology, source_versions={"profiles": "abc"}))
    write_jsonl(second, _simulate(ontology, source_versions={"profiles": "abc"}))

    assert first.read_bytes() == second.read_bytes()
    assert _simulate(ontology, _config(seed=12))[0] != _simulate(ontology)[0]


def test_metadata_carries_versions_and_no_timestamp(ontology) -> None:
    case = _simulate(ontology, source_versions={"profiles_sha256": "abc"})[0]

    assert case.metadata.simulator_version == "0.3.0"
    assert case.metadata.source_versions == {"profiles_sha256": "abc"}
    assert case.metadata.generated_at is None
    assert case.metadata.config_hash.startswith("sha256:")


def test_every_case_has_an_observed_positive(ontology) -> None:
    rare = _profile("OMIM:9", [_phenotype("HP:0001250", 0.001), _phenotype("HP:0001251", 0.001)])
    config = _config(max_redraws=3)

    cases = simulate_cases(rare, config, ontology=ontology)

    assert all(case.positive_phenotypes for case in cases)
    assert any(p.reason == FORCED_REASON for case in cases for p in case.positive_phenotypes)


def test_negatives_never_touch_present_terms_or_true_annotations(ontology) -> None:
    config = _config(presets={"hard": _preset(negatives_mean=4.0, noise_mean=2.0,
                                              positive_observation_rate=0.6,
                                              ontology_smoothing_rate=0.5)})
    cases = _simulate(ontology, config)
    closure = ConfounderIndex(PROFILES, ontology).annotation_closure(TRUE)

    origins = set()
    for case in cases:
        present = _present_ids(case)
        negatives = [p.hpo_id for p in case.negative_phenotypes]
        assert len(negatives) == len(set(negatives))
        for negative in case.negative_phenotypes:
            origins.add(negative.simulated_origin)
            assert not any(_related(ontology, negative.hpo_id, p) for p in present)
            if negative.simulated_origin == "negative_confounder":
                assert negative.hpo_id not in closure
                assert negative.reason == "confounder:OMIM:2"
    assert origins == {"negative_own_disease", "negative_confounder", "negative_not_annotation"}


def test_confounder_terms_skip_anything_the_true_disease_annotates(ontology) -> None:
    index = ConfounderIndex(PROFILES, ontology, min_information_content=0.0)

    assert index.confounders("OMIM:1")[0][0] == "OMIM:2"
    assert [term.hpo_id for term in index.candidate_terms("OMIM:1")] == [
        "HP:0000518",
        "HP:0001385",
    ]


def test_high_information_threshold_keeps_general_terms_out_of_the_index(ontology) -> None:
    index = ConfounderIndex(PROFILES, ontology, min_information_content=10.0)

    assert index.confounders("OMIM:1") == []


def test_sex_restricted_terms_follow_the_patient_sex(ontology) -> None:
    male_only_by_anchor = _profile(
        "OMIM:4",
        [_phenotype("HP:0001250", 0.9), _phenotype("HP:0000028", 0.9),
         _phenotype("HP:0000008", 0.9)],
    )
    for profile in (TRUE, male_only_by_anchor):
        for case in _simulate(ontology, profile=profile):
            present = _present_ids(case) | {p.hpo_id for p in case.negative_phenotypes}
            if case.patient.sex == "female":
                assert present.isdisjoint({"HP:0000047", "HP:0000028"})
            else:
                assert "HP:0000008" not in present


def test_sex_prior_follows_inheritance(ontology) -> None:
    limited = TRUE.model_copy(
        update={"genes": [DiseaseGene(symbol="G1", inheritance_hpo_ids=["HP:0001475"])]}
    )
    xld = TRUE.model_copy(
        update={
            "genes": [DiseaseGene(symbol="G1", inheritance_hpo_ids=["HP:0001423"])],
            "sex_bias": None,
        }
    )

    assert sex_prior_key(limited) == "male_limited"
    assert sex_prior_key(xld) == "x_linked_dominant"
    assert sex_prior_key(TRUE) == "unbiased"
    assert {case.patient.sex for case in _simulate(ontology, profile=limited)} == {"male"}
    females = [c for c in _simulate(ontology, profile=xld) if c.patient.sex == "female"]
    assert len(females) > 90


def test_ages_follow_onset_and_gate_late_onset_terms(ontology) -> None:
    cases = _simulate(ontology)

    for case in cases:
        age, onset = case.patient.age.value, case.patient.age_of_onset.value
        assert age >= onset
        assert case.patient.onset_category in {"infantile", "childhood"}
        assert 0.0767 <= onset <= 5.0
        if age < 16.0:
            assert "HP:0000505" not in _present_ids(case)
        assert "HP:0001627" in _present_ids(case) or any(
            p.source_hpo_id == "HP:0001627" for p in case.positive_phenotypes
        )


def test_unknown_frequency_uses_the_configured_prior(ontology) -> None:
    profile = _profile("OMIM:5", [_phenotype("HP:0001627", 1.0), _phenotype("HP:0001249", None)])
    config = _config(
        cases_per_disease_per_difficulty=800,
        difficulties=["easy"],
        presets={"easy": _preset()},
        frequency={"unknown_frequency": 0.2},
    )

    cases = simulate_cases(profile, config)
    rate = sum("HP:0001249" in _present_ids(case) for case in cases) / len(cases)

    assert rate == pytest.approx(0.2, abs=0.04)


def test_patient_probability_is_drawn_around_the_estimate(ontology) -> None:
    profile = _profile("OMIM:6", [_phenotype("HP:0001627", 1.0), _phenotype("HP:0001251", 0.5)])
    config = _config(difficulties=["easy"], presets={"easy": _preset()})

    probabilities = {
        p.source_probability
        for case in simulate_cases(profile, config)
        for p in case.positive_phenotypes
        if p.hpo_id == "HP:0001251"
    }

    assert len(probabilities) > 5
    assert {p.source_probability for case in simulate_cases(profile, config)
            for p in case.positive_phenotypes if p.hpo_id == "HP:0001627"} == {1.0}


def test_progression_raises_probability_with_duration(ontology) -> None:
    profile = _profile(
        "OMIM:7",
        [_phenotype("HP:0001627", 1.0), _phenotype("HP:0001251", 0.2)],
        progression="progressive",
    )

    def rate(max_boost: float) -> float:
        config = _config(
            cases_per_disease_per_difficulty=400,
            difficulties=["easy"],
            presets={"easy": _preset()},
            progression={"max_boost": max_boost, "timescale_years": 1.0},
        )
        cases = simulate_cases(profile, config)
        return sum("HP:0001251" in _present_ids(case) for case in cases) / len(cases)

    assert rate(0.8) > rate(0.0) + 0.3


def test_covariate_missingness(ontology) -> None:
    config = _config(
        missingness={"sex_unknown": 1.0, "age_unknown": 1.0, "onset_unknown": 1.0,
                     "no_negatives": 1.0},
    )

    for case in _simulate(ontology, config):
        assert case.patient.sex == "unknown"
        assert case.patient.age is None and case.patient.age_of_onset is None
        assert case.patient.onset_category == "unknown"
        assert case.negative_phenotypes == []


def test_noise_comes_only_from_weighted_vocabulary(ontology) -> None:
    config = _config(presets={"hard": _preset(noise_mean=3.0, negatives_mean=0.0)})

    noise = {p.hpo_id for case in _simulate(ontology, config) for p in case.noise_phenotypes}

    assert noise == {"HP:0000988", "HP:0001385"}


def test_negative_sources_can_be_disabled(ontology) -> None:
    config = _config(
        presets={"easy": _preset(negatives_mean=5.0)},
        difficulties=["easy"],
        negatives={"source_weights": {"not_annotation": 1.0}},
    )

    origins = {
        p.simulated_origin for case in _simulate(ontology, config) for p in case.negative_phenotypes
    }

    assert origins == {"negative_not_annotation"}


def test_negative_count_is_overdispersed_and_capped() -> None:
    import random

    from rare_disease_simulator.simulation.simulator import negative_count

    rng = random.Random(3)
    counts = [negative_count(7.0, 2.0, 15, rng) for _ in range(4000)]

    assert max(counts) == 15
    assert counts.count(0) > 100
    assert sum(counts) / len(counts) == pytest.approx(6.5, abs=0.5)
    assert negative_count(7.0, None, 0, rng) == 0


def test_default_negative_mix_is_dominated_by_the_disease_own_terms(ontology) -> None:
    rich = _profile(
        "OMIM:8",
        [
            _phenotype(hpo_id, 0.3)
            for hpo_id in TERMS
            if not any(parents == (hpo_id,) for _, parents in TERMS.values())
            and hpo_id != "HP:0000964"
        ],
        negative_phenotypes=[NegativePhenotypeAssociation(hpo_id="HP:0000964", label="Eczema")],
    )
    config = SimulationConfig(cases_per_disease_per_difficulty=100, difficulties=["easy"])
    assert config.negatives.max_per_case == 15
    index = ConfounderIndex([rich, CONFOUNDER, UNRELATED], ontology, min_information_content=0.0)

    cases = simulate_cases(rich, config, ontology=ontology, confounders=index)
    origins = [p.simulated_origin for case in cases for p in case.negative_phenotypes]

    assert origins.count("negative_own_disease") > len(origins) / 2


def test_terms_merged_by_generalization_stay_in_the_case(ontology) -> None:
    profile = _profile(
        "OMIM:7", [_phenotype("HP:0001250", 1.0), _phenotype("HP:0001251", 1.0)]
    )
    config = _config(
        difficulties=["easy"],
        cases_per_disease_per_difficulty=20,
        presets={"easy": _preset(ontology_smoothing_rate=1.0)},
    )

    for case in simulate_cases(profile, config, ontology=ontology):
        assert [p.hpo_id for p in case.positive_phenotypes] == ["HP:0000707"]
        assert [(m.hpo_id, m.reason) for m in case.missing_phenotypes] == [
            ("HP:0001251", MERGED_REASON)
        ]
