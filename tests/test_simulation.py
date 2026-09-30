from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiseaseGene,
    DiseaseProfile,
    NegativePhenotypeAssociation,
    PhenotypeAssociation,
    SexBias,
)
from rare_disease_simulator.simulation.schema import SimulationConfig
from rare_disease_simulator.simulation.simulator import NoiseTerm, simulate_cases
from tests.fixtures.readers import fixture_path, read_disease_profile


def _dump_without_timestamp(cases: list) -> list[dict]:
    dumped = []
    for case in cases:
        data = case.model_dump(mode="json")
        data["metadata"].pop("generated_at", None)
        dumped.append(data)
    return dumped


def _config(**overrides) -> SimulationConfig:
    base = {
        "cases_per_disease_per_difficulty": 5,
        "difficulties": ["easy", "medium", "hard"],
        "seed": 7,
        "reporting": {"mode": "observation"},
    }
    base.update(overrides)
    return SimulationConfig(**base)


def test_simulate_cases_emits_expected_count_and_unique_ids() -> None:
    profile = read_disease_profile()

    cases = simulate_cases(profile, _config())

    assert len(cases) == 5 * 3
    assert len({case.case_id for case in cases}) == len(cases)
    assert {case.metadata.difficulty for case in cases} == {"easy", "medium", "hard"}


def test_simulation_is_deterministic_for_same_seed() -> None:
    profile = read_disease_profile()
    config = _config()

    first = simulate_cases(profile, config)
    second = simulate_cases(profile, config)

    assert _dump_without_timestamp(first) == _dump_without_timestamp(second)


def test_different_seed_changes_output() -> None:
    profile = read_disease_profile()

    first = _dump_without_timestamp(simulate_cases(profile, _config(seed=1)))
    second = _dump_without_timestamp(simulate_cases(profile, _config(seed=2)))

    assert first != second


def test_phenotype_status_partition_is_consistent() -> None:
    profile = read_disease_profile()
    profile_ids = {phenotype.hpo_id for phenotype in profile.phenotypes}

    for case in simulate_cases(profile, _config()):
        positive_ids = {p.hpo_id for p in case.positive_phenotypes}
        negative_ids = {p.hpo_id for p in case.negative_phenotypes}

        assert positive_ids.isdisjoint(negative_ids)
        assert case.patient.sex in {"male", "female", "unknown"}
        for bucket in (
            case.missing_phenotypes,
            case.unknown_phenotypes,
        ):
            for phenotype in bucket:
                assert phenotype.hpo_id in profile_ids


def test_explicit_negatives_and_noise_and_smoothing() -> None:
    profile = DiseaseProfile(
        disease_id="ORPHA:646",
        disease_name="Niemann-Pick disease type C1",
        genes=[DiseaseGene(symbol="NPC1", association_type="causal")],
        phenotypes=[
            PhenotypeAssociation(
                hpo_id="HP:0002066",
                label="Gait ataxia",
                frequency="very_frequent",
                diagnostic_role="cardinal",
            ),
        ],
        negative_phenotypes=[
            NegativePhenotypeAssociation(hpo_id="HP:0001433", label="Hepatosplenomegaly"),
        ],
        age_of_onset=AgeOfOnset(category="childhood"),
        sex_bias=SexBias(value="female"),
    )
    ontology = HpoOntology.from_tsv(fixture_path("hpo_terms.tsv"), version="fixture-0.1")
    noise = [NoiseTerm(hpo_id="HP:0001263", label="Global developmental delay")]

    cases = simulate_cases(
        profile,
        _config(cases_per_disease_per_difficulty=40, difficulties=["hard"]),
        ontology=ontology,
        noise_vocabulary=noise,
    )

    negatives_seen = any(case.negative_phenotypes for case in cases)
    noise_seen = any(case.noise_phenotypes for case in cases)
    generalized_seen = any(
        phenotype.reason == "generalized"
        for case in cases
        for phenotype in case.positive_phenotypes
    )

    assert negatives_seen
    assert noise_seen
    assert generalized_seen
    for case in cases:
        for phenotype in case.noise_phenotypes:
            assert phenotype.hpo_id == "HP:0001263"
