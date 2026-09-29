from rare_disease_simulator.profiles.schema import (
    AgeOfOnset,
    DiseaseGene,
    DiseaseProfile,
    PhenotypeAssociation,
    Provenance,
    SourceReference,
)
from rare_disease_simulator.simulation.sampling import sample_diseases, stratum

HPOA = SourceReference(name="HPO disease annotations")


def _profile(number: int, inheritance: list[str], with_onset: bool) -> DiseaseProfile:
    return DiseaseProfile(
        disease_id=f"OMIM:{number}",
        disease_name=str(number),
        genes=[DiseaseGene(symbol="G", inheritance_hpo_ids=inheritance)],
        phenotypes=[PhenotypeAssociation(hpo_id="HP:0001250", label="Seizure")],
        age_of_onset=(
            AgeOfOnset(category="infantile", provenance=[Provenance(source=HPOA)])
            if with_onset
            else None
        ),
    )


PROFILES = (
    [_profile(n, ["HP:0000007"], True) for n in range(100)]
    + [_profile(100 + n, ["HP:0000006"], False) for n in range(50)]
    + [_profile(200 + n, ["HP:0001419"], True) for n in range(3)]
)


def test_strata_combine_inheritance_class_and_onset_source() -> None:
    assert stratum(PROFILES[0]) == "autosomal_recessive|hpoa"
    assert stratum(PROFILES[100]) == "autosomal_dominant|none"
    assert stratum(PROFILES[-1]) == "x_linked|hpoa"


def test_sampling_spreads_evenly_and_is_deterministic() -> None:
    selected, quotas = sample_diseases(PROFILES, 30, seed=1)

    assert quotas == {
        "autosomal_dominant|none": 14,
        "autosomal_recessive|hpoa": 13,
        "x_linked|hpoa": 3,
    }
    assert len(selected) == 30
    assert [p.disease_id for p in selected] == [
        p.disease_id for p in sample_diseases(PROFILES, 30, seed=1)[0]
    ]
    assert [p.disease_id for p in selected] != [
        p.disease_id for p in sample_diseases(PROFILES, 30, seed=2)[0]
    ]


def test_sampling_more_than_available_returns_everything() -> None:
    selected, _ = sample_diseases(PROFILES, 1000, seed=1)

    assert len(selected) == len(PROFILES)
