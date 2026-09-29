import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from rare_disease_simulator.cli import app
from rare_disease_simulator.exports.jsonl import write_jsonl
from rare_disease_simulator.profiles.schema import AgeOfOnset
from rare_disease_simulator.simulation.schema import Age, CasePhenotype
from rare_disease_simulator.validation.cases import validate_cases
from tests.fixtures.readers import fixture_path
from tests.test_simulation_v02 import PROFILES, TRUE, _config, _simulate, build_ontology

FIXTURES = fixture_path("hpoa")
CONFIG = ("--config", str(FIXTURES / "config.yaml"))


@pytest.fixture(scope="module")
def ontology():
    return build_ontology()


def _report(cases, ontology, **kwargs):
    return validate_cases(
        cases,
        profiles={profile.disease_id: profile for profile in PROFILES},
        ontology=ontology,
        config=_config(),
        **kwargs,
    )


def _phenotype(hpo_id: str, status: str, origin: str) -> CasePhenotype:
    return CasePhenotype(hpo_id=hpo_id, label=hpo_id, status=status, simulated_origin=origin)


def test_simulated_cases_pass_and_are_reported(ontology) -> None:
    report = _report(_simulate(ontology), ontology, noise_vocabulary={"HP:0000988", "HP:0001385"})

    assert report["violations"]["count"] == 0
    assert report["cases"] == 180
    assert report["config_matches_cases"] is True
    assert report["phenotypes_per_case"]["positive"]["min"] >= 1
    assert set(report["negatives"]["by_source"]) == {
        "negative_confounder",
        "negative_not_annotation",
        "negative_own_disease",
    }
    assert report["sex"]["by_prior"]["unbiased"]["expected_male_share"] == 0.5
    assert set(report["onset"]["by_category"]) == {"infantile", "childhood"}
    assert report["onset"]["by_category"]["infantile"]["expected_share"] == 0.5
    calibration = report["calibration"]
    assert calibration["all_terms"]["pairs"] == 6
    assert calibration["ungated_terms"]["pairs"] == 0


def test_each_invariant_violation_is_detected(ontology) -> None:
    case = _simulate(ontology)[0]
    positive = case.positive_phenotypes[0].hpo_id
    broken = case.model_copy(
        update={
            "negative_phenotypes": [
                _phenotype(positive, "negative", "negative_own_disease"),
                _phenotype("HP:0007359", "negative", "negative_confounder"),
                _phenotype("HP:0000047", "negative", "negative_own_disease"),
            ],
            "noise_phenotypes": [_phenotype("HP:0000964", "noise", "noise")],
            "patient": case.patient.model_copy(
                update={"sex": "female", "age": Age(value=1.0, unit="years"),
                        "age_of_onset": Age(value=2.0, unit="years")}
            ),
        }
    )
    empty = case.model_copy(update={"positive_phenotypes": [], "case_id": "empty"})

    report = _report([broken, empty], ontology, noise_vocabulary={"HP:0000988"})

    assert report["violations"]["by_type"] == {
        "age_before_onset": 1,
        "confounder_negative_annotated_to_disease": 1,
        "negative_related_to_present": 2,
        "no_observed_positive": 1,
        "noise_not_in_vocabulary": 1,
        "sex_restricted_term_for_other_sex": 1,
    }


def test_onset_expectation_uses_the_profile_distribution(ontology) -> None:
    adult = TRUE.model_copy(update={"age_of_onset": AgeOfOnset(category="adult")})
    cases = _simulate(ontology, profile=adult)
    report = validate_cases(cases, profiles={"OMIM:1": adult}, ontology=ontology, config=_config())

    assert report["onset"]["by_category"] == {
        "adult": {"observed_share": 1.0, "expected_share": 1.0}
    }


def _run(*args: str):
    return CliRunner().invoke(app, [*CONFIG, *args])


def test_build_simulate_validate_end_to_end(tmp_path: Path) -> None:
    profiles, cases = tmp_path / "profiles.jsonl", tmp_path / "cases.jsonl"
    hp = ("--hpo-json", str(FIXTURES / "hp_mini.json"))
    noise = ("--noise-vocabulary", str(FIXTURES / "noise_vocabulary.tsv"))

    built = _run(
        "build-profiles", "--from-hpoa", *hp,
        "--hpoa", str(FIXTURES / "phenotype_mini.hpoa"),
        "--genes-to-disease", str(FIXTURES / "genes_to_disease_mini.txt"),
        "--orphanet-ages", str(FIXTURES / "orphanet_ages_mini.xml"),
        "--omim-orpha-map", str(FIXTURES / "en_product1_mini.xml"),
        "--exclude-pmids", str(FIXTURES / "heldout_pmids.txt"),
        "--output", str(profiles),
    )
    assert built.exit_code == 0, built.output
    simulated = _run(
        "simulate", "--profiles", str(profiles), "--output", str(cases), *hp, *noise,
        "--cases-per-disease", "25",
    )
    assert simulated.exit_code == 0, simulated.output

    validated = _run("validate", "--cases", str(cases), "--profiles", str(profiles), *hp, *noise)

    assert validated.exit_code == 0, validated.output
    assert "Invariant violations: 0" in validated.output
    report = json.loads((tmp_path / "cases.validation.json").read_text(encoding="utf-8"))
    assert report["cases"] == 4 * 25 * 3
    assert report["config_source"] == str(tmp_path / "cases.summary.json")
    assert report["config_matches_cases"] is True
    assert report["inputs"]["profiles"]["sha256"]
    assert report["calibration"]["all_terms"]["calibration_error"] is not None


def test_validate_exits_non_zero_on_violations(ontology, tmp_path: Path) -> None:
    case = _simulate(ontology)[0]
    cases = tmp_path / "cases.jsonl"
    write_jsonl(cases, [case.model_copy(update={"positive_phenotypes": []})])

    result = _run("validate", "--cases", str(cases), "--hpo-json", str(FIXTURES / "hp_mini.json"))

    assert result.exit_code == 1
    assert "no_observed_positive: 1" in result.output
