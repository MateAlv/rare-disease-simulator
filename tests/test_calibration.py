import json
from pathlib import Path

import pytest

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.simulation.calibration import (
    NOISE_VOCABULARY_KEY,
    CalibrationError,
    apply_calibration,
    calibration_keys,
    load_calibration,
)
from rare_disease_simulator.simulation.schema import SimulationConfig


def _write(path: Path, data: dict) -> Path:
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_keys_cover_nested_knobs_and_skip_run_shape() -> None:
    paths = {key.path for key in calibration_keys()}

    assert {
        "presets.medium.negatives_mean",
        "missingness.sex_unknown",
        "sex.p_male.x_linked_dominant",
        "negatives.source_weights.own_gene_other_disease",
        "negatives.unfilled_slots",
        "age.onset_years.adult",
        NOISE_VOCABULARY_KEY,
    } <= paths
    assert not paths & {"seed", "difficulties", "cases_per_disease_per_difficulty"}
    assert "age.onset_years.unknown" not in paths


def test_every_leaf_of_the_config_is_a_key() -> None:
    def leaves(value, prefix=""):
        if isinstance(value, dict):
            for key, item in value.items():
                yield from leaves(item, f"{prefix}{key}.")
        else:
            yield prefix[:-1]

    dumped = SimulationConfig().model_dump(mode="json")
    for name in ("seed", "difficulties", "cases_per_disease_per_difficulty"):
        dumped.pop(name)
    paths = {key.path for key in calibration_keys()}
    assert set(leaves(dumped)) <= paths


def test_load_applies_overrides_and_records_sha(tmp_path: Path) -> None:
    noise = tmp_path / "noise.tsv"
    noise.write_text("hpo_id\tweight\nHP:0001250\t2\n", encoding="utf-8")
    path = _write(
        tmp_path / "calib.json",
        {
            "presets.medium.negatives_mean": 11,
            "missingness.sex_unknown": 0.08,
            "negatives.unfilled_slots": "redistribute",
            NOISE_VOCABULARY_KEY: "noise.tsv",
            "_provenance": {"split": "r1-v1 train"},
        },
    )

    calibration = load_calibration(path)
    config = apply_calibration(SimulationConfig(), calibration.overrides)

    assert config.presets["medium"].negatives_mean == 11
    assert config.presets["easy"].negatives_mean == 10.0
    assert config.missingness.sex_unknown == 0.08
    assert config.negatives.unfilled_slots == "redistribute"
    assert calibration.noise_vocabulary == noise
    assert calibration.sha256 == sha256_file(path)
    assert calibration.notes == {"_provenance": {"split": "r1-v1 train"}}
    assert calibration.record()["overrides"]["missingness.sex_unknown"] == 0.08


@pytest.mark.parametrize(
    ("data", "message"),
    [
        ({"negatives.count.mean.medium": 11}, "unknown calibration key"),
        ({"seed": 3}, "unknown calibration key"),
        ({"missingness.sex_unknown": 1.5}, "invalid"),
        ({"negatives.own_disease_pool": "gene"}, "invalid"),
        ({NOISE_VOCABULARY_KEY: "missing.tsv"}, "not found"),
    ],
)
def test_load_rejects_bad_files(tmp_path: Path, data: dict, message: str) -> None:
    with pytest.raises(CalibrationError, match=message):
        load_calibration(_write(tmp_path / "calib.json", data))


def test_load_rejects_non_objects(tmp_path: Path) -> None:
    path = tmp_path / "calib.json"
    path.write_text("[1, 2]", encoding="utf-8")
    with pytest.raises(CalibrationError, match="JSON object"):
        load_calibration(path)

