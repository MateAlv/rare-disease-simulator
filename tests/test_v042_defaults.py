from tests.v04_golden import golden_file_lines, golden_lines


def test_defaults_reproduce_simulator_041_cases() -> None:
    assert golden_lines() == golden_file_lines()


def test_new_knobs_at_zero_reproduce_simulator_041_cases() -> None:
    lines = golden_lines(reporting={"specialize_rate": 0.0}, noise={"related_share": 0.0})

    assert lines == golden_file_lines()


def test_v05_knobs_at_their_defaults_reproduce_simulator_041_cases() -> None:
    lines = golden_lines(
        reporting={"mode": "report_model", "specialize_rate": 0.0},
        noise={"related_share": 0.0},
        age={"duration_mean_by_onset": {}},
    )

    assert lines == golden_file_lines()
