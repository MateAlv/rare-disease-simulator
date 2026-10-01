from rare_disease_simulator.simulation.reporting import load_presentation_ages
from tests.v04_golden import (
    GOLDEN_051,
    GOLDEN_052,
    GOLDEN_INDEPENDENT,
    INDEPENDENT_SCENARIOS,
    PRESENTATION_AGES,
    SCENARIOS_051,
    SCENARIOS_052,
    golden_file_lines,
    golden_lines,
)


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


def test_defaults_reproduce_simulator_050_independent_cases() -> None:
    assert golden_lines(INDEPENDENT_SCENARIOS) == golden_file_lines(GOLDEN_INDEPENDENT)


def test_v051_knobs_at_their_defaults_reproduce_older_cases() -> None:
    knobs = {
        "reporting": {"q_scale": 1.0},
        "noise": {"related_up_levels": 2, "related_down_levels": 2, "related_max_ic": None},
    }

    assert golden_lines(INDEPENDENT_SCENARIOS, **knobs) == golden_file_lines(GOLDEN_INDEPENDENT)
    assert golden_lines(**knobs) == golden_file_lines()


def test_defaults_reproduce_simulator_051_cases() -> None:
    ages = load_presentation_ages(PRESENTATION_AGES)

    assert golden_lines(SCENARIOS_051, ages) == golden_file_lines(GOLDEN_051)


def test_defaults_reproduce_simulator_052_cases() -> None:
    assert golden_lines(SCENARIOS_052) == golden_file_lines(GOLDEN_052)
