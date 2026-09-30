import gzip

from tests.v04_golden import GOLDEN, golden_lines


def _golden() -> list[str]:
    return gzip.decompress(GOLDEN.read_bytes()).decode().splitlines()


def test_defaults_reproduce_simulator_041_cases() -> None:
    assert golden_lines() == _golden()


def test_new_knobs_at_zero_reproduce_simulator_041_cases() -> None:
    lines = golden_lines(reporting={"specialize_rate": 0.0}, noise={"related_share": 0.0})

    assert lines == _golden()
