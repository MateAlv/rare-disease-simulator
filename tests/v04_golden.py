"""Cases used to check that new knobs left at their defaults reproduce older output.

Each golden file was written by this module with the code of the version it
names, in a throwaway worktree:

- ``golden_0.4.1.jsonl.gz`` (b645c60): report-model mode;
- ``golden_0.5.0_independent.jsonl.gz`` (020e41d): independent reporting,
  proportional and related noise, per-onset durations;
- ``golden_0.5.1.jsonl.gz`` (09e5802): q_scale, one-level related noise under
  an IC cap, presentation ages;
- ``golden_0.5.2.jsonl.gz`` (40f8896): budget-normalized independent reporting;
- ``golden_0.5.3.jsonl.gz`` (6833d2d): record-scope budget normalization, presentation
  ages.

Each line is a case without the fields that name the simulator build
(``simulator_version``, ``config_hash``). Lines are compared without null
fields, as ``write_jsonl`` writes cases, so a new optional field left unset
keeps the comparison byte-for-byte meaningful.
"""

from __future__ import annotations

import gzip
import json

from rare_disease_simulator.simulation.confounders import ConfounderIndex
from rare_disease_simulator.simulation.reporting import Reporting, load_report_model
from rare_disease_simulator.simulation.simulator import NoiseTerm, simulate_gene_cases
from tests.fixtures.readers import fixture_path
from tests.test_gene_first import ALL_PROFILES, PROFILE_MAP, TARGET
from tests.test_simulation_v02 import TERMS, build_ontology
from tests.test_simulation_v02 import _config as _base_config
from tests.test_v04_simulation import CARDINAL

GOLDEN = fixture_path("v04") / "golden_0.4.1.jsonl.gz"
GOLDEN_INDEPENDENT = fixture_path("v04") / "golden_0.5.0_independent.jsonl.gz"
GOLDEN_051 = fixture_path("v04") / "golden_0.5.1.jsonl.gz"
GOLDEN_052 = fixture_path("v04") / "golden_0.5.2.jsonl.gz"
GOLDEN_053 = fixture_path("v04") / "golden_0.5.3.jsonl.gz"
PRESENTATION_AGES = fixture_path("v05") / "presentation_ages.tsv"
NOISE = [
    NoiseTerm(hpo_id, TERMS[hpo_id][0])
    for hpo_id in ("HP:0000964", "HP:0000988", "HP:0000008", "HP:0000028", "HP:0000505")
]
SCENARIOS = {
    "poisson": {"reporting": {"mode": "report_model"}},
    "budget_share": {
        "reporting": {"mode": "report_model", "noise_count": "budget_share", "noise_share": 0.4}
    },
}
INDEPENDENT_SCENARIOS = {
    "independent": {
        "reporting": {"mode": "independent", "noise_count": "proportional", "noise_share": 0.24,
                      "specialize_rate": 0.3},
        "noise": {"related_share": 0.5},
        "age": {"duration_mean_by_onset": {"infantile": 9.8}},
    },
}
SCENARIOS_051 = {
    "independent_051": {
        "reporting": {"mode": "independent", "noise_count": "proportional", "noise_share": 0.24,
                      "specialize_rate": 0.3, "q_scale": 1.5},
        "noise": {"related_share": 0.5, "related_up_levels": 1, "related_down_levels": 1,
                  "related_max_ic": 3.5},
    },
}
SCENARIOS_052 = {
    "budget_normalized_052": {
        "reporting": {"mode": "independent", "noise_count": "proportional", "noise_share": 0.24,
                      "specialize_rate": 0.3, "q_scale": 1.5, "budget_normalize": True,
                      "profile_budget_histogram": {"0": 2, "1": 3, "2": 5, "3": 4, "4": 2}},
        "noise": {"related_share": 0.5},
    },
}
SCENARIOS_053 = {
    "record_scope_053": {
        "reporting": {"mode": "independent", "noise_count": "budget_share", "noise_share": 0.24,
                      "specialize_rate": 0.3, "q_scale": 1.5, "budget_normalize": True,
                      "budget_scope": "record",
                      "profile_budget_histogram": {"0": 2, "3": 3, "4": 5, "5": 4, "6": 2}},
        "noise": {"related_share": 0.5},
        "age": {"duration_mean_by_onset": {"infantile": 9.8}},
    },
}


def golden_lines(
    scenarios: dict | None = None,
    presentation_ages: object = None,
    onset_ages: object = None,
    **extra: dict,
) -> list[str]:
    ontology = build_ontology()
    reporting = Reporting.build(
        load_report_model(fixture_path("v04") / "report_model.json"),
        ontology=ontology, profiles=ALL_PROFILES, cardinal=CARDINAL,
    )
    lines: list[str] = []
    for name, overrides in (scenarios or SCENARIOS).items():
        merged = {
            key: {**overrides.get(key, {}), **extra.get(key, {})}
            for key in set(overrides) | set(extra)
        }
        config = _base_config(
            cases_per_disease_per_difficulty=40, difficulties=["medium", "hard"], **merged
        )
        index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)
        for case in simulate_gene_cases(
            TARGET, PROFILE_MAP, config, ontology=ontology, confounders=index,
            noise_vocabulary=NOISE, reporting=reporting,
            **({"presentation_ages": presentation_ages} if presentation_ages else {}),
            **({"onset_ages": onset_ages} if onset_ages else {}),
        ):
            data = json.loads(case.model_dump_json(exclude_none=True))
            del data["metadata"]["simulator_version"], data["metadata"]["config_hash"]
            lines.append(json.dumps({"scenario": name, **data}, sort_keys=True))
    return lines


def golden_file_lines(path=GOLDEN) -> list[str]:
    """The stored golden cases, without null fields."""

    raw = gzip.decompress(path.read_bytes()).decode().splitlines()
    return [json.dumps(_without_nulls(json.loads(line)), sort_keys=True) for line in raw]


def _without_nulls(value: object) -> object:
    if isinstance(value, dict):
        return {k: _without_nulls(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_without_nulls(item) for item in value]
    return value


def _write(path, lines: list[str]) -> None:
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))


if __name__ == "__main__":
    import sys

    if sys.argv[1:] == ["independent"]:
        _write(GOLDEN_INDEPENDENT, golden_lines(INDEPENDENT_SCENARIOS))
    elif sys.argv[1:] == ["0.5.2"]:
        _write(GOLDEN_052, golden_lines(SCENARIOS_052))
    elif sys.argv[1:] == ["0.5.3"]:
        from rare_disease_simulator.simulation.reporting import load_presentation_ages

        _write(GOLDEN_053, golden_lines(SCENARIOS_053, load_presentation_ages(PRESENTATION_AGES)))
    elif sys.argv[1:] == ["0.5.1"]:
        from rare_disease_simulator.simulation.reporting import load_presentation_ages

        _write(GOLDEN_051, golden_lines(SCENARIOS_051, load_presentation_ages(PRESENTATION_AGES)))
    else:
        _write(GOLDEN, golden_lines())
