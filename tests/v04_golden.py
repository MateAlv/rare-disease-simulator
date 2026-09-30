"""Cases used to check that new knobs left at 0 reproduce simulator 0.4.1 exactly.

``tests/fixtures/v04/golden_0.4.1.jsonl.gz`` was written by this module with the
0.4.1 code (b645c60). Each line is a case without the fields that name the
simulator build (``simulator_version``, ``config_hash``).
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


def golden_lines(**extra: object) -> list[str]:
    ontology = build_ontology()
    reporting = Reporting.build(
        load_report_model(fixture_path("v04") / "report_model.json"),
        ontology=ontology, profiles=ALL_PROFILES, cardinal=CARDINAL,
    )
    lines: list[str] = []
    for name, overrides in SCENARIOS.items():
        merged = {**overrides, **extra}
        if "reporting" in extra:
            merged["reporting"] = {**overrides["reporting"], **extra["reporting"]}
        config = _base_config(
            cases_per_disease_per_difficulty=40, difficulties=["medium", "hard"], **merged
        )
        index = ConfounderIndex(ALL_PROFILES, ontology, top_n=10, min_information_content=0.0)
        for case in simulate_gene_cases(
            TARGET, PROFILE_MAP, config, ontology=ontology, confounders=index,
            noise_vocabulary=NOISE, reporting=reporting,
        ):
            data = json.loads(case.model_dump_json())
            del data["metadata"]["simulator_version"], data["metadata"]["config_hash"]
            lines.append(json.dumps({"scenario": name, **data}, sort_keys=True))
    return lines


if __name__ == "__main__":
    GOLDEN.write_bytes(gzip.compress(("\n".join(golden_lines()) + "\n").encode(), mtime=0))
