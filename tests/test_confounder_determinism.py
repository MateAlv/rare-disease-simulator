import os
import subprocess
import sys
from pathlib import Path

SCRIPT = """
import random
from rare_disease_simulator.profiles.schema import DiseaseProfile, PhenotypeAssociation
from rare_disease_simulator.simulation.confounders import ConfounderIndex

rng = random.Random(0)
pool = [f"HP:{number:07d}" for number in range(1, 60)]
profiles = [
    DiseaseProfile(
        disease_id=f"OMIM:{index}",
        disease_name=f"d{index}",
        phenotypes=[
            PhenotypeAssociation(hpo_id=hpo_id, label=hpo_id, frequency_estimate=0.5)
            for hpo_id in rng.sample(pool, rng.randint(8, 20))
        ],
    )
    for index in range(80)
]
index = ConfounderIndex(profiles, top_n=5, min_information_content=0.0)
print([[(other, score.hex()) for other, score in index.confounders(p.disease_id)]
       for p in profiles])
"""


def _run(hash_seed: str) -> str:
    source = Path(__file__).resolve().parents[1] / "src"
    env = {**os.environ, "PYTHONHASHSEED": hash_seed, "PYTHONPATH": str(source)}
    completed = subprocess.run(
        [sys.executable, "-c", SCRIPT], env=env, capture_output=True, text=True, check=True
    )
    return completed.stdout


def test_confounder_scores_do_not_depend_on_the_string_hash_seed() -> None:
    assert _run("1") == _run("2") == _run("3")
