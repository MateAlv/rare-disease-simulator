from pathlib import Path

from rare_disease_simulator.config import load_config


def test_load_mvp_config_disease_gene_pairs() -> None:
    config = load_config(Path("configs/mvp.yaml"))

    disease_to_gene = {entry.disease: entry.gene for entry in config.mvp.diseases}

    assert disease_to_gene == {
        "Niemann-Pick disease type C1": "NPC1",
        "Rett syndrome / atypical Rett syndrome": "MECP2",
        "Duchenne muscular dystrophy": "DMD",
        "Cystic fibrosis": "CFTR",
        "Neuronal ceroid lipofuscinosis 6": "CLN6",
        "Pompe disease / glycogen storage disease II": "GAA",
        "Hypophosphatasia": "ALPL",
        "Noonan syndrome": "PTPN11",
    }


def test_load_mvp_config_outputs() -> None:
    config = load_config(Path("configs/mvp.yaml"))

    assert config.exports.graphens_path == Path("outputs/graphens.json")
    assert config.simulation.difficulties == ["easy", "medium", "hard"]
    assert config.llm.provider == "dummy"

