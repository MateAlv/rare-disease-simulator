"""JSON schema export helpers for core models."""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel

from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch, TextSnippet
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.simulation.schema import SimulationConfig, SyntheticCase

CoreSchemaModel = type[BaseModel]

CORE_SCHEMA_MODELS: dict[str, CoreSchemaModel] = {
    "DiseaseProfile": DiseaseProfile,
    "DiseaseProfilePatch": DiseaseProfilePatch,
    "TextSnippet": TextSnippet,
    "SyntheticCase": SyntheticCase,
    "SimulationConfig": SimulationConfig,
}


def core_json_schemas() -> dict[str, dict[str, object]]:
    """Return JSON schemas for the core public models."""

    return {name: model.model_json_schema() for name, model in CORE_SCHEMA_MODELS.items()}


def write_core_json_schemas(output_dir: Path | str) -> list[Path]:
    """Write one JSON schema file per core model."""

    schema_dir = Path(output_dir)
    schema_dir.mkdir(parents=True, exist_ok=True)

    written_paths: list[Path] = []
    for name, schema in core_json_schemas().items():
        path = schema_dir / f"{name}.schema.json"
        path.write_text(json.dumps(schema, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        written_paths.append(path)
    return written_paths
