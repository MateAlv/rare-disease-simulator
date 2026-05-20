"""Configuration loading for the rare disease simulator."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field

DEFAULT_CONFIG_PATH = Path("configs/mvp.yaml")


class MvpDisease(BaseModel):
    """Configured disease/gene pair for the MVP."""

    disease: str
    gene: str


class MvpConfig(BaseModel):
    """MVP-specific configuration."""

    diseases: list[MvpDisease] = Field(default_factory=list)


class SourcePathsConfig(BaseModel):
    """Local source file locations."""

    hpo_dir: Path = Path("data/raw/hpo")
    orphadata_dir: Path = Path("data/raw/orphadata")
    mondo_path: Path = Path("data/raw/mondo/mondo.json")


class LlmConfig(BaseModel):
    """LLM extraction settings."""

    provider: str = "dummy"
    model: str = "dummy"
    prompt_version: str = "mvp-v0"
    extraction_schema_version: str = "mvp-v0"
    cache_dir: Path = Path("data/cache/llm_extraction")
    temperature: float = 0.0


class SimulationSettings(BaseModel):
    """Synthetic case simulation settings."""

    difficulties: list[str] = Field(default_factory=lambda: ["easy", "medium", "hard"])
    cases_per_disease_per_difficulty: int = 100
    seed: int = 42


class ExportPathsConfig(BaseModel):
    """Output artifact paths."""

    profiles_path: Path = Path("outputs/profiles.jsonl")
    profile_patches_path: Path = Path("outputs/profile_patches.jsonl")
    rich_cases_path: Path = Path("outputs/rich_cases.jsonl")
    graphens_path: Path = Path("outputs/graphens.json")
    review_dir: Path = Path("outputs/review")


class AppConfig(BaseModel):
    """Top-level application configuration."""

    model_config = ConfigDict(extra="forbid")

    mvp: MvpConfig
    sources: SourcePathsConfig = Field(default_factory=SourcePathsConfig)
    llm: LlmConfig = Field(default_factory=LlmConfig)
    simulation: SimulationSettings = Field(default_factory=SimulationSettings)
    exports: ExportPathsConfig = Field(default_factory=ExportPathsConfig)


def load_config(path: Path | str = DEFAULT_CONFIG_PATH) -> AppConfig:
    """Load and validate a YAML configuration file."""

    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as file:
        data = yaml.safe_load(file) or {}
    return AppConfig.model_validate(data)

