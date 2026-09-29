"""Difficulty-specific simulation presets.

A difficulty preset controls how an underlying true phenotype set is degraded
into an observed clinical picture: how many true positives are actually
recorded, how unobserved positives split between ``missing`` and ``unknown``,
how many "asked and absent" negatives and nonspecific noise terms a case gets,
and how often specific terms are generalized to ontology ancestors.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

Difficulty = Literal["easy", "medium", "hard"]


class DifficultyPreset(BaseModel):
    """Observation/degradation parameters for a single difficulty level."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    positive_observation_rate: float = Field(ge=0.0, le=1.0)
    cardinal_observation_boost: float = Field(ge=0.0, le=1.0)
    missing_vs_unknown_split: float = Field(ge=0.0, le=1.0)
    negatives_mean: float = Field(
        ge=0.0, description="Poisson mean of asked-and-absent terms per case."
    )
    noise_mean: float = Field(ge=0.0, description="Poisson mean of noise terms per case.")
    ontology_smoothing_rate: float = Field(ge=0.0, le=1.0)


EASY = DifficultyPreset(
    positive_observation_rate=0.92,
    cardinal_observation_boost=0.08,
    missing_vs_unknown_split=0.8,
    negatives_mean=3.0,
    noise_mean=0.0,
    ontology_smoothing_rate=0.0,
)

MEDIUM = DifficultyPreset(
    positive_observation_rate=0.75,
    cardinal_observation_boost=0.1,
    missing_vs_unknown_split=0.6,
    negatives_mean=2.0,
    noise_mean=1.0,
    ontology_smoothing_rate=0.15,
)

HARD = DifficultyPreset(
    positive_observation_rate=0.55,
    cardinal_observation_boost=0.1,
    missing_vs_unknown_split=0.5,
    negatives_mean=1.0,
    noise_mean=2.0,
    ontology_smoothing_rate=0.35,
)

DIFFICULTY_PRESETS: dict[Difficulty, DifficultyPreset] = {
    "easy": EASY,
    "medium": MEDIUM,
    "hard": HARD,
}


def get_difficulty_preset(difficulty: Difficulty) -> DifficultyPreset:
    """Return the built-in preset for a difficulty level."""

    try:
        return DIFFICULTY_PRESETS[difficulty]
    except KeyError as exc:
        raise ValueError(f"unknown difficulty: {difficulty!r}") from exc
