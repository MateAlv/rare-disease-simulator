"""Difficulty-specific simulation presets.

A difficulty preset controls how an underlying true phenotype set is degraded
into an observed clinical picture: how many true positives are actually
recorded, how unobserved positives split between ``missing`` and ``unknown``,
how often explicit negatives are asked, how much nonspecific noise is added,
and how often specific terms are generalized to ontology ancestors.
"""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_simulator.simulation.schema import Difficulty


@dataclass(frozen=True)
class DifficultyPreset:
    """Observation/degradation parameters for a single difficulty level."""

    name: Difficulty
    positive_observation_rate: float
    cardinal_observation_boost: float
    missing_vs_unknown_split: float
    known_negative_rate: float
    noise_term_rate: float
    ontology_smoothing_rate: float


EASY = DifficultyPreset(
    name="easy",
    positive_observation_rate=0.92,
    cardinal_observation_boost=0.08,
    missing_vs_unknown_split=0.8,
    known_negative_rate=0.6,
    noise_term_rate=0.0,
    ontology_smoothing_rate=0.0,
)

MEDIUM = DifficultyPreset(
    name="medium",
    positive_observation_rate=0.75,
    cardinal_observation_boost=0.1,
    missing_vs_unknown_split=0.6,
    known_negative_rate=0.35,
    noise_term_rate=0.15,
    ontology_smoothing_rate=0.15,
)

HARD = DifficultyPreset(
    name="hard",
    positive_observation_rate=0.55,
    cardinal_observation_boost=0.1,
    missing_vs_unknown_split=0.5,
    known_negative_rate=0.2,
    noise_term_rate=0.35,
    ontology_smoothing_rate=0.35,
)

DIFFICULTY_PRESETS: dict[Difficulty, DifficultyPreset] = {
    "easy": EASY,
    "medium": MEDIUM,
    "hard": HARD,
}


def get_difficulty_preset(difficulty: Difficulty) -> DifficultyPreset:
    """Return the preset for a difficulty level."""

    try:
        return DIFFICULTY_PRESETS[difficulty]
    except KeyError as exc:
        raise ValueError(f"unknown difficulty: {difficulty!r}") from exc
