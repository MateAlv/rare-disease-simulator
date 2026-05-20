"""Readers for small test-only fixture files."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch, TextSnippet

FIXTURE_DIR = Path(__file__).parent


def fixture_path(name: str) -> Path:
    """Return the absolute path for a fixture file."""

    return FIXTURE_DIR / name


def read_hpo_terms(path: Path | None = None) -> list[dict[str, str]]:
    """Read the minimal HPO terms TSV fixture."""

    terms_path = path or fixture_path("hpo_terms.tsv")
    with terms_path.open("r", encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def read_orphadata_mini(path: Path | None = None) -> dict[str, Any]:
    """Read the minimal Orphadata-like JSON fixture."""

    data_path = path or fixture_path("orphadata_mini.json")
    return json.loads(data_path.read_text(encoding="utf-8"))


def read_text_snippets(path: Path | None = None) -> list[TextSnippet]:
    """Read trusted text snippet JSONL fixtures."""

    snippets_path = path or fixture_path("text_snippets.jsonl")
    snippets: list[TextSnippet] = []
    with snippets_path.open("r", encoding="utf-8") as file:
        for line in file:
            if line.strip():
                snippets.append(TextSnippet.model_validate_json(line))
    return snippets


def read_profile_patch(path: Path | None = None) -> DiseaseProfilePatch:
    """Read a fixed DiseaseProfilePatch fixture."""

    patch_path = path or fixture_path("profile_patch.json")
    return DiseaseProfilePatch.model_validate_json(patch_path.read_text(encoding="utf-8"))


def read_disease_profile(path: Path | None = None) -> dict[str, Any]:
    """Read the expected merged disease profile fixture."""

    profile_path = path or fixture_path("disease_profile.json")
    return json.loads(profile_path.read_text(encoding="utf-8"))

