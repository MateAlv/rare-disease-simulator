"""LLM extraction orchestration."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch, TextSnippet


class LlmProvider(Protocol):
    """Provider interface for converting trusted text snippets into profile patches."""

    def extract(self, snippet: TextSnippet) -> DiseaseProfilePatch:
        """Extract a disease profile patch from a trusted text snippet."""


class DummyLlmProvider:
    """Deterministic provider for fixture-driven development and tests."""

    def __init__(self, patch: DiseaseProfilePatch) -> None:
        self._patch = patch

    @classmethod
    def from_patch_file(cls, path: Path | str) -> DummyLlmProvider:
        """Load a fixed patch from a JSON fixture file."""

        patch_path = Path(path)
        return cls(DiseaseProfilePatch.model_validate_json(patch_path.read_text(encoding="utf-8")))

    def extract(self, snippet: TextSnippet) -> DiseaseProfilePatch:
        """Return the fixed patch after checking it targets the same disease and gene."""

        if snippet.disease_id != self._patch.disease_id:
            raise ValueError(
                f"dummy patch disease_id={self._patch.disease_id!r} does not match "
                f"snippet disease_id={snippet.disease_id!r}"
            )
        if snippet.gene != self._patch.gene:
            raise ValueError(
                f"dummy patch gene={self._patch.gene!r} does not match "
                f"snippet gene={snippet.gene!r}"
            )
        return self._patch.model_copy(deep=True)
