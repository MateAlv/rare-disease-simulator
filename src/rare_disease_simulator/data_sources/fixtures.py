"""Runnable fixture loading for local vertical-slice demos."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.exports.jsonl import read_model_jsonl
from rare_disease_simulator.llm_extraction.schema import DiseaseProfilePatch, TextSnippet
from rare_disease_simulator.profiles.schema import DiseaseProfile


@dataclass(frozen=True)
class FixtureBundle:
    """Files required to run the local fixture profile-building path."""

    fixture_dir: Path
    hpo_terms_path: Path
    text_snippets_path: Path
    profile_patch_path: Path
    structured_profile_path: Path

    @classmethod
    def from_dir(cls, fixture_dir: Path | str) -> FixtureBundle:
        """Build a fixture bundle from a directory."""

        base_dir = Path(fixture_dir)
        return cls(
            fixture_dir=base_dir,
            hpo_terms_path=base_dir / "hpo_terms.tsv",
            text_snippets_path=base_dir / "text_snippets.jsonl",
            profile_patch_path=base_dir / "profile_patch.json",
            structured_profile_path=base_dir / "structured_profile.json",
        )

    def require_files(self) -> None:
        """Raise FileNotFoundError if any fixture file is missing."""

        for path in [
            self.hpo_terms_path,
            self.text_snippets_path,
            self.profile_patch_path,
            self.structured_profile_path,
        ]:
            if not path.exists():
                raise FileNotFoundError(path)


def load_fixture_hpo_ontology(bundle: FixtureBundle) -> HpoOntology:
    """Load the fixture HPO ontology."""

    return HpoOntology.from_tsv(bundle.hpo_terms_path, version="fixture-0.1")


def load_fixture_text_snippets(bundle: FixtureBundle) -> list[TextSnippet]:
    """Load fixture trusted text snippets."""

    return read_model_jsonl(bundle.text_snippets_path, TextSnippet)


def load_fixture_profile_patch(bundle: FixtureBundle) -> DiseaseProfilePatch:
    """Load the fixed fixture profile patch."""

    return DiseaseProfilePatch.model_validate_json(
        bundle.profile_patch_path.read_text(encoding="utf-8")
    )


def load_fixture_structured_profile(bundle: FixtureBundle) -> DiseaseProfile:
    """Load the structured profile seed used by the fixture merge path."""

    return DiseaseProfile.model_validate_json(
        bundle.structured_profile_path.read_text(encoding="utf-8")
    )

