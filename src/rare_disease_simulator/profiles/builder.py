"""Disease profile construction from structured sources and validated patches."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rare_disease_simulator.data_sources.fixtures import (
    FixtureBundle,
    load_fixture_hpo_ontology,
    load_fixture_profile_patch,
    load_fixture_structured_profile,
    load_fixture_text_snippets,
)
from rare_disease_simulator.exports.jsonl import write_jsonl
from rare_disease_simulator.llm_extraction.merge_patch import merge_patch_into_profile
from rare_disease_simulator.llm_extraction.validator import (
    DiseaseProfilePatchValidator,
    PatchValidationResult,
)
from rare_disease_simulator.profiles.schema import DiseaseProfile
from rare_disease_simulator.validation.hpo import HpoValidator


@dataclass(frozen=True)
class ProfileBuildResult:
    """Profiles built from a source path plus validation metadata."""

    profiles: list[DiseaseProfile]
    patch_validation: PatchValidationResult
    warnings: list[str]


def build_profiles_from_fixtures(fixture_dir: Path | str) -> ProfileBuildResult:
    """Run the first fixture vertical slice: structured profile + validated patch."""

    bundle = FixtureBundle.from_dir(fixture_dir)
    bundle.require_files()

    ontology = load_fixture_hpo_ontology(bundle)
    structured_profile = load_fixture_structured_profile(bundle)
    patch = load_fixture_profile_patch(bundle)
    source_text = _find_matching_source_text(bundle, patch.disease_id, patch.gene)

    patch_validation = DiseaseProfilePatchValidator(HpoValidator(ontology)).validate(
        patch,
        source_text=source_text,
    )
    merge_result = merge_patch_into_profile(structured_profile, patch, patch_validation)

    return ProfileBuildResult(
        profiles=[merge_result.profile],
        patch_validation=patch_validation,
        warnings=merge_result.warnings,
    )


def write_profiles_jsonl(path: Path | str, profiles: list[DiseaseProfile]) -> int:
    """Write disease profiles to JSONL."""

    return write_jsonl(path, profiles)


def _find_matching_source_text(bundle: FixtureBundle, disease_id: str, gene: str) -> str:
    snippets = load_fixture_text_snippets(bundle)
    for snippet in snippets:
        if snippet.disease_id == disease_id and snippet.gene == gene:
            return snippet.text
    raise ValueError(f"no fixture text snippet found for disease_id={disease_id!r}, gene={gene!r}")
