"""HPO validation helpers."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_simulator.data_sources.hpo import HpoOntology


@dataclass(frozen=True)
class HpoValidationIssue:
    """Validation issue for a candidate HPO ID."""

    hpo_id: str
    reason: str


class HpoValidator:
    """Validator facade around an HPO ontology cache."""

    def __init__(self, ontology: HpoOntology) -> None:
        self.ontology = ontology

    def is_valid_hpo_id(self, hpo_id: str) -> bool:
        """Return whether an HPO ID exists."""

        return self.ontology.is_valid_hpo_id(hpo_id)

    def get_label(self, hpo_id: str) -> str | None:
        """Return an HPO label for a known ID."""

        return self.ontology.get_label(hpo_id)

    def get_ancestors(self, hpo_id: str) -> tuple[str, ...]:
        """Return known ancestor IDs for a term."""

        return self.ontology.get_ancestors(hpo_id)

    def get_direct_parents(self, hpo_id: str) -> tuple[str, ...]:
        """Return direct parent IDs for a term."""

        return self.ontology.get_direct_parents(hpo_id)

    def is_phenotypic_abnormality(self, hpo_id: str) -> bool:
        """Return whether a known HPO ID is phenotypic."""

        return self.ontology.is_phenotypic_abnormality(hpo_id)

    def validate_phenotypic_hpo_ids(self, hpo_ids: list[str]) -> list[HpoValidationIssue]:
        """Return issues for unknown or non-phenotypic HPO IDs."""

        issues: list[HpoValidationIssue] = []
        for hpo_id in hpo_ids:
            if not self.is_valid_hpo_id(hpo_id):
                issues.append(HpoValidationIssue(hpo_id=hpo_id, reason="unknown_hpo_id"))
            elif not self.is_phenotypic_abnormality(hpo_id):
                issues.append(
                    HpoValidationIssue(
                        hpo_id=hpo_id,
                        reason="not_phenotypic_abnormality",
                    )
                )
        return issues
