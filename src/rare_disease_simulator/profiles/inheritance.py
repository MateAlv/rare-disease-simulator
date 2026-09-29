"""Sex bias implied by modes of inheritance, shared by profiles and the simulator."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import SexBiasValue

MENDELIAN_INHERITANCE = "HP:0034345"
AUTOSOMAL_MODES = ("HP:0000006", "HP:0000007", "HP:0032113")
MALE_ONLY_MODES = ("HP:0001419", "HP:0001450")
MALE_LIMITED = "HP:0001475"
FEMALE_LIMITED = "HP:0034344"
SEX_LIMITED_EXPRESSION = (MALE_LIMITED, FEMALE_LIMITED)


def derive_sex_bias(
    inheritance: Sequence[str], ontology: HpoOntology | None
) -> tuple[SexBiasValue | None, list[str]]:
    """Sex bias the inheritance terms imply, and the terms it rests on.

    Male-/female-limited expression -> that sex. All Mendelian modes X-linked
    recessive or Y-linked -> male. All Mendelian modes autosomal -> none. Any
    other combination (X-linked dominant or unspecified, mitochondrial, mixed
    modes) gives None: the inheritance alone does not determine it. Without an
    ontology, "is a" falls back to identity and every term counts as a mode.
    """

    def is_a(hpo_id: str, anchor: str) -> bool:
        if ontology is None:
            return hpo_id == anchor
        return ontology.is_a(hpo_id, anchor)

    def any_is_a(hpo_id: str, anchors: Iterable[str]) -> bool:
        return any(is_a(hpo_id, anchor) for anchor in anchors)

    male_limited = [hpo_id for hpo_id in inheritance if is_a(hpo_id, MALE_LIMITED)]
    female_limited = [hpo_id for hpo_id in inheritance if is_a(hpo_id, FEMALE_LIMITED)]
    modes = [
        hpo_id
        for hpo_id in inheritance
        if (
            ontology.is_a(hpo_id, MENDELIAN_INHERITANCE)
            if ontology is not None
            else hpo_id not in SEX_LIMITED_EXPRESSION
        )
    ]

    if male_limited and not female_limited:
        return "male", male_limited
    if female_limited and not male_limited:
        return "female", female_limited
    if male_limited or female_limited:
        return None, []
    if modes and all(any_is_a(hpo_id, MALE_ONLY_MODES) for hpo_id in modes):
        return "male", modes
    if modes and all(any_is_a(hpo_id, AUTOSOMAL_MODES) for hpo_id in modes):
        return "none", modes
    return None, []
