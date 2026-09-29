"""Deterministic, stratified choice of which diseases to simulate."""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Sequence

from rare_disease_simulator.profiles.schema import DiseaseProfile

INHERITANCE_CLASSES: tuple[tuple[str, str], ...] = (
    ("HP:0000006", "autosomal_dominant"),
    ("HP:0000007", "autosomal_recessive"),
    ("HP:0001417", "x_linked"),
    ("HP:0001419", "x_linked"),
    ("HP:0001423", "x_linked"),
    ("HP:0001427", "mitochondrial"),
    ("HP:0001450", "y_linked"),
)


def inheritance_class(profile: DiseaseProfile) -> str:
    """One coarse inheritance class: a single mode, ``mixed``, ``other`` or ``none``."""

    ids = {hpo_id for gene in profile.genes for hpo_id in gene.inheritance_hpo_ids}
    if not ids:
        return "none"
    classes = {label for hpo_id, label in INHERITANCE_CLASSES if hpo_id in ids}
    if len(classes) == 1:
        return classes.pop()
    return "mixed" if classes else "other"


def onset_source(profile: DiseaseProfile) -> str:
    """Source of the disease onset: ``hpoa``, ``orphanet``, ``orphanet_via_omim`` or ``none``."""

    onset = profile.age_of_onset
    if onset is None:
        return "none"
    names = [item.source.name for item in onset.provenance]
    if "Orphanet OMIM-ORPHA alignments" in names:
        return "orphanet_via_omim"
    if names and names[0].startswith("Orphanet"):
        return "orphanet"
    return "hpoa"


def stratum(profile: DiseaseProfile) -> str:
    return f"{inheritance_class(profile)}|{onset_source(profile)}"


def sample_diseases(
    profiles: Sequence[DiseaseProfile], count: int, seed: int
) -> tuple[list[DiseaseProfile], dict[str, int]]:
    """Pick ``count`` profiles spread evenly over (inheritance class, onset source) strata.

    Each stratum gets an equal share; strata smaller than their share give
    all their diseases and the remainder is spread over the others. Within a
    stratum the choice is a seeded sample of the sorted disease ids. The
    result keeps the input order.
    """

    strata: dict[str, list[DiseaseProfile]] = defaultdict(list)
    for profile in profiles:
        strata[stratum(profile)].append(profile)
    quotas = {key: 0 for key in strata}
    remaining = min(count, len(profiles))
    open_keys = sorted(strata)
    while remaining > 0 and open_keys:
        share, extra = divmod(remaining, len(open_keys))
        next_keys = []
        for position, key in enumerate(open_keys):
            capacity = len(strata[key]) - quotas[key]
            take = min(capacity, share + (1 if position < extra else 0))
            quotas[key] += take
            remaining -= take
            if len(strata[key]) > quotas[key]:
                next_keys.append(key)
        open_keys = next_keys

    rng = random.Random(seed)
    chosen: set[str] = set()
    for key in sorted(strata):
        ids = sorted(profile.disease_id for profile in strata[key])
        chosen.update(rng.sample(ids, quotas[key]))
    selected = [profile for profile in profiles if profile.disease_id in chosen]
    return selected, {key: quotas[key] for key in sorted(strata) if quotas[key]}
