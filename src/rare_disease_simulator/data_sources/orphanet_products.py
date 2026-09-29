"""Readers for Orphanet bulk XML products (Orphadata "JDBOR" files).

- ``en_product9_ages.xml``: average age of onset per disorder.
- ``en_product1.xml``: disorders with their cross-references. Only exact
  (``E``), validated OMIM references of active disorders become OMIM <-> ORPHA
  mappings; ``BTNT``/``NTBT``/``ND`` and not-yet-validated references are
  counted and ignored, so facts never move between non-equivalent concepts.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

NO_DATA_ONSET = "No data available"
EXACT_MAPPING = "E"
VALIDATED_STATUS = "Validated"
INACTIVE_FLAGS = frozenset({"Inactive", "Obsolete entity", "Deprecated entity"})


@dataclass(frozen=True)
class OrphanetOnsets:
    """Average-age-of-onset labels per ``ORPHA:<code>``."""

    onsets: dict[str, tuple[str, ...]]
    version: str | None


@dataclass(frozen=True)
class OmimOrphaMap:
    """Exact OMIM <-> ORPHA mappings from Orphanet cross-references."""

    omim_to_orpha: dict[str, tuple[str, ...]]
    orpha_to_omim: dict[str, tuple[str, ...]]
    version: str | None
    stats: dict[str, int] = field(default_factory=dict)


def read_orphanet_onsets(path: Path | str) -> OrphanetOnsets:
    """Read ``AverageAgeOfOnset`` names per disorder, dropping "No data available"."""

    onsets: dict[str, tuple[str, ...]] = {}
    version: str | None = None
    for event, element in ET.iterparse(Path(path), events=("start", "end")):
        if event == "start":
            if element.tag == "JDBOR" and version is None:
                version = element.get("date")
            continue
        if element.tag != "Disorder":
            continue
        code = (element.findtext("OrphaCode") or "").strip()
        names = [
            (onset.findtext("Name") or "").strip()
            for onset in element.iterfind("AverageAgeOfOnsetList/AverageAgeOfOnset")
        ]
        kept = tuple(sorted({name for name in names if name and name != NO_DATA_ONSET}))
        if code and kept:
            onsets[f"ORPHA:{code}"] = kept
        element.clear()
    return OrphanetOnsets(onsets=onsets, version=version)


def read_omim_orpha_map(path: Path | str) -> OmimOrphaMap:
    """Read exact (``E``), validated OMIM references of each active ORPHA disorder.

    ``stats`` counts every OMIM reference as ``omim_<relation>_<status>`` plus
    the exact references skipped because their disorder is inactive, and the
    OMIM ids with more than one exact ORPHA mapping.
    """

    omim_to_orpha: dict[str, set[str]] = {}
    orpha_to_omim: dict[str, set[str]] = {}
    stats: Counter[str] = Counter()
    version: str | None = None
    for event, element in ET.iterparse(Path(path), events=("start", "end")):
        if event == "start":
            if element.tag == "JDBOR" and version is None:
                version = element.get("date")
            continue
        if element.tag != "Disorder":
            continue
        code = (element.findtext("OrphaCode") or "").strip()
        flags = {
            (flag.findtext("Label") or "").strip()
            for flag in element.iterfind("DisorderFlagList/DisorderFlag")
        }
        for reference in element.iterfind("ExternalReferenceList/ExternalReference"):
            if (reference.findtext("Source") or "").strip() != "OMIM":
                continue
            relation = (reference.findtext("DisorderMappingRelation/Name") or "").strip()
            relation_code = relation.split(" ", 1)[0] or "blank"
            status = (reference.findtext("DisorderMappingValidationStatus/Name") or "").strip()
            validated = status == VALIDATED_STATUS
            stats[f"omim_{relation_code}_{'validated' if validated else 'not_validated'}"] += 1
            omim_code = (reference.findtext("Reference") or "").strip()
            if not code or not omim_code or relation_code != EXACT_MAPPING or not validated:
                continue
            if flags & INACTIVE_FLAGS:
                stats["exact_skipped_inactive_disorder"] += 1
                continue
            orpha_id = f"ORPHA:{code}"
            omim_id = f"OMIM:{omim_code}"
            omim_to_orpha.setdefault(omim_id, set()).add(orpha_id)
            orpha_to_omim.setdefault(orpha_id, set()).add(omim_id)
        element.clear()
    stats["omim_ids_mapped"] = len(omim_to_orpha)
    stats["omim_ids_with_several_orpha"] = sum(1 for v in omim_to_orpha.values() if len(v) > 1)
    return OmimOrphaMap(
        omim_to_orpha={key: tuple(sorted(value)) for key, value in omim_to_orpha.items()},
        orpha_to_omim={key: tuple(sorted(value)) for key, value in orpha_to_omim.items()},
        version=version,
        stats=dict(sorted(stats.items())),
    )
