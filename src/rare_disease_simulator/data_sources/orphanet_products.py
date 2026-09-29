"""Readers for Orphanet bulk XML products (Orphadata "JDBOR" files).

- ``en_product9_ages.xml``: average age of onset per disorder.
- ``en_product1.xml`` (``alignments_omim_orpha.xml``): cross-references, used
  here only for exact, validated OMIM <-> ORPHA mappings.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

NO_DATA_ONSET = "No data available"
EXACT_MAPPING_PREFIX = "E "
VALIDATED_STATUS = "Validated"


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
    """Read exact (``E``), validated OMIM references of each ORPHA disorder."""

    omim_to_orpha: dict[str, set[str]] = {}
    orpha_to_omim: dict[str, set[str]] = {}
    version: str | None = None
    for event, element in ET.iterparse(Path(path), events=("start", "end")):
        if event == "start":
            if element.tag == "JDBOR" and version is None:
                version = element.get("date")
            continue
        if element.tag != "Disorder":
            continue
        code = (element.findtext("OrphaCode") or "").strip()
        for reference in element.iterfind("ExternalReferenceList/ExternalReference"):
            if (reference.findtext("Source") or "").strip() != "OMIM":
                continue
            relation = (reference.findtext("DisorderMappingRelation/Name") or "").strip()
            status = (reference.findtext("DisorderMappingValidationStatus/Name") or "").strip()
            omim_code = (reference.findtext("Reference") or "").strip()
            if not code or not omim_code:
                continue
            if not relation.startswith(EXACT_MAPPING_PREFIX) or status != VALIDATED_STATUS:
                continue
            orpha_id = f"ORPHA:{code}"
            omim_id = f"OMIM:{omim_code}"
            omim_to_orpha.setdefault(omim_id, set()).add(orpha_id)
            orpha_to_omim.setdefault(orpha_id, set()).add(omim_id)
        element.clear()
    return OmimOrphaMap(
        omim_to_orpha={key: tuple(sorted(value)) for key, value in omim_to_orpha.items()},
        orpha_to_omim={key: tuple(sorted(value)) for key, value in orpha_to_omim.items()},
        version=version,
    )
