"""HPO disease-annotation ingestion (the structured backbone).

Parses ``phenotype.hpoa`` — the curated HPO disease-phenotype annotation file —
which already carries, per disease: HPO term, frequency, onset, sex, an explicit
NOT qualifier for negatives, inheritance, and clinical-course rows. This is the
no-LLM backbone: frequency/onset/sex/negatives come from here, not the model.

Downloads use the stable OBO PURLs and are meant to run once into ``data/raw``.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from rare_disease_simulator.data_sources.http_client import HttpClient

logger = logging.getLogger(__name__)

HPO_JSON_URL = "https://purl.obolibrary.org/obo/hp.json"
PHENOTYPE_HPOA_URL = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
GENES_TO_PHENOTYPE_URL = "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"

# HPO frequency sub-ontology terms -> our internal frequency categories.
HPO_FREQUENCY_TERMS: dict[str, str] = {
    "HP:0040280": "obligate",
    "HP:0040281": "very_frequent",
    "HP:0040282": "frequent",
    "HP:0040283": "occasional",
    "HP:0040284": "very_rare",
    "HP:0040285": "excluded",
}

# HPO onset sub-ontology terms -> our internal onset categories.
HPO_ONSET_TERMS: dict[str, str] = {
    "HP:0030674": "antenatal",
    "HP:0011460": "antenatal",
    "HP:0011461": "antenatal",
    "HP:0003577": "neonatal",
    "HP:0003623": "neonatal",
    "HP:0003593": "infantile",
    "HP:0011463": "childhood",
    "HP:0003621": "juvenile",
    "HP:0003581": "adult",
}


@dataclass(frozen=True, slots=True)
class HpoaRow:
    """One raw ``phenotype.hpoa`` line, lightly normalized."""

    line_number: int
    disease_id: str
    disease_name: str
    negated: bool
    hpo_id: str
    references: tuple[str, ...]
    evidence: str | None
    onset: str | None
    frequency: str | None
    sex: str | None
    modifier: str | None
    aspect: str


@dataclass
class HpoaPhenotype:
    """A single disease-phenotype annotation row."""

    hpo_id: str
    negative: bool
    frequency_raw: str | None
    frequency_category: str | None
    onset_raw: str | None
    onset_category: str | None
    sex: str | None
    modifier: str | None
    references: list[str] = field(default_factory=list)
    evidence: str | None = None


@dataclass
class DiseaseAnnotations:
    """All HPOA rows for one disease, split by aspect."""

    disease_id: str
    disease_name: str
    phenotypes: list[HpoaPhenotype] = field(default_factory=list)
    negative_phenotypes: list[HpoaPhenotype] = field(default_factory=list)
    inheritance_terms: list[str] = field(default_factory=list)
    clinical_course_terms: list[str] = field(default_factory=list)


def download_hpo_release(client: HttpClient, target_dir: Path | str) -> dict[str, Path]:
    """Download the HPO ontology and annotation files into ``target_dir``."""

    out_dir = Path(target_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "hp.json": HPO_JSON_URL,
        "phenotype.hpoa": PHENOTYPE_HPOA_URL,
        "genes_to_phenotype.txt": GENES_TO_PHENOTYPE_URL,
    }
    written: dict[str, Path] = {}
    for name, url in files.items():
        logger.info("downloading %s", url)
        path = out_dir / name
        path.write_text(client.get(url), encoding="utf-8")
        written[name] = path
    return written


def read_hpoa_header(path: Path | str) -> dict[str, str]:
    """Return the ``#key: value`` metadata lines at the top of a ``phenotype.hpoa``."""

    header: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as file:
        for line in file:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition(":")
            header[key.strip()] = value.strip().strip('"')
    return header


def iter_hpoa_rows(path: Path | str) -> Iterator[HpoaRow]:
    """Yield every data row of ``phenotype.hpoa`` in file order."""

    with Path(path).open("r", encoding="utf-8", newline="") as file:
        line_number = 0
        header: list[str] | None = None
        for raw_line in file:
            line_number += 1
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields, strict=False))
            yield HpoaRow(
                line_number=line_number,
                disease_id=_clean(row.get("database_id")) or "",
                disease_name=_clean(row.get("disease_name")) or "",
                negated=(_clean(row.get("qualifier")) or "").upper() == "NOT",
                hpo_id=_clean(row.get("hpo_id")) or "",
                references=tuple(
                    ref.strip() for ref in (row.get("reference") or "").split(";") if ref.strip()
                ),
                evidence=_clean(row.get("evidence")),
                onset=_clean(row.get("onset")),
                frequency=_clean(row.get("frequency")),
                sex=(_clean(row.get("sex")) or "").upper() or None,
                modifier=_clean(row.get("modifier")),
                aspect=(_clean(row.get("aspect")) or "").upper(),
            )


def read_reference_mask(path: Path | str) -> frozenset[str]:
    """Read a held-out reference list (one ``PMID:n`` or bare PMID per line)."""

    references: set[str] = set()
    with Path(path).open("r", encoding="utf-8") as file:
        for line in file:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            references.add(value if ":" in value else f"PMID:{value}")
    return frozenset(references)


def _clean(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


def parse_hpoa(path: Path | str, disease_ids: set[str]) -> dict[str, DiseaseAnnotations]:
    """Parse ``phenotype.hpoa``, keeping only rows for the requested disease ids."""

    hpoa_path = Path(path)
    wanted = {disease_id.strip() for disease_id in disease_ids}
    result: dict[str, DiseaseAnnotations] = {}

    with hpoa_path.open("r", encoding="utf-8", newline="") as file:
        rows = (line for line in file if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        for row in reader:
            disease_id = (row.get("database_id") or "").strip()
            if disease_id not in wanted:
                continue
            annotations = result.setdefault(
                disease_id,
                DiseaseAnnotations(
                    disease_id=disease_id,
                    disease_name=(row.get("disease_name") or "").strip(),
                ),
            )
            _ingest_row(annotations, row)

    return result


def _ingest_row(annotations: DiseaseAnnotations, row: dict[str, str]) -> None:
    aspect = (row.get("aspect") or "").strip().upper()
    hpo_id = (row.get("hpo_id") or "").strip()
    if not hpo_id:
        return

    if aspect == "I":
        annotations.inheritance_terms.append(hpo_id)
        return
    if aspect == "C":
        annotations.clinical_course_terms.append(hpo_id)
        return
    if aspect != "P":
        return

    negative = (row.get("qualifier") or "").strip().upper() == "NOT"
    frequency_raw = (row.get("frequency") or "").strip() or None
    onset_raw = (row.get("onset") or "").strip() or None
    phenotype = HpoaPhenotype(
        hpo_id=hpo_id,
        negative=negative,
        frequency_raw=frequency_raw,
        frequency_category=_frequency_category(frequency_raw),
        onset_raw=onset_raw,
        onset_category=HPO_ONSET_TERMS.get(onset_raw or ""),
        sex=(row.get("sex") or "").strip() or None,
        modifier=(row.get("modifier") or "").strip() or None,
        references=[ref for ref in (row.get("reference") or "").split(";") if ref.strip()],
        evidence=(row.get("evidence") or "").strip() or None,
    )
    if negative:
        annotations.negative_phenotypes.append(phenotype)
    else:
        annotations.phenotypes.append(phenotype)


def frequency_category_for_probability(probability: float) -> str:
    """Map a probability onto the HPO frequency category whose range contains it."""

    if probability >= 0.99:
        return "obligate"
    if probability >= 0.80:
        return "very_frequent"
    if probability >= 0.30:
        return "frequent"
    if probability >= 0.05:
        return "occasional"
    if probability > 0.0:
        return "very_rare"
    return "excluded"


def _frequency_category(frequency_raw: str | None) -> str | None:
    """Map an HPOA frequency value (HP term, ratio, or percent) to a category."""

    if not frequency_raw:
        return None
    if frequency_raw in HPO_FREQUENCY_TERMS:
        return HPO_FREQUENCY_TERMS[frequency_raw]

    percent = _frequency_percent(frequency_raw)
    if percent is None:
        return None
    return frequency_category_for_probability(percent)


def _frequency_percent(frequency_raw: str) -> float | None:
    value = frequency_raw.strip()
    if value.endswith("%"):
        try:
            return float(value[:-1]) / 100.0
        except ValueError:
            return None
    if "/" in value:
        numerator, _, denominator = value.partition("/")
        try:
            denom = float(denominator)
            return float(numerator) / denom if denom else None
        except ValueError:
            return None
    return None
