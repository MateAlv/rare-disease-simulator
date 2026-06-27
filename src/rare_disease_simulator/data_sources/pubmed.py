"""PubMed / PMC fetching via NCBI E-utilities.

Strategy: relevance-ranked ``esearch`` (optionally restricted to reviews) for a
disease, ``efetch`` for abstracts, then ``elink`` to discover a PMC id and
``efetch`` against PMC for full text — which only succeeds for the open-access
subset, so non-OA articles gracefully degrade to abstract-only.

Abstracts and full text are copyrighted: keep them in local raw storage, never
redistribute. The caller records this in provenance.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from xml.etree import ElementTree

from rare_disease_simulator.data_sources.http_client import HttpClient, HttpError

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class PubmedArticle:
    """A single PubMed record with optional PMC full text."""

    pmid: str
    title: str
    abstract: str
    journal: str | None = None
    year: str | None = None
    publication_types: list[str] = field(default_factory=list)
    mesh_terms: list[str] = field(default_factory=list)
    pmcid: str | None = None
    full_text: str | None = None

    @property
    def is_open_access(self) -> bool:
        return self.full_text is not None


def search_pubmed(
    client: HttpClient,
    *,
    query: str,
    retmax: int = 5,
    reviews_only: bool = False,
) -> list[str]:
    """Return relevance-ranked PMIDs for a query."""

    term = f"({query}) AND review[pt]" if reviews_only else query
    body = client.get(
        f"{EUTILS_BASE}/esearch.fcgi",
        {
            "db": "pubmed",
            "term": term,
            "retmax": str(retmax),
            "retmode": "json",
            "sort": "relevance",
        },
    )
    payload = json.loads(body)
    return list(payload.get("esearchresult", {}).get("idlist", []))


def fetch_abstracts(client: HttpClient, pmids: list[str]) -> list[PubmedArticle]:
    """Fetch and parse abstract records for a list of PMIDs."""

    if not pmids:
        return []
    body = client.get(
        f"{EUTILS_BASE}/efetch.fcgi",
        {"db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml"},
    )
    return parse_pubmed_xml(body)


def parse_pubmed_xml(xml_text: str) -> list[PubmedArticle]:
    """Parse an efetch PubMed XML payload into articles."""

    root = ElementTree.fromstring(xml_text)
    articles: list[PubmedArticle] = []
    for node in root.findall(".//PubmedArticle"):
        pmid = _text(node, ".//PMID")
        if pmid is None:
            continue
        abstract_parts = [
            _join_abstract_section(part) for part in node.findall(".//Abstract/AbstractText")
        ]
        articles.append(
            PubmedArticle(
                pmid=pmid,
                title=_text(node, ".//ArticleTitle") or "",
                abstract="\n".join(part for part in abstract_parts if part).strip(),
                journal=_text(node, ".//Journal/Title"),
                year=_text(node, ".//JournalIssue/PubDate/Year")
                or _text(node, ".//JournalIssue/PubDate/MedlineDate"),
                publication_types=[
                    pt.text.strip()
                    for pt in node.findall(".//PublicationTypeList/PublicationType")
                    if pt.text
                ],
                mesh_terms=[
                    name.text.strip()
                    for name in node.findall(".//MeshHeadingList/MeshHeading/DescriptorName")
                    if name.text
                ],
            )
        )
    return articles


def attach_pmc_full_text(client: HttpClient, article: PubmedArticle) -> PubmedArticle:
    """Best-effort: resolve a PMC id and fetch open-access full text in place."""

    try:
        pmcid = _link_pmcid(client, article.pmid)
        if pmcid is None:
            return article
        article.pmcid = pmcid
        article.full_text = _fetch_pmc_full_text(client, pmcid)
    except HttpError as exc:
        logger.warning("PMC full text unavailable for PMID %s: %s", article.pmid, exc)
    return article


def _link_pmcid(client: HttpClient, pmid: str) -> str | None:
    body = client.get(
        f"{EUTILS_BASE}/elink.fcgi",
        {"dbfrom": "pubmed", "db": "pmc", "id": pmid, "retmode": "json"},
    )
    payload = json.loads(body)
    for linkset in payload.get("linksets", []):
        for db in linkset.get("linksetdbs", []):
            links = db.get("links", [])
            if links:
                return f"PMC{links[0]}"
    return None


def _fetch_pmc_full_text(client: HttpClient, pmcid: str) -> str | None:
    numeric = pmcid.removeprefix("PMC")
    body = client.get(
        f"{EUTILS_BASE}/efetch.fcgi",
        {"db": "pmc", "id": numeric, "retmode": "xml"},
    )
    try:
        root = ElementTree.fromstring(body)
    except ElementTree.ParseError:
        return None
    text = " ".join(segment.strip() for segment in root.itertext() if segment.strip())
    return text or None


def _join_abstract_section(node: ElementTree.Element) -> str:
    label = node.get("Label")
    text = "".join(node.itertext()).strip()
    if not text:
        return ""
    return f"{label}: {text}" if label else text


def _text(node: ElementTree.Element, path: str) -> str | None:
    found = node.find(path)
    if found is None:
        return None
    text = "".join(found.itertext()).strip()
    return text or None
