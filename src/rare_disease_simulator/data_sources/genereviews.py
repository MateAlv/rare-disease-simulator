"""Best-effort GeneReviews fetching via the NCBI Bookshelf.

GeneReviews chapters are expert-authored and the highest-value text source, but
they are not offered as a clean per-disease API: we ``esearch`` the ``books``
db, resolve a Bookshelf accession via ``esummary``, then retrieve the chapter
page and reduce it to text. This is intentionally tolerant — Bookshelf markup
changes, and a miss should degrade to "no GeneReviews text", not an error.

GeneReviews content is copyrighted: store locally, never redistribute.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass

from rare_disease_simulator.data_sources.http_client import HttpClient, HttpError

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BOOKSHELF_BASE = "https://www.ncbi.nlm.nih.gov/books"

_TAG = re.compile(r"<[^>]+>")
_SCRIPT_STYLE = re.compile(r"<(script|style)\b.*?</\1>", re.IGNORECASE | re.DOTALL)
_WS = re.compile(r"[ \t\r\f\v]+")
_BLANK_LINES = re.compile(r"\n{3,}")


@dataclass
class GeneReviewsChapter:
    """A located GeneReviews chapter with extracted text."""

    accession: str
    title: str
    url: str
    text: str


def fetch_genereviews(
    client: HttpClient, *, query: str
) -> GeneReviewsChapter | None:
    """Find and fetch the best GeneReviews chapter for a query, if any."""

    accession, title = _search_bookshelf(client, query)
    if accession is None:
        logger.info("no GeneReviews chapter found for %r", query)
        return None

    url = f"{BOOKSHELF_BASE}/{accession}/"
    try:
        html = client.get(url)
    except HttpError as exc:
        logger.warning("could not fetch GeneReviews chapter %s: %s", accession, exc)
        return None

    text = html_to_text(html)
    if not text:
        return None
    return GeneReviewsChapter(accession=accession, title=title or accession, url=url, text=text)


def _search_bookshelf(client: HttpClient, query: str) -> tuple[str | None, str | None]:
    body = client.get(
        f"{EUTILS_BASE}/esearch.fcgi",
        {
            "db": "books",
            "term": f"{query} AND GeneReviews[Filter]",
            "retmax": "1",
            "retmode": "json",
        },
    )
    idlist = json.loads(body).get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return None, None

    uid = idlist[0]
    summary = client.get(
        f"{EUTILS_BASE}/esummary.fcgi",
        {"db": "books", "id": uid, "retmode": "json"},
    )
    record = json.loads(summary).get("result", {}).get(uid, {})
    accession = record.get("ackey") or record.get("accession") or record.get("bookacc")
    title = record.get("title")
    return accession, title


def html_to_text(html: str) -> str:
    """Reduce an HTML page to readable text, conservatively."""

    without_scripts = _SCRIPT_STYLE.sub(" ", html)
    without_tags = _TAG.sub(" ", without_scripts)
    unescaped = (
        without_tags.replace("&nbsp;", " ")
        .replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&#x2019;", "'")
    )
    lines = [_WS.sub(" ", line).strip() for line in unescaped.splitlines()]
    collapsed = "\n".join(line for line in lines if line)
    return _BLANK_LINES.sub("\n\n", collapsed).strip()
