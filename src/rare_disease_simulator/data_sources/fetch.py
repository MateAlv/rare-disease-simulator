"""Per-disease source-acquisition orchestrator.

Pulls everything reliably automatable for one disease and writes an auditable
bundle to ``data/raw/<slug>/``:

- ``hpo_annotations.json`` — structured backbone (frequency/onset/sex/negatives)
- ``pubmed.jsonl`` — abstracts (+ PMC-OA full text where available)
- ``genereviews.txt`` — expert chapter text, when found
- ``orphadata.json`` — raw Orphadata resources, when available
- ``text_snippets.jsonl`` — normalized TextSnippet records (the LLM input)
- ``manifest.json`` — sources, URLs, counts, timestamps, and per-source errors

Every source is wrapped so one failure never aborts the run: the goal of this
step is to gather *anything* fetchable, then review what is worth keeping.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import asdict, dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path

from rare_disease_simulator.data_sources import genereviews, hpo_annotations, orphadata, pubmed
from rare_disease_simulator.data_sources.http_client import HttpClient
from rare_disease_simulator.exports.jsonl import write_jsonl
from rare_disease_simulator.llm_extraction.schema import SourceReference, TextSnippet

logger = logging.getLogger(__name__)


@dataclass
class DiseaseQuery:
    """Identifiers used to drive source acquisition for one disease."""

    name: str
    gene: str
    orpha_id: str | None = None
    omim: list[str] = field(default_factory=list)

    @property
    def slug(self) -> str:
        return _slug(f"{self.name}-{self.gene}")

    @property
    def primary_disease_id(self) -> str:
        if self.orpha_id:
            return self.orpha_id
        if self.omim:
            return self.omim[0]
        return f"NAME:{self.slug}"

    @property
    def hpoa_disease_ids(self) -> set[str]:
        ids: set[str] = set(self.omim)
        if self.orpha_id:
            ids.add(self.orpha_id)
        return ids


@dataclass
class FetchManifest:
    """Provenance and outcome summary for one disease fetch."""

    disease_id: str
    disease_name: str
    gene: str
    retrieved_at: str
    sources: dict[str, dict[str, object]] = field(default_factory=dict)
    snippet_count: int = 0


def fetch_disease_sources(
    query: DiseaseQuery,
    *,
    client: HttpClient | None = None,
    output_root: Path | str = "data/raw",
    hpoa_path: Path | str | None = None,
    pubmed_retmax: int = 5,
    fetch_pmc_full_text: bool = True,
) -> FetchManifest:
    """Fetch all available sources for one disease and write the bundle."""

    http = client or HttpClient()
    out_dir = Path(output_root) / query.slug
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    manifest = FetchManifest(
        disease_id=query.primary_disease_id,
        disease_name=query.name,
        gene=query.gene,
        retrieved_at=datetime.now(tz=UTC).isoformat(),
    )
    snippets: list[TextSnippet] = []

    _run_hpoa(query, out_dir, hpoa_path, manifest)
    _run_pubmed(query, http, out_dir, manifest, snippets, pubmed_retmax, fetch_pmc_full_text)
    _run_genereviews(query, http, out_dir, manifest, snippets)
    _run_orphadata(query, http, out_dir, manifest)

    snippets_path = out_dir / "text_snippets.jsonl"
    manifest.snippet_count = write_jsonl(snippets_path, snippets)
    (out_dir / "manifest.json").write_text(
        json.dumps(asdict(manifest), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    logger.info("fetched %s: %d snippet(s)", query.slug, manifest.snippet_count)
    return manifest


def _run_hpoa(
    query: DiseaseQuery,
    out_dir: Path,
    hpoa_path: Path | str | None,
    manifest: FetchManifest,
) -> None:
    if hpoa_path is None or not query.hpoa_disease_ids:
        manifest.sources["hpo_annotations"] = {"status": "skipped", "reason": "no hpoa/ids"}
        return
    try:
        annotations = hpo_annotations.parse_hpoa(hpoa_path, query.hpoa_disease_ids)
        payload = {disease_id: asdict(ann) for disease_id, ann in annotations.items()}
        (out_dir / "hpo_annotations.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        phenotype_total = sum(len(ann.phenotypes) for ann in annotations.values())
        manifest.sources["hpo_annotations"] = {
            "status": "ok",
            "source": str(hpoa_path),
            "diseases_matched": len(annotations),
            "phenotypes": phenotype_total,
        }
    except (OSError, ValueError) as exc:
        logger.warning("hpoa parse failed: %s", exc)
        manifest.sources["hpo_annotations"] = {"status": "error", "error": str(exc)}


def _run_pubmed(
    query: DiseaseQuery,
    http: HttpClient,
    out_dir: Path,
    manifest: FetchManifest,
    snippets: list[TextSnippet],
    retmax: int,
    fetch_full_text: bool,
) -> None:
    try:
        search = f"{query.name} OR {query.gene}"
        pmids = pubmed.search_pubmed(http, query=search, retmax=retmax, reviews_only=True)
        if not pmids:
            pmids = pubmed.search_pubmed(http, query=search, retmax=retmax)
        articles = pubmed.fetch_abstracts(http, pmids)
        if fetch_full_text:
            articles = [pubmed.attach_pmc_full_text(http, article) for article in articles]

        write_jsonl(out_dir / "pubmed.jsonl", [_article_record(a) for a in articles])
        for article in articles:
            snippets.append(_pubmed_snippet(query, article))
        manifest.sources["pubmed"] = {
            "status": "ok",
            "query": search,
            "pmids": pmids,
            "articles": len(articles),
            "open_access": sum(1 for a in articles if a.is_open_access),
        }
    except Exception as exc:  # noqa: BLE001 - best-effort source must not abort the run
        logger.warning("pubmed fetch failed: %s", exc)
        manifest.sources["pubmed"] = {"status": "error", "error": str(exc)}


def _run_genereviews(
    query: DiseaseQuery,
    http: HttpClient,
    out_dir: Path,
    manifest: FetchManifest,
    snippets: list[TextSnippet],
) -> None:
    try:
        chapter = genereviews.fetch_genereviews(http, query=f"{query.gene} {query.name}")
        if chapter is None:
            manifest.sources["genereviews"] = {"status": "not_found"}
            return
        (out_dir / "genereviews.txt").write_text(chapter.text, encoding="utf-8")
        snippets.append(
            TextSnippet(
                snippet_id=f"{query.slug}-genereviews",
                disease_id=query.primary_disease_id,
                disease_name=query.name,
                gene=query.gene,
                source=SourceReference(
                    name="gene_reviews",
                    url=chapter.url,
                    retrieved_at=date.today(),
                    section="GeneReviews chapter",
                    license="copyright NCBI/University of Washington; do not redistribute",
                    redistribution_allowed=False,
                ),
                text=chapter.text,
            )
        )
        manifest.sources["genereviews"] = {
            "status": "ok",
            "accession": chapter.accession,
            "url": chapter.url,
            "chars": len(chapter.text),
        }
    except Exception as exc:  # noqa: BLE001 - best-effort source must not abort the run
        logger.warning("genereviews fetch failed: %s", exc)
        manifest.sources["genereviews"] = {"status": "error", "error": str(exc)}


def _run_orphadata(
    query: DiseaseQuery,
    http: HttpClient,
    out_dir: Path,
    manifest: FetchManifest,
) -> None:
    if not query.orpha_id:
        manifest.sources["orphadata"] = {"status": "skipped", "reason": "no orpha_id"}
        return
    try:
        result = orphadata.fetch_orphadata(http, query.orpha_id)
        (out_dir / "orphadata.json").write_text(
            json.dumps(asdict(result), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        manifest.sources["orphadata"] = {
            "status": "ok" if result.resources else "empty",
            "resources": sorted(result.resources),
            "errors": result.errors,
        }
    except Exception as exc:  # noqa: BLE001 - best-effort source must not abort the run
        logger.warning("orphadata fetch failed: %s", exc)
        manifest.sources["orphadata"] = {"status": "error", "error": str(exc)}


def _article_record(article: pubmed.PubmedArticle) -> dict[str, object]:
    return asdict(article)


def _pubmed_snippet(query: DiseaseQuery, article: pubmed.PubmedArticle) -> TextSnippet:
    section = "PMC open-access full text" if article.is_open_access else "PubMed abstract"
    text = article.full_text or article.abstract or article.title
    return TextSnippet(
        snippet_id=f"{query.slug}-pmid-{article.pmid}",
        disease_id=query.primary_disease_id,
        disease_name=query.name,
        gene=query.gene,
        source=SourceReference(
            name="pubmed",
            url=f"https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/",
            retrieved_at=date.today(),
            section=section,
            license="copyright respective publisher; do not redistribute",
            redistribution_allowed=False,
        ),
        text=text,
    )


def _slug(value: str) -> str:
    lowered = value.lower()
    cleaned = re.sub(r"[^a-z0-9]+", "-", lowered).strip("-")
    return cleaned or "disease"
