import json

from rare_disease_simulator.data_sources.fetch import DiseaseQuery, fetch_disease_sources
from rare_disease_simulator.data_sources.genereviews import html_to_text
from rare_disease_simulator.data_sources.orphadata import normalize_orpha_code
from rare_disease_simulator.data_sources.pubmed import parse_pubmed_xml
from tests.fixtures.readers import fixture_path

PUBMED_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle><MedlineCitation>
    <PMID>111</PMID>
    <Article>
      <ArticleTitle>Niemann-Pick type C review</ArticleTitle>
      <Abstract><AbstractText Label="BACKGROUND">Ataxia is common.</AbstractText></Abstract>
      <Journal><Title>J Test</Title>
        <JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue></Journal>
      <PublicationTypeList><PublicationType>Review</PublicationType></PublicationTypeList>
    </Article>
    <MeshHeadingList><MeshHeading><DescriptorName>Niemann-Pick Diseases</DescriptorName>
      </MeshHeading></MeshHeadingList>
  </MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


class FakeHttpClient:
    """URL/param-dispatching stand-in for HttpClient used in fetch tests."""

    def get(self, url: str, params: dict | None = None) -> str:
        params = params or {}
        if "esearch.fcgi" in url and params.get("db") == "pubmed":
            return json.dumps({"esearchresult": {"idlist": ["111"]}})
        if "efetch.fcgi" in url and params.get("db") == "pubmed":
            return PUBMED_XML
        if "elink.fcgi" in url:
            return json.dumps({"linksets": [{"linksetdbs": [{"links": []}]}]})
        if "esearch.fcgi" in url and params.get("db") == "books":
            return json.dumps({"esearchresult": {"idlist": ["222"]}})
        if "esummary.fcgi" in url and params.get("db") == "books":
            return json.dumps({"result": {"222": {"ackey": "NBK1296", "title": "NPC GeneReviews"}}})
        if "/books/" in url:
            return "<html><body><h1>Clinical</h1><p>Vertical gaze palsy occurs.</p></body></html>"
        if "api.orphadata.com" in url:
            return json.dumps({"data": {"results": "ok"}})
        raise AssertionError(f"unexpected URL: {url} {params}")


def test_parse_pubmed_xml_extracts_fields() -> None:
    articles = parse_pubmed_xml(PUBMED_XML)

    assert len(articles) == 1
    article = articles[0]
    assert article.pmid == "111"
    assert article.title == "Niemann-Pick type C review"
    assert "Ataxia is common" in article.abstract
    assert article.year == "2020"
    assert "Review" in article.publication_types
    assert article.mesh_terms == ["Niemann-Pick Diseases"]
    assert article.is_open_access is False


def test_html_to_text_strips_markup() -> None:
    text = html_to_text("<h1>Clinical</h1><p>Vertical gaze palsy &amp; ataxia.</p>")

    assert "Clinical" in text
    assert "Vertical gaze palsy & ataxia." in text
    assert "<" not in text


def test_normalize_orpha_code_variants() -> None:
    assert normalize_orpha_code("ORPHA:646") == "646"
    assert normalize_orpha_code("ORPHA646") == "646"
    assert normalize_orpha_code("646") == "646"


def test_fetch_disease_sources_writes_bundle(tmp_path) -> None:
    query = DiseaseQuery(
        name="Niemann-Pick disease type C1",
        gene="NPC1",
        orpha_id="ORPHA:646",
        omim=["OMIM:257220"],
    )

    manifest = fetch_disease_sources(
        query,
        client=FakeHttpClient(),
        output_root=tmp_path,
        hpoa_path=fixture_path("phenotype_mini.hpoa"),
    )

    bundle = tmp_path / query.slug
    assert (bundle / "manifest.json").exists()
    assert (bundle / "text_snippets.jsonl").exists()
    assert manifest.sources["hpo_annotations"]["status"] == "ok"
    assert manifest.sources["hpo_annotations"]["phenotypes"] == 3
    assert manifest.sources["pubmed"]["status"] == "ok"
    assert manifest.sources["genereviews"]["status"] == "ok"
    assert manifest.sources["orphadata"]["status"] == "ok"
    # one PubMed abstract + one GeneReviews chapter
    assert manifest.snippet_count == 2
