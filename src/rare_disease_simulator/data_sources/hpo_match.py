"""Local mention-to-HPO normalization.

The LLM extractor emits phenotype *mention strings* (e.g. "vertical gaze
palsy"), never HPO IDs. This matcher resolves a mention to an HPO term using
the local ontology's labels and synonyms, so the ontology stays owned by code
and the model cannot invent identifiers.

Matching is intentionally conservative: an exact normalized match against a
label or synonym is preferred; a fuzzy fallback is only offered above a
threshold and is surfaced for review rather than trusted blindly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from difflib import SequenceMatcher

from rare_disease_simulator.data_sources.hpo import HpoOntology

_NON_ALNUM = re.compile(r"[^a-z0-9]+")


def normalize_mention(text: str) -> str:
    """Lowercase, strip punctuation, and collapse whitespace for matching."""

    lowered = text.lower().strip()
    collapsed = _NON_ALNUM.sub(" ", lowered)
    return " ".join(collapsed.split())


@dataclass(frozen=True)
class HpoMatch:
    """A resolved mention-to-HPO match candidate."""

    hpo_id: str
    label: str
    score: float
    method: str
    matched_text: str


class HpoMatcher:
    """Resolve free-text phenotype mentions to HPO ids via labels and synonyms."""

    def __init__(self, ontology: HpoOntology, *, fuzzy_threshold: float = 0.9) -> None:
        self.ontology = ontology
        self.fuzzy_threshold = fuzzy_threshold
        self._exact_index: dict[str, str] = {}
        self._build_index()

    def _build_index(self) -> None:
        for hpo_id, term in self.ontology.terms.items():
            for surface in (term.label, *term.synonyms):
                key = normalize_mention(surface)
                # First writer wins, but a label always overrides a synonym.
                if key not in self._exact_index or surface == term.label:
                    self._exact_index[key] = hpo_id

    def match(self, mention: str) -> HpoMatch | None:
        """Return the best match for a mention, or None if nothing is confident."""

        normalized = normalize_mention(mention)
        if not normalized:
            return None

        hpo_id = self._exact_index.get(normalized)
        if hpo_id is not None:
            return HpoMatch(
                hpo_id=hpo_id,
                label=self.ontology.get_label(hpo_id) or mention,
                score=1.0,
                method="exact",
                matched_text=mention,
            )
        return self._fuzzy_match(mention, normalized)

    def _fuzzy_match(self, mention: str, normalized: str) -> HpoMatch | None:
        best_score = 0.0
        best_id: str | None = None
        for key, hpo_id in self._exact_index.items():
            score = SequenceMatcher(None, normalized, key).ratio()
            if score > best_score:
                best_score = score
                best_id = hpo_id
        if best_id is None or best_score < self.fuzzy_threshold:
            return None
        return HpoMatch(
            hpo_id=best_id,
            label=self.ontology.get_label(best_id) or mention,
            score=round(best_score, 3),
            method="fuzzy",
            matched_text=mention,
        )
