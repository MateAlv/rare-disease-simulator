"""Confounder diseases and their discriminative terms, for "asked and absent" negatives.

A clinician rules out look-alike diseases by asking about their hallmark
features. The index finds, for a true disease, the ``top_n`` most similar
diseases by IC-weighted HPO overlap and offers their terms that the true
disease does not annotate at any ontology level.

Similarity is a cosine over information-content weights: each disease is the
set of its positive terms plus their ancestors (when an ontology is given),
term IC is ``-ln(share of diseases carrying the term)``, and only terms with
IC >= ``min_information_content`` enter the inverted index, so the posting
lists stay short. Everything is computed lazily per queried disease and cached.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import DiseaseProfile, SexRestriction


@dataclass(frozen=True)
class ConfounderTerm:
    """A term a confounder disease has and the true disease never annotates."""

    hpo_id: str
    label: str
    weight: float
    confounder_id: str
    sex_restriction: SexRestriction | None


class ConfounderIndex:
    """IC-weighted inverted index over disease profiles."""

    def __init__(
        self,
        profiles: Sequence[DiseaseProfile],
        ontology: HpoOntology | None = None,
        *,
        top_n: int = 10,
        min_information_content: float = 2.0,
        unknown_frequency: float = 0.5,
    ) -> None:
        self._ontology = ontology
        self._top_n = top_n
        self._unknown_frequency = unknown_frequency
        self._profiles = {profile.disease_id: profile for profile in profiles}

        expanded: dict[str, frozenset[str]] = {}
        document_frequency: dict[str, int] = defaultdict(int)
        for disease_id in sorted(self._profiles):
            terms = self._expanded_terms(self._profiles[disease_id])
            expanded[disease_id] = terms
            for hpo_id in terms:
                document_frequency[hpo_id] += 1
        total = max(len(expanded), 1)
        self.information_content = {
            hpo_id: math.log(total / count) for hpo_id, count in document_frequency.items()
        }

        self._vectors: dict[str, dict[str, float]] = {}
        postings: dict[str, list[str]] = defaultdict(list)
        for disease_id, terms in expanded.items():
            # Sorted, so float sums (norms, cosines) do not depend on the
            # process's string-hash order and near-ties rank the same every run.
            vector = {
                hpo_id: self.information_content[hpo_id]
                for hpo_id in sorted(terms)
                if self.information_content[hpo_id] >= min_information_content
            }
            self._vectors[disease_id] = vector
            for hpo_id in vector:
                postings[hpo_id].append(disease_id)
        self._postings = dict(postings)
        self._min_information_content = min_information_content
        self._norms = {
            disease_id: math.sqrt(sum(weight * weight for weight in vector.values()))
            for disease_id, vector in self._vectors.items()
        }
        self._confounder_cache: dict[str, list[tuple[str, float]]] = {}
        self._term_cache: dict[str, list[ConfounderTerm]] = {}
        self._merged_cache: dict[tuple[str, frozenset[str]], list[ConfounderTerm]] = {}

    def confounders(self, disease_id: str) -> list[tuple[str, float]]:
        """The ``top_n`` most similar other diseases, as ``(disease_id, cosine)``."""

        cached = self._confounder_cache.get(disease_id)
        if cached is not None:
            return cached
        ranked = self._rank(self._vectors.get(disease_id, {}), frozenset({disease_id}))
        self._confounder_cache[disease_id] = ranked
        return ranked

    def _rank(
        self, vector: dict[str, float], exclude: frozenset[str]
    ) -> list[tuple[str, float]]:
        norm = math.sqrt(sum(weight * weight for weight in vector.values()))
        scores: dict[str, float] = defaultdict(float)
        for hpo_id, weight in vector.items():
            for other_id in self._postings.get(hpo_id, ()):
                if other_id not in exclude:
                    scores[other_id] += weight * self._vectors[other_id][hpo_id]
        return sorted(
            (
                (other_id, score / (norm * self._norms[other_id]))
                for other_id, score in scores.items()
                if norm and self._norms[other_id]
            ),
            key=lambda item: (-item[1], item[0]),
        )[: self._top_n]

    def candidate_terms(self, disease_id: str) -> list[ConfounderTerm]:
        """Confounder terms never annotated to the disease, weighted by similarity x frequency.

        A term is excluded when it is, or is an ancestor or descendant of, any
        positive term of the true disease.
        """

        cached = self._term_cache.get(disease_id)
        if cached is not None:
            return cached
        profile = self._profiles.get(disease_id)
        if profile is None:
            return []
        terms = self._terms(profile, self.confounders(disease_id))
        self._term_cache[disease_id] = terms
        return terms

    def candidate_terms_for_profile(
        self, profile: DiseaseProfile, exclude: frozenset[str]
    ) -> list[ConfounderTerm]:
        """Confounder terms of a profile outside the index (e.g. a merged entity).

        Its similarity uses the index's information content, and the diseases
        in ``exclude`` (the entity's own profiles) are never confounders.
        """

        key = (profile.disease_id, exclude)
        cached = self._merged_cache.get(key)
        if cached is not None:
            return cached
        vector = {
            hpo_id: self.information_content[hpo_id]
            for hpo_id in sorted(self._expanded_terms(profile))
            if self.information_content.get(hpo_id, 0.0) >= self._min_information_content
            and hpo_id in self.information_content
        }
        terms = self._terms(profile, self._rank(vector, exclude))
        self._merged_cache[key] = terms
        return terms

    def _terms(
        self, profile: DiseaseProfile, confounders: list[tuple[str, float]]
    ) -> list[ConfounderTerm]:
        blocked = self.annotation_closure(profile)
        weights: dict[str, float] = defaultdict(float)
        best: dict[str, tuple[float, str]] = {}
        labels: dict[str, str] = {}
        restrictions: dict[str, set[SexRestriction | None]] = defaultdict(set)
        for other_id, similarity in confounders:
            for phenotype in self._profiles[other_id].phenotypes:
                hpo_id = phenotype.hpo_id
                if hpo_id in blocked:
                    continue
                if self._ontology is not None and not self._ontology.is_phenotypic_abnormality(
                    hpo_id
                ):
                    continue
                frequency = (
                    phenotype.frequency_estimate
                    if phenotype.frequency_estimate is not None
                    else self._unknown_frequency
                )
                contribution = similarity * frequency
                weights[hpo_id] += contribution
                if hpo_id not in best or contribution > best[hpo_id][0]:
                    best[hpo_id] = (contribution, other_id)
                labels.setdefault(hpo_id, phenotype.label)
                restrictions[hpo_id].add(phenotype.sex_restriction)
        terms = [
            ConfounderTerm(
                hpo_id=hpo_id,
                label=labels[hpo_id],
                weight=weights[hpo_id],
                confounder_id=best[hpo_id][1],
                sex_restriction=(
                    next(iter(restrictions[hpo_id])) if len(restrictions[hpo_id]) == 1 else None
                ),
            )
            for hpo_id in sorted(weights)
            if weights[hpo_id] > 0.0
        ]
        return terms

    def annotation_closure(self, profile: DiseaseProfile) -> frozenset[str]:
        """Positive terms of a disease plus all their ancestors and descendants."""

        closure: set[str] = set()
        for phenotype in profile.phenotypes:
            closure.add(phenotype.hpo_id)
            if self._ontology is not None:
                closure |= self._ontology.get_ancestor_set(phenotype.hpo_id)
                closure |= self._ontology.get_descendant_set(phenotype.hpo_id)
        return frozenset(closure)

    def _expanded_terms(self, profile: DiseaseProfile) -> frozenset[str]:
        terms: set[str] = set()
        for phenotype in profile.phenotypes:
            terms.add(phenotype.hpo_id)
            if self._ontology is not None:
                terms |= self._ontology.get_ancestor_set(phenotype.hpo_id)
        return frozenset(terms)
