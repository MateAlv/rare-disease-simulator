"""Reporting artifacts of simulator v0.4 (ADR-0011): ``report-model-v1`` and ``cardinal-v1``.

Both are built by ``diagnostic.ar-training`` from R1-train only; nothing here
is fitted. The simulator samples a patient's true phenotype first, then the
**reported** subset: a term budget ``k`` from the model's histogram, the true
cardinal terms first, and the rest without replacement with probability
proportional to the reporting model's score.

``report-model-v1`` is a JSON object; :class:`ReportModel` is its schema and
every mismatch (a missing or unknown key, an unknown feature kind, a bad
value) is an error when the file is loaded::

    {
      "format": "report-model", "format_version": 1,
      "artifact_id": "report-model-v1",
      "link": "logistic",
      "intercept": -2.0,
      "features": [
        {"name": "frequency", "kind": "frequency", "coefficient": 0.8,
         "unknown_value": 0.5, "center": 0.0, "scale": 1.0},
        {"name": "ic", "kind": "information_content", "coefficient": 0.1,
         "source": "table", "missing_value": 0.0},
        {"name": "ancestors", "kind": "ancestor_count", "coefficient": -0.02,
         "include_self": false},
        {"name": "reportability", "kind": "log_reportability", "coefficient": 0.9,
         "pseudocount": 1.0, "normalize": true},
        {"name": "cardinal", "kind": "cardinal_flag", "coefficient": 1.2}
      ],
      "tables": {"information_content": {"HP:0001250": 3.1}},
      "reportability": {"closure": "ancestors", "patients": 5000,
                        "counts": {"HP:0001250": 812}},
      "term_budget": {"histogram": {"1": 40, "2": 55, "3": 61}}
    }

Feature values, before ``(x - center) / scale`` (defaults 0 and 1):

- ``frequency``: the term's simulation frequency (merged, then the count
  estimator); ``unknown_value`` for a term without one (null: the config's
  ``frequency.unknown_frequency``);
- ``information_content``: ``source: "table"`` reads ``tables.information_content``
  (``missing_value`` when absent); ``source: "profiles"`` is
  ``-log(share of the loaded profiles annotating the term or a descendant)``
  in base ``log_base`` (``"e"`` or ``2``), ``missing_value`` for a term no
  profile carries;
- ``ancestor_count``: the number of the term's ``hp.json`` ancestors
  (``include_self`` adds one);
- ``log_reportability``: ``ln(c + pseudocount)`` with ``c`` the term's
  R1-train patient count from ``reportability.counts`` (0 when absent),
  minus ``ln(patients + pseudocount)`` when ``normalize``;
- ``cardinal_flag``: 1 when the term is cardinal for the entity, else 0.

The score is ``1 / (1 + exp(-(intercept + sum(coefficient * value))))``.

``cardinal-v1`` is a TSV with the header ``disease_id  hpo_id  kind  source``;
a term is cardinal for an entity when any of its profile ids lists it.
"""

from __future__ import annotations

import csv
import json
import math
import random
import re
from bisect import bisect_right
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from rare_disease_simulator.data_sources.hpo import HpoOntology
from rare_disease_simulator.profiles.schema import DiseaseProfile

REPORT_MODEL_FORMAT = "report-model"
REPORT_MODEL_FORMAT_VERSION = 1
CARDINAL_COLUMNS = ("disease_id", "hpo_id", "kind", "source")
CARDINAL_KINDS = frozenset({"diagnostic_criterion", "pathognomonic", "cardinal_proxy"})
_HPO_ID_RE = re.compile(r"^HP:\d{7}$")
_DISEASE_ID_RE = re.compile(r"^(OMIM|ORPHA|DECIPHER):\d+$")


class ReportingArtifactError(ValueError):
    """A report-model or cardinal file that does not match its format."""


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class _FeatureBase(_Strict):
    name: str = Field(min_length=1)
    coefficient: float
    center: float = 0.0
    scale: float = Field(default=1.0, gt=0.0)
    description: str | None = None


class FrequencyFeature(_FeatureBase):
    kind: Literal["frequency"]
    unknown_value: float | None = Field(default=None, ge=0.0, le=1.0)


class InformationContentFeature(_FeatureBase):
    kind: Literal["information_content"]
    source: Literal["table", "profiles"]
    missing_value: float
    log_base: Literal["e", "2"] | None = None


class AncestorCountFeature(_FeatureBase):
    kind: Literal["ancestor_count"]
    include_self: bool


class LogReportabilityFeature(_FeatureBase):
    kind: Literal["log_reportability"]
    pseudocount: float = Field(gt=0.0)
    normalize: bool


class CardinalFlagFeature(_FeatureBase):
    kind: Literal["cardinal_flag"]


Feature = Annotated[
    FrequencyFeature
    | InformationContentFeature
    | AncestorCountFeature
    | LogReportabilityFeature
    | CardinalFlagFeature,
    Field(discriminator="kind"),
]


class Reportability(_Strict):
    closure: Literal["ancestors"]
    patients: int = Field(gt=0)
    counts: dict[str, int]

    @field_validator("counts")
    @classmethod
    def _check_counts(cls, value: dict[str, int]) -> dict[str, int]:
        bad = [key for key, count in value.items() if not _HPO_ID_RE.match(key) or count < 0]
        if bad:
            raise ValueError(f"counts need HPO ids and non-negative counts: {bad[:5]}")
        return value


class TermBudget(_Strict):
    histogram: dict[int, int]
    description: str | None = None

    @field_validator("histogram")
    @classmethod
    def _check_histogram(cls, value: dict[int, int]) -> dict[int, int]:
        if any(k < 0 or count < 0 for k, count in value.items()):
            raise ValueError("histogram keys and counts must be non-negative")
        if sum(count for k, count in value.items() if k >= 1) <= 0:
            raise ValueError("histogram has no case with a budget of at least 1")
        return value


class ReportModel(_Strict):
    """Schema of ``report-model-v1``."""

    format: Literal["report-model"]
    format_version: Literal[1]
    artifact_id: str = Field(min_length=1)
    link: Literal["logistic"]
    intercept: float
    features: list[Feature] = Field(min_length=1)
    tables: dict[str, dict[str, float]] = Field(default_factory=dict)
    reportability: Reportability | None = None
    term_budget: TermBudget
    provenance: dict[str, Any] | None = None
    notes: str | None = None

    @field_validator("features")
    @classmethod
    def _unique(cls, value: list[Any]) -> list[Any]:
        names = [feature.name for feature in value]
        kinds = [feature.kind for feature in value]
        if len(set(names)) != len(names) or len(set(kinds)) != len(kinds):
            raise ValueError("feature names and kinds must be unique")
        return value


def load_report_model(path: Path | str) -> ReportModel:
    """Read ``report-model-v1``; any schema mismatch raises :class:`ReportingArtifactError`."""

    model_path = Path(path)
    try:
        data = json.loads(model_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportingArtifactError(f"{model_path}: cannot read JSON: {exc}") from exc
    try:
        model = ReportModel.model_validate(data)
    except ValidationError as exc:
        raise ReportingArtifactError(
            f"{model_path}: does not match {REPORT_MODEL_FORMAT} "
            f"v{REPORT_MODEL_FORMAT_VERSION}:\n{exc}"
        ) from exc
    kinds = {feature.kind for feature in model.features}
    if "log_reportability" in kinds and model.reportability is None:
        raise ReportingArtifactError(f"{model_path}: log_reportability needs 'reportability'")
    for feature in model.features:
        if isinstance(feature, InformationContentFeature):
            if feature.source == "table" and "information_content" not in model.tables:
                raise ReportingArtifactError(
                    f"{model_path}: information_content from a table needs "
                    "tables.information_content"
                )
            if feature.source == "profiles" and feature.log_base is None:
                raise ReportingArtifactError(
                    f"{model_path}: information_content from profiles needs log_base"
                )
    return model


@dataclass(frozen=True)
class CardinalIndex:
    """``cardinal-v1``: cardinal terms per disease id, with their kinds."""

    terms: Mapping[str, frozenset[str]]
    kinds: Counter[str] = field(default_factory=Counter)
    unresolved_terms: int = 0

    def terms_for(self, disease_ids: Iterable[str]) -> frozenset[str]:
        """Terms cardinal for an entity: listed under any of its profile ids."""

        result: set[str] = set()
        for disease_id in disease_ids:
            result |= self.terms.get(disease_id, frozenset())
        return frozenset(result)

    def summary(self) -> dict[str, Any]:
        return {
            "diseases": len(self.terms),
            "pairs": sum(len(terms) for terms in self.terms.values()),
            "rows_by_kind": dict(sorted(self.kinds.items())),
            "unresolved_terms": self.unresolved_terms,
        }


def load_cardinal(path: Path | str, ontology: HpoOntology | None = None) -> CardinalIndex:
    """Read ``cardinal-v1``; a wrong header, kind or id raises :class:`ReportingArtifactError`.

    With an ontology, alternative and obsolete ids resolve to the current term;
    rows whose term does not resolve are dropped and counted.
    """

    cardinal_path = Path(path)
    terms: dict[str, set[str]] = {}
    kinds: Counter[str] = Counter()
    unresolved = 0
    with cardinal_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(header) != CARDINAL_COLUMNS:
            raise ReportingArtifactError(
                f"{cardinal_path}: header must be {list(CARDINAL_COLUMNS)} "
                f"(tab-separated), got {header!r}"
            )
        for line_number, row in enumerate(reader, start=2):
            if not row or not any(row):
                continue
            if len(row) != len(CARDINAL_COLUMNS):
                raise ReportingArtifactError(
                    f"{cardinal_path}:{line_number}: expected {len(CARDINAL_COLUMNS)} columns"
                )
            disease_id, hpo_id, kind, source = row
            if not _DISEASE_ID_RE.match(disease_id):
                raise ReportingArtifactError(
                    f"{cardinal_path}:{line_number}: bad disease id {disease_id!r}"
                )
            if not _HPO_ID_RE.match(hpo_id):
                raise ReportingArtifactError(
                    f"{cardinal_path}:{line_number}: bad HPO id {hpo_id!r}"
                )
            if kind not in CARDINAL_KINDS:
                raise ReportingArtifactError(
                    f"{cardinal_path}:{line_number}: unknown kind {kind!r} "
                    f"(expected one of {sorted(CARDINAL_KINDS)})"
                )
            if not source:
                raise ReportingArtifactError(f"{cardinal_path}:{line_number}: empty source")
            if ontology is not None:
                resolved = ontology.resolve(hpo_id).hpo_id
                if resolved is None:
                    unresolved += 1
                    continue
                hpo_id = resolved
            terms.setdefault(disease_id, set()).add(hpo_id)
            kinds[kind] += 1
    return CardinalIndex(
        terms={disease_id: frozenset(ids) for disease_id, ids in sorted(terms.items())},
        kinds=kinds,
        unresolved_terms=unresolved,
    )


class TermScorer:
    """The reporting model's score of a term, computed as the model's features define."""

    def __init__(
        self,
        model: ReportModel,
        *,
        ontology: HpoOntology | None,
        profiles: Sequence[DiseaseProfile] = (),
        unknown_frequency: float = 0.5,
    ) -> None:
        self.model = model
        self._ontology = ontology
        self._unknown_frequency = unknown_frequency
        kinds = {feature.kind for feature in model.features}
        if "ancestor_count" in kinds and ontology is None:
            raise ReportingArtifactError("the report model's ancestor_count feature needs hp.json")
        self._profile_ic: dict[str, float] = {}
        for feature in model.features:
            if isinstance(feature, InformationContentFeature) and feature.source == "profiles":
                self._profile_ic = _profile_information_content(
                    profiles, ontology, math.e if feature.log_base == "e" else 2.0
                )
        self._cache: dict[tuple[str, float | None, bool], float] = {}

    def values(self, hpo_id: str, frequency: float | None, cardinal: bool) -> dict[str, float]:
        """Raw feature values (before centering and scaling), by feature name."""

        values: dict[str, float] = {}
        for feature in self.model.features:
            values[feature.name] = self._value(feature, hpo_id, frequency, cardinal)
        return values

    def score(self, hpo_id: str, frequency: float | None, cardinal: bool) -> float:
        key = (hpo_id, frequency, cardinal)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        linear = self.model.intercept
        for feature in self.model.features:
            value = self._value(feature, hpo_id, frequency, cardinal)
            linear += feature.coefficient * (value - feature.center) / feature.scale
        score = _logistic(linear)
        self._cache[key] = score
        return score

    def _value(self, feature: Any, hpo_id: str, frequency: float | None, cardinal: bool) -> float:
        if isinstance(feature, FrequencyFeature):
            if frequency is not None:
                return frequency
            if feature.unknown_value is not None:
                return feature.unknown_value
            return self._unknown_frequency
        if isinstance(feature, InformationContentFeature):
            table = (
                self.model.tables["information_content"]
                if feature.source == "table"
                else self._profile_ic
            )
            return table.get(hpo_id, feature.missing_value)
        if isinstance(feature, AncestorCountFeature):
            assert self._ontology is not None
            return len(self._ontology.get_ancestor_set(hpo_id)) + (
                1 if feature.include_self else 0
            )
        if isinstance(feature, LogReportabilityFeature):
            assert self.model.reportability is not None
            count = self.model.reportability.counts.get(hpo_id, 0)
            value = math.log(count + feature.pseudocount)
            if feature.normalize:
                value -= math.log(self.model.reportability.patients + feature.pseudocount)
            return value
        if isinstance(feature, CardinalFlagFeature):
            return 1.0 if cardinal else 0.0
        raise ReportingArtifactError(f"unsupported feature {feature!r}")


class BudgetSampler:
    """Draws a term budget from the model's histogram, over budgets of at least 1.

    A case must show at least one positive, so the histogram's 0 bin (cases
    that reported no present term) is dropped and the rest renormalised.
    """

    def __init__(self, histogram: Mapping[int, int]) -> None:
        items = sorted((k, count) for k, count in histogram.items() if k >= 1 and count > 0)
        self.budgets = [k for k, _ in items]
        total = sum(count for _, count in items)
        self.probabilities = {k: count / total for k, count in items}
        cumulative = 0.0
        self._cumulative: list[float] = []
        for _, count in items:
            cumulative += count / total
            self._cumulative.append(cumulative)
        self.dropped_zero_cases = histogram.get(0, 0)

    def draw(self, rng: random.Random) -> int:
        index = bisect_right(self._cumulative, rng.random() * self._cumulative[-1])
        return self.budgets[min(index, len(self.budgets) - 1)]

    def mean(self) -> float:
        return sum(k * p for k, p in self.probabilities.items())


@dataclass
class Reporting:
    """Everything report-model emission needs: the model, its scorer, budgets and cardinals."""

    model: ReportModel
    scorer: TermScorer
    budget: BudgetSampler
    cardinal: CardinalIndex | None = None

    @classmethod
    def build(
        cls,
        model: ReportModel,
        *,
        ontology: HpoOntology | None,
        profiles: Sequence[DiseaseProfile] = (),
        cardinal: CardinalIndex | None = None,
        unknown_frequency: float = 0.5,
    ) -> Reporting:
        return cls(
            model=model,
            scorer=TermScorer(
                model, ontology=ontology, profiles=profiles, unknown_frequency=unknown_frequency
            ),
            budget=BudgetSampler(model.term_budget.histogram),
            cardinal=cardinal,
        )

    def cardinal_terms(self, disease_ids: Iterable[str]) -> frozenset[str]:
        return self.cardinal.terms_for(disease_ids) if self.cardinal is not None else frozenset()

    def summary(self) -> dict[str, Any]:
        return {
            "report_model": self.model.artifact_id,
            "features": [
                {"name": f.name, "kind": f.kind, "coefficient": f.coefficient}
                for f in self.model.features
            ],
            "intercept": self.model.intercept,
            "budget_probabilities": {
                str(k): round(p, 6) for k, p in self.budget.probabilities.items()
            },
            "budget_mean": round(self.budget.mean(), 4),
            "budget_zero_cases_dropped": self.budget.dropped_zero_cases,
            "cardinal": self.cardinal.summary() if self.cardinal is not None else None,
        }


def _profile_information_content(
    profiles: Sequence[DiseaseProfile], ontology: HpoOntology | None, base: float
) -> dict[str, float]:
    counts: Counter[str] = Counter()
    for profile in profiles:
        terms: set[str] = set()
        for phenotype in profile.phenotypes:
            terms.add(phenotype.hpo_id)
            if ontology is not None:
                terms |= ontology.get_ancestor_set(phenotype.hpo_id)
        counts.update(terms)
    total = len(profiles)
    return {
        hpo_id: -math.log(count / total, base) for hpo_id, count in sorted(counts.items())
    }


def _logistic(value: float) -> float:
    if value >= 0:
        return 1.0 / (1.0 + math.exp(-value))
    exp = math.exp(value)
    return exp / (1.0 + exp)
