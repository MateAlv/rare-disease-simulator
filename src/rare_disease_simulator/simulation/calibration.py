"""Calibration files: override simulator knobs by their config paths.

A calibration file is a flat JSON object produced outside this repo (the
training repo fits it on R1-train only). Each key is a dotted path into
:class:`SimulationConfig`, e.g. ``presets.medium.negatives_mean`` or
``missingness.sex_unknown``; the accepted paths are generated from the schema
by :func:`calibration_keys`. Unknown keys are errors. Two kinds of extra keys
are allowed:

- ``noise.vocabulary_path``: a noise vocabulary TSV (``hpo_id``, optional
  ``label`` and ``weight``), relative to the calibration file;
- keys starting with ``_`` (e.g. ``_provenance``): notes, recorded but unused.

Run-shape fields (``seed``, ``difficulties``,
``cases_per_disease_per_difficulty``) are set on the command line, not here.
"""

from __future__ import annotations

import json
import types
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Union, get_args, get_origin

from pydantic import BaseModel, ValidationError

from rare_disease_simulator.build_info import sha256_file
from rare_disease_simulator.simulation.schema import SimulationConfig

NOISE_VOCABULARY_KEY = "noise.vocabulary_path"
RUN_SHAPE_FIELDS = frozenset({"seed", "difficulties", "cases_per_disease_per_difficulty"})


class CalibrationError(ValueError):
    """A calibration file that cannot be applied."""


@dataclass(frozen=True)
class CalibrationKey:
    """One knob a calibration file may set."""

    path: str
    type: str
    default: Any
    description: str | None = None


@dataclass
class Calibration:
    """A parsed calibration file."""

    path: Path
    sha256: str
    overrides: dict[str, Any]
    noise_vocabulary: Path | None = None
    notes: dict[str, Any] = field(default_factory=dict)

    def record(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "sha256": self.sha256,
            "overrides": self.overrides,
            "noise_vocabulary": str(self.noise_vocabulary) if self.noise_vocabulary else None,
            "notes": self.notes,
        }


def calibration_keys() -> list[CalibrationKey]:
    """Every settable path of :class:`SimulationConfig`, with type and default."""

    keys: list[CalibrationKey] = []
    _walk(SimulationConfig, SimulationConfig().model_dump(mode="json"), "", keys)
    keys.append(
        CalibrationKey(
            path=NOISE_VOCABULARY_KEY,
            type="path",
            default=None,
            description="Noise vocabulary TSV (hpo_id, label, weight), relative to the file.",
        )
    )
    return keys


def format_calibration_keys(keys: list[CalibrationKey] | None = None) -> str:
    """Markdown table of the accepted keys (the README's list is this output)."""

    lines = ["| Key | Type | Default |", "| --- | --- | --- |"]
    for key in keys or calibration_keys():
        default = json.dumps(key.default)
        lines.append(f"| `{key.path}` | {key.type} | `{default}` |")
    return "\n".join(lines)


def load_calibration(path: Path | str) -> Calibration:
    """Read and check a calibration file; unknown keys raise :class:`CalibrationError`."""

    calibration_path = Path(path)
    try:
        data = json.loads(calibration_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise CalibrationError(f"{calibration_path}: invalid JSON: {exc.msg}") from exc
    if not isinstance(data, dict):
        raise CalibrationError(f"{calibration_path}: expected a JSON object of dotted keys")

    known = {key.path for key in calibration_keys()}
    unknown = sorted(key for key in data if key not in known and not key.startswith("_"))
    if unknown:
        raise CalibrationError(
            f"{calibration_path}: unknown calibration key(s): {', '.join(unknown)} "
            "(run `rare-disease-simulator calibration-keys` for the accepted list)"
        )

    noise_vocabulary = None
    if NOISE_VOCABULARY_KEY in data:
        value = data[NOISE_VOCABULARY_KEY]
        if not isinstance(value, str) or not value:
            raise CalibrationError(f"{calibration_path}: {NOISE_VOCABULARY_KEY} must be a path")
        noise_vocabulary = Path(value)
        if not noise_vocabulary.is_absolute():
            noise_vocabulary = calibration_path.parent / noise_vocabulary
        if not noise_vocabulary.is_file():
            raise CalibrationError(f"{calibration_path}: noise vocabulary not found: {value}")

    overrides = {
        key: value
        for key, value in sorted(data.items())
        if not key.startswith("_") and key != NOISE_VOCABULARY_KEY
    }
    calibration = Calibration(
        path=calibration_path,
        sha256=sha256_file(calibration_path),
        overrides=overrides,
        noise_vocabulary=noise_vocabulary,
        notes={key: value for key, value in sorted(data.items()) if key.startswith("_")},
    )
    apply_calibration(SimulationConfig(), calibration.overrides)
    return calibration


def apply_calibration(
    config: SimulationConfig, overrides: Mapping[str, Any]
) -> SimulationConfig:
    """Return ``config`` with each dotted-path override set, re-validated."""

    data = config.model_dump(mode="json")
    for path, value in overrides.items():
        target = data
        parts = path.split(".")
        for part in parts[:-1]:
            if not isinstance(target.get(part), dict):
                raise CalibrationError(f"cannot set {path}: {part} is not a section")
            target = target[part]
        target[parts[-1]] = value
    try:
        return SimulationConfig.model_validate(data)
    except ValidationError as exc:
        raise CalibrationError(f"calibration values are invalid: {exc}") from exc


def _walk(model: type[BaseModel], defaults: Any, prefix: str, keys: list[CalibrationKey]) -> None:
    for name, info in model.model_fields.items():
        if not prefix and name in RUN_SHAPE_FIELDS:
            continue
        path = f"{prefix}{name}"
        _add(info.annotation, defaults[name], path, info.description, keys)


def _add(
    annotation: Any, default: Any, path: str, description: str | None, keys: list[CalibrationKey]
) -> None:
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        _walk(annotation, default, f"{path}.", keys)
        return
    if get_origin(annotation) is dict:
        key_type, value_type = get_args(annotation)
        if get_origin(key_type) is Literal:
            # An "unknown" onset category has no window and is never drawn.
            for key in (key for key in get_args(key_type) if key != "unknown"):
                _add(value_type, default.get(key), f"{path}.{key}", description, keys)
            return
    keys.append(
        CalibrationKey(
            path=path, type=_type_name(annotation), default=default, description=description
        )
    )


def _type_name(annotation: Any) -> str:
    origin = get_origin(annotation)
    if origin in (Union, types.UnionType):
        return " or ".join(_type_name(arg) for arg in get_args(annotation))
    if origin is Literal:
        return " / ".join(json.dumps(arg) for arg in get_args(annotation))
    if origin is tuple:
        return "[" + ", ".join(_type_name(arg) for arg in get_args(annotation)) + "]"
    if origin is list:
        return f"list of {_type_name(get_args(annotation)[0])}"
    if annotation is type(None):
        return "null"
    return getattr(annotation, "__name__", str(annotation))
