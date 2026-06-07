"""JSONL readers and writers for simulator artifacts."""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path
from typing import Any, TypeVar

from pydantic import BaseModel, ValidationError

ModelT = TypeVar("ModelT", bound=BaseModel)
JsonRecord = BaseModel | Mapping[str, Any]


class JsonlError(ValueError):
    """Base error for invalid JSONL input."""


class JsonlDecodeError(JsonlError):
    """JSONL decoding error that includes file and line context."""

    def __init__(self, path: Path, line_number: int, message: str) -> None:
        super().__init__(f"{path}:{line_number}: {message}")
        self.path = path
        self.line_number = line_number


class JsonlValidationError(JsonlError):
    """Pydantic validation error that includes file and line context."""

    def __init__(self, path: Path, line_number: int, error: ValidationError) -> None:
        super().__init__(f"{path}:{line_number}: {error}")
        self.path = path
        self.line_number = line_number
        self.error = error


def _record_to_json(record: JsonRecord) -> str:
    if isinstance(record, BaseModel):
        return record.model_dump_json(exclude_none=True)
    return json.dumps(record, ensure_ascii=False, sort_keys=True)


def iter_jsonl(path: Path | str) -> Iterator[dict[str, Any]]:
    """Yield raw JSON objects from a JSONL file, skipping blank lines."""

    jsonl_path = Path(path)
    with jsonl_path.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                value = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise JsonlDecodeError(jsonl_path, line_number, exc.msg) from exc
            if not isinstance(value, dict):
                raise JsonlDecodeError(jsonl_path, line_number, "expected a JSON object")
            yield value


def read_jsonl(path: Path | str) -> list[dict[str, Any]]:
    """Read all raw JSON objects from a JSONL file."""

    return list(iter_jsonl(path))


def iter_model_jsonl(path: Path | str, model: type[ModelT]) -> Iterator[ModelT]:
    """Yield Pydantic models from a JSONL file."""

    jsonl_path = Path(path)
    for line_number, record in _iter_jsonl_with_line_numbers(jsonl_path):
        try:
            yield model.model_validate(record)
        except ValidationError as exc:
            raise JsonlValidationError(jsonl_path, line_number, exc) from exc


def read_model_jsonl(path: Path | str, model: type[ModelT]) -> list[ModelT]:
    """Read all records from a JSONL file into Pydantic models."""

    return list(iter_model_jsonl(path, model))


def write_jsonl(path: Path | str, records: Iterable[JsonRecord], *, append: bool = False) -> int:
    """Write records to JSONL and return the number of records written."""

    jsonl_path = Path(path)
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    mode = "a" if append else "w"
    with jsonl_path.open(mode, encoding="utf-8") as file:
        for record in records:
            file.write(_record_to_json(record))
            file.write("\n")
            count += 1
    return count


def append_jsonl(path: Path | str, records: Iterable[JsonRecord]) -> int:
    """Append records to a JSONL file and return the number of records written."""

    return write_jsonl(path, records, append=True)


def _iter_jsonl_with_line_numbers(path: Path) -> Iterator[tuple[int, dict[str, Any]]]:
    with path.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                value = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise JsonlDecodeError(path, line_number, exc.msg) from exc
            if not isinstance(value, dict):
                raise JsonlDecodeError(path, line_number, "expected a JSON object")
            yield line_number, value
