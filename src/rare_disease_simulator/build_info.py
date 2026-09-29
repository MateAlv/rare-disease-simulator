"""Reproducibility helpers: file checksums and the simulator's git revision."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

_CHUNK_SIZE = 1 << 20


def sha256_file(path: Path | str) -> str:
    """Return the hex SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as file:
        for chunk in iter(lambda: file.read(_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_revision(repo_dir: Path | str | None = None) -> dict[str, object]:
    """Return ``{"sha": ..., "dirty": ...}`` for the repository holding this package."""

    cwd = Path(repo_dir) if repo_dir is not None else Path(__file__).resolve().parent
    try:
        sha = _git(cwd, "rev-parse", "HEAD")
        status = _git(cwd, "status", "--porcelain", "--untracked-files=no")
    except (OSError, subprocess.CalledProcessError):
        return {"sha": None, "dirty": None}
    return {"sha": sha, "dirty": bool(status)}


def _git(cwd: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    )
    return completed.stdout.strip()
