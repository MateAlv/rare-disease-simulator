"""Best-effort Orphadata API client.

The Orphanet/Orphadata REST API exposes per-ORPHAcode resources (cross
references, phenotype associations, natural history, etc.). Exact endpoint
paths have changed across Orphadata revisions, so this client is intentionally
configuration-driven: it tries a set of named resource templates for an ORPHA
code, stores whatever JSON comes back, and records failures instead of raising.

This lets the first live run reveal which endpoints actually resolve; correct
any 404s by editing ``DEFAULT_RESOURCES`` once verified.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field

from rare_disease_simulator.data_sources.http_client import HttpClient, HttpError

logger = logging.getLogger(__name__)

ORPHADATA_BASE = "https://api.orphadata.com"

# resource name -> path template (``{code}`` is the numeric ORPHA code).
# Best-effort defaults; verify against the live API on first run.
DEFAULT_RESOURCES: dict[str, str] = {
    "cross_references": "/rd-cross-referencing/orphacodes/{code}",
    "phenotypes": "/rd-phenotypes/orphacodes/{code}",
    "natural_history": "/rd-natural-history/orphacodes/{code}",
    "epidemiology": "/rd-epidemiology/orphacodes/{code}",
    "classification": "/rd-classification/orphacodes/{code}",
}


@dataclass
class OrphadataResult:
    """Raw Orphadata resources fetched for one ORPHA code."""

    orpha_code: str
    resources: dict[str, object] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)


def normalize_orpha_code(orpha_id: str) -> str:
    """Return the bare numeric code from ``ORPHA:646`` / ``ORPHA646`` / ``646``."""

    return orpha_id.upper().replace("ORPHA:", "").replace("ORPHA", "").strip()


def fetch_orphadata(
    client: HttpClient,
    orpha_id: str,
    *,
    base_url: str = ORPHADATA_BASE,
    resources: dict[str, str] | None = None,
) -> OrphadataResult:
    """Fetch all configured Orphadata resources for an ORPHA code, tolerantly."""

    code = normalize_orpha_code(orpha_id)
    templates = resources or DEFAULT_RESOURCES
    result = OrphadataResult(orpha_code=code)

    for name, template in templates.items():
        url = f"{base_url}{template.format(code=code)}"
        try:
            body = client.get(url, {"Accept": "application/json"})
            result.resources[name] = json.loads(body)
        except HttpError as exc:
            logger.warning("orphadata %s failed: %s", name, exc)
            result.errors[name] = str(exc)
        except json.JSONDecodeError as exc:
            logger.warning("orphadata %s returned non-JSON: %s", name, exc)
            result.errors[name] = f"non-JSON response: {exc}"

    return result
