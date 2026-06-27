"""Minimal dependency-free HTTP client for source fetchers.

Uses the standard library only. Adds the etiquette every public biomedical API
expects: a descriptive User-Agent, a global request rate limit, bounded
retries with backoff, and explicit timeouts. All fetchers go through this so
politeness and error handling live in one place.
"""

from __future__ import annotations

import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

DEFAULT_USER_AGENT = (
    "rare-disease-simulator/0.1 (research; mailto:mateoalvarez2012@gmail.com)"
)


class HttpError(RuntimeError):
    """Raised when a request fails after exhausting retries."""

    def __init__(self, url: str, status: int | None, message: str) -> None:
        super().__init__(f"{url} -> {status}: {message}")
        self.url = url
        self.status = status


@dataclass
class HttpClient:
    """Rate-limited, retrying HTTP GET client."""

    min_interval_seconds: float = 0.34  # ~3 req/s, NCBI's no-key limit
    timeout_seconds: float = 30.0
    max_retries: int = 3
    backoff_seconds: float = 1.0
    user_agent: str = DEFAULT_USER_AGENT
    extra_params: dict[str, str] = field(default_factory=dict)
    _last_request_at: float = field(default=0.0, init=False, repr=False)

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        wait = self.min_interval_seconds - elapsed
        if wait > 0:
            time.sleep(wait)
        self._last_request_at = time.monotonic()

    def get(self, url: str, params: dict[str, str] | None = None) -> str:
        """GET a URL with shared params merged in; return the decoded body."""

        merged = {**self.extra_params, **(params or {})}
        full_url = url
        if merged:
            full_url = f"{url}?{urllib.parse.urlencode(merged)}"

        last_error: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            self._throttle()
            request = urllib.request.Request(full_url, headers={"User-Agent": self.user_agent})
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    charset = response.headers.get_content_charset() or "utf-8"
                    return response.read().decode(charset, errors="replace")
            except urllib.error.HTTPError as exc:
                last_error = exc
                if exc.code in {429, 500, 502, 503, 504} and attempt < self.max_retries:
                    delay = self.backoff_seconds * attempt
                    logger.warning("HTTP %s on %s; retrying in %.1fs", exc.code, full_url, delay)
                    time.sleep(delay)
                    continue
                raise HttpError(full_url, exc.code, exc.reason or "http error") from exc
            except urllib.error.URLError as exc:
                last_error = exc
                if attempt < self.max_retries:
                    delay = self.backoff_seconds * attempt
                    logger.warning("network error on %s: %s; retry in %.1fs", full_url, exc, delay)
                    time.sleep(delay)
                    continue
                raise HttpError(full_url, None, str(exc.reason)) from exc

        raise HttpError(full_url, None, str(last_error))
