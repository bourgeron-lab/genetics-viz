"""Upstream ghfc-ngs documentation, fetched at runtime.

The pipeline's own reference for the parameters file lives in its repository, at
``documentation/params.md`` of bourgeron-lab/ghfc-ngs. The cohort page shows it
beside a cohort's actual ``params.yml`` so a key and its meaning can be read
together.

It is fetched rather than vendored: the repository is public, so no credentials
are involved, and the document is under active development — a copy committed
here would go stale silently. The cost is a network dependency, so a failure is
a degraded state that the panel reports, never an exception reaching the page.

Two things about rendering it are worth knowing before changing anything:

- The document is read **outside its own repository**, so its repo-relative
  links have to be rewritten to absolute GitHub URLs or they 404.
- Its 73 in-page anchors work only because ``markdown2``'s ``header-ids`` extra
  happens to generate the same slugs GitHub does. That was verified against
  every anchor the document uses; do not drop it from :data:`MARKDOWN_EXTRAS`.
"""

from __future__ import annotations

import logging
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Optional

logger = logging.getLogger(__name__)

#: The document itself, raw.
PARAMS_DOC_RAW = (
    "https://raw.githubusercontent.com/bourgeron-lab/ghfc-ngs"
    "/main/documentation/params.md"
)
#: Where to send someone who wants the rendered original.
PARAMS_DOC_PAGE = (
    "https://github.com/bourgeron-lab/ghfc-ngs/blob/main/documentation/params.md"
)

_BLOB = "https://github.com/bourgeron-lab/ghfc-ngs/blob/main/"
_TREE = "https://github.com/bourgeron-lab/ghfc-ngs/tree/main/"

#: markdown2 extras for the rendered document. ``tables`` and
#: ``fenced-code-blocks`` are NiceGUI's defaults, restated because this document
#: is mostly tables (18 of them) and would be unreadable without them.
#: ``header-ids`` is what gives the in-page anchors something to land on.
MARKDOWN_EXTRAS = ["fenced-code-blocks", "tables", "header-ids"]

#: A document under active development, so not cached for long.
_TTL = timedelta(hours=1)
#: Failures are cached too, briefly, so a dead network is not re-dialled on
#: every tab open.
_FAILURE_TTL = timedelta(minutes=1)

_TIMEOUT_SECONDS = 10
_USER_AGENT = "genetics-viz"

#: Markdown links: ``[text](target)``.
_LINK_RE = re.compile(r"\[([^\]]*)\]\(([^)\s]+)\)")
#: Leading ``../`` hops. params.md sits in documentation/, so ``../x`` is the
#: repository root's ``x``.
_UPWARD_RE = re.compile(r"^(?:\.\./)+")
#: Inline code in link text.
_BACKTICK_RE = re.compile(r"`([^`]+)`")


@dataclass(frozen=True)
class FetchedDoc:
    """The document, or why it is not here."""

    #: Link-rewritten markdown, ready for ``ui.markdown``. None on failure.
    markdown: Optional[str]
    #: Human-readable reason, shown in place of the content. None on success.
    error: Optional[str]
    fetched_at: datetime

    @property
    def ok(self) -> bool:
        return self.markdown is not None


def _absolutise_target(target: str) -> Optional[str]:
    """Return the absolute GitHub URL for a repo-relative link target.

    None when the target needs no rewriting (an in-page anchor, an already
    absolute URL, a mail link).
    """
    if target.startswith(("#", "http://", "https://", "mailto:")):
        return None

    path = _UPWARD_RE.sub("", target)
    if not path:
        return None

    # GitHub serves a directory under /tree/ and a file under /blob/; using the
    # wrong one 404s. The only directory link in the document is
    # ``../params_example/``, recognisable by its trailing slash — taken before
    # any fragment, so ``README.md#anchor`` is still a file.
    base = _TREE if path.split("#", 1)[0].endswith("/") else _BLOB
    return f"{base}{path}"


def absolutise_links(markdown: str) -> str:
    """Rewrite the document's repo-relative links to absolute GitHub URLs.

    Rewritten links become raw HTML anchors opening in a new tab, which keeps
    the whole job in the markdown and needs no JavaScript. In-page anchors are
    left alone: they are resolved inside the panel.
    """

    def replace(match: re.Match[str]) -> str:
        text, target = match.group(1), match.group(2)
        url = _absolutise_target(target)
        if url is None:
            return match.group(0)
        # markdown2 does not process markdown inside an inline HTML span, so
        # backticked link text would otherwise render with literal backticks.
        label = _BACKTICK_RE.sub(r"<code>\1</code>", text)
        return f'<a href="{url}" target="_blank" rel="noopener">{label}</a>'

    return _LINK_RE.sub(replace, markdown)


#: The single cached result. One document, so one slot.
_cached: Optional[FetchedDoc] = None


def _is_fresh(doc: FetchedDoc) -> bool:
    ttl = _TTL if doc.ok else _FAILURE_TTL
    return datetime.now() - doc.fetched_at < ttl


def fetch_params_doc(*, force: bool = False) -> FetchedDoc:
    """Fetch the parameters reference, cached.

    Blocking: this is a network read, so call it through
    ``asyncio.to_thread`` — NiceGUI's event loop is single-threaded and a stalled
    fetch would freeze every session.
    """
    global _cached

    if not force and _cached is not None and _is_fresh(_cached):
        return _cached

    try:
        request = urllib.request.Request(
            PARAMS_DOC_RAW, headers={"User-Agent": _USER_AGENT}
        )
        with urllib.request.urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
            raw = response.read().decode("utf-8")
        doc = FetchedDoc(
            markdown=absolutise_links(raw), error=None, fetched_at=datetime.now()
        )
    except (urllib.error.URLError, OSError, UnicodeDecodeError, ValueError) as e:
        # Every failure mode is the same as far as the page is concerned: no
        # document, and a sentence explaining why.
        logger.warning("Could not fetch %s: %s", PARAMS_DOC_RAW, e)
        doc = FetchedDoc(
            markdown=None,
            error=f"Could not reach github.com: {e}",
            fetched_at=datetime.now(),
        )

    _cached = doc
    return doc


def clear_doc_cache() -> None:
    """Drop the cached document."""
    global _cached
    _cached = None
