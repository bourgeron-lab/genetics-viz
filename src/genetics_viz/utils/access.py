"""What each role may fetch from a data directory's ``/data-N`` route.

Curators and administrators get every file: they curate from the raw tables
and the somalier report. Readers and students get only what the pages
themselves load in the browser: igv.js's alignment and coverage tracks, and
the REViewer images. Everything a reader may take away goes through the
role-aware Download buttons, which build their files in Python; students do
not get those either.

The rule is on paths, not roles' good faith: without it, a reader could fetch
``validations/``, ``diagnostics/`` or the cohort tables (with relationships
the reader wording hides) by typing a URL.
"""

from __future__ import annotations

import re

#: Roles that may fetch any file of a data directory they can open.
FULL_ACCESS_ROLES = frozenset({"curator", "administrator"})

#: Roles that may use the Download buttons (TSV exports, Search's TSV).
DOWNLOAD_ROLES = frozenset({"reader", "curator", "administrator"})

#: Files the pages load in the browser, relative to the data directory. The
#: sample directory may be sharded (``samples/<s1>/<s2>/<id>/``).
_VIEWER_FILES = re.compile(
    r"^samples/(?:[^/]+/)+(?:"
    # igv.js: alignments, by1000 / VAF bedgraphs and their tabix indexes
    r"sequences/[^/]+\.(?:cram|crai|bam|bai|bedgraph|bedgraph\.gz|bedgraph\.gz\.tbi)"
    # igv.js: the normalised coverage bigWig
    r"|coverage/[^/]+\.(?:bw|bigwig)"
    # REViewer images of flagged STR genotypes
    r"|str/expansionhunter/reviewer/[^/]+\.(?:svg|png)"
    r")$"
)


def is_traversal(path: str) -> bool:
    """Whether ``path`` has a ``..`` or ``.`` segment (or is absolute)."""
    return path.startswith("/") or any(p in ("..", ".") for p in path.split("/"))


def is_viewer_file(path: str) -> bool:
    """Whether a page loads ``path`` itself (an igv.js track, a REViewer image)."""
    return not is_traversal(path) and bool(_VIEWER_FILES.match(path))


def may_fetch(role: str, path: str) -> bool:
    """Whether ``role`` may fetch ``path`` (relative to the data directory).

    An unknown role is treated as the most restricted one.
    """
    if is_traversal(path):
        return False
    return role in FULL_ACCESS_ROLES or is_viewer_file(path)


def may_download(role: str) -> bool:
    """Whether ``role`` gets the Download buttons."""
    return role in DOWNLOAD_ROLES
