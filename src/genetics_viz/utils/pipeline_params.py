"""The ghfc-ngs workflow parameters file of a cohort.

Every cohort directory may hold ``<NAME>.params.yml``, the parameter file the
workflow was launched with. Its ``steps`` array is the list of pipeline steps
requested for the cohort, and it is the only source of that list until the
pipeline writes a run record (see :mod:`genetics_viz.utils.cohort_state`).

Functions here take paths rather than a :class:`~genetics_viz.models.Cohort`,
so the module stays free of a models import -- the same reason given in
``utils/ancestry.py`` for its loosely typed cohort argument.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

logger = logging.getLogger(__name__)

#: Preferred first. The workflow writes ``.yml``; tolerating the other spelling
#: costs one extra stat and saves a cohort from silently losing both new tabs.
_PARAMS_SUFFIXES = (".params.yml", ".params.yaml")

#: ``(path, mtime, size)`` -> parsed mapping. Keyed on the stat rather than the
#: path alone so an edited file is re-read without any invalidation plumbing.
_params_cache: Dict[Tuple[str, float, int], Dict[str, Any]] = {}


def find_params_file(cohort_dir: Path, cohort_name: str) -> Optional[Path]:
    """Return the cohort's workflow parameters file, or None when absent."""
    for suffix in _PARAMS_SUFFIXES:
        candidate = cohort_dir / f"{cohort_name}{suffix}"
        if candidate.is_file():
            return candidate
    return None


def read_params_text(path: Path) -> str:
    """Return the file verbatim, comments included.

    The Parameters tab shows this rather than a re-serialised mapping: the real
    files carry substantial explanation in comments (the ``ancestry_panel_name``
    invalidation warning, the depth/GQ threshold measurements) that parsing
    throws away.
    """
    try:
        return path.read_text()
    except OSError as e:
        logger.warning("Could not read params file %s: %s", path, e)
        return ""


def load_params(path: Path) -> Dict[str, Any]:
    """Parse the parameters file, cached on its stat. Returns {} on failure."""
    try:
        stat = path.stat()
    except OSError as e:
        logger.warning("Could not stat params file %s: %s", path, e)
        return {}

    key = (str(path), stat.st_mtime, stat.st_size)
    cached = _params_cache.get(key)
    if cached is not None:
        return cached

    try:
        with open(path, "r") as f:
            data = yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError) as e:
        # A malformed params file is a degraded state, not an error: the rest of
        # the cohort page still works without it.
        logger.warning("Could not parse params file %s: %s", path, e)
        data = {}

    if not isinstance(data, dict):
        logger.warning("Params file %s is not a mapping", path)
        data = {}

    _params_cache[key] = data
    return data


def get_requested_steps(params: Dict[str, Any]) -> List[str]:
    """Return the ``steps`` array as a list of strings.

    Coerced defensively: the key may be missing, a scalar, or carry non-string
    entries, and none of those should break the page that renders it.
    """
    raw = params.get("steps")
    if isinstance(raw, str):
        return [raw]
    if not isinstance(raw, (list, tuple)):
        return []
    return [str(step) for step in raw if step is not None]


def clear_params_cache() -> None:
    """Drop every cached parameters mapping."""
    _params_cache.clear()
