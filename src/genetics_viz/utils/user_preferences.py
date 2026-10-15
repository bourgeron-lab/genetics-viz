"""Per-user preferences, stored server-side next to the application config.

Today this holds one thing: how a user wants gene badges coloured from
GeneTrek genesets (:class:`GenePrefs`). The file is
``<config stem>.preferences.yaml`` beside the config file (see
``config_model.get_preferences_path``)::

    default:            # set by an administrator; absent -> packaged default
      tier_colors: ["#c62828", "#ef6c00", "#f9a825", "#1565c0"]
      sets: {hc_ndd: 1, constraint: 3, autism_combined_v7: show}
    users:
      alice:
        tier_colors: [...]
        sets: {...}

It is deliberately not the config file. ``save_config`` rewrites that file
whole from the in-memory singleton, truncating before it locks; profile saves
mixed into it would race with admin saves. Here every write takes a lock,
re-reads the file, changes one entry and replaces the file atomically.

Set keys are GeneTrek names or aliases (see ``GenesetInfo.key``). A key
GeneTrek no longer knows is kept, so a list that comes back is not lost, and
ignored when scoring.
"""

from __future__ import annotations

import fcntl
import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import yaml

from genetics_viz.config_model import get_preferences_path
from genetics_viz.utils.auth import get_current_user, is_admin
from genetics_viz.utils.genetrek import GENETREK_CONFIG

logger = logging.getLogger(__name__)

#: Number of colour tiers. The profile page and the scorer both assume four.
N_TIERS = 4
SHOW = "show"

#: A set's assignment: a tier number (1 = best) or "show" (tooltip only).
Assignment = Union[int, str]


def _clean_assignment(value: Any) -> Optional[Assignment]:
    if value == SHOW:
        return SHOW
    try:
        tier = int(value)
    except (TypeError, ValueError):
        return None
    return tier if 1 <= tier <= N_TIERS else None


def _clean_colors(colors: Any) -> Tuple[str, ...]:
    packaged = tuple(GENETREK_CONFIG.get("tier_colors") or ())
    given = [str(c) for c in (colors or []) if isinstance(c, str) and c.startswith("#")]
    return tuple(given[:N_TIERS]) + packaged[len(given[:N_TIERS]) : N_TIERS]


@dataclass(frozen=True)
class GenePrefs:
    """How one user colours and lists genesets."""

    tier_colors: Tuple[str, ...] = field(default_factory=lambda: _clean_colors(None))
    #: Set key -> tier (1..N_TIERS) or "show". Unlisted sets are hidden.
    sets: Dict[str, Assignment] = field(default_factory=dict)
    #: Where these came from: "user", "admin default" or "packaged default".
    source: str = "user"

    @classmethod
    def from_dict(cls, raw: Any, source: str = "user") -> "GenePrefs":
        raw = raw if isinstance(raw, dict) else {}
        sets: Dict[str, Assignment] = {}
        for key, value in (raw.get("sets") or {}).items():
            cleaned = _clean_assignment(value)
            if cleaned is not None:
                sets[str(key)] = cleaned
        return cls(_clean_colors(raw.get("tier_colors")), sets, source)

    def to_dict(self) -> Dict[str, Any]:
        return {"tier_colors": list(self.tier_colors), "sets": dict(self.sets)}

    @property
    def revision(self) -> Tuple:
        """Hashable identity, for caches keyed on the preferences."""
        return (self.tier_colors, tuple(sorted(self.sets.items(), key=str)))


def packaged_default() -> GenePrefs:
    return GenePrefs.from_dict(
        {
            "tier_colors": GENETREK_CONFIG.get("tier_colors"),
            "sets": GENETREK_CONFIG.get("default_sets"),
        },
        source="packaged default",
    )


# ---------------------------------------------------------------------------
# File access
# ---------------------------------------------------------------------------

#: (path, mtime_ns, size) -> parsed file. One slot: there is one file.
_cache: Optional[Tuple[Tuple[str, int, int], Dict[str, Any]]] = None


def clear_preferences_cache() -> None:
    global _cache
    _cache = None


def _read_file(path: Path) -> Dict[str, Any]:
    global _cache
    try:
        st = path.stat()
    except FileNotFoundError:
        return {}
    key = (str(path), st.st_mtime_ns, st.st_size)
    if _cache is not None and _cache[0] == key:
        return _cache[1]
    try:
        raw = yaml.safe_load(path.read_text()) or {}
    except (OSError, yaml.YAMLError) as e:
        logger.warning("Ignoring unreadable preferences file %s: %s", path, e)
        raw = {}
    if not isinstance(raw, dict):
        raw = {}
    _cache = (key, raw)
    return raw


def _update_file(change) -> None:
    """Apply ``change(data)`` to the file under an exclusive lock, atomically."""
    path = get_preferences_path()
    if path is None:
        raise RuntimeError("No application config loaded")
    lock_path = path.with_name(path.name + ".lock")
    with open(lock_path, "a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            clear_preferences_cache()
            data = dict(_read_file(path))
            change(data)
            fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
            try:
                with os.fdopen(fd, "w") as fh:
                    yaml.safe_dump(data, fh, default_flow_style=False, sort_keys=False)
                os.replace(tmp, path)
            except BaseException:
                Path(tmp).unlink(missing_ok=True)
                raise
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            clear_preferences_cache()


def _read() -> Dict[str, Any]:
    path = get_preferences_path()
    return _read_file(path) if path is not None else {}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_default_gene_prefs() -> GenePrefs:
    """The admin-set default, else the packaged one."""
    raw = _read().get("default")
    if isinstance(raw, dict) and raw.get("sets"):
        return GenePrefs.from_dict(raw, source="admin default")
    return packaged_default()


def has_own_gene_prefs(username: str) -> bool:
    users = _read().get("users") or {}
    return isinstance(users, dict) and isinstance(users.get(username), dict)


def get_gene_prefs(username: Optional[str]) -> GenePrefs:
    """A user's preferences, else the default."""
    users = _read().get("users") or {}
    raw = users.get(username) if isinstance(users, dict) and username else None
    if isinstance(raw, dict):
        return GenePrefs.from_dict(raw, source="user")
    return get_default_gene_prefs()


def save_gene_prefs(username: str, prefs: GenePrefs) -> None:
    """Save the current user's own preferences.

    Raises ``PermissionError`` for anyone else's: the profile page only ever
    edits its own user, and this guard makes that true at the backend too.
    """
    if username != get_current_user():
        raise PermissionError("Users can only save their own preferences")

    def change(data: Dict[str, Any]) -> None:
        users = data.get("users") if isinstance(data.get("users"), dict) else {}
        users[username] = prefs.to_dict()
        data["users"] = users

    _update_file(change)


def reset_gene_prefs(username: str) -> None:
    """Drop a user's own preferences, so they follow the default again."""
    if username != get_current_user() and not is_admin():
        raise PermissionError("Users can only reset their own preferences")
    _drop_user(username)


def _drop_user(username: str) -> None:
    def change(data: Dict[str, Any]) -> None:
        users = data.get("users")
        if isinstance(users, dict):
            users.pop(username, None)

    _update_file(change)


def remove_user_preferences(username: str) -> None:
    """Forget a removed user's preferences. Called by the admin users page."""
    if not is_admin():
        raise PermissionError("Only administrators can remove users")
    if has_own_gene_prefs(username):
        _drop_user(username)


def save_default_gene_prefs(prefs: GenePrefs) -> None:
    """Make ``prefs`` the default for every user without their own."""
    if not is_admin():
        raise PermissionError("Only administrators can set the default")

    def change(data: Dict[str, Any]) -> None:
        data["default"] = prefs.to_dict()

    _update_file(change)
