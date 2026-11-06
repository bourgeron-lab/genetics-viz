"""Data store utilities — multi-directory registry with per-user selection."""

import os
from pathlib import Path
from typing import Iterable, Optional

from genetics_viz.config_model import DataDirectoryConfig
from genetics_viz.models import DataStore
from genetics_viz.utils.sharding import clear_sharding_cache

# ---------------------------------------------------------------------------
# Multi-store registry
# ---------------------------------------------------------------------------

_data_stores: dict[str, DataStore] = {}  # path_str -> DataStore
_static_prefix_map: dict[str, str] = {}  # path_str -> "/data-0", …
_default_data_dir: str = ""
_dir_descriptions: dict[str, str] = {}  # path_str -> description
#: Next ``/data-N`` index. Never reused: a static route cannot be unmounted, so
#: a prefix freed by removing a directory stays bound to that directory (it
#: answers 404 once the directory is out of the registry).
_next_prefix_index: int = 0


def init_all_data_stores(configs: list[DataDirectoryConfig]) -> None:
    """Initialise one DataStore per configured directory and assign static prefixes."""
    global _data_stores, _static_prefix_map, _default_data_dir, _dir_descriptions
    global _next_prefix_index

    _data_stores = {}
    _static_prefix_map = {}
    _dir_descriptions = {}
    _default_data_dir = ""

    for i, dc in enumerate(configs):
        path_str = str(dc.path)
        store = DataStore(data_dir=Path(dc.path))
        try:
            store.load()
            print(f"  Loaded data directory: {dc.path} ({len(store.cohorts)} cohorts)")
        except FileNotFoundError as e:
            print(f"  Warning: {dc.path}: {e}")
        _data_stores[path_str] = store
        _static_prefix_map[path_str] = f"/data-{i}"
        _dir_descriptions[path_str] = dc.description
        if dc.default and not _default_data_dir:
            _default_data_dir = path_str

    _next_prefix_index = len(configs)

    # Fallback default: first directory
    if not _default_data_dir and _data_stores:
        _default_data_dir = next(iter(_data_stores))


# ---------------------------------------------------------------------------
# Per-user accessors (read app.storage.user in request context)
# ---------------------------------------------------------------------------


#: What ``app.storage.user`` raises outside a request (the poller, a thread
#: started without the context, the unit tests).
_NO_CONTEXT = (RuntimeError, AssertionError, KeyError)


def _norm(path: str) -> str:
    return os.path.normpath(path)


def dirs_for(granted: Optional[Iterable[str]]) -> list[str]:
    """The registered directories among ``granted`` (``None``: all), in order."""
    if granted is None:
        return list(_data_stores)
    wanted = {_norm(g) for g in granted}
    return [p for p in _data_stores if _norm(p) in wanted]


def user_data_dirs() -> list[str]:
    """The registered data directories the current user may open, in order.

    Every directory without a user context (the poller, the gene index), or
    for a user granted none in particular.
    """
    from genetics_viz.utils.auth import granted_data_dirs

    try:
        granted = granted_data_dirs()
    except _NO_CONTEXT:
        granted = None
    return dirs_for(granted)


def current_data_dir() -> Optional[str]:
    """The user's selected directory if they may open it, else their first one.

    The default directory comes first when the user may open it. ``None`` when
    the user may open none.
    """
    allowed = user_data_dirs()
    try:
        from nicegui import app

        chosen = app.storage.user.get("data_dir")
    except _NO_CONTEXT:
        chosen = None
    if chosen in allowed:
        return chosen
    if _default_data_dir in allowed:
        return _default_data_dir
    return allowed[0] if allowed else None


def get_data_store() -> DataStore:
    """Return the DataStore for the current user's selected data directory.

    Falls back to the user's first directory (the default one when allowed)
    if the stored path is no longer valid or no longer theirs.
    """
    current = current_data_dir()
    if current is not None:
        return _data_stores[current]
    if not _data_stores:
        raise RuntimeError("No data stores initialised")
    raise RuntimeError("No data directory is open to this account")


def get_data_store_or_none() -> Optional[DataStore]:
    """Return the current DataStore or None."""
    try:
        return get_data_store()
    except RuntimeError:
        return None


def get_static_prefix() -> str:
    """Return the static-file URL prefix for the current user's data directory."""
    current = current_data_dir()
    if current is not None:
        return _static_prefix_map[current]
    return "/data"


def dir_label(path_str: str) -> str:
    """A directory's name, with its description: ``wgs — Whole genomes``."""
    name = Path(path_str).name
    desc = _dir_descriptions.get(path_str, "")
    return f"{name} — {desc}" if desc else name


def get_data_dir_options() -> list[dict]:
    """Return option dicts for the data-directory dropdown: the user's own.

    Each dict has ``value`` (path string) and ``label`` (directory basename).
    """
    return [{"value": p, "label": dir_label(p)} for p in user_data_dirs()]


def all_data_dir_options() -> list[dict]:
    """Every registered directory as an option (the admin pages)."""
    return [{"value": p, "label": dir_label(p)} for p in _data_stores]


def get_default_data_dir_path() -> str:
    """Return the default data directory path string."""
    return _default_data_dir


def get_static_prefix_map() -> dict[str, str]:
    """Return the full {path_str: prefix} mapping (for app.py startup)."""
    return _static_prefix_map


def mount_data_files(prefix: str, path_str: str) -> None:
    """Serve a data directory's files under ``prefix``, to its users only.

    Not ``app.add_static_files``: its route checks nothing, and the directory
    holds patient genomes. Starlette's ``StaticFiles`` still does the serving,
    so Range requests (igv.js reading CRAMs and tabix files) work as before.
    Anonymous requests get 401, not a redirect to the login page, since they
    come from igv.js, ``<img>`` tags and links rather than page navigation. A
    user not granted the directory, or a reader or student asking for anything
    but a file the pages load themselves, gets 403.
    """
    from fastapi import HTTPException, Request
    from fastapi.responses import Response
    from nicegui import app as nicegui_app
    from starlette.staticfiles import StaticFiles

    from genetics_viz.utils.access import is_traversal, may_fetch
    from genetics_viz.utils.auth import get_current_role, is_authenticated

    # check_dir=False: a missing directory is a degraded state, as in the registry.
    files = StaticFiles(directory=path_str, check_dir=False)

    @nicegui_app.get(prefix + "/{path:path}")
    async def data_file(request: Request, path: str = "") -> Response:
        if not is_authenticated():
            raise HTTPException(status_code=401)
        if _static_prefix_map.get(path_str) != prefix or is_traversal(path):
            raise HTTPException(status_code=404)
        # Only the directories granted to the user, and by role only the
        # files a page of theirs loads (utils/access.py).
        if path_str not in user_data_dirs() or not may_fetch(get_current_role(), path):
            raise HTTPException(status_code=403)
        response = await files.get_response(path, request.scope)
        # private: the F5 and any other shared cache must not keep patient data.
        response.headers["Cache-Control"] = "private, max-age=3600"
        return response


def get_all_data_stores() -> dict[str, DataStore]:
    """Return the full {path_str: DataStore} mapping (read-only access for the poller)."""
    return _data_stores


# ---------------------------------------------------------------------------
# Hot-add / remove (for admin page)
# ---------------------------------------------------------------------------


def add_data_store(path: str, description: str = "") -> str:
    """Add and load a new data directory at runtime. Returns its static prefix."""
    global _next_prefix_index
    path_str = str(path)
    if path_str in _data_stores:
        return _static_prefix_map[path_str]

    prefix = f"/data-{_next_prefix_index}"
    _next_prefix_index += 1

    store = DataStore(data_dir=Path(path))
    try:
        store.load()
    except FileNotFoundError:
        pass

    _data_stores[path_str] = store
    _static_prefix_map[path_str] = prefix
    _dir_descriptions[path_str] = description
    clear_sharding_cache()

    mount_data_files(prefix, path_str)

    return prefix


def remove_data_store(path: str) -> None:
    """Remove a data directory from the registry."""
    path_str = str(path)
    _data_stores.pop(path_str, None)
    _static_prefix_map.pop(path_str, None)
    _dir_descriptions.pop(path_str, None)
    clear_sharding_cache()
