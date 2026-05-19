"""Data store utilities — multi-directory registry with per-user selection."""

from pathlib import Path
from typing import Optional

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


def init_all_data_stores(configs: list[DataDirectoryConfig]) -> None:
    """Initialise one DataStore per configured directory and assign static prefixes."""
    global _data_stores, _static_prefix_map, _default_data_dir, _dir_descriptions

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

    # Fallback default: first directory
    if not _default_data_dir and _data_stores:
        _default_data_dir = next(iter(_data_stores))


# ---------------------------------------------------------------------------
# Per-user accessors (read app.storage.user in request context)
# ---------------------------------------------------------------------------


def get_data_store() -> DataStore:
    """Return the DataStore for the current user's selected data directory.

    Falls back to the default directory if no user context or if the stored
    path is no longer valid.
    """
    try:
        from nicegui import app

        data_dir = app.storage.user.get("data_dir")
        if data_dir and data_dir in _data_stores:
            return _data_stores[data_dir]
    except Exception:
        pass
    # Fallback to default
    if _default_data_dir and _default_data_dir in _data_stores:
        return _data_stores[_default_data_dir]
    if _data_stores:
        return next(iter(_data_stores.values()))
    raise RuntimeError("No data stores initialised")


def get_data_store_or_none() -> Optional[DataStore]:
    """Return the current DataStore or None."""
    try:
        return get_data_store()
    except RuntimeError:
        return None


def get_static_prefix() -> str:
    """Return the static-file URL prefix for the current user's data directory."""
    try:
        from nicegui import app

        data_dir = app.storage.user.get("data_dir")
        if data_dir and data_dir in _static_prefix_map:
            return _static_prefix_map[data_dir]
    except Exception:
        pass
    if _default_data_dir and _default_data_dir in _static_prefix_map:
        return _static_prefix_map[_default_data_dir]
    if _static_prefix_map:
        return next(iter(_static_prefix_map.values()))
    return "/data"


def get_data_dir_options() -> list[dict]:
    """Return option dicts for the data-directory dropdown.

    Each dict has ``value`` (path string) and ``label`` (directory basename).
    """
    options = []
    for path_str in _data_stores:
        name = Path(path_str).name
        desc = _dir_descriptions.get(path_str, "")
        label = f"{name} — {desc}" if desc else name
        options.append({"value": path_str, "label": label})
    return options


def get_default_data_dir_path() -> str:
    """Return the default data directory path string."""
    return _default_data_dir


def get_static_prefix_map() -> dict[str, str]:
    """Return the full {path_str: prefix} mapping (for app.py startup)."""
    return _static_prefix_map


def get_all_data_stores() -> dict[str, DataStore]:
    """Return the full {path_str: DataStore} mapping (read-only access for the poller)."""
    return _data_stores


# ---------------------------------------------------------------------------
# Hot-add / remove (for admin page)
# ---------------------------------------------------------------------------


def add_data_store(path: str, description: str = "") -> str:
    """Add and load a new data directory at runtime. Returns its static prefix."""
    path_str = str(path)
    if path_str in _data_stores:
        return _static_prefix_map[path_str]

    idx = len(_static_prefix_map)
    prefix = f"/data-{idx}"

    store = DataStore(data_dir=Path(path))
    try:
        store.load()
    except FileNotFoundError:
        pass

    _data_stores[path_str] = store
    _static_prefix_map[path_str] = prefix
    _dir_descriptions[path_str] = description
    clear_sharding_cache()

    # Register static files
    from nicegui import app as nicegui_app

    nicegui_app.add_static_files(prefix, path_str)

    return prefix


def remove_data_store(path: str) -> None:
    """Remove a data directory from the registry."""
    path_str = str(path)
    _data_stores.pop(path_str, None)
    _static_prefix_map.pop(path_str, None)
    _dir_descriptions.pop(path_str, None)
    clear_sharding_cache()
