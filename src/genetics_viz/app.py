"""
NiceGUI web application for genetics-viz.

This is the main entry point for the application. All pages are now modular:
- genetics_viz.pages.cohort: Home, Cohort, Family, and Variant pages
- genetics_viz.pages.validation: Validation Statistics, File, and All pages
"""

import os
from pathlib import Path
from typing import Optional

from nicegui import ui

from genetics_viz.models import DataStore

# Import pages to register routes - this triggers all @ui.page decorators
from genetics_viz.pages import cohort, validation  # noqa: F401
from genetics_viz.utils import data as data_module

# Legacy: Keep local reference for backward compatibility
_data_store: Optional[DataStore] = None


def run_app(
    data_dir: Path,
    host: str = "127.0.0.1",
    port: int = 8080,
    reload: bool = False,
) -> None:
    """
    Initialize and run the NiceGUI application.

    Args:
        data_dir: Path to the data directory containing cohort data
        host: Host address to bind the server to
        port: Port to run the server on
        reload: Enable auto-reload for development
    """
    global _data_store

    # Store config in environment for reload mode
    if reload:
        os.environ["GENETICS_VIZ_DATA_DIR"] = str(data_dir)
        os.environ["GENETICS_VIZ_HOST"] = host
        os.environ["GENETICS_VIZ_PORT"] = str(port)
        os.environ["GENETICS_VIZ_RELOAD"] = "1"

    # Initialize the data store - both local and in shared module
    _data_store = DataStore(data_dir=data_dir)
    data_module._data_store = _data_store

    try:
        _data_store.load()
        print(f"Loaded {len(_data_store.cohorts)} cohorts")
    except FileNotFoundError as e:
        print(f"Warning: {e}")
        print("The application will start but no cohorts will be available.")

    # Configure NiceGUI
    ui.run(
        host=host,
        port=port,
        reload=reload,
        title="Genetics-Viz",
        favicon="🧬",
        dark=False,
    )


# Auto-initialize and run when module is reloaded (for reload mode)
if os.environ.get("GENETICS_VIZ_RELOAD") == "1" and _data_store is None:
    data_dir_env = os.environ.get("GENETICS_VIZ_DATA_DIR")
    if data_dir_env:
        _data_store = DataStore(data_dir=Path(data_dir_env))
        data_module._data_store = _data_store
        try:
            _data_store.load()
            print(f"[Reload] Loaded {len(_data_store.cohorts)} cohorts")
        except FileNotFoundError as e:
            print(f"[Reload] Warning: {e}")

        # Call ui.run() at module level for reload mode
        host = os.environ.get("GENETICS_VIZ_HOST", "127.0.0.1")
        port = int(os.environ.get("GENETICS_VIZ_PORT", "8080"))
        ui.run(
            host=host,
            port=port,
            reload=True,
            title="Genetics-Viz",
            favicon="🧬",
            dark=False,
        )
