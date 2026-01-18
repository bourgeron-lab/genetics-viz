"""Data store utilities."""

from pathlib import Path
from typing import Optional

from genetics_viz.models import DataStore

# Global data store - will be initialized when the app starts
_data_store: Optional[DataStore] = None


def get_data_store() -> DataStore:
    """Get the global data store instance."""
    if _data_store is None:
        raise RuntimeError("Data store not initialized")
    return _data_store


def get_data_store_or_none() -> Optional[DataStore]:
    """Get the global data store instance or None if not initialized."""
    return _data_store


def init_data_store(data_dir: Path) -> DataStore:
    """Initialize the global data store."""
    global _data_store
    _data_store = DataStore(data_dir=data_dir)
    return _data_store
