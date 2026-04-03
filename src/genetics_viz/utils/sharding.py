"""Two-level sharding utilities for sample and family directories.

WGS data directories may use a sharding scheme to avoid thousands of entries
in a single directory. Given an entity ID:

1. Strip separators (-, ., _)
2. shard1 = last character (uppercased)
3. shard2 = second-to-last character (uppercased)

Path: samples/<shard1>/<shard2>/<original_id>/
"""

import os
from pathlib import Path

_sharding_cache: dict[tuple[str, str], bool] = {}


def compute_shard_prefix(entity_id: str) -> tuple[str, str]:
    """Compute shard keys from an entity ID.

    Returns (shard1, shard2) where shard1 is the last character and shard2 is
    the second-to-last character of the stripped (no separators) ID, uppercased.
    """
    stripped = entity_id.replace("-", "").replace(".", "").replace("_", "")
    return stripped[-1].upper(), stripped[-2].upper()


def _is_sharded(base_dir: Path, entity_type: str) -> bool:
    """Check if a directory uses sharded layout (cached).

    Returns True if ALL immediate children of base_dir/entity_type are
    single-character directories.
    """
    cache_key = (str(base_dir), entity_type)
    if cache_key in _sharding_cache:
        return _sharding_cache[cache_key]

    entity_dir = base_dir / entity_type
    if not entity_dir.is_dir():
        _sharding_cache[cache_key] = False
        return False

    has_children = False
    for entry in os.scandir(entity_dir):
        if not entry.is_dir():
            continue
        has_children = True
        if len(entry.name) != 1:
            _sharding_cache[cache_key] = False
            return False

    _sharding_cache[cache_key] = has_children
    return has_children


def get_entity_path(data_dir: Path, entity_type: str, entity_id: str) -> Path:
    """Return the full filesystem path for an entity (sample or family).

    Uses a try-sharded-first strategy: checks the sharded path first,
    falls back to flat. This handles hybrid directories where both
    shard buckets and direct entity folders coexist.
    """
    shard1, shard2 = compute_shard_prefix(entity_id)
    sharded_path = data_dir / entity_type / shard1 / shard2 / entity_id
    if sharded_path.is_dir():
        return sharded_path
    return data_dir / entity_type / entity_id


def get_entity_url_segment(data_dir: Path, entity_type: str, entity_id: str) -> str:
    """Return the relative URL path segment for an entity.

    E.g. 'samples/J/Z/C000EZJ' (sharded) or 'samples/C000EZJ' (flat).
    Mirrors get_entity_path logic: tries sharded first, falls back to flat.
    """
    shard1, shard2 = compute_shard_prefix(entity_id)
    sharded_path = data_dir / entity_type / shard1 / shard2 / entity_id
    if sharded_path.is_dir():
        return f"{entity_type}/{shard1}/{shard2}/{entity_id}"
    return f"{entity_type}/{entity_id}"


def get_sample_path(data_dir: Path, sample_id: str) -> Path:
    """Return the filesystem path for a sample directory."""
    return get_entity_path(data_dir, "samples", sample_id)


def get_family_path(data_dir: Path, family_id: str) -> Path:
    """Return the filesystem path for a family directory."""
    return get_entity_path(data_dir, "families", family_id)


def get_sample_url(data_dir: Path, sample_id: str) -> str:
    """Return the relative URL segment for a sample directory."""
    return get_entity_url_segment(data_dir, "samples", sample_id)


def get_family_url(data_dir: Path, family_id: str) -> str:
    """Return the relative URL segment for a family directory."""
    return get_entity_url_segment(data_dir, "families", family_id)


def clear_sharding_cache() -> None:
    """Clear the sharding detection cache."""
    _sharding_cache.clear()
