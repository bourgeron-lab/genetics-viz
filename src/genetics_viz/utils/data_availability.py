"""Data availability checking for samples and families.

Provides functions to check which data files and directories exist
for a given sample or family, used by the home page search and
standalone family page.
"""

from pathlib import Path
from typing import Any

from genetics_viz.utils.sharding import get_family_path, get_sample_path


def check_sample_availability(data_dir: Path, sample_id: str) -> dict[str, Any]:
    """Check which data files exist for a sample.

    Returns a dict with boolean flags for each data type.
    Blocking I/O — callers should wrap in asyncio.to_thread().
    """
    sample_path = get_sample_path(data_dir, sample_id)
    seq_dir = sample_path / "sequences"

    return {
        "exists": sample_path.is_dir(),
        "cram": (seq_dir / f"{sample_id}.GRCh38_GIABv3.cram").exists(),
        "cram_index": (seq_dir / f"{sample_id}.GRCh38_GIABv3.cram.crai").exists(),
        "bedgraph": (seq_dir / f"{sample_id}.by1000.bedgraph.gz").exists(),
        "bedgraph_index": (seq_dir / f"{sample_id}.by1000.bedgraph.gz.tbi").exists(),
        "vaf_bedgraph": (seq_dir / f"{sample_id}.vaf.bedgraph.gz").exists(),
        "deepvariant": (sample_path / "deepvariant").is_dir(),
        "svs": (sample_path / "svs").is_dir(),
    }


def check_family_availability(
    data_dir: Path, family_id: str, sample_ids: list[str] | None = None
) -> dict[str, Any]:
    """Check which data files/directories exist for a family.

    Returns a dict with family-level and per-sample availability.
    Blocking I/O — callers should wrap in asyncio.to_thread().
    """
    family_path = get_family_path(data_dir, family_id)

    result: dict[str, Any] = {
        "exists": family_path.is_dir(),
        "pedigree": (family_path / f"{family_id}.pedigree.tsv").exists(),
        "vcfs": (family_path / "vcfs").is_dir(),
        "wombat": (family_path / "wombat").is_dir(),
        "wisecondorx": (family_path / "svs" / "wisecondorx").is_dir(),
        "extractor": (family_path / "extractor").is_dir(),
    }

    if sample_ids:
        result["samples"] = {
            sid: check_sample_availability(data_dir, sid) for sid in sample_ids
        }

    return result
