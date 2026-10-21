"""The diagnostic files under ``<data_dir>/diagnostics/``, one per variant kind.

The pages that read every kind (Diagnostic All and Statistics, the cohort and
home summaries, the family Diagnostics panel) loop over this list, so a new
kind is one entry here, as for :mod:`genetics_viz.utils.validation_files`.
"""

from pathlib import Path
from typing import List, Tuple

from .validation_files import TYPE_COLORS

__all__ = [
    "DIAGNOSTIC_KINDS",
    "DIAGNOSTIC_TYPES",
    "TYPE_COLORS",
    "diagnostic_path",
    "diagnostic_paths",
]

#: (type label shown in the pages, file name), in display order.
DIAGNOSTIC_KINDS: List[Tuple[str, str]] = [
    ("SNV", "snvs.tsv"),
    ("SV", "svs.tsv"),
    ("STR", "strs.tsv"),
    ("UPD", "upds.tsv"),
]

DIAGNOSTIC_TYPES: List[str] = [label for label, _ in DIAGNOSTIC_KINDS]


def diagnostic_path(data_dir: Path, dtype: str) -> Path:
    """``<data_dir>/diagnostics/<file>`` for a type label (``SNV``, ``SV``, ``STR``, ``UPD``)."""
    name = dict(DIAGNOSTIC_KINDS)[dtype]
    return data_dir / "diagnostics" / name


def diagnostic_paths(data_dir: Path) -> List[Tuple[str, Path]]:
    """Every kind's file, existing or not, in display order."""
    return [(dtype, diagnostic_path(data_dir, dtype)) for dtype in DIAGNOSTIC_TYPES]
