"""The curation files under ``<data_dir>/validations/``, one per variant kind.

The validation pages read every kind; listing them here keeps a new kind a
one-line change instead of a hunt through hand-copied ``snvs.tsv``/``svs.tsv``
branches.
"""

from pathlib import Path
from typing import Dict, List, Tuple

#: (type label shown in the pages, file name), in display order.
VALIDATION_KINDS: List[Tuple[str, str]] = [
    ("SNV", "snvs.tsv"),
    ("SV", "svs.tsv"),
    ("STR", "strs.tsv"),
]

VALIDATION_TYPES: List[str] = [label for label, _ in VALIDATION_KINDS]

#: Chart colour per type.
TYPE_COLORS: Dict[str, str] = {"SNV": "#3b82f6", "SV": "#8b5cf6", "STR": "#0d9488"}


def validation_path(data_dir: Path, vtype: str) -> Path:
    """``<data_dir>/validations/<file>`` for a type label (``SNV``, ``SV``, ``STR``)."""
    name = dict(VALIDATION_KINDS)[vtype]
    return data_dir / "validations" / name


def validation_paths(data_dir: Path) -> List[Tuple[str, Path]]:
    """Every kind's file, existing or not, in display order."""
    return [(vtype, validation_path(data_dir, vtype)) for vtype in VALIDATION_TYPES]
