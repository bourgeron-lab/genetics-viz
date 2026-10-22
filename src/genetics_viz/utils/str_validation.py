"""Curator validation of STR genotypes: ``<data_dir>/validations/strs.tsv``.

Same layout and vocabulary as the SNV and SV files, so ``validation_loader``,
the badges and the validation pages read it unchanged:

- ``Variant`` is ``<chrom>:<pos>-<end>:<LOCUS>``. The coordinates make it
  readable on the validation pages; the locus ID tells apart the loci sharing
  a region (ARX_1 / ARX_2). It does not name the catalog, so a validation
  survives a new catalog label for the same locus.
- ``Validation`` is ``present`` (expansion confirmed), ``absent`` (artefact)
  or ``uncertain``.
- ``CuratedRepeats`` holds the repeat count a curator settled on (by
  REViewer, a long-read or a PCR), when they give one.

Append-only under an exclusive ``flock``, like the other curation files.
"""

import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from .auth import can_write
from .tsv import append_tsv_row
from .validation_files import validation_path

HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "CuratedRepeats", "Comment", "Ignore", "Timestamp",
]  # fmt: skip

VALIDATIONS = ["present", "absent", "uncertain"]

INHERITANCES = [
    "unknown", "de novo", "paternal", "maternal", "not paternal",
    "not maternal", "either", "homozygous",
]  # fmt: skip


def str_validation_file(data_dir: Path) -> Path:
    return validation_path(data_dir, "STR")


def str_variant_key(chrom: Any, pos: Any, end: Any, locus: str) -> str:
    """``chr4:3074876-3074933:HTT``."""
    return f"{chrom}:{pos}-{end}:{locus}"


def save_str_validation(
    data_dir: Path,
    *,
    family_id: str,
    variant_key: str,
    sample: str,
    user: str,
    validation: str,
    inheritance: str = "unknown",
    curated_repeats: str = "",
    comment: str = "",
    timestamp: Optional[str] = None,
) -> None:
    """Append one validation. Refused for viewers and for unknown values.

    The role check sits here as well as on the form: the save button is only
    hidden from viewers, which a crafted event would not respect.
    """
    if not can_write():
        raise PermissionError("Only curators and administrators can save validations")
    if validation not in VALIDATIONS:
        raise ValueError(f"Unknown validation: {validation!r}")
    if inheritance not in INHERITANCES:
        raise ValueError(f"Unknown inheritance: {inheritance!r}")
    curated_repeats = (curated_repeats or "").strip()
    if curated_repeats and not curated_repeats.isdigit():
        raise ValueError("Curated repeats must be a whole number of units")
    if not user:
        raise ValueError("No user to record the validation under")

    fields = [
        family_id, variant_key, sample, user, inheritance, validation,
        curated_repeats, comment, "0", timestamp or datetime.now().isoformat(),
    ]  # fmt: skip
    append_tsv_row(str_validation_file(data_dir), HEADER, fields)


def read_str_validations(
    data_dir: Path, family_id: str, variant_key: Optional[str] = None
) -> List[Dict[str, str]]:
    """A family's STR validations, oldest first, optionally at one locus."""
    path = str_validation_file(data_dir)
    if not path.exists():
        return []
    with open(path, "r") as f:
        return [
            row
            for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
            if row.get("FID") == family_id
            and (variant_key is None or row.get("Variant") == variant_key)
        ]
