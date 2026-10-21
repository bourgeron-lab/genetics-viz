"""What the family Diagnostics panel says about an STR or UPD diagnostic.

A diagnostic stores its record's key, which is unreadable on its own for these
kinds (``chr19:45770204-45770264:DMPK``, ``chr8:UPD``). The panel shows the
kind and what was found instead:

- an STR reads ``STR`` ``~209 repeats``, after its gene badge. The count is a
  curator's ``CuratedRepeats`` when the member's standing validation records
  one, else ExpansionHunter's longer allele (``~``) from the family genotype
  table of the newest catalog;
- a UPD reads ``UPD`` ``chr8 · maternal``, the origin from its standing
  validation.

The panel draws the kind as a tag; :func:`describe_str` and
:func:`describe_upd` return the text after it and the tooltip.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import polars as pl

from .str_eh import (
    eh_path,
    find_eh_labels,
    get_family_eh_dir,
    longest_allele,
    read_eh_table_cached,
)
from .str_validation import read_str_validations
from .upd import (
    chrom_of_key,
    latest_validation,
    origin_label,
    read_upd_validations,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RepeatCount:
    """One member's repeat count at a locus."""

    count: int
    #: A curator's ``CuratedRepeats``, not ExpansionHunter's estimate.
    curated: bool
    #: ``"5/209"``, ExpansionHunter's repeat counts, for the tooltip.
    genotype: str = ""
    motif: str = ""

    @property
    def text(self) -> str:
        return f"{self.count}" if self.curated else f"~{self.count}"


def str_locus(key: str) -> str:
    """``DMPK`` for ``chr19:45770204-45770264:DMPK``."""
    return key.rsplit(":", 1)[-1]


def _family_genotypes(data_dir: Path, family_id: str) -> Optional[pl.DataFrame]:
    directory = get_family_eh_dir(data_dir, family_id)
    labels = find_eh_labels(directory, family_id)
    if not labels:
        return None
    try:
        return read_eh_table_cached(
            eh_path(directory, family_id, labels[0], "genotypes")
        )
    except Exception:
        logger.warning(
            "Failed to read the STR genotypes of %s", family_id, exc_info=True
        )
        return None


def str_repeat_counts(
    data_dir: Path, family_id: str, key: str, samples: List[str]
) -> Dict[str, RepeatCount]:
    """Each member's repeat count at an STR key; members with none are left out."""
    validations = read_str_validations(data_dir, family_id, key)
    fam = _family_genotypes(data_dir, family_id)
    rows: Dict[str, Dict[str, Any]] = {}
    if fam is not None and "locus" in fam.columns:
        for row in fam.filter(
            (pl.col("locus") == str_locus(key)) & pl.col("sample").is_in(samples)
        ).to_dicts():
            rows[row["sample"]] = row

    out: Dict[str, RepeatCount] = {}
    for sample in samples:
        row = rows.get(sample, {})
        genotype = "/".join(
            str(row[k]) for k in ("allele1", "allele2") if row.get(k) is not None
        )
        motif = row.get("motif") or ""
        curated = (latest_validation(validations, sample) or {}).get("CuratedRepeats")
        if curated and curated.strip().isdigit():
            out[sample] = RepeatCount(int(curated), True, genotype, motif)
            continue
        count = longest_allele(row) if row else None
        if count is not None:
            out[sample] = RepeatCount(count, False, genotype, motif)
    return out


def describe_str(
    data_dir: Path, family_id: str, entry: Dict[str, Any]
) -> Tuple[str, str]:
    """``("~209 repeats", tooltip)`` for a merged diagnostic row; ``""`` with no count."""
    key = entry.get("Variant", "")
    samples = [s.strip() for s in (entry.get("Sample") or "").split(",") if s.strip()]
    counts = str_repeat_counts(data_dir, family_id, key, samples)
    if not counts:
        text = ""
    elif len(samples) == 1:
        text = f"{counts[samples[0]].text} repeats"
    else:
        text = (
            " · ".join(f"{s} {counts[s].text}" for s in samples if s in counts)
            + " repeats"
        )
    if any(c.curated for c in counts.values()):
        text += " (curated)"

    lines = [f"Short tandem repeat at {key.rsplit(':', 1)[0]} ({str_locus(key)})"]
    motif = next((c.motif for c in counts.values() if c.motif), "")
    if motif:
        lines.append(f"Motif: {motif}")
    for sample in samples:
        c = counts.get(sample)
        if c is None:
            lines.append(f"{sample}: no genotype in the family table")
            continue
        detail = "curated" if c.curated else "ExpansionHunter, longer allele"
        if c.genotype:
            detail += f"; genotype {c.genotype}"
        lines.append(f"{sample}: {c.text} repeats ({detail})")
    if entry.get("Impact"):
        lines.append(f"Status: {entry['Impact']}")
    return text, "\n".join(lines)


def describe_upd(
    data_dir: Path, family_id: str, entry: Dict[str, Any]
) -> Tuple[str, str]:
    """``("chr8 · maternal", tooltip)`` for a merged diagnostic row."""
    key = entry.get("Variant", "")
    chrom = chrom_of_key(key)
    samples = [s.strip() for s in (entry.get("Sample") or "").split(",") if s.strip()]
    validations = read_upd_validations(data_dir, family_id, key)
    origins = []
    lines = [f"Uniparental disomy of chromosome {chrom}"]
    for sample in samples:
        latest = latest_validation(validations, sample)
        if latest is None:
            lines.append(f"{sample}: not validated")
            continue
        origin = origin_label(latest.get("Inheritance") or "")
        if origin:
            origins.append(origin)
        lines.append(
            f"{sample}: validated {latest.get('Validation', '')}"
            + (f", {origin}" if origin else "")
        )
    if entry.get("Impact"):
        lines.append(f"Found by: {entry['Impact']}")
    text = f"chr{chrom}"
    if origins:
        text += " · " + "/".join(dict.fromkeys(origins))
    return text, "\n".join(lines)
