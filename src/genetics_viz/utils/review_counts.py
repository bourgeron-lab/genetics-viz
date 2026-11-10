"""What the family page's tab badges count: content a curator has not looked at yet.

Each count covers the whole family, not the member selection, like the header
QC badge, and is blocking but pure: the page runs it off the event loop and
draws the number on the tab (``components/tab_badge.py``).

- **Sample QC**: unacknowledged high and medium flags of the family's own QC
  files. ``UPD_CANDIDATE`` flags are left to the UPD tab when there is one.
- **Wombat**: de novo calls (configs matching ``de_?novo``), one per variant
  and carrier, after the default exclusions (LCR, gnomAD-filtered, likely
  MNV), of HIGH or MODERATE impact, with no validation in ``snvs.tsv``. When
  Sample QC holds the de novo table back (a wrong relationship), the table is
  not read: the badge says why instead of counting.
- **STR**: genotypes of the newest catalog the pipeline flags (severity high
  or medium), one per locus and member, with no validation in ``strs.tsv``.
- **UPD**: the candidates no curator validated present or absent
  (:func:`genetics_viz.utils.upd.count_to_review`).
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)

from . import upd
from .sample_qc import (
    FlagRow,
    count_by_severity,
    gates_denovo,
    load_family_qc,
    severity_summary,
)
from .sample_qc_acks import load_acks, split_acknowledged, unexplained_denovo_risk
from .sharding import get_family_path
from .str_eh import eh_path, find_eh_labels, flagged_expr, get_family_eh_dir
from .str_eh import read_eh_table_cached
from .str_validation import curated_statuses, str_validation_file
from .validation_status import UNVALIDATED, row_status
from .variant_rows import gnomad_filtered, in_lcr, likely_mnv, split_terms
from .vep import get_consequence_impact
from .wombat import (
    DENOVO_CONFIG,
    aggregate_by_variant,
    find_wombat_files,
    load_config_frame,
    variant_key,
)

logger = logging.getLogger(__name__)

#: Impacts a de novo call needs to count, most severe first.
DENOVO_IMPACTS = ("HIGH", "MODERATE")


@dataclass
class ReviewCount:
    """One tab's badge: how many, how severe the worst, and why."""

    n: int = 0
    #: ``high`` or ``medium``: the badge's colour.
    worst: str = "medium"
    tooltip: str = ""
    #: Nothing was counted because the content is held back; the badge says so.
    held: bool = False

    @property
    def shown(self) -> bool:
        return self.n > 0 or self.held


def qc_to_review(flags: Iterable[FlagRow], upd_tab: bool) -> ReviewCount:
    """The unacknowledged high and medium flags among ``flags``.

    ``flags`` are the family's active (unacknowledged) flags. With ``upd_tab``
    a ``UPD_CANDIDATE`` flag is left out: the UPD tab counts it.
    """
    flags = list(flags)
    kept = [f for f in flags if not (upd_tab and upd.is_upd_flag(f))]
    counts = count_by_severity(kept)
    n = counts.get("high", 0) + counts.get("medium", 0)
    if not n:
        return ReviewCount()
    tooltip = f"Sample QC: {severity_summary(counts)} not acknowledged"
    left_out = count_by_severity([f for f in flags if f not in kept])
    if left_out.get("high", 0) + left_out.get("medium", 0):
        tooltip += " (UPD candidates are counted on the UPD tab)"
    return ReviewCount(
        n=n, worst="high" if counts.get("high") else "medium", tooltip=tooltip
    )


def _impact(row: Dict) -> str:
    """The most severe of ``DENOVO_IMPACTS`` a call has, else ``""``.

    ``VEP_IMPACT`` when the table has it (comma-joined over transcripts after
    :func:`aggregate_by_variant`), else the impacts of its consequence terms.
    """
    raw = row.get("VEP_IMPACT")
    if raw:
        impacts = {p.strip().upper() for p in str(raw).split(",")}
    else:
        impacts = {
            get_consequence_impact(t) for t in split_terms(row.get("VEP_Consequence"))
        }
    return next((i for i in DENOVO_IMPACTS if i in impacts), "")


def denovo_to_review(data_dir: Path, family_id: str) -> ReviewCount:
    """The family's de novo calls of HIGH or MODERATE impact never validated."""
    wombat_dir = get_family_path(data_dir, family_id) / "wombat"
    if not wombat_dir.is_dir():
        return ReviewCount()
    files = [
        wf
        for wf in find_wombat_files(wombat_dir, family_id)
        if DENOVO_CONFIG.search(wf["wombat_config"])
    ]
    if not files:
        return ReviewCount()
    if gates_denovo(unexplained_denovo_risk(data_dir, family_id)):
        return ReviewCount(
            held=True,
            worst="high",
            tooltip="De novo calls held back: Sample QC finds a wrong relationship"
            " in this family",
        )
    validations = load_validation_map(data_dir / "validations" / "snvs.tsv", family_id)
    found: Dict[Tuple[str, str], str] = {}
    for wf in files:
        df = load_config_frame(wf["file_path"])
        if df is None:
            continue
        for row in aggregate_by_variant(df).to_dicts():
            if in_lcr(row) or gnomad_filtered(row) or likely_mnv(row):
                continue
            impact = _impact(row)
            if not impact:
                continue
            key, sample = variant_key(row), row.get("sample") or ""
            status: Dict = {}
            add_validation_status_to_row(status, validations, key, sample)
            if row_status(status) != UNVALIDATED:
                continue
            if found.get((key, sample)) != "HIGH":
                found[(key, sample)] = impact
    if not found:
        return ReviewCount()
    n_high = sum(1 for i in found.values() if i == "HIGH")
    return ReviewCount(
        n=len(found),
        worst="high" if n_high else "medium",
        tooltip=f"{len(found)} de novo call(s) of HIGH or MODERATE impact not yet"
        f" validated ({n_high} HIGH), after the LCR, gnomAD and MNV exclusions",
    )


def str_to_review(data_dir: Path, family_id: str) -> ReviewCount:
    """The flagged genotypes of the newest catalog no curator validated."""
    directory = get_family_eh_dir(data_dir, family_id)
    labels = find_eh_labels(directory, family_id)
    if not labels:
        return ReviewCount()
    fam = read_eh_table_cached(eh_path(directory, family_id, labels[0], "genotypes"))
    if fam is None or len(fam) == 0:
        return ReviewCount()
    flagged = fam.filter(flagged_expr(fam))
    if len(flagged) == 0:
        return ReviewCount()
    statuses = curated_statuses(
        load_validation_map(str_validation_file(data_dir), family_id),
        flagged,
        add_validation_status_to_row,
        "Validation",
    )
    open_rows = [
        r
        for r in flagged.to_dicts()
        if (r.get("locus"), r.get("sample")) not in statuses
    ]
    if not open_rows:
        return ReviewCount()
    pairs = {(r.get("locus"), r.get("sample")) for r in open_rows}
    high = {
        (r.get("locus"), r.get("sample"))
        for r in open_rows
        if r.get("severity") == "high"
    }
    return ReviewCount(
        n=len(pairs),
        worst="high" if high else "medium",
        tooltip=f"{len(pairs)} flagged STR genotype(s) not yet validated"
        f" ({len(high)} in the pathogenic range)",
    )


def upd_to_review(data_dir: Path, family_id: str) -> ReviewCount:
    """The UPD candidates no curator validated present or absent."""
    n = upd.count_to_review(upd.family_candidates(data_dir, family_id))
    if not n:
        return ReviewCount()
    return ReviewCount(
        n=n, worst="medium", tooltip=f"{n} UPD candidate(s) no curator settled"
    )


def qc_family_to_review(data_dir: Path, family_id: str, upd_tab: bool) -> ReviewCount:
    """:func:`qc_to_review` of the family's own QC files, acknowledgements applied."""
    try:
        qc = load_family_qc(data_dir, family_id)
    except Exception:
        logger.warning("Failed to read sample QC of %s", family_id, exc_info=True)
        return ReviewCount()
    if qc is None:
        return ReviewCount()
    active, _ = split_acknowledged(qc.flags, load_acks(data_dir))
    return qc_to_review(active, upd_tab)
