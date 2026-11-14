"""Wombat tab data: the family's result files, their table rows and columns.

No UI here except reading the user's gene scorer: the tab renders what these
functions return.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import polars as pl

from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.models import Family
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    genomic_sort_key,
    get_column_group,
    get_column_sorting,
    get_column_type,
    get_display_label,
    reorder_columns_by_group,
)
from genetics_viz.utils.gene_scoring import (
    build_gene_badges,
    get_gene_scorer,
    split_genes,
)
from genetics_viz.utils.genetrek import get_snapshot
from genetics_viz.utils.inheritance import annotate_moi, family_pedigree
from genetics_viz.utils.variant_rows import add_score_badges, add_vep_badges
from genetics_viz.utils.wombat import (  # noqa: F401 - re-exported for the tab
    DENOVO_CONFIG,
    find_wombat_files,
    load_config_frame,
    variant_key,
    wombat_file_pattern,
)
from genetics_viz.utils.vep import get_highest_priority_consequence

logger = logging.getLogger(__name__)


def add_curation_status(
    rows: List[Dict[str, Any]], data_dir: Path, family_id: str
) -> None:
    """(Re)set every row's Validation and Diagnostic from the curation files."""
    validation_map = load_validation_map(
        data_dir / "validations" / "snvs.tsv", family_id
    )
    diagnostic_map = load_diagnostic_map(
        data_dir / "diagnostics" / "snvs.tsv", family_id
    )
    for row in rows:
        key, sample_id = variant_key(row), row.get("sample", "")
        add_validation_status_to_row(row, validation_map, key, sample_id)
        add_diagnostic_status_to_row(row, diagnostic_map, key, sample_id)


def build_rows(
    df: pl.DataFrame,
    data_dir: Path,
    family_id: str,
    family: Optional[Family],
    gene_scorer: Any = None,
) -> Tuple[List[Dict[str, Any]], Set[str], Set[str]]:
    """The table rows of one config, with their badges, curation status and MOI.

    ``df`` is the config's table after ``wombat.aggregate_by_variant``.
    Returns the rows, then the VEP consequences and ClinVar terms the colour
    maps lack (for ``filters.notify_unknown_terms``). Without ``gene_scorer``
    it reads the user's, so call it on the event loop; with one taken there
    first, it can run in a thread.
    """
    rows = df.to_dicts()
    unknown_consequences: Set[str] = set()
    unknown_clinvar: Set[str] = set()
    # Gene badges, coloured by the user's GeneTrek tiers
    if gene_scorer is None:
        gene_scorer = get_gene_scorer()
    for row in rows:
        row["Variant"] = variant_key(row)
        add_vep_badges(row, unknown_consequences, unknown_clinvar)
        # Consequence sort priority for custom sorting
        row["_consequence_priority"] = get_highest_priority_consequence(
            row.get("VEP_Consequence", "")
        )
        row["GeneBadges"] = build_gene_badges(
            split_genes(row.get("VEP_SYMBOL", "")), gene_scorer
        )
        row["VEP_Gene_badges"] = build_gene_badges(
            split_genes(row.get("VEP_Gene", ""), ","),
            gene_scorer,
            cap=None,
            sort=False,
        )
    add_curation_status(rows, data_dir, family_id)
    for row in rows:
        add_score_badges(row)
    # MOI: does each variant fit its gene's mode of inheritance?
    annotate_moi(
        rows,
        get_snapshot(),
        pedigree=family_pedigree(family.samples if family else []),
    )
    return rows, unknown_consequences, unknown_clinvar


#: The cohort and all-cohort columns, from the gene index's variant counts.
COUNT_COLUMNS = ["cohort_carriers", "cohort_freq", "seen_families", "seen_cohorts"]


def add_variant_counts(
    rows: List[Dict[str, Any]],
    counts: Optional[pl.DataFrame],
    data_dir: str,
    cohort: str,
    cohort_size: int,
    cohort_indexed: bool,
) -> None:
    """Set each row's :data:`COUNT_COLUMNS` from ``gene_index.variant_counts``.

    ``counts`` holds the rows' variants in every directory the user may open;
    ``None`` when the index is not loaded yet, which leaves the cells empty.
    Carriers are individuals of the family's own cohort, the frequency over
    ``cohort_size`` (its pedigree); both stay empty when the cohort has no
    indexed table. Families count every cohort, this family included: it
    carries the variant, its row says so, even when its cohort is not indexed.
    """
    if counts is None:
        return
    found: Dict[str, List[Dict[str, Any]]] = {}
    for c in counts.iter_rows(named=True):
        found.setdefault(c["variant"], []).append(c)
    for row in rows:
        here = found.get(row.get("Variant") or "", [])
        own = [c for c in here if c["data_dir"] == data_dir and c["cohort"] == cohort]
        others = [c for c in here if c not in own]
        if cohort_indexed:
            carriers = sum(c["carriers"] for c in own)
            row["cohort_carriers"] = carriers
            row["cohort_freq"] = (
                round(carriers / cohort_size, 4) if cohort_size else None
            )
            own_families = max(sum(c["families"] for c in own), 1)
        else:
            row["cohort_carriers"] = row["cohort_freq"] = None
            own_families = 1
        by_cohort = {cohort: own_families}
        for c in others:
            by_cohort[c["cohort"]] = by_cohort.get(c["cohort"], 0) + c["families"]
        row["seen_families"] = sum(by_cohort.values())
        row["seen_cohorts"] = ", ".join(
            f"{name} ({n})" for name, n in sorted(by_cohort.items())
        )


def wombat_columns(df_columns: List[str]) -> List[str]:
    """Every column the table can show: Variant and MOI first, the curation last.

    ``df_columns`` are the aggregated table's. The coordinates are folded into
    Variant; ``n_grouped`` and the ``_`` columns (the picked HGVS, for
    diagnostics) are not shown. The variant counts come before the curation.
    """
    return reorder_columns_by_group(
        ["Variant", "MOI"]
        + [
            col
            for col in df_columns
            if col not in ["#CHROM", "POS", "REF", "ALT", "n_grouped"]
            and not col.startswith("_")
            and col not in COUNT_COLUMNS
        ]
        + COUNT_COLUMNS
        + ["Validation", "Diagnostic"]
    )


#: Columns drawn as badge lists or gene badges: (cellType, badgesField, sortField).
_BADGE_COLUMNS: Dict[str, Tuple[str, str, Optional[str]]] = {
    "VEP_Consequence": ("badge_list", "ConsequenceBadges", "_consequence_priority"),
    "VEP_CLIN_SIG": ("badge_list", "ClinVarBadges", None),
    "VEP_SYMBOL": ("gene_badge", "GeneBadges", None),
    "VEP_Gene": ("gene_badge", "VEP_Gene_badges", None),
    "MOI": ("badge_list", "MoiBadges", "_moi_rank"),
    # The denovo step's evidence (De novo sub-tab)
    "dnm_strict": ("badge_list", "DnmStrictBadges", None),
    "dnm_flags": ("badge_list", "DnmFlagBadges", "_dnm_n_flags"),
    "dnm_origin": ("badge_list", "DnmOriginBadges", None),
}


def table_columns(all_columns: List[str]) -> List[Dict[str, Any]]:
    """The DataTable column definitions: a View action, then every column."""
    cols: List[Dict[str, Any]] = [
        {
            "id": "actions",
            "header": "",
            "cellType": "action",
            "actionName": "view_variant",
            "actionIcon": "visibility",
            "actionColor": "#1976d2",
            "actionTooltip": "View in IGV",
            "sortable": False,
        }
    ]
    for col in all_columns:
        col_def: Dict[str, Any] = {
            "id": col,
            "header": get_display_label(col),
            "group": get_column_group(col),
            "sorting": get_column_sorting(col),
            "sortable": True,
        }
        if col == "Validation":
            col_def["cellType"] = "validation"
        elif col == "seen_families":
            # The count links to the variant page: every carrier, every cohort
            col_def.update(cellType="link", href="/variant/{Variant}")
        elif col == "Diagnostic":
            col_def["cellType"] = "diagnostic"
        elif col in _BADGE_COLUMNS:
            cell_type, badges_field, sort_field = _BADGE_COLUMNS[col]
            col_def["cellType"] = cell_type
            col_def["badgesField"] = badges_field
            if sort_field:
                col_def["sortField"] = sort_field
        elif get_column_type(col) in ("int", "float"):
            col_def["cellType"] = "number"
        else:
            col_def["cellType"] = "score_badge"
        apply_width_constraints(col_def, col)
        cols.append(col_def)
    return cols


def presort_rows(
    rows: List[Dict[str, Any]],
    sorting: List[Dict[str, Any]],
    columns: List[Dict[str, Any]],
) -> None:
    """Sort ``rows`` in place as the table's saved sorting would.

    ``DataTable`` does not reorder its initial rows to match
    ``initial_sorting``, so a refresh would otherwise show them unsorted.
    Missing values sort last.
    """
    if not sorting:
        return
    col_id = sorting[0]["id"]
    desc = sorting[0].get("desc", False)
    col_def = next((c for c in columns if c.get("id") == col_id), {})
    sort_field = col_def.get("sortField", col_id)
    sort_type = col_def.get("sorting", "")
    if sort_type == "genomic":
        rows.sort(
            key=lambda r: (
                r.get(sort_field) is None,
                genomic_sort_key(r.get(sort_field, "")),
            ),
            reverse=desc,
        )
    elif sort_type == "numerical":

        def _num_key(r):
            v = r.get(sort_field)
            if v is None:
                return (True, 0.0)
            try:
                return (False, float(v))
            except (ValueError, TypeError):
                return (True, 0.0)

        rows.sort(key=_num_key, reverse=desc)
    else:
        rows.sort(
            key=lambda r: (r.get(sort_field) is None, r.get(sort_field, "")),
            reverse=desc,
        )


def count_rows(path: Path) -> int:
    """Data lines of a table, read as bytes: a count without parsing it."""
    try:
        with open(path, "rb") as f:
            lines = sum(1 for line in f if line.strip())
    except OSError:
        return 0
    return max(lines - 1, 0)
