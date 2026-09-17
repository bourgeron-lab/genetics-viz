"""Polygenic scores tab component for the family page.

The pipeline writes z-scores wide - one row per sample, one column per score,
around a hundred of them. A hundred columns is not a usable table, so the frame
is transposed here: one row per trait, one column per family member.
"""

import logging
from typing import Any, Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.ancestry import (
    NULL_VALUES,
    find_ancestry_files,
    get_ancestry_dir,
    get_pgs_zscore_thresholds,
    probe_pgs_data,
)
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

__all__ = ["probe_pgs_data", "render_pgs_tab"]

#: The score name as the pipeline writes it, e.g. ``brain_biton2020/ICV`` or
#: ``social_behaviour_saint-pourcain2025/PSP/earlychildhood``. Kept whole rather
#: than split on "/", so the column reads as the trait it identifies.
_TRAIT_COL = "Trait"

#: The longest trait name in the current catalogue runs to 54 characters.
_TRAIT_MIN_WIDTH = 380


def _to_float(value: Any) -> Optional[float]:
    """Coerce a cell to a float, returning None when it is not numeric.

    A score column that is entirely ``nan`` in the source file comes back from
    polars as an all-null string column, so the values reaching here are not
    always numbers.
    """
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if result != result else result  # drop NaN


def render_pgs_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the Polygenic Scores tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    ancestry_dir = get_ancestry_dir(store.data_dir, family_id)
    if not ancestry_dir.exists():
        ui.label(f"No ancestry directory found at: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    files = find_ancestry_files(store.data_dir, family_id)
    if files is None or files.pgs_zscore_path is None:
        ui.label(f"No PGS z-scores found in: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    try:
        df = read_tsv_or_none(
            files.pgs_zscore_path,
            infer_schema_length=10000,
            null_values=NULL_VALUES,
        )
    except Exception as e:
        ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
        logger.warning("Failed to read %s", files.pgs_zscore_path, exc_info=True)
        return

    if df is None or len(df) == 0:
        ui.label("No PGS z-scores in file").classes("text-gray-500 italic")
        return

    if "IID" not in df.columns:
        ui.label("PGS z-score file has no IID column").classes("text-red-500 mt-4")
        return

    score_names = [c for c in df.columns if c != "IID"]
    if not score_names:
        ui.label("PGS z-score file has no score columns").classes(
            "text-gray-500 italic"
        )
        return

    # Transpose: sample_id -> {score_name: value}
    by_sample: Dict[str, Dict[str, Any]] = {
        row["IID"]: row for row in df.to_dicts() if row.get("IID")
    }
    thresholds = get_pgs_zscore_thresholds()

    with ui.row().classes("items-center gap-2 mb-2"):
        ui.label(f"{len(score_names)} polygenic scores").classes(
            "text-lg font-semibold text-blue-700"
        )
        ui.label(files.label).classes("text-xs text-gray-500")

    @ui.refreshable
    def render_pgs_table() -> None:
        members = [iid for iid in by_sample if iid in selected_members["value"]]
        if not members:
            ui.label("No members selected").classes("text-gray-500 italic")
            return

        all_rows: List[Dict[str, Any]] = []
        for score in score_names:
            row: Dict[str, Any] = {_TRAIT_COL: score}
            for iid in members:
                row[iid] = _to_float(by_sample[iid].get(score))
            all_rows.append(row)

        columns: List[Dict[str, Any]] = [
            {
                "id": _TRAIT_COL,
                "header": _TRAIT_COL,
                "sortable": True,
                "minWidth": _TRAIT_MIN_WIDTH,
                "filter": {"type": "text", "placeholder": "Filter traits..."},
            },
        ]
        for iid in members:
            columns.append(
                {
                    "id": iid,
                    "header": iid,
                    "sortable": True,
                    "sorting": "numerical",
                    "cellType": "color_scale",
                    "thresholds": thresholds,
                    "minWidth": 100,
                }
            )

        dt_ref: Dict[str, Any] = {"dt": None}

        def on_filter(e: Dict[str, Any]) -> None:
            filters = e.get("filters", {})
            filtered = all_rows

            trait_text = (filters.get(_TRAIT_COL) or "").strip().lower()
            if trait_text:
                filtered = [r for r in filtered if trait_text in r[_TRAIT_COL].lower()]

            if dt_ref["dt"]:
                dt_ref["dt"].update_data(filtered)

        with ui.card().classes("w-full"):
            dt_ref["dt"] = DataTable(
                columns=columns,
                rows=all_rows,
                row_key=_TRAIT_COL,
                pagination={"rowsPerPage": 25},
                sticky_first_column=True,
                on_filter=on_filter,
            )

    render_pgs_table()
    data_table_refreshers.append(render_pgs_table.refresh)
