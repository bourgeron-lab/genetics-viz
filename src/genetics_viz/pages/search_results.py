"""Search page results: the table, its presets, Save/Download and row actions.

Rows arrive built and filtered (``utils/search_rows.py``); this module shows
them and keeps them in step with curation saved from the dialogs.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.search_stats import show_stats_dialog
from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    get_column_group,
    get_column_sorting,
    get_column_type,
    get_display_label,
)
from genetics_viz.utils.export import is_export_column, rows_to_tsv
from genetics_viz.utils.genetrek import search_url
from genetics_viz.utils.search_rows import (
    CurationMaps,
    refresh_curation,
    sort_rows,
    sv_variant_key,
    validation_file_columns,
)
from genetics_viz.utils.view_presets import VIEW_PRESETS
from genetics_viz.utils.wisecondorx import WISECONDORX_CONFIG, build_call_colors

#: Cell renderers by column; any other column is a number or a score badge.
_CELL_TYPES: Dict[str, Dict[str, Any]] = {
    "Variant": {
        "cellType": "curated_locus",
        "curatedField": "IsCurated",
        "tooltipField": "_curated_tooltip",
    },
    "Validation": {"cellType": "validation"},
    "Diagnostic": {"cellType": "diagnostic"},
    "VEP_Consequence": {"cellType": "badge_list", "badgesField": "ConsequenceBadges"},
    "VEP_CLIN_SIG": {"cellType": "badge_list", "badgesField": "ClinVarBadges"},
    "gene": {"cellType": "gene_badge", "badgesField": "GeneBadges"},
    "VEP_SYMBOL": {"cellType": "gene_badge", "badgesField": "GeneBadges"},
    "VEP_Gene": {"cellType": "gene_badge", "badgesField": "VEP_Gene_badges"},
    "MOI": {
        "cellType": "badge_list",
        "badgesField": "MoiBadges",
        "sortField": "_moi_rank",
    },
    "FID": {"cellType": "link", "href": "/cohort/{_cohort_name}/family/{FID}"},
}


@dataclass
class SearchResults:
    """One search's rows and the table's view state across refreshes."""

    rows: List[Dict[str, Any]]
    all_columns: List[str]
    selected_cols: Dict[str, List[str]]
    selected_preset: Dict[str, str]
    has_wcx: bool
    maps: CurationMaps
    cohort_name: str
    sample_to_family: Dict[str, str]
    data_dir: Path
    page_client: Any
    #: Table state for persistence across refreshes
    table_state: Dict[str, Any] = field(
        default_factory=lambda: {"sorting": [], "page": 0}
    )


def result_column_defs(all_columns: List[str]) -> List[Dict[str, Any]]:
    """The table's columns: the action, then every result column."""
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
        if col in _CELL_TYPES:
            col_def.update(_CELL_TYPES[col])
        elif col == "wisecondorX":
            col_def["cellType"] = "cnv_call"
            col_def["callColors"] = build_call_colors()
            col_def["tooltipField"] = "_wisecondorx_tooltip"
        elif get_column_type(col) in ("int", "float"):
            col_def["cellType"] = "number"
        else:
            col_def["cellType"] = "score_badge"
        apply_width_constraints(col_def, col)
        cols.append(col_def)
    return cols


def render_results(res: SearchResults) -> None:
    """Build the results table in the current container."""

    @ui.refreshable
    def render_results_table():
        # All per-component filters (impact, exclude, validation, calls)
        # are already applied — use the rows directly
        rows = res.rows.copy()
        # An SV's region (curated when curated) opens GeneTrek
        # on the genes it overlaps.
        for r in rows:
            if r.get("_source_type") == "wisecondorx":
                r["_genetrek_url"] = search_url([r["Variant"]])

        columns = result_column_defs(res.all_columns)

        # Reference to the DataTable for column visibility updates
        search_dt: Dict[str, Any] = {"ref": None}

        def _apply_col_visibility():
            if search_dt["ref"]:
                visible = ["actions"] + list(res.selected_cols["value"])
                search_dt["ref"].set_column_visibility(visible)

        _toolbar(res, rows, _apply_col_visibility)

        # CNV calls legend (only when WCX data is present)
        if res.has_wcx:
            _cnv_legend()

        def _refresh_all_rows():
            """Reload validation + diagnostic data for all rows."""
            refresh_curation(res.rows, res.maps.load())
            with res.page_client:
                render_results_table.refresh()

        # Restore table state (sorting / page) across refreshes
        saved_sorting = res.table_state.get("sorting", [])
        if saved_sorting:
            col_id = saved_sorting[0]["id"]
            col_def = next((c for c in columns if c.get("id") == col_id), {})
            sort_rows(
                rows,
                col_def.get("sortField", col_id),
                col_def.get("sorting", ""),
                saved_sorting[0].get("desc", False),
            )

        # Create table with all columns, initial visibility from preset
        search_dt["ref"] = DataTable(
            columns=columns,
            rows=rows,
            row_key="Variant",
            pagination={"rowsPerPage": 50},
            visible_columns=["actions"] + list(res.selected_cols["value"]),
            on_row_action=lambda e: _on_row_action(res, e, _refresh_all_rows),
            initial_sorting=saved_sorting,
            initial_page=res.table_state.get("page", 0),
            state_holder=res.table_state,
        )

    render_results_table()


def _toolbar(res: SearchResults, rows: List[Dict[str, Any]], apply_visibility) -> None:
    """Count, preset, Columns, Stats and Save above the table."""
    with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
        ui.label(f"Results ({len(rows)} rows)").classes(
            "text-lg font-semibold text-blue-700"
        )

        # Preset dropdown
        preset_select = ui.select(
            options={p["name"]: p["name"] for p in VIEW_PRESETS},
            value=res.selected_preset["name"],
            label="Preset",
        ).classes("w-48")

        ui.space()  # Push column selector to the right

        # Column selector dialog
        col_dialog, _sync_col_selector = build_column_selector(
            all_columns=res.all_columns,
            selected_cols=res.selected_cols,
            on_visibility_change=apply_visibility,
            presets=VIEW_PRESETS,
        )
        ui.button(
            "Columns",
            icon="view_column",
            on_click=col_dialog.open,
        ).props("outline color=blue size=sm")

        # --- Stats button + dialog ---
        ui.button(
            "Stats",
            icon="bar_chart",
            on_click=lambda _e=None, r=rows: show_stats_dialog(r),
        ).props("outline color=blue size=sm")

        _save_menu(res, rows)

    # Preset change handler
    def on_preset_change(e):
        """Handle preset selection change."""
        preset = next((p for p in VIEW_PRESETS if p["name"] == e.value), None)
        if not preset:
            return
        # Filter columns to only those available in the data
        res.selected_cols["value"] = [
            col for col in preset.get("columns", []) if col in res.all_columns
        ]
        res.selected_preset["name"] = e.value
        apply_visibility()
        _sync_col_selector()

    # Connect preset change handler
    preset_select.on_value_change(on_preset_change)


def _save_menu(res: SearchResults, rows: List[Dict[str, Any]]) -> None:
    """Save as a validation file (curators) or download the shown columns."""

    def _save_validation_file(filename: str):
        """Write TSV to to_validate directory."""
        save_cols = validation_file_columns(
            res.selected_preset["name"], res.all_columns, VIEW_PRESETS
        )
        tsv_content = rows_to_tsv(rows, save_cols)
        to_validate_dir = res.data_dir / "to_validate"
        to_validate_dir.mkdir(parents=True, exist_ok=True)
        file_path = to_validate_dir / filename
        if file_path.exists():
            ui.notify(f"File '{filename}' already exists!", type="warning")
            return
        file_path.write_text(tsv_content)
        ui.notify(
            f"Saved {filename} to to_validate/",
            type="positive",
            position="top",
        )

    def _download_tsv(filename: str):
        """Trigger browser download of TSV with displayed columns."""
        download_cols = [
            col for col in res.selected_cols["value"] if is_export_column(col)
        ]
        tsv_content = rows_to_tsv(rows, download_cols)
        ui.download(
            tsv_content.encode("utf-8"),
            filename,
            media_type="text/tab-separated-values",
        )

    save_dlg, save_input = _build_save_dialog(
        "Save as validation file", "Save", "save", _save_validation_file
    )
    dl_dlg, dl_input = _build_save_dialog(
        "Download TSV", "Download", "download", _download_tsv
    )

    with ui.button("Save", icon="save").props("outline color=blue size=sm"):
        with ui.menu().classes("min-w-[250px]"):
            if can_write():
                ui.menu_item(
                    "Save as validation file",
                    on_click=lambda: _open_dialog(save_dlg, save_input),
                )
            ui.menu_item(
                "Download TSV",
                on_click=lambda: _open_dialog(dl_dlg, dl_input),
            )


def _build_save_dialog(title, confirm_label, confirm_icon, on_confirm):
    """Pre-build a save/download dialog. Returns (dialog, name_input)."""
    with ui.dialog() as dlg, ui.card().classes("p-4 min-w-[400px]"):
        ui.label(title).classes("text-lg font-semibold mb-2")

        with ui.row().classes("items-center w-full gap-1"):
            name_input = (
                ui.input(label="Filename", value="")
                .props("outlined dense")
                .classes("flex-grow")
            )
            ui.label(".tsv").classes("text-gray-500 text-sm mt-2")

        with ui.row().classes("justify-end gap-2 w-full mt-4"):
            ui.button("Cancel", on_click=dlg.close).props("flat")

            def _on_confirm():
                raw_name = name_input.value.strip()
                if not raw_name:
                    ui.notify("Please enter a filename", type="warning")
                    return
                if not raw_name.endswith(".tsv"):
                    raw_name += ".tsv"
                dlg.close()
                on_confirm(raw_name)

            ui.button(confirm_label, icon=confirm_icon, on_click=_on_confirm).props(
                "unelevated color=blue"
            )
    return dlg, name_input


def _open_dialog(dlg, name_input):
    """Reset filename, open dialog, auto-focus and select text."""
    name_input.value = f"genetics-viz.{datetime.now().strftime('%Y-%m-%d.%H-%M')}"
    dlg.open()
    ui.run_javascript(
        f"setTimeout(function() {{"
        f'  var el = getElement({name_input.id}).$el.querySelector("input");'
        f"  if(el){{el.focus();el.select();}}"
        f"}}, 200);"
    )


def _cnv_legend() -> None:
    with ui.row().classes("gap-4 items-center"):
        ui.label("CNV calls:").classes("text-sm font-semibold")
        for _cfg_key in [
            "robust_loss",
            "permissive_loss",
            "robust_gain",
            "permissive_gain",
        ]:
            cfg = WISECONDORX_CONFIG[_cfg_key]
            _op = "≤" if "loss" in _cfg_key else "≥"
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {cfg["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{cfg["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2{_op}{cfg['ratio_threshold']} & Z{_op}{cfg['zscore_threshold']})"
                ).classes("text-xs text-gray-600")


def _sv_gene_and_impact(row_data: Dict[str, Any]) -> tuple[str, str]:
    """An SV row's genes and GAIN/LOSS, as its diagnostic records them."""
    gene_str = row_data.get("gene", "")
    gene_syms = []
    if gene_str and gene_str != "-":
        for gp in str(gene_str).split(","):
            sym = gp.split(":")[0].strip()
            if sym:
                gene_syms.append(sym)
    gene = ", ".join(gene_syms) if gene_syms else ""
    sv_call = str(row_data.get("wisecondorX", row_data.get("call", ""))).upper()
    impact = "GAIN" if "GAIN" in sv_call else ("LOSS" if "LOSS" in sv_call else sv_call)
    return gene, impact


def _on_row_action(res: SearchResults, e: Dict[str, Any], refresh_all_rows) -> None:
    """Handle row action click: a diagnostic, or the SNV/SV dialog."""
    action = e.get("action", "")
    row_data = e.get("row", {})
    variant_str = row_data.get("Variant", "")
    sample_id = row_data.get("sample", "")
    source_type = row_data.get("_source_type", "wombat")

    # Get family from sample
    family_id = res.sample_to_family.get(sample_id)

    if not family_id:
        ui.notify(
            f"Could not find family for sample {sample_id}",
            type="warning",
        )
        return

    if action == "open_diagnostic":
        if source_type == "wisecondorx":
            # SV diagnostic
            gene, impact = _sv_gene_and_impact(row_data)
            show_diagnostic_dialog(
                family_id=family_id,
                variant_key=sv_variant_key(row_data),
                gene=gene,
                impact=impact,
                sample=sample_id,
                variant_type="sv",
                diagnostic_file=res.maps.sv_diagnostic_file,
                on_save_callback=lambda s: refresh_all_rows(),
            )
        else:
            # SNV diagnostic
            show_diagnostic_dialog(
                family_id=family_id,
                variant_key=variant_str,
                gene=row_data.get("VEP_SYMBOL", ""),
                impact=row_data.get("VEP_Consequence", ""),
                sample=sample_id,
                variant_type="snv",
                diagnostic_file=res.maps.snv_diagnostic_file,
                on_save_callback=lambda s: refresh_all_rows(),
            )
        return

    try:
        if source_type == "wisecondorx":
            # SV dialog
            locus_for_dialog = row_data.get("_original_locus", variant_str)
            parts = locus_for_dialog.split(":")
            range_parts = parts[1].split("-") if len(parts) == 2 else []
            if len(range_parts) == 2:
                show_sv_dialog(
                    cohort_name=res.cohort_name,
                    family_id=family_id,
                    chrom=parts[0],
                    start=range_parts[0],
                    end=range_parts[1],
                    sample=sample_id,
                    sv_data={
                        **row_data,
                        "call": row_data.get("wisecondorX", row_data.get("call", "")),
                    },
                    on_validation_saved=refresh_all_rows,
                )
            else:
                ui.notify("Invalid SV format", type="warning")
        else:
            # Wombat variant dialog
            parts = variant_str.split(":")
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
                show_variant_dialog(
                    cohort_name=res.cohort_name,
                    family_id=family_id,
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    sample=sample_id,
                    variant_data=dict(row_data),
                    on_save_callback=lambda s: refresh_all_rows(),
                )
            else:
                ui.notify(
                    "Invalid variant format. Expected chr:pos:ref:alt",
                    type="warning",
                )
    except Exception as ex:
        ui.notify(f"Error parsing variant: {ex}", type="warning")
