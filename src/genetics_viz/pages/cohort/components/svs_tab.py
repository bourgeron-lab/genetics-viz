"""SVs tab component for family page."""

from pathlib import Path
from typing import Any, Callable, Dict, List

from nicegui import ui

from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.utils.sharding import get_family_path
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    get_column_group,
    get_column_sorting,
    get_column_type,
    get_display_label,
    reorder_columns_by_group,
)
from genetics_viz.utils.cytobands import CHROM_ORDER, VALIDATION_COLORS
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.wisecondorx import (
    WISECONDORX_CONFIG,
    build_call_colors,
    build_color_thresholds,
    infer_sv_type,
    parse_wisecondorx_bed_for_display,
)


_GENE_BADGE_COLUMNS = {
    "genic_symbol",
    "genic_ensg",
    "exonic_symbol",
    "exonic_ensg",
    "VEP_Gene",
}


def render_svs_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the SVs tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    svs_dir = get_family_path(store.data_dir, family_id) / "svs"

    if not svs_dir.exists():
        ui.label(f"No SVs directory found at: {svs_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Create subtabs for different SV callers
    with ui.tabs().classes("w-full") as svs_subtabs:
        wisecondorx_tab = ui.tab("WisecondorX")

    with ui.tab_panels(svs_subtabs, value=wisecondorx_tab).classes("w-full"):
        # WisecondorX subtab
        with ui.tab_panel(wisecondorx_tab):
            render_wisecondorx_subtab(
                store=store,
                family_id=family_id,
                svs_dir=svs_dir,
                selected_members=selected_members,
                data_table_refreshers=data_table_refreshers,
                cohort_name=cohort_name,
            )


def render_wisecondorx_subtab(
    store: Any,
    family_id: str,
    svs_dir: Path,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    cohort_name: str,
) -> None:
    """Render the WisecondorX subtab content.

    Args:
        store: DataStore instance
        family_id: Family ID
        svs_dir: Path to SVs directory
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
        cohort_name: Cohort name
    """
    wisecondorx_dir = svs_dir / "wisecondorx"
    aberrations_file = wisecondorx_dir / f"{family_id}_aberrations.annotated.bed"

    if not wisecondorx_dir.exists():
        ui.label(f"No WisecondorX directory found at: {wisecondorx_dir}").classes(
            "text-gray-500 italic"
        )
        return

    if not aberrations_file.exists():
        ui.label(f"No aberrations file found at: {aberrations_file}").classes(
            "text-gray-500 italic"
        )
        return

    with ui.card().classes("w-full p-4"):
        ui.label("WisecondorX Aberrations").classes(
            "text-lg font-semibold text-blue-700 mb-2"
        )
        with ui.row().classes("gap-4"):
            ui.label("File Path:").classes("font-semibold")
            ui.label(str(aberrations_file)).classes("text-sm text-gray-600 font-mono")

        # Gene badge legend
        with ui.row().classes("gap-4 mt-3 items-center"):
            ui.label("Gene badges:").classes("text-sm font-semibold")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #ffffff; color: black; border: 2px solid black; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(exonic - black border)").classes("text-xs text-gray-600")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #ffffff; color: black; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(genic - no border)").classes("text-xs text-gray-600")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #8b0000; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(color indicates geneset importance)").classes(
                    "text-xs text-gray-600"
                )

        # CNV call legend
        with ui.row().classes("gap-4 mt-2 items-center"):
            ui.label("CNV calls:").classes("text-sm font-semibold")

            robust_loss = WISECONDORX_CONFIG["robust_loss"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {robust_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_loss["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≤{robust_loss['ratio_threshold']} & Z≤{robust_loss['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {permissive_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_loss["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≤{permissive_loss['ratio_threshold']} & Z≤{permissive_loss['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            robust_gain = WISECONDORX_CONFIG["robust_gain"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {robust_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_gain["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≥{robust_gain['ratio_threshold']} & Z≥{robust_gain['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            permissive_gain = WISECONDORX_CONFIG["permissive_gain"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {permissive_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_gain["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≥{permissive_gain['ratio_threshold']} & Z≥{permissive_gain['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

    # Display BED file content in a table
    try:
        df = parse_wisecondorx_bed_for_display(aberrations_file)
        if df is None or len(df) == 0:
            ui.label("File is empty").classes("text-gray-500 italic")
            return

        # Convert to list of dicts for NiceGUI table
        all_rows = df.to_dicts()

        # Store original chr:start-end for each row before any modifications
        for row in all_rows:
            if "chr:start-end" in row:
                row["_original_locus"] = row["chr:start-end"]

        # Function to reload and apply validation + diagnostic data
        def reload_validations():
            """Reload validation and diagnostic data and update all rows."""
            # Load validation data from svs.tsv
            validation_file = store.data_dir / "validations" / "svs.tsv"
            validation_map = load_validation_map(validation_file, family_id)
            # Load diagnostic data from svs.tsv
            diagnostic_file = store.data_dir / "diagnostics" / "svs.tsv"
            diagnostic_map = load_diagnostic_map(diagnostic_file, family_id)

            # Add Validation status to each row
            for row in all_rows:
                # Reset to original locus and clear curated flags
                if "_original_locus" in row:
                    row["chr:start-end"] = row["_original_locus"]
                row["IsCurated"] = False
                row["_curated_tooltip"] = ""
                row.pop("OriginalLocus", None)

                # For SVs, variant_key is the original chr:start-end format
                variant_key = row.get("_original_locus", row.get("chr:start-end", ""))

                sample_id = row.get("sample", "")

                # Construct variant key in the format stored in svs.tsv
                # Format: chr:start-end:type (e.g., chr1:1000-2000:del or chr1:1000-2000:dup)
                # First, determine type from call
                sv_call = row.get("call", "")
                if (
                    "GAIN" in str(sv_call).upper()
                    or "gain" in str(sv_call).lower()
                    or "Gain" in str(sv_call)
                ):
                    sv_type = "dup"
                elif (
                    "LOSS" in str(sv_call).upper()
                    or "loss" in str(sv_call).lower()
                    or "Loss" in str(sv_call)
                ):
                    sv_type = "del"
                else:
                    # Try to infer from ratio
                    ratio = row.get("ratio", 0)
                    try:
                        ratio_val = float(ratio) if ratio else 0
                        sv_type = "dup" if ratio_val > 0 else "del"
                    except (ValueError, TypeError):
                        sv_type = "del"  # Default to deletion

                # Construct the full variant key
                full_variant_key = f"{variant_key}:{sv_type}"

                add_validation_status_to_row(
                    row, validation_map, full_variant_key, sample_id
                )
                add_diagnostic_status_to_row(
                    row, diagnostic_map, full_variant_key, sample_id
                )

                # Check if there are "present" validations with curated boundaries
                # If so, update the chr:start-end display to show curated values
                # Store original coordinates separately for dialog opening
                map_key = (full_variant_key, sample_id)
                if map_key in validation_map:
                    validations = validation_map[map_key]
                    # Find present validations with curated boundaries (not ignored)
                    present_with_curated = [
                        v
                        for v in validations
                        if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
                    ]
                    if present_with_curated:
                        # Sort by timestamp (most recent first)
                        present_with_curated.sort(key=lambda v: v[6], reverse=True)
                        most_recent = present_with_curated[0]
                        curated_start = most_recent[4]
                        curated_end = most_recent[5]

                        # Parse original chr:start-end
                        parts = variant_key.split(":")
                        if len(parts) == 2:
                            chrom = parts[0]
                            range_parts = parts[1].split("-")
                            if len(range_parts) == 2:
                                orig_start = range_parts[0]
                                orig_end = range_parts[1]

                                # Store original coordinates for dialog opening
                                row["OriginalLocus"] = variant_key

                                # Use curated values if provided, otherwise keep original
                                new_start = (
                                    curated_start if curated_start else orig_start
                                )
                                new_end = curated_end if curated_end else orig_end

                                # Update the display value and mark as curated
                                row["chr:start-end"] = f"{chrom}:{new_start}-{new_end}"
                                row["IsCurated"] = True
                                row["_curated_tooltip"] = (
                                    f"Original: {variant_key}\n"
                                    f"Curated: {row['chr:start-end']}"
                                )
                                # Recompute svlen from curated boundaries
                                try:
                                    row["svlen"] = int(new_end) - int(new_start)
                                except (ValueError, TypeError):
                                    pass

            # Add gene badge information for all rows
            gene_scorer = get_gene_scorer()
            for row in all_rows:
                # Process main gene column
                gene_str = row.get("gene", "")
                if gene_str and gene_str != "-":
                    # Parse gene string format: "SYMBOL:type,SYMBOL2:type"
                    gene_badges = []
                    for gene_part in str(gene_str).split(","):
                        if ":" in gene_part:
                            symbol = gene_part.split(":")[0].strip()
                            gene_type = gene_part.split(":")[1].strip()
                        else:
                            symbol = gene_part.strip()
                            gene_type = ""

                        if symbol:
                            score, _ = gene_scorer.get_gene_score_and_sets(symbol)
                            color = gene_scorer.get_gene_color(symbol)
                            tooltip = gene_scorer.get_gene_tooltip(symbol)
                            is_exonic = gene_type.strip() == "exonic"
                            gene_badges.append(
                                {
                                    "label": symbol,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "type": gene_type,
                                    "score": score,
                                    "borderColor": "black" if is_exonic else None,
                                }
                            )

                    # Sort by score (descending)
                    gene_badges.sort(key=lambda x: x["score"], reverse=True)

                    # Limit to first 6 genes and add "+X genes" indicator if needed
                    total_genes = len(gene_badges)
                    if total_genes > 6:
                        gene_badges = gene_badges[:6]
                        # Add a "+X genes" badge
                        remaining_count = total_genes - 6
                        gene_badges.append(
                            {
                                "label": f"+{remaining_count} genes",
                                "color": "#9e9e9e",  # grey color
                                "tooltip": f"{remaining_count} more genes",
                                "type": "",
                            }
                        )

                    row["GeneBadges"] = gene_badges
                else:
                    row["GeneBadges"] = []

                # Process genic_symbol, exonic_symbol columns
                for col_name in ["genic_symbol", "exonic_symbol"]:
                    col_value = row.get(col_name, "")
                    if col_value and col_value != "-":
                        badges = []
                        symbols = [
                            s.strip() for s in str(col_value).split(",") if s.strip()
                        ]
                        for symbol in symbols:
                            color = gene_scorer.get_gene_color(symbol)
                            tooltip = gene_scorer.get_gene_tooltip(symbol)
                            is_exonic = col_name == "exonic_symbol"
                            badges.append(
                                {
                                    "label": symbol,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "isExonic": is_exonic,
                                    "borderColor": "black" if is_exonic else None,
                                }
                            )
                        row[f"{col_name}_badges"] = badges
                    else:
                        row[f"{col_name}_badges"] = []

                # Process genic_ensg, exonic_ensg, VEP_Gene columns (ENSG IDs)
                for col_name in ["genic_ensg", "exonic_ensg", "VEP_Gene"]:
                    col_value = row.get(col_name, "")
                    if col_value and col_value != "-":
                        badges = []
                        ensgs = [
                            e.strip() for e in str(col_value).split(",") if e.strip()
                        ]
                        for ensg in ensgs:
                            color = gene_scorer.get_gene_color(ensg)
                            tooltip = gene_scorer.get_gene_tooltip(ensg)
                            is_exonic = col_name == "exonic_ensg"
                            badges.append(
                                {
                                    "label": ensg,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "isExonic": is_exonic,
                                    "borderColor": "black" if is_exonic else None,
                                }
                            )
                        row[f"{col_name}_badges"] = badges
                    else:
                        row[f"{col_name}_badges"] = []

        # Initial load of validations
        reload_validations()

        # Get all columns (add Validation column), group same-group columns together
        all_columns = reorder_columns_by_group(
            list(df.columns) + ["Validation", "Diagnostic"]
        )

        # All columns visible by default except gene ID and symbol columns and type
        unchecked_columns = {
            "genic_ensg",
            "exonic_ensg",
            "genic_symbol",
            "exonic_symbol",
            "type",
        }
        selected_cols = {
            "value": [col for col in all_columns if col not in unchecked_columns]
        }

        # Define all possible call values
        all_call_values = [
            "Robust LOSS",
            "Robust GAIN",
            "Permissive LOSS",
            "Permissive Gain",
            "Below threshold",
        ]
        # Default: all selected except "Below threshold"
        selected_calls = {
            "value": [call for call in all_call_values if call != "Below threshold"]
        }

        # Create a container for the data table
        data_container = ui.column().classes("w-full")

        # Capture the client context for use in callbacks
        from nicegui import context

        page_client = context.client

        with data_container:

            @ui.refreshable
            def render_data_table():
                # Filter rows by selected members if 'sample' column exists
                if "sample" in df.columns:
                    rows = [
                        r
                        for r in all_rows
                        if r.get("sample") in selected_members["value"]
                    ]
                else:
                    rows = all_rows

                # Filter rows by selected call values if 'call' column exists
                if "call" in df.columns:
                    rows = [r for r in rows if r.get("call") in selected_calls["value"]]

                ratio_thresholds = build_color_thresholds("ratio")
                zscore_thresholds = build_color_thresholds("zscore")
                call_colors = build_call_colors()

                def make_columns(visible_cols):
                    cols = [
                        {
                            "id": "actions",
                            "header": "",
                            "cellType": "action",
                            "actionName": "view_sv",
                            "actionIcon": "visibility",
                            "actionColor": "#1976d2",
                            "actionTooltip": "View in IGV",
                            "sortable": False,
                        }
                    ]
                    for col in visible_cols:
                        col_def: Dict[str, Any] = {
                            "id": col,
                            "header": get_display_label(col),
                            "group": get_column_group(col),
                            "sorting": get_column_sorting(col),
                            "sortable": True,
                        }
                        if col == "chr:start-end":
                            col_def["cellType"] = "curated_locus"
                            col_def["curatedField"] = "IsCurated"
                            col_def["tooltipField"] = "_curated_tooltip"
                        elif col == "Validation":
                            col_def["cellType"] = "validation"
                        elif col == "Diagnostic":
                            col_def["cellType"] = "diagnostic"
                        elif col == "call":
                            col_def["cellType"] = "cnv_call"
                            col_def["callColors"] = call_colors
                        elif col == "ratio":
                            col_def["cellType"] = "color_scale"
                            col_def["thresholds"] = ratio_thresholds
                        elif col == "zscore":
                            col_def["cellType"] = "color_scale"
                            col_def["thresholds"] = zscore_thresholds
                        elif col == "gene":
                            col_def["cellType"] = "gene_badge"
                            col_def["badgesField"] = "GeneBadges"
                        elif col in _GENE_BADGE_COLUMNS:
                            col_def["cellType"] = "gene_badge"
                            col_def["badgesField"] = f"{col}_badges"
                        else:
                            col_type = get_column_type(col)
                            if col_type in ("int", "float"):
                                col_def["cellType"] = "number"
                        apply_width_constraints(col_def, col)
                        cols.append(col_def)
                    return cols

                with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
                    ui.label(f"Data ({len(rows)} rows)").classes(
                        "text-lg font-semibold text-blue-700"
                    )

                    ui.space()  # Push buttons to the right

                    # Column selector
                    with ui.button("Columns", icon="view_column").props(
                        "outline color=blue size=sm"
                    ):
                        with ui.menu():
                            ui.label("Show/Hide Columns:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2"):
                                    checkboxes: Dict[str, Any] = {}

                                    def select_all():
                                        selected_cols["value"] = list(all_columns)
                                        render_data_table.refresh()

                                    def select_none():
                                        selected_cols["value"] = []
                                        render_data_table.refresh()

                                    ui.button("All", on_click=select_all).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button("None", on_click=select_none).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")

                                ui.separator()

                                for col in all_columns:
                                    checkboxes[col] = ui.checkbox(
                                        col,
                                        value=col in selected_cols["value"],
                                        on_change=lambda e, c=col: handle_col_change(
                                            c, e.value
                                        ),
                                    ).classes("text-sm")

                    # Call filter
                    with ui.button("Filter Call", icon="filter_list").props(
                        "outline color=blue size=sm"
                    ):
                        with ui.menu():
                            ui.label("Filter by Call:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2"):
                                    call_checkboxes: Dict[str, Any] = {}

                                    def select_all_calls():
                                        selected_calls["value"] = list(all_call_values)
                                        update_call_filter()

                                    def select_none_calls():
                                        selected_calls["value"] = []
                                        update_call_filter()

                                    ui.button("All", on_click=select_all_calls).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button("None", on_click=select_none_calls).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")

                                ui.separator()

                                for call_value in all_call_values:
                                    call_checkboxes[call_value] = ui.checkbox(
                                        call_value,
                                        value=call_value in selected_calls["value"],
                                        on_change=lambda e,
                                        c=call_value: handle_call_change(c, e.value),
                                    ).classes("text-sm")

                    # Stats button
                    def show_stats_dialog(current_rows=rows):
                        from collections import Counter

                        with (
                            ui.dialog().props("full-width") as stats_dialog,
                            ui.card().classes("w-full"),
                        ):
                            with ui.column().classes("w-full p-4"):
                                # Header
                                with ui.row().classes(
                                    "items-center justify-between w-full mb-2"
                                ):
                                    with ui.row().classes("items-center gap-3"):
                                        ui.label("SV Statistics").classes(
                                            "text-xl font-bold text-blue-900"
                                        )
                                        ui.label(f"{len(current_rows)} SVs").classes(
                                            "text-sm text-gray-500"
                                        )
                                    ui.button(
                                        icon="close",
                                        on_click=lambda: stats_dialog.close(),
                                    ).props("flat round")

                                # --- Chromosome distribution (stacked bar) ---
                                chrom_counts: Dict[str, Dict[str, int]] = {
                                    c: {} for c in CHROM_ORDER
                                }
                                for r in current_rows:
                                    locus = r.get("chr:start-end", "")
                                    if ":" in locus:
                                        chrom = locus.split(":")[0].replace("chr", "")
                                        status = r.get("Validation", "") or "TODO"
                                        if chrom in chrom_counts:
                                            chrom_counts[chrom][status] = (
                                                chrom_counts[chrom].get(status, 0) + 1
                                            )

                                all_statuses = sorted(
                                    {s for cc in chrom_counts.values() for s in cc}
                                )

                                ui.label("SVs per Chromosome").classes(
                                    "text-lg font-semibold text-gray-800 mt-2"
                                )

                                stacked_series = [
                                    {
                                        "name": status,
                                        "type": "bar",
                                        "stack": "total",
                                        "data": [
                                            chrom_counts[c].get(status, 0)
                                            for c in CHROM_ORDER
                                        ],
                                        "itemStyle": {
                                            "color": VALIDATION_COLORS.get(
                                                status, "#94a3b8"
                                            )
                                        },
                                    }
                                    for status in all_statuses
                                ]
                                ui.echart(
                                    {
                                        "tooltip": {
                                            "trigger": "axis",
                                            "axisPointer": {"type": "shadow"},
                                        },
                                        "legend": {"data": all_statuses, "top": 0},
                                        "grid": {"top": 30},
                                        "xAxis": {
                                            "type": "category",
                                            "data": list(CHROM_ORDER),
                                            "name": "Chromosome",
                                        },
                                        "yAxis": {"type": "value", "name": "Count"},
                                        "series": stacked_series,
                                    }
                                ).classes("w-full h-64")

                                # --- Pie charts row ---
                                with ui.row().classes("w-full gap-4 flex-wrap mt-4"):
                                    # Call type pie chart
                                    with ui.column().classes("flex-1 min-w-[400px]"):
                                        ui.label("Call Type Distribution").classes(
                                            "text-lg font-semibold text-gray-800"
                                        )
                                        call_counts = Counter(
                                            r.get("call", "N/A") for r in current_rows
                                        )
                                        call_colors_map = build_call_colors()
                                        call_data = [
                                            {
                                                "name": call,
                                                "value": count,
                                                "itemStyle": {
                                                    "color": call_colors_map.get(
                                                        call, "#94a3b8"
                                                    )
                                                },
                                            }
                                            for call, count in call_counts.most_common()
                                        ]
                                        ui.echart(
                                            {
                                                "tooltip": {"trigger": "item"},
                                                "series": [
                                                    {
                                                        "type": "pie",
                                                        "radius": "70%",
                                                        "data": call_data,
                                                        "label": {
                                                            "formatter": "{b}: {c} ({d}%)"
                                                        },
                                                    }
                                                ],
                                            }
                                        ).classes("w-full h-64")

                                    # Validation status pie chart
                                    with ui.column().classes("flex-1 min-w-[400px]"):
                                        ui.label("Validation Status").classes(
                                            "text-lg font-semibold text-gray-800"
                                        )
                                        val_counts = Counter(
                                            r.get("Validation", "") or "TODO"
                                            for r in current_rows
                                        )
                                        val_data = [
                                            {
                                                "name": status,
                                                "value": count,
                                                "itemStyle": {
                                                    "color": VALIDATION_COLORS.get(
                                                        status, "#94a3b8"
                                                    )
                                                },
                                            }
                                            for status, count in val_counts.most_common()
                                        ]
                                        ui.echart(
                                            {
                                                "tooltip": {"trigger": "item"},
                                                "series": [
                                                    {
                                                        "type": "pie",
                                                        "radius": "70%",
                                                        "data": val_data,
                                                        "label": {
                                                            "formatter": "{b}: {c} ({d}%)"
                                                        },
                                                    }
                                                ],
                                            }
                                        ).classes("w-full h-64")

                        stats_dialog.open()

                    ui.button(
                        "Stats",
                        icon="bar_chart",
                        on_click=show_stats_dialog,
                    ).props("outline color=blue size=sm")

                def _on_save_refresh():
                    """Refresh callback for validation/diagnostic saves."""
                    reload_validations()
                    with page_client:
                        for rfn in data_table_refreshers:
                            rfn()

                def on_row_action(e):
                    action = e.get("action", "")
                    row_data = e.get("row", {})
                    locus = row_data.get("chr:start-end", "")
                    sample_id = row_data.get("sample", "")

                    if not locus or not sample_id:
                        ui.notify(
                            "Missing locus or sample information",
                            type="warning",
                        )
                        return

                    if action == "open_diagnostic":
                        # Extract gene symbols from gene column
                        gene_str = row_data.get("gene", "")
                        gene_symbols = []
                        if gene_str and gene_str != "-":
                            for gp in str(gene_str).split(","):
                                sym = gp.split(":")[0].strip()
                                if sym:
                                    gene_symbols.append(sym)
                        gene = ", ".join(gene_symbols) if gene_symbols else ""

                        # Impact from call type
                        sv_call = str(row_data.get("call", "")).upper()
                        if "GAIN" in sv_call:
                            impact = "GAIN"
                        elif "LOSS" in sv_call:
                            impact = "LOSS"
                        else:
                            impact = sv_call

                        # Use original locus + sv_type for variant key
                        # Must match the full_variant_key format used in
                        # reload_validations() and sv_dialog.py
                        orig_locus = row_data.get("_original_locus", locus)
                        sv_type = infer_sv_type(row_data)
                        full_variant_key = f"{orig_locus}:{sv_type}"
                        diag_file = store.data_dir / "diagnostics" / "svs.tsv"
                        show_diagnostic_dialog(
                            family_id=family_id,
                            variant_key=full_variant_key,
                            gene=gene,
                            impact=impact,
                            sample=sample_id,
                            variant_type="sv",
                            diagnostic_file=diag_file,
                            on_save_callback=lambda s: _on_save_refresh(),
                        )
                    else:
                        # Default: open SV validation dialog
                        try:
                            locus_to_parse = row_data.get("OriginalLocus", locus)
                            parts = locus_to_parse.split(":")
                            if len(parts) == 2:
                                chrom = parts[0]
                                range_parts = parts[1].split("-")
                                if len(range_parts) == 2:
                                    start = range_parts[0]
                                    end = range_parts[1]

                                    show_sv_dialog(
                                        cohort_name=cohort_name,
                                        family_id=family_id,
                                        chrom=chrom,
                                        start=start,
                                        end=end,
                                        sample=sample_id,
                                        sv_data=row_data,
                                        on_validation_saved=_on_save_refresh,
                                    )
                                else:
                                    ui.notify(
                                        "Invalid locus format. Expected chr:start-end",
                                        type="warning",
                                    )
                            else:
                                ui.notify(
                                    "Invalid locus format. Expected chr:start-end",
                                    type="warning",
                                )
                        except Exception as ex:
                            ui.notify(
                                f"Error parsing locus: {ex}",
                                type="warning",
                            )

                with ui.card().classes("w-full"):
                    DataTable(
                        columns=make_columns(selected_cols["value"]),
                        rows=rows,
                        pagination={"rowsPerPage": 10},
                        on_row_action=on_row_action,
                    )

                def handle_col_change(col_name, is_checked):
                    if is_checked and col_name not in selected_cols["value"]:
                        selected_cols["value"].append(col_name)
                    elif not is_checked and col_name in selected_cols["value"]:
                        selected_cols["value"].remove(col_name)

                    # Reorder to match all_columns order
                    selected_cols["value"] = [
                        col for col in all_columns if col in selected_cols["value"]
                    ]

                    render_data_table.refresh()

                def handle_call_change(call_value, is_checked):
                    if is_checked and call_value not in selected_calls["value"]:
                        selected_calls["value"].append(call_value)
                    elif not is_checked and call_value in selected_calls["value"]:
                        selected_calls["value"].remove(call_value)
                    update_call_filter()

                def update_call_filter():
                    render_data_table.refresh()

            data_table_refreshers.append(render_data_table.refresh)
            render_data_table()

    except Exception as e:
        ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
