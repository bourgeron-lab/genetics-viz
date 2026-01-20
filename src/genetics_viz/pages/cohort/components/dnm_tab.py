"""DNM (de novo mutations) tab component for family page."""

import glob
from pathlib import Path
from typing import Any, Callable, Dict, List, Set

import polars as pl
from nicegui import ui

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog

# Table slot template for DNM data with view button and validation icons
DNM_TABLE_SLOT = r"""
    <q-tr :props="props">
        <q-td key="actions" :props="props">
            <q-btn 
                flat 
                dense 
                size="sm" 
                icon="visibility" 
                color="blue"
                @click="$parent.$emit('view_variant', props.row)"
            >
                <q-tooltip>View in IGV</q-tooltip>
            </q-btn>
        </q-td>
        <q-td v-for="col in props.cols.filter(c => c.name !== 'actions')" :key="col.name" :props="props">
            <template v-if="col.name === 'Validation'">
                <span v-if="col.value === 'present' || col.value === 'in phase MNV'" style="display: flex; align-items: center; gap: 4px;">
                    <q-icon name="check_circle" color="green" size="sm">
                        <q-tooltip>Validated as {{ col.value }}</q-tooltip>
                    </q-icon>
                    <span v-if="props.row.ValidationInheritance === 'de novo'" style="font-weight: bold;">dnm</span>
                    <span v-else-if="props.row.ValidationInheritance === 'homozygous'" style="font-weight: bold;">hom</span>
                    <span v-if="col.value === 'in phase MNV'" style="font-size: 0.75em; color: #666;">MNV</span>
                </span>
                <q-icon v-else-if="col.value === 'absent'" name="cancel" color="red" size="sm">
                    <q-tooltip>Validated as absent</q-tooltip>
                </q-icon>
                <q-icon v-else-if="col.value === 'uncertain' || col.value === 'different'" name="help" color="orange" size="sm">
                    <q-tooltip>Validation uncertain or different</q-tooltip>
                </q-icon>
                <q-icon v-else-if="col.value === 'conflicting'" name="bolt" color="amber-9" size="sm">
                    <q-tooltip>Conflicting validations</q-tooltip>
                </q-icon>
            </template>
            <template v-else>
                {{ col.value }}
            </template>
        </q-td>
    </q-tr>
"""

# LoF and coding variant categories
LOF_BASE_NAMES: Set[str] = {
    "transcript_ablation",
    "splice_acceptor",
    "splice_donor",
    "stop_gained",
    "frameshift",
    "stop_lost",
    "start_lost",
}

CODING_BASE_NAMES: Set[str] = {
    "missense",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering",
    "coding_sequence",
    "synonymous",
    "UTR5",
    "UTR3",
    "3_prime_utr",
    "5_prime_utr",
}


def impact_matches_category(impact: str, base_names: Set[str]) -> bool:
    """Check if impact matches any base name, with or without _variant suffix."""
    for base in base_names:
        if impact == base or impact == f"{base}_variant":
            return True
    return False


def get_dnm_display_label(col: str) -> str:
    """Get display label for DNM column."""
    if col == "chr:pos:ref:alt":
        return "Variant"
    elif col.startswith("VEP_"):
        if col == "VEP_CLIN_SIG":
            return "ClinVar"
        return col[4:]
    elif col == "fafmax_faf95_max_genomes":
        return "gnomAD 4.1 WGS"
    elif col == "nhomalt_genomes":
        return "gnomAD 4.1 nhomalt WGS"
    elif col == "genomes_filters":
        return "gnomAD 4.1 WGS filter"
    return col


def render_dnm_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the SNVs dnm tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    vcfs_dir = store.data_dir / "families" / family_id / "vcfs"

    if not vcfs_dir.exists():
        ui.label(f"No vcfs directory found at: {vcfs_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Find all matching DNM files
    pattern_str = str(vcfs_dir / f"{family_id}.rare.*.annotated.dnm.tsv")
    dnm_files = glob.glob(pattern_str)

    if len(dnm_files) == 0:
        ui.label(
            f"No DNM files found matching pattern: {family_id}.rare.*.annotated.dnm.tsv"
        ).classes("text-gray-500 italic")
        return

    if len(dnm_files) > 1:
        ui.label(
            "Multiple DNM files found. Please ensure only one file matches the pattern:"
        ).classes("text-orange-600 font-semibold mb-2")
        for f in dnm_files:
            ui.label(f"• {Path(f).name}").classes("text-sm text-gray-600 ml-4")
        return

    # Load the single DNM file
    dnm_file = Path(dnm_files[0])

    try:
        df = pl.read_csv(
            dnm_file,
            separator="\t",
            infer_schema_length=100,
        )

        all_rows = df.to_dicts()

        # Ensure required columns exist
        if "chr:pos:ref:alt" not in df.columns:
            ui.label("Error: 'chr:pos:ref:alt' column not found in DNM file").classes(
                "text-red-500"
            )
            return

        if "sample_id" not in df.columns:
            ui.label("Error: 'sample_id' column not found in DNM file").classes(
                "text-red-500"
            )
            return

        # Parse variant key into components and extract impact priority
        for row in all_rows:
            variant_key = row.get("chr:pos:ref:alt", "")
            parts = variant_key.split(":")
            if len(parts) == 4:
                row["#CHROM"] = parts[0]
                row["POS"] = parts[1]
                row["REF"] = parts[2]
                row["ALT"] = parts[3]
            else:
                row["#CHROM"] = ""
                row["POS"] = ""
                row["REF"] = ""
                row["ALT"] = ""

            # Extract impact priority from highest_impact
            if "highest_impact" in row and row["highest_impact"]:
                impact_val = str(row["highest_impact"])
                if len(impact_val) > 3 and impact_val[2] == "_":
                    row["Impact_priority"] = impact_val[:2]
                    row["highest_impact"] = impact_val[3:]
                else:
                    row["Impact_priority"] = "99"
            else:
                row["Impact_priority"] = "99"

        # Sort by Impact_priority
        all_rows.sort(key=lambda r: r.get("Impact_priority", "99"))

        # Load validation data
        validation_file = store.data_dir / "validations" / "snvs.tsv"
        validation_map = load_validation_map(validation_file, family_id)

        # Add validation status
        for row in all_rows:
            variant_key = row.get("chr:pos:ref:alt", "")
            sample_id = row.get("sample_id", "")
            add_validation_status_to_row(row, validation_map, variant_key, sample_id)

        # All available columns
        all_columns = (
            ["chr:pos:ref:alt"]
            + [col for col in df.columns if col != "chr:pos:ref:alt"]
            + ["Impact_priority", "Validation"]
        )

        # Default columns
        default_visible = [
            "chr:pos:ref:alt",
            "sample_id",
            "genomes_filter",
            "genomes_filters",
            "fafmax_faf95_max_genomes",
            "nhomalt_genomes",
            "LCR",
            "gene",
            "highest_impact",
            "Validation",
        ]
        initial_selected = [col for col in default_visible if col in all_columns]

        selected_cols = {"value": initial_selected}

        # Get all unique highest_impact values
        all_impacts = sorted(
            set(
                str(row.get("highest_impact", ""))
                for row in all_rows
                if row.get("highest_impact")
            )
        )

        # Filter state
        filter_exclude_lcr = {"value": False}
        filter_exclude_gnomad = {"value": False}
        filter_selected_impacts: Dict[str, List[str]] = {"value": list(all_impacts)}

        # Create container for the data table
        data_container = ui.column().classes("w-full")

        with data_container:
            # Filter panel
            with ui.card().classes("w-full mb-4 p-4"):
                ui.label("Filters").classes("text-lg font-semibold text-blue-700 mb-2")
                with ui.row().classes("gap-4 items-start"):
                    ui.checkbox(
                        "Exclude LCR",
                        value=filter_exclude_lcr["value"],
                        on_change=lambda e: (
                            filter_exclude_lcr.update({"value": e.value}),
                            render_dnm_table.refresh(),
                        ),
                    )

                    ui.checkbox(
                        "Exclude gnomAD 4.1 WGS filter",
                        value=filter_exclude_gnomad["value"],
                        on_change=lambda e: (
                            filter_exclude_gnomad.update({"value": e.value}),
                            render_dnm_table.refresh(),
                        ),
                    )

                    # Impact filter menu
                    with ui.button("Filter by Impact", icon="filter_list").props(
                        "outline color=blue"
                    ):
                        with ui.menu():
                            ui.label("Select Impact Types:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2 flex-wrap"):
                                    impact_checkboxes: Dict[str, Any] = {}

                                    def select_all_impacts():
                                        filter_selected_impacts["value"] = list(
                                            all_impacts
                                        )
                                        for cb in impact_checkboxes.values():
                                            cb.value = True
                                        render_dnm_table.refresh()

                                    def select_none_impacts():
                                        filter_selected_impacts["value"] = []
                                        for cb in impact_checkboxes.values():
                                            cb.value = False
                                        render_dnm_table.refresh()

                                    def select_lof():
                                        selected = [
                                            i
                                            for i in all_impacts
                                            if impact_matches_category(
                                                i, LOF_BASE_NAMES
                                            )
                                        ]
                                        filter_selected_impacts["value"] = selected
                                        for impact, cb in impact_checkboxes.items():
                                            cb.value = impact in selected
                                        render_dnm_table.refresh()

                                    def select_coding():
                                        selected = [
                                            i
                                            for i in all_impacts
                                            if impact_matches_category(
                                                i, CODING_BASE_NAMES
                                            )
                                            or impact_matches_category(
                                                i, LOF_BASE_NAMES
                                            )
                                        ]
                                        filter_selected_impacts["value"] = selected
                                        for impact, cb in impact_checkboxes.items():
                                            cb.value = impact in selected
                                        render_dnm_table.refresh()

                                    ui.button("All", on_click=select_all_impacts).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button(
                                        "None", on_click=select_none_impacts
                                    ).props("size=sm flat dense").classes("text-xs")
                                    ui.button("LoF", on_click=select_lof).props(
                                        "size=sm flat dense color=orange"
                                    ).classes("text-xs")
                                    ui.button("Coding", on_click=select_coding).props(
                                        "size=sm flat dense color=purple"
                                    ).classes("text-xs")

                                ui.separator()

                                for impact in all_impacts:

                                    def handle_impact_change(imp, val):
                                        if val:
                                            if (
                                                imp
                                                not in filter_selected_impacts["value"]
                                            ):
                                                filter_selected_impacts["value"].append(
                                                    imp
                                                )
                                        else:
                                            if imp in filter_selected_impacts["value"]:
                                                filter_selected_impacts["value"].remove(
                                                    imp
                                                )
                                        render_dnm_table.refresh()

                                    impact_checkboxes[impact] = ui.checkbox(
                                        impact,
                                        value=impact
                                        in filter_selected_impacts["value"],
                                        on_change=lambda e,
                                        i=impact: handle_impact_change(i, e.value),
                                    ).classes("text-sm")

            # Capture the client context for use in callbacks
            from nicegui import context

            page_client = context.client

            @ui.refreshable
            def render_dnm_table():
                # Filter rows by selected members
                rows = [
                    r
                    for r in all_rows
                    if r.get("sample_id") in selected_members["value"]
                ]

                total_before_filters = len(rows)

                # Apply filters
                if filter_exclude_lcr["value"]:
                    rows = [r for r in rows if r.get("LCR") != "LCR"]

                if filter_exclude_gnomad["value"]:
                    rows = [r for r in rows if not r.get("genomes_filters")]

                if filter_selected_impacts["value"]:
                    rows = [
                        r
                        for r in rows
                        if str(r.get("highest_impact", ""))
                        in filter_selected_impacts["value"]
                    ]

                filtered_out = total_before_filters - len(rows)

                def make_columns(visible_cols):
                    cols = [
                        {
                            "name": "actions",
                            "label": "",
                            "field": "actions",
                            "sortable": False,
                            "align": "center",
                        }
                    ]
                    cols.extend(
                        [
                            {
                                "name": col,
                                "label": get_dnm_display_label(col),
                                "field": col,
                                "sortable": True,
                                "align": "left",
                            }
                            for col in visible_cols
                        ]
                    )
                    return cols

                with ui.row().classes("items-center gap-4 mt-4 mb-2"):
                    label_text = f"Data ({len(rows)} rows"
                    if filtered_out > 0:
                        label_text += f", {filtered_out} filtered out)"
                    else:
                        label_text += ")"
                    ui.label(label_text).classes("text-lg font-semibold text-blue-700")

                    # Column selector
                    with ui.button("Select Columns", icon="view_column").props(
                        "outline color=blue"
                    ):
                        with ui.menu():
                            ui.label("Show/Hide Columns:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2"):
                                    col_checkboxes: Dict[str, Any] = {}

                                    def col_select_all():
                                        selected_cols["value"] = list(all_columns)
                                        update_table()

                                    def col_select_none():
                                        selected_cols["value"] = []
                                        update_table()

                                    ui.button("All", on_click=col_select_all).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button("None", on_click=col_select_none).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")

                                ui.separator()

                                for col in all_columns:
                                    col_checkboxes[col] = ui.checkbox(
                                        col,
                                        value=col in selected_cols["value"],
                                        on_change=lambda e, c=col: handle_col_change(
                                            c, e.value
                                        ),
                                    ).classes("text-sm")

                with ui.card().classes("w-full"):
                    dnm_table = (
                        ui.table(
                            columns=make_columns(selected_cols["value"]),
                            rows=rows,
                            pagination={"rowsPerPage": 10},
                        )
                        .classes("w-full")
                        .props("dense flat")
                    )

                    dnm_table.add_slot("body", DNM_TABLE_SLOT)

                    def on_view_dnm_variant(e):
                        row_data = e.args
                        variant_key = row_data.get("chr:pos:ref:alt", "")
                        parts = variant_key.split(":")

                        if len(parts) == 4:
                            chrom, pos, ref, alt = parts
                            sample_val = row_data.get("sample_id", "")

                            # Callback to update the Validation column in the table
                            def on_save(validation_status: str):
                                # Reload validation data from file
                                validation_map = load_validation_map(
                                    validation_file, family_id
                                )
                                # Re-add validation status to all rows
                                for row in all_rows:
                                    variant_key = row.get("chr:pos:ref:alt", "")
                                    sample_id = row.get("sample_id", "")
                                    add_validation_status_to_row(
                                        row, validation_map, variant_key, sample_id
                                    )
                                # Refresh the table using the captured client context
                                with page_client:
                                    ui.timer(0.1, render_dnm_table.refresh, once=True)

                            # Show dialog
                            show_variant_dialog(
                                cohort_name=cohort_name,
                                family_id=family_id,
                                chrom=chrom,
                                pos=pos,
                                ref=ref,
                                alt=alt,
                                sample=sample_val,
                                variant_data=row_data,
                                on_save_callback=on_save,
                            )

                    dnm_table.on("view_variant", on_view_dnm_variant)

                def handle_col_change(col_name, is_checked):
                    if is_checked and col_name not in selected_cols["value"]:
                        selected_cols["value"].append(col_name)
                    elif not is_checked and col_name in selected_cols["value"]:
                        selected_cols["value"].remove(col_name)
                    update_table()

                def update_table():
                    render_dnm_table.refresh()

            render_dnm_table()

            def refresh_dnm_table():
                render_dnm_table.refresh()

            data_table_refreshers.append(refresh_dnm_table)

    except Exception as e:
        ui.label(f"Error reading DNM file: {e}").classes("text-red-500 mt-4")
