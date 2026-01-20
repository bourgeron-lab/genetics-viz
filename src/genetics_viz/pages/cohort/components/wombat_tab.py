"""Wombat tab component for family page."""

import re
from typing import Any, Callable, Dict, List

import polars as pl
from nicegui import ui

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog

# Table slot template for wombat data with view button and validation icons
WOMBAT_TABLE_SLOT = r"""
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


def get_wombat_display_label(col: str) -> str:
    """Get display label for wombat column, removing VEP_ prefix and renaming gnomAD columns."""
    if col == "fafmax_faf95_max_genomes":
        return "gnomAD 4.1 WGS"
    elif col == "nhomalt_genomes":
        return "gnomAD 4.1 nhomalt WGS"
    elif col == "VEP_CLIN_SIG":
        return "ClinVar"
    elif col.startswith("VEP_"):
        return col[4:]  # Remove VEP_ prefix
    else:
        return col


def render_wombat_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the Wombat tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    wombat_dir = store.data_dir / "families" / family_id / "wombat"

    if not wombat_dir.exists():
        ui.label(f"No wombat directory found at: {wombat_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Parse wombat TSV files
    pattern = re.compile(
        rf"{re.escape(family_id)}\.rare\.([^.]+)\.annotated\.(.+?)\.tsv$"
    )

    wombat_files = []
    for tsv_file in wombat_dir.glob("*.tsv"):
        match = pattern.match(tsv_file.name)
        if match:
            vep_config = match.group(1)
            wombat_config = match.group(2)
            wombat_files.append(
                {
                    "file_path": tsv_file,
                    "vep_config": vep_config,
                    "wombat_config": wombat_config,
                }
            )

    if not wombat_files:
        ui.label(
            f"No wombat TSV files found matching pattern in: {wombat_dir}"
        ).classes("text-gray-500 italic")
        return

    # Create dictionaries to store data for each wombat config
    wombat_data: Dict[str, Dict[str, Any]] = {}

    # Create subtabs for each wombat config
    with ui.tabs().classes("w-full") as wombat_subtabs:
        subtab_refs = {}
        for wf in wombat_files:
            subtab_refs[wf["wombat_config"]] = ui.tab(wf["wombat_config"])

    with ui.tab_panels(wombat_subtabs, value=list(subtab_refs.values())[0]).classes(
        "w-full"
    ):
        for wf in wombat_files:
            with ui.tab_panel(subtab_refs[wf["wombat_config"]]):
                config_name = wf["wombat_config"]

                with ui.card().classes("w-full p-4"):
                    ui.label(f"Wombat Configuration: {wf['wombat_config']}").classes(
                        "text-lg font-semibold text-blue-700 mb-2"
                    )
                    with ui.row().classes("gap-4 mb-4"):
                        ui.label("VEP Config:").classes("font-semibold")
                        ui.badge(wf["vep_config"]).props("color=indigo")
                    with ui.row().classes("gap-4"):
                        ui.label("File Path:").classes("font-semibold")
                        ui.label(str(wf["file_path"])).classes(
                            "text-sm text-gray-600 font-mono"
                        )

                # Display TSV content in a table
                try:
                    df = pl.read_csv(
                        wf["file_path"],
                        separator="\t",
                        infer_schema_length=100,
                    )

                    # Convert to list of dicts for NiceGUI table
                    all_rows = df.to_dicts()

                    # Store in the wombat_data dict keyed by config
                    wombat_data[config_name] = {
                        "df": df,
                        "all_rows": all_rows,
                    }

                    # Load validation data from snvs.tsv
                    validation_file = store.data_dir / "validations" / "snvs.tsv"
                    validation_map = load_validation_map(validation_file, family_id)

                    # Add concatenated Variant column and Validation status to each row
                    for row in all_rows:
                        chrom = row.get("#CHROM", "")
                        pos = row.get("POS", "")
                        ref = row.get("REF", "")
                        alt = row.get("ALT", "")
                        sample_id = row.get("sample", "")
                        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                        row["Variant"] = variant_key

                        add_validation_status_to_row(
                            row, validation_map, variant_key, sample_id
                        )

                    # All available columns (add Variant and Validation columns)
                    all_columns = (
                        ["Variant"]
                        + [
                            col
                            for col in df.columns
                            if col not in ["#CHROM", "POS", "REF", "ALT"]
                        ]
                        + ["Validation"]
                    )

                    wombat_data[config_name]["all_columns"] = all_columns

                    # Default columns to display
                    default_visible = [
                        "Variant",
                        "VEP_Consequence",
                        "VEP_SYMBOL",
                        "VEP_CLIN_SIG",
                        "fafmax_faf95_max_genomes",
                        "nhomalt_genomes",
                        "sample",
                        "sample_gt",
                        "father_gt",
                        "mother_gt",
                        "Validation",
                    ]
                    initial_selected = [
                        col for col in default_visible if col in all_columns
                    ]

                    selected_cols = {"value": initial_selected}
                    wombat_data[config_name]["selected_cols"] = selected_cols

                    # Create a container for the data table
                    data_container = ui.column().classes("w-full")

                    # Capture the client context for use in callbacks
                    from nicegui import context

                    page_client = context.client

                    with data_container:

                        @ui.refreshable
                        def render_data_table(cfg=config_name):
                            data = wombat_data[cfg]
                            df_local = data["df"]
                            all_rows_local = data["all_rows"]
                            all_columns_local = data["all_columns"]
                            selected_cols_local = data["selected_cols"]

                            # Filter rows by selected members
                            if "sample" in df_local.columns:
                                rows = [
                                    r
                                    for r in all_rows_local
                                    if r.get("sample") in selected_members["value"]
                                ]
                            else:
                                rows = all_rows_local

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
                                            "label": get_wombat_display_label(col),
                                            "field": col,
                                            "sortable": True,
                                            "align": "left",
                                        }
                                        for col in visible_cols
                                    ]
                                )
                                return cols

                            with ui.row().classes("items-center gap-4 mt-4 mb-2"):
                                ui.label(f"Data ({len(rows)} rows)").classes(
                                    "text-lg font-semibold text-blue-700"
                                )

                                # Column selector
                                with ui.button(
                                    "Select Columns", icon="view_column"
                                ).props("outline color=blue"):
                                    with ui.menu():
                                        ui.label("Show/Hide Columns:").classes(
                                            "px-4 py-2 font-semibold text-sm"
                                        )
                                        ui.separator()

                                        with ui.column().classes("p-2"):
                                            with ui.row().classes("gap-2 mb-2"):
                                                checkboxes: Dict[str, Any] = {}

                                                def select_all():
                                                    selected_cols_local["value"] = list(
                                                        all_columns_local
                                                    )
                                                    update_table()

                                                def select_none():
                                                    selected_cols_local["value"] = []
                                                    update_table()

                                                ui.button(
                                                    "All", on_click=select_all
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )
                                                ui.button(
                                                    "None", on_click=select_none
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )

                                            ui.separator()

                                            for col in all_columns_local:
                                                checkboxes[col] = ui.checkbox(
                                                    col,
                                                    value=col
                                                    in selected_cols_local["value"],
                                                    on_change=lambda e,
                                                    c=col: handle_col_change(
                                                        c, e.value
                                                    ),
                                                ).classes("text-sm")

                            with ui.card().classes("w-full"):
                                data_table = (
                                    ui.table(
                                        columns=make_columns(
                                            selected_cols_local["value"]
                                        ),
                                        rows=rows,
                                        pagination={"rowsPerPage": 10},
                                    )
                                    .classes("w-full")
                                    .props("dense flat")
                                )

                                data_table.add_slot("body", WOMBAT_TABLE_SLOT)

                                def on_view_variant(e):
                                    row_data = e.args
                                    chrom = row_data.get("#CHROM", "")
                                    pos = row_data.get("POS", "")
                                    ref = row_data.get("REF", "")
                                    alt = row_data.get("ALT", "")
                                    sample_val = row_data.get("sample", "")

                                    # Callback to update the Validation column in the table
                                    def on_save(validation_status: str):
                                        # Reload validation data from file
                                        validation_file = (
                                            store.data_dir / "validations" / "snvs.tsv"
                                        )
                                        validation_map = load_validation_map(
                                            validation_file, family_id
                                        )
                                        # Re-add validation status to all rows
                                        for row in all_rows:
                                            chrom = row.get("#CHROM", "")
                                            pos = row.get("POS", "")
                                            ref = row.get("REF", "")
                                            alt = row.get("ALT", "")
                                            sample_id = row.get("sample", "")
                                            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                                            add_validation_status_to_row(
                                                row,
                                                validation_map,
                                                variant_key,
                                                sample_id,
                                            )
                                        # Refresh the table using the captured client context
                                        with page_client:
                                            ui.timer(
                                                0.1,
                                                render_data_table.refresh,
                                                once=True,
                                            )

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

                                data_table.on("view_variant", on_view_variant)

                            def handle_col_change(col_name, is_checked):
                                if (
                                    is_checked
                                    and col_name not in selected_cols_local["value"]
                                ):
                                    selected_cols_local["value"].append(col_name)
                                elif (
                                    not is_checked
                                    and col_name in selected_cols_local["value"]
                                ):
                                    selected_cols_local["value"].remove(col_name)
                                update_table()

                            def update_table():
                                visible = [
                                    c
                                    for c in all_columns_local
                                    if c in selected_cols_local["value"]
                                ]
                                data_table.columns = make_columns(visible)
                                data_table.update()

                                for col, checkbox in checkboxes.items():
                                    checkbox.value = col in selected_cols_local["value"]

                        data_table_refreshers.append(render_data_table.refresh)
                        render_data_table()

                except Exception as e:
                    ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
