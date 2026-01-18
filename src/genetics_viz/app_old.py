"""
NiceGUI web application for genetics-viz.
"""

import os
import re
from pathlib import Path
from typing import Optional

import polars as pl
from nicegui import ui

# Import shared components
from genetics_viz.components.header import create_header
from genetics_viz.models import DataStore

# Import pages to register routes
from genetics_viz.pages import validation  # noqa: F401
from genetics_viz.utils import data as data_module

# Import shared utilities - these are used by the modular pages
from genetics_viz.utils.data import get_data_store

# Legacy: Keep local reference for backward compatibility with pages still in this file
_data_store: Optional[DataStore] = None


@ui.page("/")
def home_page() -> None:
    """Render the home/welcome page."""
    create_header()

    with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
        ui.label("Welcome to Genetics-Viz").classes(
            "text-3xl font-bold mb-2 text-blue-900"
        )
        ui.label("Select a cohort to explore").classes("text-lg text-gray-600 mb-6")

        try:
            store = get_data_store()

            if not store.cohorts:
                with ui.card().classes("w-full p-6 bg-yellow-50"):
                    ui.label("⚠️ No cohorts found").classes(
                        "text-xl font-semibold text-yellow-800"
                    )
                    ui.label(
                        f"No valid cohorts were found in: {store.cohorts_dir}"
                    ).classes("text-gray-600")
                    ui.label(
                        "Make sure each cohort directory contains a .pedigree.tsv file."
                    ).classes("text-gray-500 text-sm")
                return

            # Display cohorts as cards
            ui.label("Available Cohorts").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.row().classes("w-full flex-wrap gap-4"):
                for cohort in sorted(store.cohorts.values(), key=lambda c: c.name):
                    cohort_name = cohort.name  # Capture name for lambda
                    with (
                        ui.card()
                        .classes(
                            "cursor-pointer hover:shadow-lg transition-shadow w-80 border-l-4 border-blue-500"
                        )
                        .on(
                            "click",
                            lambda _, n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                        )
                    ):
                        with ui.card_section():
                            ui.label(cohort.name).classes(
                                "text-xl font-bold text-blue-700"
                            )

                        with ui.card_section():
                            with ui.row().classes("gap-6"):
                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_families)).classes(
                                        "text-3xl font-bold text-blue-600"
                                    )
                                    ui.label("Families").classes(
                                        "text-sm text-gray-500"
                                    )

                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_samples)).classes(
                                        "text-3xl font-bold text-green-600"
                                    )
                                    ui.label("Samples").classes("text-sm text-gray-500")

                        with ui.card_section().classes("bg-gray-50"):
                            ui.label(f"📄 {cohort.pedigree_file.name}").classes(
                                "text-xs text-gray-400"
                            )

        except RuntimeError as e:
            ui.label(f"Error: {e}").classes("text-red-500")


@ui.page("/cohort/{cohort_name}")
def cohort_page(cohort_name: str) -> None:
    """Render the cohort detail page."""
    create_header()

    try:
        store = get_data_store()
        cohort = store.get_cohort(cohort_name)

        if cohort is None:
            with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
                ui.label(f"Cohort not found: {cohort_name}").classes(
                    "text-xl text-red-500"
                )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.label(cohort_name).classes("font-semibold")

            # Cohort header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")

            # Main content: two-panel layout
            with ui.row().classes("w-full gap-6 flex flex-row"):
                # Left panel: Families table
                with ui.column().classes("flex-1 min-w-0"):
                    ui.label("Families").classes(
                        "text-xl font-semibold mb-2 text-blue-800"
                    )

                    families_data = cohort.get_families_summary()

                    # Create table with selection and pagination
                    families_table = ui.table(
                        columns=[
                            {
                                "name": "family_id",
                                "label": "Family ID",
                                "field": "Family ID",
                                "sortable": True,
                            },
                            {
                                "name": "members",
                                "label": "Members",
                                "field": "Members",
                                "sortable": True,
                            },
                        ],
                        rows=families_data,
                        row_key="Family ID",
                        selection="single",
                        pagination={"rowsPerPage": 10},
                    ).classes("w-full")
                    families_table.props("dense")

                    # Add custom slot for Family ID column to make it a link
                    families_table.add_slot(
                        "body-cell-family_id",
                        r"""
                            <q-td :props="props">
                                <a :href="'/cohort/"""
                        + cohort_name
                        + r"""/family/' + props.row['Family ID']" 
                                   class="text-blue-600 hover:text-blue-800 underline cursor-pointer">
                                    {{ props.row['Family ID'] }}
                                </a>
                            </q-td>
                        """,
                    )

                # Right panel: Family members (shown when family selected)
                with ui.column().classes("flex-1 min-w-0"):
                    members_label = ui.label("Select a family to view members").classes(
                        "text-xl font-semibold mb-2 text-gray-400"
                    )
                    members_container = ui.column().classes("w-full")

                def on_family_select(e) -> None:
                    """Handle family selection."""
                    members_container.clear()

                    # Access selection from the table's selected property
                    selection = families_table.selected

                    if not selection:
                        members_label.text = "Select a family to view members"
                        members_label.classes(
                            remove="text-blue-800", add="text-gray-400"
                        )
                        return

                    selected_row = selection[0]
                    family_id = selected_row.get("Family ID")

                    members_label.text = f"Members of Family: {family_id}"
                    members_label.classes(remove="text-gray-400", add="text-blue-800")

                    members_data = cohort.get_family_members(family_id)

                    with members_container:
                        ui.table(
                            columns=[
                                {
                                    "name": "sample_id",
                                    "label": "Sample ID",
                                    "field": "Sample ID",
                                    "sortable": True,
                                },
                                {
                                    "name": "father",
                                    "label": "Father",
                                    "field": "Father",
                                },
                                {
                                    "name": "mother",
                                    "label": "Mother",
                                    "field": "Mother",
                                },
                                {"name": "sex", "label": "Sex", "field": "Sex"},
                                {
                                    "name": "phenotype",
                                    "label": "Phenotype",
                                    "field": "Phenotype",
                                },
                            ],
                            rows=members_data,
                            row_key="Sample ID",
                        ).classes("w-full").props("dense")

                families_table.on("selection", on_family_select)

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")


@ui.page("/cohort/{cohort_name}/family/{family_id}")
def family_page(cohort_name: str, family_id: str) -> None:
    """Render the family detail page."""
    create_header()

    try:
        store = get_data_store()
        cohort = store.get_cohort(cohort_name)

        if cohort is None:
            with ui.column().classes("w-full px-6 py-6"):
                ui.label(f"Cohort not found: {cohort_name}").classes(
                    "text-xl text-red-500"
                )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        family = cohort.families.get(family_id)
        if family is None:
            with ui.column().classes("w-full px-6 py-6"):
                ui.label(f"Family not found: {family_id}").classes(
                    "text-xl text-red-500"
                )
                ui.button(
                    "← Back to Cohort",
                    on_click=lambda: ui.navigate.to(f"/cohort/{cohort_name}"),
                )
            return

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.link(cohort_name, f"/cohort/{cohort_name}").classes(
                    "text-blue-600 hover:text-blue-800"
                )
                ui.label("/").classes("text-gray-400")
                ui.label(family_id).classes("font-semibold")

            # Family header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"👨‍👩‍👧‍👦 Family: {family_id}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{family.num_samples} members").props("color=blue")
                ui.badge(f"{family.num_founders} founders").props("color=teal")

            # Family members table
            # ui.label("Family Members").classes(
            #     "text-2xl font-semibold mb-4 text-blue-800"
            # )

            members_data = cohort.get_family_members(family_id)

            # Track selected members for filtering (default: all selected)
            selected_members = {"value": [m["Sample ID"] for m in members_data]}
            member_checkboxes = {}

            # Store refresh functions for all data tables
            data_table_refreshers = []

            with ui.card().classes("w-full"):
                # Member selection checkboxes
                with ui.column().classes("p-4 bg-blue-50"):
                    with ui.row().classes("items-center gap-2 mb-2"):
                        ui.label("Select Members to Display:").classes(
                            "font-semibold text-blue-800"
                        )

                        def select_all_members():
                            selected_members["value"] = [
                                m["Sample ID"] for m in members_data
                            ]
                            for cb in member_checkboxes.values():
                                cb.value = True
                            # Refresh all data tables
                            for refresher in data_table_refreshers:
                                refresher()

                        def select_none_members():
                            selected_members["value"] = []
                            for cb in member_checkboxes.values():
                                cb.value = False
                            # Refresh all data tables
                            for refresher in data_table_refreshers:
                                refresher()

                        ui.button("All", on_click=select_all_members).props(
                            "size=sm flat dense"
                        ).classes("text-xs")
                        ui.button("None", on_click=select_none_members).props(
                            "size=sm flat dense"
                        ).classes("text-xs")

                    # Create HTML table with checkboxes and member info
                    table_html = """
                    <table class="w-full text-sm">
                        <thead class="bg-blue-100">
                            <tr>
                                <th class="px-3 py-2 text-left font-semibold">Select</th>
                                <th class="px-3 py-2 text-left font-semibold"></th>
                                <th class="px-3 py-2 text-left font-semibold">Sample ID</th>
                                <th class="px-3 py-2 text-left font-semibold">Father</th>
                                <th class="px-3 py-2 text-left font-semibold">Mother</th>
                                <th class="px-3 py-2 text-left font-semibold">Sex</th>
                                <th class="px-3 py-2 text-left font-semibold">Phenotype</th>
                            </tr>
                        </thead>
                        <tbody>
                    """

                    for idx, member in enumerate(members_data):
                        sample_id = member["Sample ID"]
                        bg_class = "bg-white" if idx % 2 == 0 else "bg-gray-50"
                        table_html += f'''
                            <tr class="{bg_class} border-b border-gray-200">
                                <td class="px-3 py-2" id="checkbox-cell-{idx}"></td>
                                <td class="px-3 py-2" id="only-button-cell-{idx}"></td>
                                <td class="px-3 py-2 font-medium">{sample_id}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Father", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Mother", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Sex", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Phenotype", "-")}</td>
                            </tr>
                        '''

                    table_html += """
                        </tbody>
                    </table>
                    """

                    ui.html(table_html, sanitize=False)

                    # Create checkboxes and insert them into the table cells
                    for idx, member in enumerate(members_data):
                        sample_id = member["Sample ID"]

                        def make_change_handler(sid):
                            def handler(e):
                                if e.value and sid not in selected_members["value"]:
                                    selected_members["value"].append(sid)
                                elif not e.value and sid in selected_members["value"]:
                                    selected_members["value"].remove(sid)
                                # Refresh all data tables
                                for refresher in data_table_refreshers:
                                    refresher()

                            return handler

                        with ui.element().classes(f"checkbox-cell-{idx}"):
                            member_checkboxes[sample_id] = ui.checkbox(
                                "",
                                value=True,
                                on_change=make_change_handler(sample_id),
                            )

                        # Create "Only" button for this member
                        def make_only_handler(sid):
                            def handler():
                                # Unselect all members first
                                selected_members["value"] = [sid]
                                # Update all checkboxes
                                for s_id, checkbox in member_checkboxes.items():
                                    checkbox.value = s_id == sid
                                # Refresh all data tables
                                for refresher in data_table_refreshers:
                                    refresher()

                            return handler

                        with ui.element().classes(f"only-button-cell-{idx}"):
                            ui.button(
                                "only", on_click=make_only_handler(sample_id)
                            ).props("size=xs flat dense color=blue").classes("text-xs")

                    # Move checkboxes and only buttons into table cells using JavaScript
                    ui.run_javascript(f"""
                        for (let i = 0; i < {len(members_data)}; i++) {{
                            const checkbox = document.querySelector('.checkbox-cell-' + i);
                            const checkboxCell = document.getElementById('checkbox-cell-' + i);
                            if (checkbox && checkboxCell) {{
                                checkboxCell.appendChild(checkbox);
                            }}
                            
                            const onlyButton = document.querySelector('.only-button-cell-' + i);
                            const onlyButtonCell = document.getElementById('only-button-cell-' + i);
                            if (onlyButton && onlyButtonCell) {{
                                onlyButtonCell.appendChild(onlyButton);
                            }}
                        }}
                    """)

            # Analysis tabs section
            # ui.label("Analysis").classes(
            #     "text-2xl font-semibold mb-4 mt-8 text-blue-800"
            # )

            with ui.tabs().classes("w-full") as tabs:
                wombat_tab = ui.tab("Wombat")
                snvs_tab = ui.tab("SNVs dnm")

            with ui.tab_panels(tabs, value=wombat_tab).classes("w-full"):
                # Wombat tab panel
                with ui.tab_panel(wombat_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):
                    # Look for wombat files
                    wombat_dir = store.data_dir / "families" / family_id / "wombat"

                    if not wombat_dir.exists():
                        ui.label(f"No wombat directory found at: {wombat_dir}").classes(
                            "text-gray-500 italic"
                        )
                    else:
                        # Parse wombat TSV files
                        # Pattern: {family_id}.rare.{vep_config}.annotated.{wombat_config}.tsv
                        # family_id can contain: -, _, .
                        # vep_config can contain: -, _
                        # wombat_config can contain: -, _, .
                        # Exclude: {family_id}.rare.{vep_config}.annotated.tsv.gz
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
                        else:
                            # Create dictionaries to store data for each wombat config
                            wombat_data = {}  # config_name -> {"df": df, "rows": rows, "columns": columns, etc}

                            # Create subtabs for each wombat config
                            with ui.tabs().classes("w-full") as wombat_subtabs:
                                subtab_refs = {}
                                for wf in wombat_files:
                                    subtab_refs[wf["wombat_config"]] = ui.tab(
                                        wf["wombat_config"]
                                    )

                            with ui.tab_panels(
                                wombat_subtabs, value=list(subtab_refs.values())[0]
                            ).classes("w-full"):
                                for wf in wombat_files:
                                    with ui.tab_panel(subtab_refs[wf["wombat_config"]]):
                                        config_name = wf[
                                            "wombat_config"
                                        ]  # Capture for this iteration

                                        with ui.card().classes("w-full p-4"):
                                            ui.label(
                                                f"Wombat Configuration: {wf['wombat_config']}"
                                            ).classes(
                                                "text-lg font-semibold text-blue-700 mb-2"
                                            )
                                            with ui.row().classes("gap-4 mb-4"):
                                                ui.label("VEP Config:").classes(
                                                    "font-semibold"
                                                )
                                                ui.badge(wf["vep_config"]).props(
                                                    "color=indigo"
                                                )
                                            with ui.row().classes("gap-4"):
                                                ui.label("File Path:").classes(
                                                    "font-semibold"
                                                )
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
                                            validation_file = (
                                                store.data_dir
                                                / "validations"
                                                / "snvs.tsv"
                                            )
                                            validation_map = {}  # (variant_key, sample_id) -> list of (validation_status, inheritance)

                                            if validation_file.exists():
                                                import csv

                                                with open(validation_file, "r") as f:
                                                    reader = csv.DictReader(
                                                        f, delimiter="\t"
                                                    )
                                                    for row in reader:
                                                        fid = row.get("FID")
                                                        variant_key = row.get("Variant")
                                                        sample_id = row.get("Sample")
                                                        validation_status = row.get(
                                                            "Validation"
                                                        )
                                                        inheritance = row.get(
                                                            "Inheritance"
                                                        )
                                                        # Filter by family_id and check required fields
                                                        if (
                                                            fid == family_id
                                                            and variant_key
                                                            and sample_id
                                                        ):
                                                            map_key = (
                                                                variant_key,
                                                                sample_id,
                                                            )
                                                            if (
                                                                map_key
                                                                not in validation_map
                                                            ):
                                                                validation_map[
                                                                    map_key
                                                                ] = []
                                                            validation_map[
                                                                map_key
                                                            ].append(
                                                                (
                                                                    validation_status,
                                                                    inheritance,
                                                                )
                                                            )

                                            # Add concatenated Variant column and Validation status to each row
                                            for row in all_rows:
                                                chrom = row.get("#CHROM", "")
                                                pos = row.get("POS", "")
                                                ref = row.get("REF", "")
                                                alt = row.get("ALT", "")
                                                sample_id = row.get("sample", "")
                                                variant_key = (
                                                    f"{chrom}:{pos}:{ref}:{alt}"
                                                )
                                                row["Variant"] = variant_key

                                                # Determine validation status
                                                map_key = (variant_key, sample_id)
                                                if map_key in validation_map:
                                                    validations = validation_map[
                                                        map_key
                                                    ]
                                                    validation_statuses = [
                                                        v[0] for v in validations
                                                    ]
                                                    unique_validations = set(
                                                        validation_statuses
                                                    )

                                                    if len(unique_validations) > 1:
                                                        # Conflicting validations
                                                        row["Validation"] = (
                                                            "conflicting"
                                                        )
                                                        row["ValidationInheritance"] = (
                                                            ""
                                                        )
                                                    elif (
                                                        "present" in unique_validations
                                                    ):
                                                        row["Validation"] = "present"
                                                        # Check if any validation is de novo
                                                        is_de_novo = any(
                                                            v[1] == "de novo"
                                                            for v in validations
                                                            if v[0] == "present"
                                                        )
                                                        row["ValidationInheritance"] = (
                                                            "de novo"
                                                            if is_de_novo
                                                            else ""
                                                        )
                                                    elif "absent" in unique_validations:
                                                        row["Validation"] = "absent"
                                                        row["ValidationInheritance"] = (
                                                            ""
                                                        )
                                                    else:
                                                        # uncertain or different
                                                        row["Validation"] = "uncertain"
                                                        row["ValidationInheritance"] = (
                                                            ""
                                                        )
                                                else:
                                                    row["Validation"] = ""
                                                    row["ValidationInheritance"] = ""

                                            # All available columns (add Variant and Validation columns)
                                            all_columns = (
                                                ["Variant"]
                                                + [
                                                    col
                                                    for col in df.columns
                                                    if col
                                                    not in [
                                                        "#CHROM",
                                                        "POS",
                                                        "REF",
                                                        "ALT",
                                                    ]
                                                ]
                                                + ["Validation"]
                                            )

                                            # Store columns in the data dict
                                            wombat_data[config_name]["all_columns"] = (
                                                all_columns
                                            )

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
                                            # Only keep columns that actually exist in the file
                                            initial_selected = [
                                                col
                                                for col in default_visible
                                                if col in all_columns
                                            ]

                                            # Store selected columns
                                            selected_cols = {"value": initial_selected}
                                            wombat_data[config_name][
                                                "selected_cols"
                                            ] = selected_cols

                                            # Create a container for the data table that can be refreshed
                                            data_container = ui.column().classes(
                                                "w-full"
                                            )

                                            with data_container:

                                                @ui.refreshable
                                                def render_data_table(
                                                    cfg=config_name,
                                                ):  # Capture config_name as default arg
                                                    # Get data for this config
                                                    data = wombat_data[cfg]
                                                    df_local = data["df"]
                                                    all_rows_local = data["all_rows"]
                                                    all_columns_local = data[
                                                        "all_columns"
                                                    ]
                                                    selected_cols_local = data[
                                                        "selected_cols"
                                                    ]

                                                    # Filter rows by selected members (if 'sample' column exists)
                                                    if "sample" in df_local.columns:
                                                        rows = [
                                                            r
                                                            for r in all_rows_local
                                                            if r.get("sample")
                                                            in selected_members["value"]
                                                        ]
                                                    else:
                                                        rows = all_rows_local

                                                    # Create columns definition
                                                    def make_columns(visible_cols):
                                                        def get_display_label(col):
                                                            """Get display label for column, removing VEP_ prefix and renaming gnomAD columns."""
                                                            if (
                                                                col
                                                                == "fafmax_faf95_max_genomes"
                                                            ):
                                                                return "gnomAD 4.1 WGS"
                                                            elif (
                                                                col == "nhomalt_genomes"
                                                            ):
                                                                return "gnomAD 4.1 nhomalt WGS"
                                                            elif col == "VEP_CLIN_SIG":
                                                                return "ClinVar"
                                                            elif col.startswith("VEP_"):
                                                                return col[
                                                                    4:
                                                                ]  # Remove VEP_ prefix
                                                            else:
                                                                return col

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
                                                                    "label": get_display_label(
                                                                        col
                                                                    ),
                                                                    "field": col,
                                                                    "sortable": True,
                                                                    "align": "left",
                                                                }
                                                                for col in visible_cols
                                                            ]
                                                        )
                                                        return cols

                                                    with ui.row().classes(
                                                        "items-center gap-4 mt-4 mb-2"
                                                    ):
                                                        ui.label(
                                                            f"Data ({len(rows)} rows)"
                                                        ).classes(
                                                            "text-lg font-semibold text-blue-700"
                                                        )

                                                        # Column selector
                                                        with ui.button(
                                                            "Select Columns",
                                                            icon="view_column",
                                                        ).props("outline color=blue"):
                                                            with ui.menu() as col_menu:
                                                                ui.label(
                                                                    "Show/Hide Columns:"
                                                                ).classes(
                                                                    "px-4 py-2 font-semibold text-sm"
                                                                )
                                                                ui.separator()

                                                                with (
                                                                    ui.column().classes(
                                                                        "p-2"
                                                                    )
                                                                ):
                                                                    # Select All / Deselect All buttons
                                                                    with ui.row().classes(
                                                                        "gap-2 mb-2"
                                                                    ):

                                                                        def select_all():
                                                                            selected_cols_local[
                                                                                "value"
                                                                            ] = list(
                                                                                all_columns_local
                                                                            )
                                                                            update_table()

                                                                        def select_none():
                                                                            selected_cols_local[
                                                                                "value"
                                                                            ] = []
                                                                            update_table()

                                                                        ui.button(
                                                                            "All",
                                                                            on_click=select_all,
                                                                        ).props(
                                                                            "size=sm flat dense"
                                                                        ).classes(
                                                                            "text-xs"
                                                                        )
                                                                        ui.button(
                                                                            "None",
                                                                            on_click=select_none,
                                                                        ).props(
                                                                            "size=sm flat dense"
                                                                        ).classes(
                                                                            "text-xs"
                                                                        )

                                                                    ui.separator()

                                                                    # Checkboxes for each column
                                                                    checkboxes = {}
                                                                    for col in all_columns_local:
                                                                        checkboxes[
                                                                            col
                                                                        ] = ui.checkbox(
                                                                            col,
                                                                            value=col
                                                                            in selected_cols_local[
                                                                                "value"
                                                                            ],
                                                                            on_change=lambda e,
                                                                            c=col: handle_col_change(
                                                                                c,
                                                                                e.value,
                                                                            ),
                                                                        ).classes(
                                                                            "text-sm"
                                                                        )

                                                    with ui.card().classes(
                                                        "w-full"
                                                    ) as table_card:
                                                        data_table = (
                                                            ui.table(
                                                                columns=make_columns(
                                                                    selected_cols_local[
                                                                        "value"
                                                                    ]
                                                                ),
                                                                rows=rows,
                                                                pagination={
                                                                    "rowsPerPage": 10
                                                                },
                                                            )
                                                            .classes("w-full")
                                                            .props("dense flat")
                                                        )

                                                        # Add view button and validation icons to rows
                                                        data_table.add_slot(
                                                            "body",
                                                            r"""
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
                                                                            <span v-if="col.value === 'present'" style="display: flex; align-items: center; gap: 4px;">
                                                                                <q-icon name="check_circle" color="green" size="sm">
                                                                                    <q-tooltip>Validated as present</q-tooltip>
                                                                                </q-icon>
                                                                                <span v-if="props.row.ValidationInheritance === 'de novo'" style="font-weight: bold;">dnm</span>
                                                                            </span>
                                                                            <q-icon v-else-if="col.value === 'absent'" name="cancel" color="red" size="sm">
                                                                                <q-tooltip>Validated as absent</q-tooltip>
                                                                            </q-icon>
                                                                            <q-icon v-else-if="col.value === 'uncertain'" name="help" color="orange" size="sm">
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
                                                            """,
                                                        )

                                                        # Handle view button click
                                                        import json
                                                        import urllib.parse

                                                        def on_view_variant(e):
                                                            row_data = e.args
                                                            # Get from original columns if available
                                                            chrom = row_data.get(
                                                                "#CHROM", ""
                                                            )
                                                            pos = row_data.get(
                                                                "POS", ""
                                                            )
                                                            ref = row_data.get(
                                                                "REF", ""
                                                            )
                                                            alt = row_data.get(
                                                                "ALT", ""
                                                            )
                                                            sample_val = row_data.get(
                                                                "sample", ""
                                                            )

                                                            # Encode all variant data as URL parameter
                                                            data_encoded = (
                                                                urllib.parse.quote(
                                                                    json.dumps(row_data)
                                                                )
                                                            )

                                                            url = f"/variant/{cohort_name}/{family_id}?chrom={chrom}&pos={pos}&ref={ref}&alt={alt}&sample={sample_val}&data={data_encoded}"
                                                            ui.navigate.to(url)

                                                        data_table.on(
                                                            "view_variant",
                                                            on_view_variant,
                                                        )

                                                    def handle_col_change(
                                                        col_name, is_checked
                                                    ):
                                                        if (
                                                            is_checked
                                                            and col_name
                                                            not in selected_cols_local[
                                                                "value"
                                                            ]
                                                        ):
                                                            selected_cols_local[
                                                                "value"
                                                            ].append(col_name)
                                                        elif (
                                                            not is_checked
                                                            and col_name
                                                            in selected_cols_local[
                                                                "value"
                                                            ]
                                                        ):
                                                            selected_cols_local[
                                                                "value"
                                                            ].remove(col_name)
                                                        update_table()

                                                    def update_table():
                                                        # Preserve column order from original
                                                        visible = [
                                                            c
                                                            for c in all_columns_local
                                                            if c
                                                            in selected_cols_local[
                                                                "value"
                                                            ]
                                                        ]
                                                        data_table.columns = (
                                                            make_columns(visible)
                                                        )
                                                        data_table.update()

                                                        # Update checkboxes
                                                        for (
                                                            col,
                                                            checkbox,
                                                        ) in checkboxes.items():
                                                            checkbox.value = (
                                                                col
                                                                in selected_cols_local[
                                                                    "value"
                                                                ]
                                                            )

                                                # Store the refresh function so checkboxes can trigger it
                                                data_table_refreshers.append(
                                                    render_data_table.refresh
                                                )

                                                # Initial call to render
                                                render_data_table()

                                        except Exception as e:
                                            ui.label(
                                                f"Error reading file: {e}"
                                            ).classes("text-red-500 mt-4")

                # SNVs dnm tab panel
                with ui.tab_panel(snvs_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):
                    # Look for DNM files
                    vcfs_dir = store.data_dir / "families" / family_id / "vcfs"

                    if not vcfs_dir.exists():
                        ui.label(f"No vcfs directory found at: {vcfs_dir}").classes(
                            "text-gray-500 italic"
                        )
                    else:
                        # Find all matching DNM files
                        import glob

                        pattern = str(
                            vcfs_dir / f"{family_id}.rare.*.annotated.dnm.tsv"
                        )
                        dnm_files = glob.glob(pattern)

                        if len(dnm_files) == 0:
                            ui.label(
                                f"No DNM files found matching pattern: {family_id}.rare.*.annotated.dnm.tsv"
                            ).classes("text-gray-500 italic")
                        elif len(dnm_files) > 1:
                            ui.label(
                                "Multiple DNM files found. Please ensure only one file matches the pattern:"
                            ).classes("text-orange-600 font-semibold mb-2")
                            for f in dnm_files:
                                ui.label(f"• {Path(f).name}").classes(
                                    "text-sm text-gray-600 ml-4"
                                )
                        else:
                            # Load the single DNM file
                            dnm_file = Path(dnm_files[0])

                            try:
                                df = pl.read_csv(
                                    dnm_file,
                                    separator="\t",
                                    infer_schema_length=100,
                                )

                                # Convert to list of dicts
                                all_rows = df.to_dicts()

                                # Ensure chr:pos:ref:alt and sample_id columns exist
                                if "chr:pos:ref:alt" not in df.columns:
                                    ui.label(
                                        "Error: 'chr:pos:ref:alt' column not found in DNM file"
                                    ).classes("text-red-500")
                                elif "sample_id" not in df.columns:
                                    ui.label(
                                        "Error: 'sample_id' column not found in DNM file"
                                    ).classes("text-red-500")
                                else:
                                    # Parse variant key into components for each row
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

                                        # Extract "XX" prefix from highest_impact into Impact_priority
                                        if (
                                            "highest_impact" in row
                                            and row["highest_impact"]
                                        ):
                                            impact_val = str(row["highest_impact"])
                                            if (
                                                len(impact_val) > 3
                                                and impact_val[2] == "_"
                                            ):
                                                # Extract priority (first 2 digits)
                                                row["Impact_priority"] = impact_val[:2]
                                                # Keep full impact name with digits removed
                                                row["highest_impact"] = impact_val[3:]
                                            else:
                                                # No priority prefix, set to empty or high value
                                                row["Impact_priority"] = "99"
                                        else:
                                            row["Impact_priority"] = "99"

                                    # Sort by Impact_priority (ascending)
                                    all_rows.sort(
                                        key=lambda r: r.get("Impact_priority", "99")
                                    )

                                    # Load validation data
                                    validation_file = (
                                        store.data_dir / "validations" / "snvs.tsv"
                                    )
                                    validation_map = {}

                                    if validation_file.exists():
                                        import csv

                                        with open(validation_file, "r") as f:
                                            reader = csv.DictReader(f, delimiter="\t")
                                            for vrow in reader:
                                                vfid = vrow.get("FID")
                                                vkey = vrow.get("Variant")
                                                vsample = vrow.get("Sample")
                                                vstatus = vrow.get("Validation")
                                                vinherit = vrow.get("Inheritance")
                                                # Filter by family_id and check required fields
                                                if (
                                                    vfid == family_id
                                                    and vkey
                                                    and vsample
                                                ):
                                                    map_key = (vkey, vsample)
                                                    if map_key not in validation_map:
                                                        validation_map[map_key] = []
                                                    validation_map[map_key].append(
                                                        (vstatus, vinherit)
                                                    )

                                    # Add validation status
                                    for row in all_rows:
                                        variant_key = row.get("chr:pos:ref:alt", "")
                                        sample_id = row.get("sample_id", "")
                                        map_key = (variant_key, sample_id)

                                        if map_key in validation_map:
                                            validations = validation_map[map_key]
                                            validation_statuses = [
                                                v[0] for v in validations
                                            ]
                                            unique_validations = set(
                                                validation_statuses
                                            )

                                            if len(unique_validations) > 1:
                                                row["Validation"] = "conflicting"
                                                row["ValidationInheritance"] = ""
                                            elif "present" in unique_validations:
                                                row["Validation"] = "present"
                                                is_de_novo = any(
                                                    v[1] == "de novo"
                                                    for v in validations
                                                    if v[0] == "present"
                                                )
                                                row["ValidationInheritance"] = (
                                                    "de novo" if is_de_novo else ""
                                                )
                                            elif "absent" in unique_validations:
                                                row["Validation"] = "absent"
                                                row["ValidationInheritance"] = ""
                                            else:
                                                row["Validation"] = "uncertain"
                                                row["ValidationInheritance"] = ""
                                        else:
                                            row["Validation"] = ""
                                            row["ValidationInheritance"] = ""

                                    # Filter by selected members
                                    filtered_rows = [
                                        r
                                        for r in all_rows
                                        if r.get("sample_id")
                                        in selected_members["value"]
                                    ]

                                    # All available columns
                                    all_columns = (
                                        ["chr:pos:ref:alt"]
                                        + [
                                            col
                                            for col in df.columns
                                            if col != "chr:pos:ref:alt"
                                        ]
                                        + ["Impact_priority", "Validation"]
                                    )

                                    # Default columns to display
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
                                    # Only keep columns that actually exist in the file
                                    initial_selected = [
                                        col
                                        for col in default_visible
                                        if col in all_columns
                                    ]

                                    # Store selected columns
                                    selected_cols = {"value": initial_selected}

                                    # Get all unique highest_impact values
                                    all_impacts = sorted(
                                        set(
                                            str(row.get("highest_impact", ""))
                                            for row in all_rows
                                            if row.get("highest_impact")
                                        )
                                    )

                                    # Define LoF and coding variant categories (base names)
                                    lof_base_names = {
                                        "transcript_ablation",
                                        "splice_acceptor",
                                        "splice_donor",
                                        "stop_gained",
                                        "frameshift",
                                        "stop_lost",
                                        "start_lost",
                                    }
                                    coding_base_names = {
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

                                    # Helper function to check if impact matches a category
                                    def impact_matches_category(impact, base_names):
                                        """Check if impact matches any base name, with or without _variant suffix"""
                                        for base in base_names:
                                            if (
                                                impact == base
                                                or impact == f"{base}_variant"
                                            ):
                                                return True
                                        return False

                                    # Filter state (all impacts selected by default)
                                    filter_exclude_lcr = {"value": False}
                                    filter_exclude_gnomad = {"value": False}
                                    filter_selected_impacts = {
                                        "value": list(all_impacts)
                                    }

                                    # Create a container for the data table that can be refreshed
                                    data_container = ui.column().classes("w-full")

                                    with data_container:
                                        # Filter panel
                                        with ui.card().classes("w-full mb-4 p-4"):
                                            ui.label("Filters").classes(
                                                "text-lg font-semibold text-blue-700 mb-2"
                                            )
                                            with ui.row().classes("gap-4 items-start"):
                                                # LCR filter
                                                ui.checkbox(
                                                    "Exclude LCR",
                                                    value=filter_exclude_lcr["value"],
                                                    on_change=lambda e: (
                                                        filter_exclude_lcr.update(
                                                            {"value": e.value}
                                                        ),
                                                        render_dnm_table.refresh(),
                                                    ),
                                                )

                                                # gnomAD filter
                                                ui.checkbox(
                                                    "Exclude gnomAD 4.1 WGS filter",
                                                    value=filter_exclude_gnomad[
                                                        "value"
                                                    ],
                                                    on_change=lambda e: (
                                                        filter_exclude_gnomad.update(
                                                            {"value": e.value}
                                                        ),
                                                        render_dnm_table.refresh(),
                                                    ),
                                                )

                                                # Impact filter menu
                                                with ui.button(
                                                    "Filter by Impact",
                                                    icon="filter_list",
                                                ).props("outline color=blue"):
                                                    with ui.menu():
                                                        ui.label(
                                                            "Select Impact Types:"
                                                        ).classes(
                                                            "px-4 py-2 font-semibold text-sm"
                                                        )
                                                        ui.separator()

                                                        with ui.column().classes("p-2"):
                                                            # Quick select buttons
                                                            with ui.row().classes(
                                                                "gap-2 mb-2 flex-wrap"
                                                            ):
                                                                # Store checkbox references
                                                                impact_checkboxes = {}

                                                                def select_all_impacts():
                                                                    filter_selected_impacts[
                                                                        "value"
                                                                    ] = list(
                                                                        all_impacts
                                                                    )
                                                                    # Update checkboxes
                                                                    for cb in impact_checkboxes.values():
                                                                        cb.value = True
                                                                    render_dnm_table.refresh()

                                                                def select_none_impacts():
                                                                    filter_selected_impacts[
                                                                        "value"
                                                                    ] = []
                                                                    # Update checkboxes
                                                                    for cb in impact_checkboxes.values():
                                                                        cb.value = False
                                                                    render_dnm_table.refresh()

                                                                def select_lof():
                                                                    selected = [
                                                                        i
                                                                        for i in all_impacts
                                                                        if impact_matches_category(
                                                                            i,
                                                                            lof_base_names,
                                                                        )
                                                                    ]
                                                                    filter_selected_impacts[
                                                                        "value"
                                                                    ] = selected
                                                                    # Update checkboxes
                                                                    for (
                                                                        impact,
                                                                        cb,
                                                                    ) in impact_checkboxes.items():
                                                                        cb.value = (
                                                                            impact
                                                                            in selected
                                                                        )
                                                                    render_dnm_table.refresh()

                                                                def select_coding():
                                                                    # Include both coding and LoF variants (LoF are coding too)
                                                                    selected = [
                                                                        i
                                                                        for i in all_impacts
                                                                        if impact_matches_category(
                                                                            i,
                                                                            coding_base_names,
                                                                        )
                                                                        or impact_matches_category(
                                                                            i,
                                                                            lof_base_names,
                                                                        )
                                                                    ]
                                                                    filter_selected_impacts[
                                                                        "value"
                                                                    ] = selected
                                                                    # Update checkboxes
                                                                    for (
                                                                        impact,
                                                                        cb,
                                                                    ) in impact_checkboxes.items():
                                                                        cb.value = (
                                                                            impact
                                                                            in selected
                                                                        )
                                                                    render_dnm_table.refresh()

                                                                ui.button(
                                                                    "All",
                                                                    on_click=select_all_impacts,
                                                                ).props(
                                                                    "size=sm flat dense"
                                                                ).classes("text-xs")
                                                                ui.button(
                                                                    "None",
                                                                    on_click=select_none_impacts,
                                                                ).props(
                                                                    "size=sm flat dense"
                                                                ).classes("text-xs")
                                                                ui.button(
                                                                    "LoF",
                                                                    on_click=select_lof,
                                                                ).props(
                                                                    "size=sm flat dense color=orange"
                                                                ).classes("text-xs")
                                                                ui.button(
                                                                    "Coding",
                                                                    on_click=select_coding,
                                                                ).props(
                                                                    "size=sm flat dense color=purple"
                                                                ).classes("text-xs")

                                                            ui.separator()

                                                            # Checkboxes for each impact
                                                            for impact in all_impacts:

                                                                def handle_impact_change(
                                                                    imp, val
                                                                ):
                                                                    if val:
                                                                        if (
                                                                            imp
                                                                            not in filter_selected_impacts[
                                                                                "value"
                                                                            ]
                                                                        ):
                                                                            filter_selected_impacts[
                                                                                "value"
                                                                            ].append(
                                                                                imp
                                                                            )
                                                                    else:
                                                                        if (
                                                                            imp
                                                                            in filter_selected_impacts[
                                                                                "value"
                                                                            ]
                                                                        ):
                                                                            filter_selected_impacts[
                                                                                "value"
                                                                            ].remove(
                                                                                imp
                                                                            )
                                                                    render_dnm_table.refresh()

                                                                impact_checkboxes[
                                                                    impact
                                                                ] = ui.checkbox(
                                                                    impact,
                                                                    value=impact
                                                                    in filter_selected_impacts[
                                                                        "value"
                                                                    ],
                                                                    on_change=lambda e,
                                                                    i=impact: handle_impact_change(
                                                                        i, e.value
                                                                    ),
                                                                ).classes("text-sm")

                                        @ui.refreshable
                                        def render_dnm_table():
                                            # Filter rows by selected members
                                            rows = [
                                                r
                                                for r in all_rows
                                                if r.get("sample_id")
                                                in selected_members["value"]
                                            ]

                                            # Track total for filtering statistics
                                            total_before_filters = len(rows)

                                            # Apply filters
                                            if filter_exclude_lcr["value"]:
                                                rows = [
                                                    r
                                                    for r in rows
                                                    if r.get("LCR") != "LCR"
                                                ]

                                            if filter_exclude_gnomad["value"]:
                                                rows = [
                                                    r
                                                    for r in rows
                                                    if not r.get("genomes_filters")
                                                ]

                                            if filter_selected_impacts["value"]:
                                                rows = [
                                                    r
                                                    for r in rows
                                                    if str(r.get("highest_impact", ""))
                                                    in filter_selected_impacts["value"]
                                                ]

                                            # Calculate filtered out count
                                            filtered_out = total_before_filters - len(
                                                rows
                                            )

                                            # Helper function for column display labels
                                            def get_dnm_display_label(col):
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

                                            # Create columns definition
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
                                                            "label": get_dnm_display_label(
                                                                col
                                                            ),
                                                            "field": col,
                                                            "sortable": True,
                                                            "align": "left",
                                                        }
                                                        for col in visible_cols
                                                    ]
                                                )
                                                return cols

                                            with ui.row().classes(
                                                "items-center gap-4 mt-4 mb-2"
                                            ):
                                                ui.label(
                                                    f"Data ({len(rows)} rows"
                                                    + (
                                                        f", {filtered_out} filtered out)"
                                                        if filtered_out > 0
                                                        else ")"
                                                    )
                                                ).classes(
                                                    "text-lg font-semibold text-blue-700"
                                                )

                                                # Column selector
                                                with ui.button(
                                                    "Select Columns",
                                                    icon="view_column",
                                                ).props("outline color=blue"):
                                                    with ui.menu() as col_menu:
                                                        ui.label(
                                                            "Show/Hide Columns:"
                                                        ).classes(
                                                            "px-4 py-2 font-semibold text-sm"
                                                        )
                                                        ui.separator()

                                                        with ui.column().classes("p-2"):
                                                            # Select All / Deselect All buttons
                                                            with ui.row().classes(
                                                                "gap-2 mb-2"
                                                            ):

                                                                def select_all():
                                                                    selected_cols[
                                                                        "value"
                                                                    ] = list(
                                                                        all_columns
                                                                    )
                                                                    update_table()

                                                                def select_none():
                                                                    selected_cols[
                                                                        "value"
                                                                    ] = []
                                                                    update_table()

                                                                ui.button(
                                                                    "All",
                                                                    on_click=select_all,
                                                                ).props(
                                                                    "size=sm flat dense"
                                                                ).classes("text-xs")
                                                                ui.button(
                                                                    "None",
                                                                    on_click=select_none,
                                                                ).props(
                                                                    "size=sm flat dense"
                                                                ).classes("text-xs")

                                                            ui.separator()

                                                            # Checkboxes for each column
                                                            checkboxes = {}
                                                            for col in all_columns:
                                                                checkboxes[col] = (
                                                                    ui.checkbox(
                                                                        col,
                                                                        value=col
                                                                        in selected_cols[
                                                                            "value"
                                                                        ],
                                                                        on_change=lambda e,
                                                                        c=col: handle_col_change(
                                                                            c,
                                                                            e.value,
                                                                        ),
                                                                    ).classes("text-sm")
                                                                )

                                            with ui.card().classes(
                                                "w-full"
                                            ) as table_card:
                                                dnm_table = (
                                                    ui.table(
                                                        columns=make_columns(
                                                            selected_cols["value"]
                                                        ),
                                                        rows=rows,
                                                        pagination={"rowsPerPage": 10},
                                                    )
                                                    .classes("w-full")
                                                    .props("dense flat")
                                                )

                                                # Add view button and validation icons
                                                dnm_table.add_slot(
                                                    "body",
                                                    r"""
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
                                                                    <span v-if="col.value === 'present'" style="display: flex; align-items: center; gap: 4px;">
                                                                        <q-icon name="check_circle" color="green" size="sm">
                                                                            <q-tooltip>Validated as present</q-tooltip>
                                                                        </q-icon>
                                                                        <span v-if="props.row.ValidationInheritance === 'de novo'" style="font-weight: bold;">dnm</span>
                                                                    </span>
                                                                    <q-icon v-else-if="col.value === 'absent'" name="cancel" color="red" size="sm">
                                                                        <q-tooltip>Validated as absent</q-tooltip>
                                                                    </q-icon>
                                                                    <q-icon v-else-if="col.value === 'uncertain'" name="help" color="orange" size="sm">
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
                                                    """,
                                                )

                                                # Handle view button click
                                                import json
                                                import urllib.parse

                                                def on_view_dnm_variant(e):
                                                    row_data = e.args
                                                    variant_key = row_data.get(
                                                        "chr:pos:ref:alt", ""
                                                    )
                                                    parts = variant_key.split(":")

                                                    if len(parts) == 4:
                                                        chrom, pos, ref, alt = parts
                                                        sample_val = row_data.get(
                                                            "sample_id", ""
                                                        )

                                                        # Encode all variant data as URL parameter
                                                        data_encoded = (
                                                            urllib.parse.quote(
                                                                json.dumps(row_data)
                                                            )
                                                        )

                                                        url = f"/variant/{cohort_name}/{family_id}?chrom={chrom}&pos={pos}&ref={ref}&alt={alt}&sample={sample_val}&data={data_encoded}"
                                                        ui.navigate.to(url)

                                                dnm_table.on(
                                                    "view_variant", on_view_dnm_variant
                                                )

                                        def handle_col_change(col_name, is_checked):
                                            if (
                                                is_checked
                                                and col_name
                                                not in selected_cols["value"]
                                            ):
                                                selected_cols["value"].append(col_name)
                                            elif (
                                                not is_checked
                                                and col_name in selected_cols["value"]
                                            ):
                                                selected_cols["value"].remove(col_name)
                                            update_table()

                                        def update_table():
                                            render_dnm_table.refresh()

                                        # Initial render
                                        render_dnm_table()

                                    # Add to refreshers so table updates when member selection changes
                                    def refresh_dnm_table():
                                        render_dnm_table.refresh()

                                    data_table_refreshers.append(refresh_dnm_table)

                            except Exception as e:
                                ui.label(f"Error reading DNM file: {e}").classes(
                                    "text-red-500 mt-4"
                                )

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )


@ui.page("/variant/{cohort_name}/{family_id}")
def variant_page(
    cohort_name: str,
    family_id: str,
    chrom: str = "",
    pos: str = "",
    ref: str = "",
    alt: str = "",
    sample: str = "",
    data: str = "",
) -> None:
    """Render the variant visualization page with IGV.js."""
    import json
    import urllib.parse

    from nicegui import app as nicegui_app

    create_header()

    # Store the page client context for JavaScript execution
    from nicegui import context

    page_client = context.client

    # Get all variant data from query params (URL encoded JSON)
    try:
        variant_data = json.loads(urllib.parse.unquote(data)) if data else {}
    except Exception:
        variant_data = {}

    try:
        store = get_data_store()

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.link(cohort_name, f"/cohort/{cohort_name}").classes(
                    "text-blue-600 hover:text-blue-800"
                )
                ui.label("/").classes("text-gray-400")
                ui.link(family_id, f"/cohort/{cohort_name}/family/{family_id}").classes(
                    "text-blue-600 hover:text-blue-800"
                )
                ui.label("/").classes("text-gray-400")
                ui.label("Variant").classes("font-semibold")

            # Variant title
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
            ui.label(f"🧬 Variant: {variant_key}").classes(
                "text-3xl font-bold text-blue-900 mb-2"
            )

            # Sample section with add menu
            with ui.row().classes("items-center gap-4 mb-2"):
                ui.label(f"Sample: {sample}").classes("text-lg text-gray-600")

                # Track additional samples
                additional_samples = {"value": []}

                # Get family members to enable add parents/family options
                cohort = store.get_cohort(cohort_name)
                family_members = []
                sample_parents = {"father": None, "mother": None}

                if cohort:
                    members_data = cohort.get_family_members(family_id)
                    family_members = [m["Sample ID"] for m in members_data]
                    # Find current sample's parents
                    for member in members_data:
                        if member["Sample ID"] == sample:
                            sample_parents["father"] = member.get("Father")
                            sample_parents["mother"] = member.get("Mother")
                            break

                    # Automatically add parents if available
                    for parent_type, parent_id in sample_parents.items():
                        if (
                            parent_id
                            and parent_id != "-"
                            and parent_id != "0"
                            and parent_id != sample
                        ):
                            sample_cram = (
                                store.data_dir
                                / f"samples/{parent_id}/sequences/{parent_id}.GRCh38_GIABv3.cram"
                            )
                            if sample_cram.exists():
                                additional_samples["value"].append(parent_id)

                # Menu to add samples
                with ui.button("Add Samples", icon="add").props(
                    "outline color=blue size=sm"
                ):
                    with ui.menu():
                        ui.menu_item("Add Parents", on_click=lambda: add_parents())
                        ui.menu_item("Add Family", on_click=lambda: add_family())
                        ui.separator()
                        with ui.row().classes("items-center gap-2 px-4 py-2"):
                            barcode_input = (
                                ui.input("Barcode").classes("flex-grow").props("dense")
                            )
                            ui.button(
                                "Add",
                                icon="add",
                                on_click=lambda: add_sample(barcode_input.value),
                            ).props("flat dense size=sm")

            # Display additional samples
            additional_samples_container = ui.column().classes("gap-1 mb-4")

            with additional_samples_container:
                pass  # Will be populated dynamically

            def get_relationship_label(sample_id):
                """Get relationship label for a sample."""
                if sample_id == sample:
                    return "(carrier)"

                # Check if it's a parent
                if sample_id == sample_parents["father"]:
                    return "(father)"
                if sample_id == sample_parents["mother"]:
                    return "(mother)"

                # Check if it's a sibling (shares both parents)
                if cohort and sample_parents["father"] and sample_parents["mother"]:
                    members_data = cohort.get_family_members(family_id)
                    for member in members_data:
                        if member["Sample ID"] == sample_id:
                            member_father = member.get("Father")
                            member_mother = member.get("Mother")
                            if (
                                member_father == sample_parents["father"]
                                and member_mother == sample_parents["mother"]
                                and member_father
                                and member_mother
                                and member_father != "-"
                                and member_father != "0"
                                and member_mother != "-"
                                and member_mother != "0"
                            ):
                                return "(sibling)"
                            break

                return ""

            def refresh_additional_samples():
                additional_samples_container.clear()
                with additional_samples_container:
                    if additional_samples["value"]:
                        ui.label("Additional Samples:").classes(
                            "text-sm font-semibold text-gray-700"
                        )
                        for add_sample_id in additional_samples["value"]:
                            with ui.row().classes("items-center gap-2"):
                                label_text = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                                ui.label(label_text).classes("text-sm text-gray-600")
                                ui.button(
                                    icon="delete",
                                    on_click=lambda sid=add_sample_id: remove_sample(
                                        sid
                                    ),
                                ).props("flat dense size=xs color=red")

            # Show parents in the additional samples display if they were auto-added
            refresh_additional_samples()

            def add_sample(sample_id):
                if (
                    sample_id
                    and sample_id not in additional_samples["value"]
                    and sample_id != sample
                ):
                    # Check if sample exists
                    sample_cram = (
                        store.data_dir
                        / f"samples/{sample_id}/sequences/{sample_id}.GRCh38_GIABv3.cram"
                    )
                    if sample_cram.exists():
                        additional_samples["value"].append(sample_id)
                        refresh_additional_samples()
                        refresh_igv()
                    else:
                        ui.notify(
                            f"CRAM file not found for sample: {sample_id}",
                            type="warning",
                        )

            def add_parents():
                added = []
                for parent_type, parent_id in sample_parents.items():
                    if parent_id and parent_id != "-" and parent_id != "0":
                        if (
                            parent_id not in additional_samples["value"]
                            and parent_id != sample
                        ):
                            sample_cram = (
                                store.data_dir
                                / f"samples/{parent_id}/sequences/{parent_id}.GRCh38_GIABv3.cram"
                            )
                            if sample_cram.exists():
                                additional_samples["value"].append(parent_id)
                                added.append(parent_id)
                if added:
                    refresh_additional_samples()
                    refresh_igv()
                    ui.notify(f"Added parents: {', '.join(added)}", type="positive")
                else:
                    ui.notify("No parents to add or files not found", type="warning")

            def add_family():
                added = []
                for member_id in family_members:
                    if (
                        member_id not in additional_samples["value"]
                        and member_id != sample
                    ):
                        sample_cram = (
                            store.data_dir
                            / f"samples/{member_id}/sequences/{member_id}.GRCh38_GIABv3.cram"
                        )
                        if sample_cram.exists():
                            additional_samples["value"].append(member_id)
                            added.append(member_id)
                if added:
                    refresh_additional_samples()
                    refresh_igv()
                    ui.notify(f"Added {len(added)} family members", type="positive")
                else:
                    ui.notify("No additional family members to add", type="warning")

            def remove_sample(sample_id):
                if sample_id in additional_samples["value"]:
                    additional_samples["value"].remove(sample_id)
                    refresh_additional_samples()
                    refresh_igv()

            # Panel 1: Variant details
            ui.label("Variant Details").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.card().classes("w-full mb-6"):
                with ui.column().classes("p-4 gap-4"):
                    # Primary fields - displayed prominently
                    primary_fields = [
                        "VEP_CANONICAL",
                        "VEP_Consequence",
                        "VEP_SYMBOL",
                        "VEP_HGVSp",
                        "VEP_LoF",
                        "fafmax_faf95_max_genomes",
                    ]

                    # Fields to exclude from secondary section (already shown in title or primary)
                    excluded_fields = primary_fields + ["#CHROM", "POS", "REF", "ALT"]

                    with ui.row().classes("gap-6 flex-wrap items-center"):
                        for field in primary_fields:
                            if field in variant_data:
                                with ui.column().classes("gap-0"):
                                    ui.label(field).classes(
                                        "text-xs font-semibold text-gray-500"
                                    )
                                    ui.label(
                                        str(variant_data[field])
                                        if variant_data[field] is not None
                                        else "-"
                                    ).classes("text-base text-gray-900 font-medium")

                    # Secondary fields - collapsible
                    other_fields = {
                        k: v
                        for k, v in variant_data.items()
                        if k not in excluded_fields
                    }

                    if other_fields:
                        ui.separator()

                        show_more = {"value": False}

                        def toggle_more():
                            show_more["value"] = not show_more["value"]
                            more_button.text = (
                                "See less ▲" if show_more["value"] else "See more ▼"
                            )
                            details_container.set_visibility(show_more["value"])

                        more_button = (
                            ui.button("See more ▼", on_click=toggle_more)
                            .props("flat dense")
                            .classes("text-sm text-blue-600")
                        )

                        with ui.column().classes("gap-2 mt-2") as details_container:
                            with ui.element("div").classes("grid grid-cols-4 gap-4"):
                                for key, value in other_fields.items():
                                    with ui.column().classes("gap-0"):
                                        ui.label(key).classes(
                                            "text-xs font-semibold text-gray-500"
                                        )
                                        ui.label(
                                            str(value) if value is not None else "-"
                                        ).classes("text-sm text-gray-800 break-all")

                        details_container.set_visibility(False)

            # Panel 2: IGV.js viewer
            ui.label("Sequencing Data Viewer").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            # Build CRAM file paths
            cram_path = f"samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram"
            cram_index_path = (
                f"samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram.crai"
            )

            cram_full_path = store.data_dir / cram_path

            if not cram_full_path.exists():
                with ui.card().classes("w-full p-6 bg-yellow-50"):
                    ui.label("⚠️ CRAM file not found").classes(
                        "text-xl font-semibold text-yellow-800"
                    )
                    ui.label(f"Expected path: {cram_full_path}").classes(
                        "text-gray-600 text-sm font-mono"
                    )
            else:
                with ui.card().classes("w-full"):
                    # Create a container for IGV with dynamic height
                    def calculate_igv_height():
                        """Calculate height based on number of tracks (300px per track + 200px for controls)."""
                        num_tracks = 1 + len(additional_samples["value"])
                        return max(500, num_tracks * 300 + 200)

                    igv_container = (
                        ui.element("div")
                        .classes("w-full")
                        .style(f"height: {calculate_igv_height()}px;")
                    )
                    igv_id = f"igv-{id(igv_container)}"
                    igv_container._props["id"] = igv_id

                    # Serve the CRAM files via NiceGUI's static route
                    # Register the data directory as a static path
                    nicegui_app.add_static_files("/data", str(store.data_dir))

                    # IGV.js configuration
                    locus = f"{chrom}:{int(pos) - 100}-{int(pos) + 100}" if pos else ""

                    # Load IGV.js and initialize
                    ui.add_head_html("""
                        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
                    """)

                    # Store browser instance globally for refresh
                    browser_var = f"igvBrowser_{igv_id.replace('-', '_')}"

                    def build_igv_config():
                        """Build IGV config with all samples."""
                        tracks = []
                        # Add main sample
                        main_label = (
                            f"{sample} {get_relationship_label(sample)}".strip()
                        )
                        tracks.append(
                            {
                                "name": main_label,
                                "type": "alignment",
                                "format": "cram",
                                "url": f"/data/samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram",
                                "indexURL": f"/data/samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram.crai",
                                "height": 300,
                                "colorBy": "strand",
                            }
                        )
                        # Add additional samples
                        for add_sample_id in additional_samples["value"]:
                            track_label = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                            tracks.append(
                                {
                                    "name": track_label,
                                    "type": "alignment",
                                    "format": "cram",
                                    "url": f"/data/samples/{add_sample_id}/sequences/{add_sample_id}.GRCh38_GIABv3.cram",
                                    "indexURL": f"/data/samples/{add_sample_id}/sequences/{add_sample_id}.GRCh38_GIABv3.cram.crai",
                                    "height": 300,
                                    "colorBy": "strand",
                                }
                            )
                        return {
                            "genome": "hg38",
                            "locus": locus,
                            "tracks": tracks,
                        }

                    def refresh_igv():
                        """Refresh IGV browser with updated tracks."""
                        config = build_igv_config()
                        # Update container height (may fail if called from deleted context)
                        try:
                            new_height = calculate_igv_height()
                            igv_container.style(f"height: {new_height}px;")
                        except RuntimeError:
                            pass  # Container update failed, but JS update will still work
                        # Update IGV tracks using page client context
                        page_client.run_javascript(
                            f"""
                            if (window.{browser_var}) {{
                                window.{browser_var}.removeAllTracks();
                                const tracks = {json.dumps(config["tracks"])};
                                for (const track of tracks) {{
                                    window.{browser_var}.loadTrack(track);
                                }}
                            }}
                        """,
                            timeout=5.0,
                        )

                    # Initialize IGV after page load
                    igv_config = build_igv_config()
                    ui.run_javascript(f'''
                        setTimeout(function() {{
                            var igvDiv = document.getElementById("{igv_id}");
                            if (igvDiv && typeof igv !== 'undefined') {{
                                igv.createBrowser(igvDiv, {json.dumps(igv_config)})
                                    .then(function(browser) {{
                                        window.{browser_var} = browser;
                                        console.log("IGV browser created successfully");
                                    }})
                                    .catch(function(error) {{
                                        console.error("Error creating IGV browser:", error);
                                    }});
                            }} else {{
                                console.error("IGV container not found or igv not loaded");
                            }}
                        }}, 500);
                    ''')

                    # Show file info
                    with ui.row().classes("items-center gap-4 p-4 bg-gray-50"):
                        ui.label("CRAM:").classes("font-semibold text-sm")
                        ui.label(cram_path).classes("text-xs text-gray-600 font-mono")

            # Panel 3: Variant Validation
            ui.label("Variant Validation").classes(
                "text-2xl font-semibold mb-4 text-blue-800 mt-6"
            )

            with ui.card().classes("w-full p-6"):
                with ui.column().classes("gap-4"):
                    # Get system username as default
                    import getpass

                    default_user = getpass.getuser()

                    # Single row with all controls
                    with ui.row().classes("items-center gap-4 w-full"):
                        ui.label("User:").classes("font-semibold")
                        user_input = (
                            ui.input("Username").props("outlined dense").classes("w-48")
                        )
                        user_input.value = default_user

                        ui.label("Inheritance:").classes("font-semibold ml-4")
                        inheritance_select = (
                            ui.select(
                                [
                                    "unknown",
                                    "de novo",
                                    "paternal",
                                    "maternal",
                                    "either",
                                ],
                                value="unknown",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                        ui.label("Validation:").classes("font-semibold ml-4")
                        validation_select = (
                            ui.select(
                                ["uncertain", "present", "absent", "different"],
                                value="uncertain",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                        ui.button(
                            "Save Validation",
                            icon="save",
                            on_click=lambda: save_validation(),
                        ).props("color=blue").classes("ml-4")

                    # Container for validation history table
                    validation_history_container = ui.column().classes("w-full mt-4")

                    def load_validations():
                        """Load and display existing validations for this variant."""
                        validation_file = store.data_dir / "validations" / "snvs.tsv"

                        validation_history_container.clear()

                        if not validation_file.exists():
                            with validation_history_container:
                                ui.label("No validations recorded yet").classes(
                                    "text-gray-500 text-sm italic"
                                )
                            return

                        try:
                            # Read TSV file
                            import csv

                            validations = []
                            variant_key = f"{chrom}:{pos}:{ref}:{alt}"

                            with open(validation_file, "r") as f:
                                reader = csv.DictReader(f, delimiter="\t")
                                for row in reader:
                                    # Filter for current family, variant and sample
                                    if (
                                        row.get("FID") == family_id
                                        and row.get("Variant") == variant_key
                                        and row.get("Sample") == sample
                                    ):
                                        validations.append(row)

                            with validation_history_container:
                                if validations:
                                    ui.label("Previous Validations:").classes(
                                        "font-semibold text-gray-700 mb-2"
                                    )

                                    def delete_validation(timestamp):
                                        """Delete a specific validation entry from snvs.tsv."""
                                        import fcntl

                                        try:
                                            # Read all lines
                                            with open(validation_file, "r") as f:
                                                lines = f.readlines()

                                            # Filter out the validation to delete
                                            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                                            filtered_lines = []
                                            header_line = None

                                            for i, line in enumerate(lines):
                                                if i == 0:
                                                    # Keep header
                                                    header_line = line
                                                    filtered_lines.append(line)
                                                else:
                                                    # Parse line to check if it matches
                                                    parts = line.strip().split("\t")
                                                    if len(parts) >= 7:
                                                        line_fid = parts[0]
                                                        line_variant = parts[1]
                                                        line_timestamp = (
                                                            parts[6]
                                                            if len(parts) > 6
                                                            else ""
                                                        )
                                                        # Keep line if it doesn't match the one to delete
                                                        if not (
                                                            line_fid == family_id
                                                            and line_variant
                                                            == variant_key
                                                            and line_timestamp
                                                            == timestamp
                                                        ):
                                                            filtered_lines.append(line)

                                            # Write back with file locking
                                            with open(validation_file, "w") as f:
                                                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                                                try:
                                                    f.writelines(filtered_lines)
                                                    f.flush()
                                                finally:
                                                    fcntl.flock(
                                                        f.fileno(), fcntl.LOCK_UN
                                                    )

                                            ui.notify(
                                                "Validation deleted successfully",
                                                type="positive",
                                            )
                                            load_validations()  # Refresh the display
                                        except Exception as e:
                                            ui.notify(
                                                f"Error deleting validation: {str(e)}",
                                                type="negative",
                                            )

                                    # Create table data with delete action column
                                    columns = [
                                        {
                                            "name": "sample",
                                            "label": "Sample",
                                            "field": "Sample",
                                            "align": "left",
                                        },
                                        {
                                            "name": "user",
                                            "label": "User",
                                            "field": "User",
                                            "align": "left",
                                        },
                                        {
                                            "name": "inheritance",
                                            "label": "Inheritance",
                                            "field": "Inheritance",
                                            "align": "left",
                                        },
                                        {
                                            "name": "validation",
                                            "label": "Validation",
                                            "field": "Validation",
                                            "align": "left",
                                        },
                                        {
                                            "name": "timestamp",
                                            "label": "Timestamp",
                                            "field": "Timestamp",
                                            "align": "left",
                                        },
                                        {
                                            "name": "actions",
                                            "label": "",
                                            "field": "actions",
                                            "align": "center",
                                        },
                                    ]

                                    validation_table = (
                                        ui.table(
                                            columns=columns,
                                            rows=validations,
                                            row_key="Timestamp",
                                        )
                                        .classes("w-full")
                                        .props("dense flat")
                                    )

                                    # Add delete button slot
                                    validation_table.add_slot(
                                        "body",
                                        r"""
                                            <q-tr :props="props">
                                                <q-td v-for="col in props.cols.filter(c => c.name !== 'actions')" :key="col.name" :props="props">
                                                    {{ col.value }}
                                                </q-td>
                                                <q-td key="actions" :props="props">
                                                    <q-btn 
                                                        flat 
                                                        dense 
                                                        size="xs" 
                                                        icon="delete" 
                                                        color="red"
                                                        @click="$parent.$emit('delete_validation', props.row.Timestamp)"
                                                    >
                                                        <q-tooltip>Delete this validation</q-tooltip>
                                                    </q-btn>
                                                </q-td>
                                            </q-tr>
                                        """,
                                    )

                                    # Handle delete button click
                                    validation_table.on(
                                        "delete_validation",
                                        lambda e: delete_validation(e.args),
                                    )
                                else:
                                    ui.label(
                                        "No validations recorded for this variant yet"
                                    ).classes("text-gray-500 text-sm italic")
                        except Exception as e:
                            with validation_history_container:
                                ui.label(
                                    f"Error loading validations: {str(e)}"
                                ).classes("text-red-500 text-sm")

                    def save_validation():
                        """Save validation to TSV file with file locking."""
                        import fcntl
                        from datetime import datetime

                        # Prepare validation directory
                        validation_dir = store.data_dir / "validations"
                        validation_dir.mkdir(parents=True, exist_ok=True)

                        validation_file = validation_dir / "snvs.tsv"

                        # Prepare data row
                        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                        row_data = [
                            family_id,
                            variant_key,
                            sample,
                            user_input.value or default_user,
                            inheritance_select.value,
                            validation_select.value,
                            timestamp,
                        ]

                        # Write with file locking
                        try:
                            # Check if file exists to determine if we need header
                            file_exists = validation_file.exists()

                            with open(validation_file, "a") as f:
                                # Acquire exclusive lock
                                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                                try:
                                    # Write header if new file
                                    if not file_exists:
                                        header = [
                                            "FID",
                                            "Variant",
                                            "Sample",
                                            "User",
                                            "Inheritance",
                                            "Validation",
                                            "Timestamp",
                                        ]
                                        f.write("\t".join(header) + "\n")

                                    # Write data row
                                    f.write("\t".join(row_data) + "\n")
                                    f.flush()
                                finally:
                                    # Release lock
                                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

                            ui.notify(
                                f"Validation saved successfully at {timestamp}",
                                type="positive",
                            )
                            # Reload validations to show the new entry
                            load_validations()
                        except Exception as e:
                            ui.notify(
                                f"Error saving validation: {str(e)}", type="negative"
                            )

                    # Load existing validations on page load
                    load_validations()

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )


# Validation pages are now in genetics_viz.pages.validation module


def run_app(
    data_dir: Path,
    host: str = "127.0.0.1",
    port: int = 8080,
    reload: bool = False,
) -> None:
    """
    Initialize and run the NiceGUI application.

    Args:
        data_dir: Path to the data directory containing cohort data
        host: Host address to bind the server to
        port: Port to run the server on
        reload: Enable auto-reload for development
    """
    global _data_store

    # Store config in environment for reload mode
    if reload:
        os.environ["GENETICS_VIZ_DATA_DIR"] = str(data_dir)
        os.environ["GENETICS_VIZ_HOST"] = host
        os.environ["GENETICS_VIZ_PORT"] = str(port)
        os.environ["GENETICS_VIZ_RELOAD"] = "1"

    # Initialize the data store - both local and in shared module
    _data_store = DataStore(data_dir=data_dir)
    data_module._data_store = _data_store

    try:
        _data_store.load()
        print(f"Loaded {len(_data_store.cohorts)} cohorts")
    except FileNotFoundError as e:
        print(f"Warning: {e}")
        print("The application will start but no cohorts will be available.")

    # Configure NiceGUI
    ui.run(
        host=host,
        port=port,
        reload=reload,
        title="Genetics-Viz",
        favicon="🧬",
        dark=False,
    )


# Auto-initialize and run when module is reloaded (for reload mode)
if os.environ.get("GENETICS_VIZ_RELOAD") == "1" and _data_store is None:
    data_dir_env = os.environ.get("GENETICS_VIZ_DATA_DIR")
    if data_dir_env:
        _data_store = DataStore(data_dir=Path(data_dir_env))
        data_module._data_store = _data_store
        try:
            _data_store.load()
            print(f"[Reload] Loaded {len(_data_store.cohorts)} cohorts")
        except FileNotFoundError as e:
            print(f"[Reload] Warning: {e}")

        # Call ui.run() at module level for reload mode
        host = os.environ.get("GENETICS_VIZ_HOST", "127.0.0.1")
        port = int(os.environ.get("GENETICS_VIZ_PORT", "8080"))
        ui.run(
            host=host,
            port=port,
            reload=True,
            title="Genetics-Viz",
            favicon="🧬",
            dark=False,
        )
