"""Validation file page - displays a specific to_validate file."""

import csv
from typing import Any, Dict, List

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.tables import VALIDATION_TABLE_SLOT
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.data import get_data_store


def _load_validation_map(validation_file_path) -> Dict[tuple, List[tuple]]:
    """Load validation data from snvs.tsv into a lookup map.

    Returns:
        Dictionary mapping (fid, variant_key, sample_id) to list of (validation_status, inheritance, comment, ignore)
    """
    validation_map: Dict[tuple, List[tuple]] = {}

    if validation_file_path.exists():
        with open(validation_file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for vrow in reader:
                fid = vrow.get("FID")
                variant_key = vrow.get("Variant")
                sample_id = vrow.get("Sample")
                validation_status = vrow.get("Validation")
                inheritance = vrow.get("Inheritance")
                comment = vrow.get("Comment", "")
                ignore = vrow.get("Ignore", "0")

                # Only include non-ignored validations
                if fid and variant_key and sample_id and ignore != "1":
                    map_key = (fid, variant_key, sample_id)
                    if map_key not in validation_map:
                        validation_map[map_key] = []
                    validation_map[map_key].append(
                        (validation_status, inheritance, comment, ignore)
                    )

    return validation_map


def _add_validation_status_to_rows(
    file_data: List[Dict[str, Any]],
    validation_map: Dict[tuple, List[tuple]],
    fid_col: str,
    variant_col: str,
    sample_col: str,
) -> None:
    """Add Validation status to each row based on validation map."""
    for row in file_data:
        fid = row.get(fid_col, "")
        variant = row.get(variant_col, "")
        sample = row.get(sample_col, "")

        map_key = (fid, variant, sample)
        if map_key in validation_map:
            validations = validation_map[map_key]
            validation_statuses = [v[0] for v in validations]
            # Normalize "in phase MNV" to "present" for conflict detection
            normalized_statuses = [
                "present" if s == "in phase MNV" else s for s in validation_statuses
            ]
            unique_validations = set(normalized_statuses)

            if len(unique_validations) > 1:
                row["Validation"] = "conflicting"
                row["ValidationInheritance"] = ""
            elif "present" in unique_validations:
                # Check if any is specifically "in phase MNV"
                if "in phase MNV" in validation_statuses:
                    row["Validation"] = "in phase MNV"
                else:
                    row["Validation"] = "present"
                # Check inheritance - prioritize de novo, then homozygous
                is_de_novo = any(
                    v[1] == "de novo"
                    for v in validations
                    if v[0] in ("present", "in phase MNV")
                )
                is_homozygous = any(
                    v[1] == "homozygous"
                    for v in validations
                    if v[0] in ("present", "in phase MNV")
                )
                if is_de_novo:
                    row["ValidationInheritance"] = "de novo"
                elif is_homozygous:
                    row["ValidationInheritance"] = "homozygous"
                else:
                    row["ValidationInheritance"] = ""
            elif "absent" in unique_validations:
                row["Validation"] = "absent"
                row["ValidationInheritance"] = ""
            else:
                row["Validation"] = "uncertain"
                row["ValidationInheritance"] = ""
        else:
            row["Validation"] = ""
            row["ValidationInheritance"] = ""


@ui.page("/validation/file/{filename}")
def validation_file_page(filename: str) -> None:
    """Render a specific to_validate file."""
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()
        to_validate_dir = store.data_dir / "to_validate"
        file_path = to_validate_dir / f"{filename}.tsv"

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"🔍 Validating: {filename}").classes(
                    "text-3xl font-bold text-blue-900"
                )

            if not file_path.exists():
                ui.label(f"File not found: {filename}.tsv").classes(
                    "text-red-500 text-lg"
                )
                return

            # Read TSV file
            file_data: List[Dict[str, Any]] = []
            headers: List[str] = []
            with open(file_path, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                headers = list(reader.fieldnames or [])
                for row in reader:
                    file_data.append(dict(row))

            if not file_data:
                ui.label("No data in selected file").classes("text-gray-500 italic")
                return

            # Check for required columns (case-insensitive)
            headers_lower = {h.lower(): h for h in headers}
            has_variant = "variant" in headers_lower
            has_sample = "sample" in headers_lower
            has_fid = "fid" in headers_lower

            # Get actual column names from file
            variant_col = headers_lower.get("variant", "Variant")
            sample_col = headers_lower.get("sample", "Sample")
            fid_col = headers_lower.get("fid", "FID")

            if not (has_variant and has_sample and has_fid):
                missing = []
                if not has_variant:
                    missing.append("Variant")
                if not has_sample:
                    missing.append("Sample")
                if not has_fid:
                    missing.append("FID")
                ui.label(
                    f"⚠️ Warning: Missing required columns: {', '.join(missing)}"
                ).classes("text-orange-600 text-sm mb-2")

            # Load validation data
            validation_file = store.data_dir / "validations" / "snvs.tsv"
            validation_map = _load_validation_map(validation_file)

            # Add Validation status to each row
            _add_validation_status_to_rows(
                file_data, validation_map, fid_col, variant_col, sample_col
            )

            # Filter state - all statuses selected by default
            all_validation_statuses = [
                "present",
                "absent",
                "uncertain",
                "conflicting",
                "TODO",
            ]
            filter_validations: Dict[str, List[str]] = {
                "value": list(all_validation_statuses)
            }

            # Create filter menu first (before table)
            create_validation_filter_menu(
                all_statuses=all_validation_statuses,
                filter_state=filter_validations,
                on_change=lambda: refresh_table(),
            )

            # Table container
            table_container = ui.column().classes("w-full")

            # Capture the client context for use in callbacks
            from nicegui import context

            page_client = context.client

            @ui.refreshable
            def refresh_table():
                """Refresh the table with current filters."""
                table_container.clear()

                # Apply filters
                filtered_data = file_data.copy()
                if filter_validations["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Validation", "") in filter_validations["value"]
                        or (
                            "TODO" in filter_validations["value"]
                            and not row.get("Validation")
                        )
                    ]

                with table_container:
                    # Show count
                    if filter_validations["value"] != all_validation_statuses:
                        ui.label(
                            f"Showing {len(filtered_data)} of {len(file_data)} variants"
                        ).classes("text-sm text-gray-600 mb-2")
                    else:
                        ui.label(f"{len(filtered_data)} variants to validate").classes(
                            "text-sm text-gray-600 mb-2"
                        )

                    # Prepare columns for table
                    columns: List[Dict[str, Any]] = [
                        {"name": "actions", "label": "", "field": "actions"}
                    ]
                    for header in headers:
                        columns.append(
                            {
                                "name": header,
                                "label": header,
                                "field": header,
                                "sortable": True,
                                "align": "left",
                            }
                        )
                    # Add Validation column
                    columns.append(
                        {
                            "name": "Validation",
                            "label": "Validation",
                            "field": "Validation",
                            "sortable": True,
                            "align": "left",
                        }
                    )

                    # Create table with pagination
                    validation_table = (
                        ui.table(
                            columns=columns,
                            rows=filtered_data,
                            row_key=variant_col if has_variant else "Variant",
                            pagination={"rowsPerPage": 50},
                        )
                        .classes("w-full")
                        .props("dense flat")
                    )

                    # Add custom slot for view button and validation icons
                    validation_table.add_slot("body", VALIDATION_TABLE_SLOT)

                    # Handle view button click
                    def on_view_variant(e):
                        row_data = e.args
                        family_id = row_data.get(fid_col, "")
                        variant_str = row_data.get(variant_col, "")
                        sample_id = row_data.get(sample_col, "")

                        try:
                            parts = variant_str.split(":")
                            if len(parts) == 4:
                                chrom, pos, ref, alt = parts

                                # Find the cohort from family_id
                                cohort_name = None
                                for c_name, cohort in store.cohorts.items():
                                    if family_id in cohort.families:
                                        cohort_name = c_name
                                        break

                                if not cohort_name:
                                    ui.notify(
                                        f"Could not find cohort for family {family_id}",
                                        type="warning",
                                    )
                                    return

                                # Create variant data dict
                                variant_data = dict(row_data)

                                # Callback to update the Validation column in the table
                                def on_save(validation_status: str):
                                    # Reload validation data from file
                                    validation_map = _load_validation_map(
                                        validation_file
                                    )
                                    # Re-add validation status to rows
                                    _add_validation_status_to_rows(
                                        file_data,
                                        validation_map,
                                        fid_col,
                                        variant_col,
                                        sample_col,
                                    )
                                    # Refresh the table display using the captured client context
                                    with page_client:
                                        ui.timer(0.1, refresh_table, once=True)

                                # Show dialog
                                show_variant_dialog(
                                    cohort_name=cohort_name,
                                    family_id=family_id,
                                    chrom=chrom,
                                    pos=pos,
                                    ref=ref,
                                    alt=alt,
                                    sample=sample_id,
                                    variant_data=variant_data,
                                    on_save_callback=on_save,
                                )
                            else:
                                ui.notify(
                                    "Invalid variant format. Expected chr:pos:ref:alt",
                                    type="warning",
                                )
                        except Exception as ex:
                            ui.notify(f"Error parsing variant: {ex}", type="warning")

                    validation_table.on("view_variant", on_view_variant)

            # Initial render
            refresh_table()

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
