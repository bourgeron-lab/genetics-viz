"""Validation file page - displays a specific to_validate file."""

import asyncio
import csv
import re
import shutil
from pathlib import Path
from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.search_stats import show_stats_dialog
from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    genomic_sort_key,
    get_column_group,
    get_column_sorting,
    get_column_type,
    get_display_label,
)
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.score_colors import get_score_color
from genetics_viz.utils.validation_badges import build_validation_badge
from genetics_viz.utils.wisecondorx import build_call_colors, infer_sv_type
from genetics_viz.utils.clinvar import (
    format_clinvar_display,
    get_clinvar_color,
)
from genetics_viz.utils.vep import (
    format_consequence_display,
    get_consequence_color,
)
from genetics_viz.utils.view_presets import VIEW_PRESETS

_SV_PATTERN = re.compile(r"^chr[^:]+:\d+-\d+$")


def _read_tsv_file(file_path: Path) -> tuple:
    """Read TSV file and return (headers, rows) tuple.

    Blocking I/O — intended to be called via asyncio.to_thread().
    """
    file_data: List[Dict[str, Any]] = []
    headers: List[str] = []
    with open(file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = list(reader.fieldnames or [])
        for row in reader:
            file_data.append(dict(row))
    return headers, file_data


def _detect_sv_format(file_data: List[Dict[str, Any]], variant_col: str) -> bool:
    """Check if any variant in the file matches the SV format (chr:start-end)."""
    return any(_SV_PATTERN.match(row.get(variant_col, "")) for row in file_data)


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
    is_sv: bool = False,
) -> None:
    """Add Validation status to each row based on validation map."""
    for row in file_data:
        fid = row.get(fid_col, "")
        variant = row.get(variant_col, "")
        sample = row.get(sample_col, "")

        if is_sv and _SV_PATTERN.match(variant):
            sv_type = infer_sv_type(row)
            lookup_variant = f"{variant}:{sv_type}"
        else:
            lookup_variant = variant

        map_key = (fid, lookup_variant, sample)
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
                # Check inheritance - prioritize de novo, then homozygous,
                # then first non-empty inheritance from present validations
                present = [
                    v for v in validations if v[0] in ("present", "in phase MNV")
                ]
                inh_values = [v[1] for v in present if v[1]]
                if "de novo" in inh_values:
                    row["ValidationInheritance"] = "de novo"
                elif "homozygous" in inh_values:
                    row["ValidationInheritance"] = "homozygous"
                elif inh_values:
                    row["ValidationInheritance"] = inh_values[0]
                else:
                    row["ValidationInheritance"] = ""
            elif "absent" in unique_validations:
                row["Validation"] = "absent"
                row["ValidationInheritance"] = ""
            else:
                row["Validation"] = "uncertain"
                row["ValidationInheritance"] = ""

            row["Validation_badge"] = build_validation_badge(
                row["Validation"], row["ValidationInheritance"], validations
            )
        else:
            row["Validation"] = ""
            row["ValidationInheritance"] = ""
            row["Validation_badge"] = None


@ui.page("/validation/file/{filename}")
async def validation_file_page(filename: str) -> None:
    """Render a specific to_validate file."""
    if redirect := check_auth():
        return redirect
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()
        to_validate_dir = store.data_dir / "to_validate"
        file_path = to_validate_dir / f"{filename}.tsv"

        with ui.column().classes("w-full px-6 py-6"):
            # Title row with trash button on the right
            with ui.row().classes("items-center w-full mb-6"):
                ui.label(f"🔍 Validating: {filename}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.space()

                # Trash file button with confirmation dialog
                trash_dialog = ui.dialog()
                with trash_dialog, ui.card().classes("p-4"):
                    ui.label("Trash this file?").classes("text-lg font-semibold mb-2")
                    ui.label(
                        f"{filename}.tsv will be moved to the trash folder."
                    ).classes("text-sm text-gray-600 mb-4")
                    with ui.row().classes("justify-end gap-2 w-full"):
                        ui.button("Cancel", on_click=trash_dialog.close).props("flat")

                        def _confirm_trash():
                            trash_dir = to_validate_dir / "trash"
                            trash_dir.mkdir(parents=True, exist_ok=True)
                            dest = trash_dir / f"{filename}.tsv"
                            counter = 1
                            while dest.exists():
                                dest = trash_dir / f"{filename}_{counter}.tsv"
                                counter += 1
                            shutil.move(str(file_path), str(dest))
                            trash_dialog.close()
                            ui.notify(
                                f"{filename}.tsv moved to trash",
                                type="positive",
                                position="top",
                            )
                            ui.navigate.to("/validation/all")

                        ui.button("Trash", on_click=_confirm_trash).props(
                            "color=red unelevated"
                        ).classes("text-white")

                ui.button(icon="delete", on_click=trash_dialog.open).props(
                    "flat round color=red"
                ).tooltip("Move file to trash")

            if not file_path.exists():
                ui.label(f"File not found: {filename}.tsv").classes(
                    "text-red-500 text-lg"
                )
                return

            # Read TSV file (offloaded to thread to avoid blocking event loop)
            headers, file_data = await asyncio.to_thread(_read_tsv_file, file_path)

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

            # Detect SV vs SNV format
            is_sv_format = (
                _detect_sv_format(file_data, variant_col) if has_variant else False
            )

            # Load validation data (offloaded to thread)
            if is_sv_format:
                validation_file = store.data_dir / "validations" / "svs.tsv"
                diagnostic_file = store.data_dir / "diagnostics" / "svs.tsv"
            else:
                validation_file = store.data_dir / "validations" / "snvs.tsv"
                diagnostic_file = store.data_dir / "diagnostics" / "snvs.tsv"
            validation_map = await asyncio.to_thread(
                _load_validation_map, validation_file
            )
            diagnostic_map = load_diagnostic_map(diagnostic_file)

            # Add Validation status to each row
            _add_validation_status_to_rows(
                file_data,
                validation_map,
                fid_col,
                variant_col,
                sample_col,
                is_sv=is_sv_format,
            )

            # Add Diagnostic status to each row
            for row in file_data:
                variant = row.get(variant_col, "")
                sample = row.get(sample_col, "")
                if is_sv_format and _SV_PATTERN.match(variant):
                    sv_type = infer_sv_type(row)
                    diag_variant_key = f"{variant}:{sv_type}"
                else:
                    diag_variant_key = variant
                add_diagnostic_status_to_row(
                    row, diagnostic_map, diag_variant_key, sample
                )

            # Build badge data for gene, consequence, and ClinVar columns
            gene_scorer = get_gene_scorer()
            for row in file_data:
                for col_name in headers:
                    col_lower = col_name.lower()

                    # Gene/symbol columns → gene_badge
                    if "symbol" in col_lower or "gene" in col_lower:
                        value = row.get(col_name, "")
                        if value and value != "-":
                            genes = [
                                g.strip() for g in str(value).split(",") if g.strip()
                            ]
                            row[f"{col_name}_badges"] = [
                                {
                                    "label": gene,
                                    "color": gene_scorer.get_gene_color(gene),
                                    "tooltip": gene_scorer.get_gene_tooltip(gene),
                                }
                                for gene in genes
                            ]
                        else:
                            row[f"{col_name}_badges"] = []

                    # Consequence / impact columns → badge_list
                    elif "impact" in col_lower or col_name == "VEP_Consequence":
                        cons_str = row.get(col_name, "")
                        if cons_str:
                            terms = [
                                t.strip() for t in str(cons_str).split("&") if t.strip()
                            ]
                            badges = []
                            seen_badges: set = set()
                            for term in terms:
                                label = format_consequence_display(term)
                                color = get_consequence_color(term)
                                key = (label, color)
                                if key not in seen_badges:
                                    seen_badges.add(key)
                                    badges.append({"label": label, "color": color})
                            row[f"{col_name}_badges"] = badges
                        else:
                            row[f"{col_name}_badges"] = []

                    # ClinVar column → badge_list
                    elif col_name == "VEP_CLIN_SIG":
                        clin_str = row.get(col_name, "")
                        if clin_str:
                            sigs = [
                                s.strip()
                                for part in str(clin_str).split(",")
                                for s in part.split("&")
                                if s.strip() and s.strip() != "."
                            ]
                            badges = []
                            seen_badges = set()
                            for sig in sigs:
                                label = format_clinvar_display(sig)
                                color = get_clinvar_color(sig)
                                key = (label, color)
                                if key not in seen_badges:
                                    seen_badges.add(key)
                                    badges.append({"label": label, "color": color})
                            row[f"{col_name}_badges"] = badges
                        else:
                            row[f"{col_name}_badges"] = []

                # Continuous score badges (after per-column processing)
                for col_name, value_str in list(row.items()):
                    if value_str and value_str != ".":
                        try:
                            value = float(value_str)
                            badge_info = get_score_color(col_name, value)
                            if badge_info:
                                row[f"{col_name}_badge"] = {
                                    "label": f"{value:.3f}",
                                    "color": badge_info["color"],
                                    "tooltip": (
                                        f"{col_name}: {value:.3f}"
                                        f" ({badge_info['label']})"
                                    ),
                                }
                        except (ValueError, TypeError):
                            pass

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

            # Column visibility state — all columns visible by default
            all_columns = list(headers) + ["Validation", "Diagnostic"]
            selected_cols: Dict[str, Any] = {"value": list(all_columns)}
            dt_ref: Dict[str, Any] = {"ref": None}
            table_state: Dict[str, Any] = {"sorting": [], "page": 0}

            def _apply_col_visibility():
                if dt_ref["ref"]:
                    visible = ["actions"] + list(selected_cols["value"])
                    dt_ref["ref"].set_column_visibility(visible)

            # Column selector dialog
            col_dialog, _sync_col_selector = build_column_selector(
                all_columns=all_columns,
                selected_cols=selected_cols,
                on_visibility_change=_apply_col_visibility,
                presets=VIEW_PRESETS,
            )

            # Toolbar row: validation filter + Columns + Stats buttons
            with ui.row().classes("items-center gap-2 w-full"):
                create_validation_filter_menu(
                    all_statuses=all_validation_statuses,
                    filter_state=filter_validations,
                    on_change=lambda: refresh_table(),
                )
                ui.space()
                ui.button(
                    "Columns",
                    icon="view_column",
                    on_click=col_dialog.open,
                ).props("outline color=blue size=sm")
                ui.button(
                    "Stats",
                    icon="bar_chart",
                    on_click=lambda: show_stats_dialog(file_data),
                ).props("outline color=blue size=sm")

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

                    # Prepare columns for table (same pattern as search/wombat pages)
                    columns: List[Dict[str, Any]] = [
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
                        col_lower = col.lower()
                        if col == "Validation":
                            col_def["cellType"] = "validation"
                        elif col == "Diagnostic":
                            col_def["cellType"] = "diagnostic"
                        elif col == "Variant":
                            col_def["sorting"] = "genomic"
                        elif (
                            "symbol" in col_lower
                            or col_lower == "gene"
                            or "gene" in col_lower
                        ):
                            col_def["cellType"] = "gene_badge"
                            col_def["badgesField"] = f"{col}_badges"
                        elif "impact" in col_lower or col == "VEP_Consequence":
                            col_def["cellType"] = "badge_list"
                            col_def["badgesField"] = f"{col}_badges"
                        elif col == "VEP_CLIN_SIG":
                            col_def["cellType"] = "badge_list"
                            col_def["badgesField"] = f"{col}_badges"
                        elif col == "wisecondorX":
                            col_def["cellType"] = "cnv_call"
                            col_def["callColors"] = build_call_colors()
                        else:
                            col_type = get_column_type(col)
                            if col_type in ("int", "float"):
                                col_def["cellType"] = "number"
                            else:
                                col_def["cellType"] = "score_badge"
                        apply_width_constraints(col_def, col)
                        columns.append(col_def)

                    def _refresh_file_data():
                        """Reload validation + diagnostic data and refresh table."""
                        updated_map = _load_validation_map(validation_file)
                        _add_validation_status_to_rows(
                            file_data,
                            updated_map,
                            fid_col,
                            variant_col,
                            sample_col,
                            is_sv=is_sv_format,
                        )
                        diag_map = load_diagnostic_map(diagnostic_file)
                        for row in file_data:
                            v = row.get(variant_col, "")
                            s = row.get(sample_col, "")
                            if is_sv_format and _SV_PATTERN.match(v):
                                dvk = f"{v}:{infer_sv_type(row)}"
                            else:
                                dvk = v
                            add_diagnostic_status_to_row(row, diag_map, dvk, s)
                        with page_client:
                            refresh_table()

                    # Handle row action click
                    def on_row_action(e):
                        action = e.get("action", "")
                        row_data = e.get("row", {})
                        family_id = row_data.get(fid_col, "")
                        variant_str = row_data.get(variant_col, "")
                        sample_id = row_data.get(sample_col, "")

                        if action == "open_diagnostic":
                            # Determine gene and impact from row data
                            gene = row_data.get(
                                "VEP_SYMBOL",
                                row_data.get("Gene", ""),
                            )
                            impact = row_data.get(
                                "VEP_Consequence",
                                row_data.get("Impact", ""),
                            )
                            if is_sv_format and _SV_PATTERN.match(variant_str):
                                vtype = "sv"
                                # Parse gene symbols from SV gene column
                                gene_str = row_data.get("gene", gene)
                                if gene_str and ":" in gene_str:
                                    syms = []
                                    for gp in str(gene_str).split(","):
                                        sym = gp.split(":")[0].strip()
                                        if sym:
                                            syms.append(sym)
                                    gene = ", ".join(syms) if syms else gene
                                sv_call = str(
                                    row_data.get(
                                        "wisecondorX",
                                        row_data.get("call", ""),
                                    )
                                ).upper()
                                if "GAIN" in sv_call:
                                    impact = "GAIN"
                                elif "LOSS" in sv_call:
                                    impact = "LOSS"
                            else:
                                vtype = "snv"

                            show_diagnostic_dialog(
                                family_id=family_id,
                                variant_key=variant_str,
                                gene=gene,
                                impact=impact,
                                sample=sample_id,
                                variant_type=vtype,
                                diagnostic_file=diagnostic_file,
                                on_save_callback=lambda s: _refresh_file_data(),
                            )
                            return

                        try:
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

                            if is_sv_format and _SV_PATTERN.match(variant_str):
                                # SV variant: chr:start-end
                                parts = variant_str.split(":")
                                chrom = parts[0]
                                range_parts = parts[1].split("-")
                                start, end = range_parts[0], range_parts[1]

                                sv_data = dict(row_data)
                                if "wisecondorX" in sv_data and "call" not in sv_data:
                                    sv_data["call"] = sv_data["wisecondorX"]

                                show_sv_dialog(
                                    cohort_name=cohort_name,
                                    family_id=family_id,
                                    chrom=chrom,
                                    start=start,
                                    end=end,
                                    sample=sample_id,
                                    sv_data=sv_data,
                                    on_validation_saved=_refresh_file_data,
                                )
                            else:
                                # SNV variant: chr:pos:ref:alt
                                parts = variant_str.split(":")
                                if len(parts) == 4:
                                    chrom, pos, ref, alt = parts
                                    show_variant_dialog(
                                        cohort_name=cohort_name,
                                        family_id=family_id,
                                        chrom=chrom,
                                        pos=pos,
                                        ref=ref,
                                        alt=alt,
                                        sample=sample_id,
                                        variant_data=dict(row_data),
                                        on_save_callback=lambda s: _refresh_file_data(),
                                    )
                                else:
                                    ui.notify(
                                        "Invalid variant format. Expected chr:pos:ref:alt",
                                        type="warning",
                                    )
                        except Exception as ex:
                            ui.notify(
                                f"Error parsing variant: {ex}",
                                type="warning",
                            )

                    # Restore saved sorting
                    saved_sorting = table_state.get("sorting", [])
                    if saved_sorting:
                        col_id = saved_sorting[0]["id"]
                        desc = saved_sorting[0].get("desc", False)
                        col_def = next(
                            (c for c in columns if c.get("id") == col_id),
                            {},
                        )
                        sort_field = col_def.get("sortField", col_id)
                        sort_type = col_def.get("sorting", "")
                        if sort_type == "genomic":
                            filtered_data.sort(
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

                            filtered_data.sort(key=_num_key, reverse=desc)
                        else:
                            filtered_data.sort(
                                key=lambda r: (
                                    r.get(sort_field) is None,
                                    r.get(sort_field, ""),
                                ),
                                reverse=desc,
                            )

                    dt_ref["ref"] = DataTable(
                        columns=columns,
                        rows=filtered_data,
                        row_key=variant_col if has_variant else "Variant",
                        pagination={"rowsPerPage": 50},
                        visible_columns=["actions"] + list(selected_cols["value"]),
                        on_row_action=on_row_action,
                        initial_sorting=saved_sorting,
                        initial_page=table_state.get("page", 0),
                        state_holder=table_state,
                    )

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
