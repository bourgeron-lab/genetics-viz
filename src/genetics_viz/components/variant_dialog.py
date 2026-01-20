"""Shared component for variant visualization in a dialog with IGV.js."""

import csv
import fcntl
import getpass
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.utils.data import get_data_store

# Path to the validation guide markdown file
VALIDATION_GUIDE_PATH = (
    Path(__file__).parent.parent.parent.parent
    / "documentation"
    / "snvs_validation_guide.md"
)


def _load_validation_guide() -> str:
    """Load the validation guide markdown content."""
    if VALIDATION_GUIDE_PATH.exists():
        return VALIDATION_GUIDE_PATH.read_text()
    return "Validation guide not found."


def _update_ignore_status(
    validation_file: Path,
    family_id: str,
    variant_key: str,
    sample: str,
    timestamp: str,
    ignore_value: str,
) -> bool:
    """Update the Ignore status for a specific validation row.

    Args:
        validation_file: Path to the snvs.tsv file
        family_id: Family ID to match
        variant_key: Variant key to match
        sample: Sample ID to match
        timestamp: Timestamp to match (unique identifier)
        ignore_value: New ignore value ("0" or "1")

    Returns:
        True if update was successful, False otherwise
    """
    if not validation_file.exists():
        return False

    # Read all rows
    rows = []
    fieldnames = []
    with open(validation_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for row in reader:
            rows.append(row)

    # Find and update the matching row
    updated = False
    for row in rows:
        if (
            row.get("FID") == family_id
            and row.get("Variant") == variant_key
            and row.get("Sample") == sample
            and row.get("Timestamp") == timestamp
        ):
            row["Ignore"] = ignore_value
            updated = True
            break

    if not updated:
        return False

    # Write back with file locking
    with open(validation_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)

    return True


def show_variant_dialog(
    cohort_name: str,
    family_id: str,
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    sample: str,
    variant_data: Dict[str, Any],
    on_save_callback: Optional[Callable[[str], None]] = None,
) -> None:
    """Show variant validation dialog with IGV viewer.

    Args:
        cohort_name: Cohort name
        family_id: Family ID
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
        sample: Sample ID
        variant_data: Additional variant data to display
        on_save_callback: Optional callback to call after saving validation with the validation status
    """
    store = get_data_store()
    variant_key = f"{chrom}:{pos}:{ref}:{alt}"

    with (
        ui.dialog().props("maximized") as dialog,
        ui.card().classes("w-full h-full"),
    ):
        with ui.column().classes("w-full h-full p-6"):
            # Header with close button
            with ui.row().classes("items-center justify-between w-full mb-4"):
                ui.label(f"🧬 {variant_key} - {sample}").classes(
                    "text-2xl font-bold text-blue-900"
                )
                ui.button(icon="close", on_click=lambda: dialog.close()).props(
                    "flat round"
                )

            # Get family members
            cohort = store.get_cohort(cohort_name)
            family_members: List[str] = []
            sample_parents: Dict[str, Optional[str]] = {
                "father": None,
                "mother": None,
            }

            if cohort:
                members_data = cohort.get_family_members(family_id)
                family_members = [m["Sample ID"] for m in members_data]
                # Find current sample's parents
                for member in members_data:
                    if member["Sample ID"] == sample:
                        sample_parents["father"] = member.get("Father")
                        sample_parents["mother"] = member.get("Mother")
                        break

            # Track additional samples
            additional_samples: Dict[str, List[str]] = {"value": []}

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

            def get_relationship_label(sample_id: str) -> str:
                """Get relationship label for a sample."""
                if sample_id == sample:
                    return "(carrier)"

                if sample_id == sample_parents["father"]:
                    return "(father)"
                if sample_id == sample_parents["mother"]:
                    return "(mother)"

                # Check if it's a sibling
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

            # Additional samples section with add menu
            with ui.row().classes("items-center gap-4 mb-2"):
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

                                def make_remove_handler(sid: str):
                                    return lambda: remove_sample(sid)

                                ui.button(
                                    icon="delete",
                                    on_click=make_remove_handler(add_sample_id),
                                ).props("flat dense size=xs color=red")

            def add_sample(sample_id: str):
                if (
                    sample_id
                    and sample_id not in additional_samples["value"]
                    and sample_id != sample
                ):
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

            def remove_sample(sample_id: str):
                if sample_id in additional_samples["value"]:
                    additional_samples["value"].remove(sample_id)
                    refresh_additional_samples()
                    refresh_igv()

            refresh_additional_samples()

            # Variant details card (collapsible, above IGV)
            ui.label("Variant Details").classes("text-xl font-semibold mb-2")

            with ui.card().classes("w-full mb-4"):
                with ui.column().classes("p-4 gap-4"):
                    primary_fields = [
                        "VEP_CANONICAL",
                        "VEP_Consequence",
                        "VEP_SYMBOL",
                        "VEP_HGVSp",
                        "VEP_LoF",
                        "fafmax_faf95_max_genomes",
                    ]

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

            # Build IGV tracks function
            def build_igv_tracks():
                tracks = []

                # Main sample track
                sample_cram = (
                    store.data_dir
                    / f"samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram"
                )
                if sample_cram.exists():
                    main_label = f"{sample} {get_relationship_label(sample)}".strip()
                    tracks.append(
                        {
                            "name": main_label,
                            "type": "alignment",
                            "format": "cram",
                            "url": f"/data/samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram",
                            "indexURL": f"/data/samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram.crai",
                            "height": 250,
                            "displayMode": "SQUISHED",
                        }
                    )

                # Additional samples tracks
                for add_sample_id in additional_samples["value"]:
                    add_sample_cram = (
                        store.data_dir
                        / f"samples/{add_sample_id}/sequences/{add_sample_id}.GRCh38_GIABv3.cram"
                    )
                    if add_sample_cram.exists():
                        track_label = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                        tracks.append(
                            {
                                "name": track_label,
                                "type": "alignment",
                                "format": "cram",
                                "url": f"/data/samples/{add_sample_id}/sequences/{add_sample_id}.GRCh38_GIABv3.cram",
                                "indexURL": f"/data/samples/{add_sample_id}/sequences/{add_sample_id}.GRCh38_GIABv3.cram.crai",
                                "height": 250,
                                "displayMode": "SQUISHED",
                            }
                        )

                return tracks

            # Calculate dynamic height: 200px base + (tracks * 250px per track) + 50px buffer
            def calculate_igv_height():
                num_tracks = 1 + len(additional_samples["value"])
                return 200 + (num_tracks * 250) + 50

            # IGV.js viewer - wrapped in card with dynamic height
            with ui.card().classes("w-full mb-6"):
                igv_container = (
                    ui.element("div")
                    .classes("w-full")
                    .style(f"height: {calculate_igv_height()}px;")
                )
                igv_id = f"igv-{id(igv_container)}"
                igv_container._props["id"] = igv_id

            browser_var = f"igvBrowser_{igv_id.replace('-', '_')}"

            # Calculate locus with padding
            try:
                pos_int = int(pos)
                padding = 50
                locus = f"{chrom}:{max(1, pos_int - padding)}-{pos_int + padding}"
            except Exception:
                locus = f"{chrom}:{pos}"

            def build_igv_config():
                return {
                    "genome": "hg38",
                    "locus": locus,
                    "tracks": build_igv_tracks(),
                }

            def refresh_igv():
                """Refresh IGV tracks when samples are added/removed."""
                config = build_igv_config()
                try:
                    new_height = calculate_igv_height()
                    igv_container.style(f"height: {new_height}px;")
                except RuntimeError:
                    pass
                ui.run_javascript(
                    f"""
                    if (window.{browser_var}) {{
                        window.{browser_var}.removeAllTracks();
                        const tracks = {json.dumps(config["tracks"])};
                        for (const track of tracks) {{
                            window.{browser_var}.loadTrack(track);
                        }}
                    }}
                """
                )

            igv_config = build_igv_config()

            # Validation section
            validation_file = store.data_dir / "validations" / "snvs.tsv"

            # Validation form - positioned first
            with ui.row().classes("items-center gap-2 mt-4 mb-2"):
                ui.label("Variant Validation").classes("text-xl font-semibold")

                # Info button for validation guide
                def show_validation_guide():
                    with (
                        ui.dialog() as guide_dialog,
                        ui.card().classes("w-full max-w-3xl"),
                    ):
                        with ui.row().classes(
                            "items-center justify-between w-full mb-4"
                        ):
                            ui.label("Validation Guide").classes("text-xl font-bold")
                            ui.button(
                                icon="close", on_click=lambda: guide_dialog.close()
                            ).props("flat round")
                        with ui.scroll_area().classes("w-full h-96"):
                            ui.markdown(_load_validation_guide())
                    guide_dialog.open()

                ui.button(icon="info", on_click=show_validation_guide).props(
                    "flat round color=blue"
                ).tooltip("Validation instructions")

            with ui.card().classes("w-full p-4 mb-4"):
                with ui.column().classes("gap-4"):
                    default_user = getpass.getuser()

                    with ui.row().classes("items-center gap-4 w-full flex-wrap"):
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
                                    "not paternal",
                                    "not maternal",
                                    "either",
                                    "homozygous",
                                ],
                                value="unknown",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                        ui.label("Validation:").classes("font-semibold ml-4")
                        validation_select = (
                            ui.select(
                                [
                                    "present",
                                    "absent",
                                    "uncertain",
                                    "different",
                                    "in phase MNV",
                                ],
                                value="present",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                    with ui.row().classes("items-center gap-4 w-full"):
                        ui.label("Comment:").classes("font-semibold")
                        comment_input = (
                            ui.input("Optional comment")
                            .props("outlined dense")
                            .classes("flex-grow")
                        )

                        ui.button(
                            "Save Validation",
                            icon="save",
                            on_click=lambda: save_validation(),
                        ).props("color=blue")

                    def save_validation():
                        """Save a validation."""
                        user = user_input.value.strip()
                        validation_status = validation_select.value
                        inheritance = inheritance_select.value or ""
                        comment = (
                            comment_input.value.strip() if comment_input.value else ""
                        )

                        if not user:
                            ui.notify("Please enter a username", type="warning")
                            return

                        timestamp = datetime.now().isoformat()

                        try:
                            # Ensure directory exists
                            validation_file.parent.mkdir(parents=True, exist_ok=True)

                            # Check if file exists
                            file_exists = validation_file.exists()

                            # Append validation
                            with open(validation_file, "a") as f:
                                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                                try:
                                    # Only write header if file is new
                                    if not file_exists:
                                        f.write(
                                            "FID\tVariant\tSample\tUser\tInheritance\tValidation\tComment\tIgnore\tTimestamp\n"
                                        )
                                    f.write(
                                        f"{family_id}\t{variant_key}\t{sample}\t{user}\t{inheritance}\t{validation_status}\t{comment}\t0\t{timestamp}\n"
                                    )
                                finally:
                                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

                            ui.notify(
                                f"Validation saved: {validation_status}",
                                type="positive",
                            )

                            # Reload validation history
                            load_validation_history()

                            # Close dialog and call callback
                            dialog.close()
                            if on_save_callback and validation_status:
                                on_save_callback(validation_status)

                        except Exception as e:
                            ui.notify(f"Error saving validation: {e}", type="negative")
                            import traceback

                            print(traceback.format_exc())

            # Validation history container - positioned after form
            validation_history_container = ui.column().classes("w-full mb-4")

            def load_validation_history():
                """Load and display validation history."""
                validations = []

                if validation_file.exists():
                    with open(validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            if (
                                row.get("FID") == family_id
                                and row.get("Variant") == variant_key
                                and row.get("Sample") == sample
                            ):
                                validations.append(row)

                validation_history_container.clear()

                with validation_history_container:
                    ui.label("Previous validations:").classes("font-semibold mb-2")
                    if not validations:
                        ui.label("No validations recorded yet").classes(
                            "text-gray-500 text-sm italic"
                        )
                    else:
                        with ui.card().classes("w-full p-2"):
                            with ui.column().classes("gap-2"):
                                for validation in validations:
                                    val_status = validation.get("Validation", "")
                                    inheritance = validation.get("Inheritance", "")
                                    user = validation.get("User", "")
                                    timestamp = validation.get("Timestamp", "")
                                    comment = validation.get("Comment", "")
                                    is_ignored = validation.get("Ignore", "0") == "1"

                                    # Format timestamp
                                    try:
                                        dt = datetime.fromisoformat(timestamp)
                                        formatted_time = dt.strftime(
                                            "%Y-%m-%d %H:%M:%S"
                                        )
                                    except Exception:
                                        formatted_time = timestamp

                                    # Color based on status
                                    if (
                                        val_status == "present"
                                        or val_status == "in phase MNV"
                                    ):
                                        color = "green"
                                        icon = "check_circle"
                                    elif val_status == "absent":
                                        color = "red"
                                        icon = "cancel"
                                    else:
                                        color = "orange"
                                        icon = "help"

                                    # Apply styling for ignored rows
                                    row_classes = "items-center gap-2 w-full"
                                    if is_ignored:
                                        row_classes += " opacity-50"

                                    with ui.row().classes(row_classes):
                                        ui.icon(icon, color=color).classes("text-sm")
                                        label_text = val_status
                                        if inheritance:
                                            label_text += f" ({inheritance})"
                                        ui.label(label_text).classes("font-semibold")
                                        ui.label(f"- {user}").classes(
                                            "text-sm text-gray-600"
                                        )
                                        ui.label(formatted_time).classes(
                                            "text-xs text-gray-500"
                                        )

                                        # Show comment if present
                                        if comment:
                                            ui.label(f'"{comment}"').classes(
                                                "text-xs text-gray-500 italic"
                                            )

                                        # Spacer to push switch to the right
                                        ui.space()

                                        # Ignore switch
                                        def make_ignore_handler(
                                            ts: str, current_ignored: bool
                                        ):
                                            def handler(e):
                                                new_value = "1" if e.value else "0"
                                                success = _update_ignore_status(
                                                    validation_file,
                                                    family_id,
                                                    variant_key,
                                                    sample,
                                                    ts,
                                                    new_value,
                                                )
                                                if success:
                                                    action = (
                                                        "ignored"
                                                        if e.value
                                                        else "restored"
                                                    )
                                                    ui.notify(
                                                        f"Validation {action}",
                                                        type="info",
                                                    )
                                                    load_validation_history()
                                                    # Trigger refresh callback to update table
                                                    if on_save_callback:
                                                        on_save_callback(
                                                            "ignored"
                                                            if e.value
                                                            else "restored"
                                                        )
                                                else:
                                                    ui.notify(
                                                        "Failed to update ignore status",
                                                        type="negative",
                                                    )

                                            return handler

                                        ignore_switch = ui.switch(
                                            "Ignore",
                                            value=is_ignored,
                                            on_change=make_ignore_handler(
                                                timestamp, is_ignored
                                            ),
                                        ).classes("text-xs")
                                        if is_ignored:
                                            ignore_switch.props("color=grey")

            # Display validation history below the form
            load_validation_history()

        dialog.open()

        # Initialize IGV after dialog is open
        def init_igv():
            ui.run_javascript(
                f"""
                var igvDiv = document.getElementById("{igv_id}");
                console.log("Trying to initialize IGV, div found:", igvDiv !== null);
                console.log("IGV library loaded:", typeof igv !== 'undefined');
                if (igvDiv && typeof igv !== 'undefined') {{
                    igv.createBrowser(igvDiv, {json.dumps(igv_config)})
                        .then(function(browser) {{
                            window.{browser_var} = browser;
                            console.log("IGV browser created in dialog");
                        }})
                        .catch(function(error) {{
                            console.error("Error creating IGV browser:", error);
                        }});
                }} else {{
                    console.error("IGV container not found or igv not loaded");
                }}
            """
            )

        ui.timer(0.5, init_igv, once=True)
