"""Variant visualization page with IGV.js viewer."""

import csv
import fcntl
import getpass
import json
import traceback
import urllib.parse
from datetime import datetime
from typing import Any, Dict, List, Optional

from nicegui import app as nicegui_app
from nicegui import context, ui

from genetics_viz.components.header import create_header
from genetics_viz.utils.data import get_data_store


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
    create_header()

    # Store the page client context for JavaScript execution
    page_client = context.client

    # Get all variant data from query params (URL encoded JSON)
    try:
        variant_data: Dict[str, Any] = (
            json.loads(urllib.parse.unquote(data)) if data else {}
        )
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
                additional_samples: Dict[str, List[str]] = {"value": []}

                # Get family members to enable add parents/family options
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

            refresh_additional_samples()

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

            # Panel 1: Variant details
            ui.label("Variant Details").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.card().classes("w-full mb-6"):
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

            # Panel 2: IGV.js viewer
            ui.label("Sequencing Data Viewer").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            cram_path = f"samples/{sample}/sequences/{sample}.GRCh38_GIABv3.cram"
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

                    def calculate_igv_height():
                        num_tracks = 1 + len(additional_samples["value"])
                        return max(500, num_tracks * 300 + 200)

                    igv_container = (
                        ui.element("div")
                        .classes("w-full")
                        .style(f"height: {calculate_igv_height()}px;")
                    )
                    igv_id = f"igv-{id(igv_container)}"
                    igv_container._props["id"] = igv_id

                    nicegui_app.add_static_files("/data", str(store.data_dir))

                    locus = f"{chrom}:{int(pos) - 100}-{int(pos) + 100}" if pos else ""

                    ui.add_head_html("""
                        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
                    """)

                    browser_var = f"igvBrowser_{igv_id.replace('-', '_')}"

                    def build_igv_config():
                        tracks = []
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
                        config = build_igv_config()
                        try:
                            new_height = calculate_igv_height()
                            igv_container.style(f"height: {new_height}px;")
                        except RuntimeError:
                            pass
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

                    with ui.row().classes("items-center gap-4 p-4 bg-gray-50"):
                        ui.label("CRAM:").classes("font-semibold text-sm")
                        ui.label(cram_path).classes("text-xs text-gray-600 font-mono")

            # Panel 3: Variant Validation
            ui.label("Variant Validation").classes(
                "text-2xl font-semibold mb-4 text-blue-800 mt-6"
            )

            with ui.card().classes("w-full p-6"):
                with ui.column().classes("gap-4"):
                    default_user = getpass.getuser()

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

                    validation_history_container = ui.column().classes("w-full mt-4")

                    def load_validations():
                        validation_file = store.data_dir / "validations" / "snvs.tsv"

                        validation_history_container.clear()

                        if not validation_file.exists():
                            with validation_history_container:
                                ui.label("No validations recorded yet").classes(
                                    "text-gray-500 text-sm italic"
                                )
                            return

                        try:
                            validations = []
                            with open(validation_file, "r") as f:
                                reader = csv.DictReader(f, delimiter="\t")
                                for row in reader:
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

                                    def delete_validation(timestamp: str):
                                        try:
                                            with open(validation_file, "r") as f:
                                                lines = f.readlines()

                                            filtered_lines = []
                                            for i, line in enumerate(lines):
                                                if i == 0:
                                                    filtered_lines.append(line)
                                                else:
                                                    parts = line.strip().split("\t")
                                                    if len(parts) >= 7:
                                                        line_fid = parts[0]
                                                        line_variant = parts[1]
                                                        line_timestamp = (
                                                            parts[6]
                                                            if len(parts) > 6
                                                            else ""
                                                        )
                                                        if not (
                                                            line_fid == family_id
                                                            and line_variant
                                                            == variant_key
                                                            and line_timestamp
                                                            == timestamp
                                                        ):
                                                            filtered_lines.append(line)

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
                                            load_validations()
                                        except Exception as e:
                                            ui.notify(
                                                f"Error deleting validation: {str(e)}",
                                                type="negative",
                                            )

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
                        validation_dir = store.data_dir / "validations"
                        validation_dir.mkdir(parents=True, exist_ok=True)

                        validation_file = validation_dir / "snvs.tsv"

                        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        row_data = [
                            family_id,
                            variant_key,
                            sample,
                            user_input.value or default_user,
                            inheritance_select.value,
                            validation_select.value,
                            timestamp,
                        ]

                        try:
                            file_exists = validation_file.exists()

                            with open(validation_file, "a") as f:
                                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                                try:
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

                                    f.write("\t".join(row_data) + "\n")
                                    f.flush()
                                finally:
                                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

                            ui.notify(
                                f"Validation saved successfully at {timestamp}",
                                type="positive",
                            )
                            load_validations()
                        except Exception as e:
                            ui.notify(
                                f"Error saving validation: {str(e)}", type="negative"
                            )

                    load_validations()

    except Exception as e:
        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
