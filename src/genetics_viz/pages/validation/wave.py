"""Wave validation page for individual sample - IGV viewer with bedgraph."""

import getpass
import json
from datetime import datetime

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.waves_loader import (
    get_wave_score_color,
    load_waves_validations,
    save_wave_validation,
)
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/wave/{sample_id}")
def wave_validation_page(sample_id: str) -> None:
    """Render the wave validation page with IGV.js viewer for a sample's bedgraph."""
    create_header()

    try:
        store = get_data_store()
        validation_file = store.data_dir / "validations" / "waves.tsv"

        # Check if bedgraph exists
        bedgraph_path = (
            store.data_dir
            / "samples"
            / sample_id
            / "sequences"
            / f"{sample_id}.by1000.bedgraph.gz"
        )

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.link("Waves Validation", "/validation/waves").classes(
                    "text-blue-600 hover:text-blue-800"
                )
                ui.label("/").classes("text-gray-400")
                ui.label(sample_id).classes("font-semibold")

            # Sample title
            ui.label(f"🌊 Wave Validation: {sample_id}").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            if not bedgraph_path.exists():
                with ui.card().classes("w-full p-6 bg-red-50"):
                    ui.label("⚠️ Bedgraph file not found").classes(
                        "text-xl font-semibold text-red-800"
                    )
                    ui.label(f"Expected path: {bedgraph_path}").classes(
                        "text-gray-600 text-sm font-mono"
                    )
                return

            # Panel 1: IGV.js viewer
            ui.label("Coverage Visualization").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.card().classes("w-full"):
                igv_container = (
                    ui.element("div").classes("w-full").style("height: 600px;")
                )
                igv_id = f"igv-{id(igv_container)}"
                igv_container._props["id"] = igv_id

                # Serve data files
                nicegui_app.add_static_files("/data", str(store.data_dir))

                # Add IGV.js library
                ui.add_head_html("""
                    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
                """)

                browser_var = f"igvBrowser_{igv_id.replace('-', '_')}"

                # Build IGV config with bedgraph track
                bedgraph_url = f"/data/samples/{sample_id}/sequences/{sample_id}.by1000.bedgraph.gz"
                bedgraph_index_url = bedgraph_url + ".tbi"

                igv_config = {
                    "genome": "hg38",
                    "locus": "chr1",  # Default view
                    "tracks": [
                        {
                            "name": f"{sample_id} Coverage",
                            "type": "wig",
                            "format": "bedgraph",
                            "url": bedgraph_url,
                            "indexURL": bedgraph_index_url,
                            "height": 400,
                            "color": "rgb(0, 0, 150)",
                            "autoscale": False,
                            "min": 0,
                            "max": 80,
                        }
                    ],
                }

                # Initialize IGV
                ui.run_javascript(
                    f"""
                    setTimeout(function() {{
                        var igvDiv = document.getElementById("{igv_id}");
                        if (igvDiv && typeof igv !== 'undefined') {{
                            igv.createBrowser(igvDiv, {json.dumps(igv_config)})
                                .then(function(browser) {{
                                    window.{browser_var} = browser;
                                    console.log("IGV browser created successfully for wave validation");
                                }})
                                .catch(function(error) {{
                                    console.error("Error creating IGV browser:", error);
                                }});
                        }} else {{
                            console.error("IGV container not found or igv not loaded");
                        }}
                    }}, 500);
                """
                )

                with ui.row().classes("items-center gap-4 p-4 bg-gray-50"):
                    ui.label("Bedgraph:").classes("font-semibold text-sm")
                    ui.label(
                        f"samples/{sample_id}/sequences/{sample_id}.by1000.bedgraph.gz"
                    ).classes("text-xs text-gray-600 font-mono")

            # Panel 2: Wave Validation Form
            ui.label("Add Validation").classes(
                "text-2xl font-semibold mb-4 text-blue-800 mt-6"
            )

            validation_history_container = ui.column().classes("w-full mb-4")

            def load_validation_history():
                """Load and display validation history."""
                waves_map = load_waves_validations(validation_file)
                validation_history_container.clear()

                waves = waves_map.get(sample_id, [])

                with validation_history_container:
                    if not waves:
                        ui.label("No validations recorded yet").classes(
                            "text-gray-500 text-sm italic"
                        )
                    else:
                        with ui.card().classes("w-full p-4"):
                            with ui.row().classes("gap-2 flex-wrap"):
                                for wave in waves:
                                    wave_label = {
                                        0: "good",
                                        1: "low wave",
                                        2: "medium wave",
                                        3: "high wave",
                                    }.get(wave, f"unknown ({wave})")

                                    color = get_wave_score_color(wave)
                                    ui.badge(wave_label, color=color).classes(
                                        "text-sm px-3 py-1"
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

                        ui.label("Wave Score:").classes("font-semibold ml-4")
                        wave_select = (
                            ui.select(
                                options={
                                    0: "0 - good",
                                    1: "1 - low wave",
                                    2: "2 - medium wave",
                                    3: "3 - high wave",
                                },
                                value=0,
                            )
                            .props("outlined dense")
                            .classes("w-56")
                        )

                        ui.button(
                            "Save Validation",
                            icon="save",
                            on_click=lambda: save_validation(),
                        ).props("color=blue").classes("ml-4")

                    def save_validation():
                        """Save a wave validation."""
                        user = user_input.value.strip()
                        wave = wave_select.value

                        if not user:
                            ui.notify("Please enter a username", type="warning")
                            return

                        if wave is None or not isinstance(wave, int):
                            ui.notify(
                                "Please select a valid wave score", type="warning"
                            )
                            return

                        # Generate timestamp
                        timestamp = datetime.now().isoformat()

                        try:
                            save_wave_validation(
                                validation_file,
                                sample_id,
                                user,
                                wave,
                                timestamp,
                            )

                            ui.notify(
                                f"Validation saved: {sample_id} - Wave {wave}",
                                type="positive",
                            )

                            # Reload validation history
                            load_validation_history()

                        except Exception as e:
                            ui.notify(f"Error saving validation: {e}", type="negative")
                            import traceback

                            print(traceback.format_exc())

            # Panel 3: Previous Validations
            ui.label("Previous Validations").classes(
                "text-2xl font-semibold mb-4 text-blue-800 mt-6"
            )

            load_validation_history()

    except Exception as e:
        ui.label(f"Error loading wave validation page: {e}").classes(
            "text-red-600 text-lg"
        )
        import traceback

        ui.label(traceback.format_exc()).classes("text-xs text-gray-600 font-mono")
