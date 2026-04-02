"""Standalone sample visualization dialog with IGV.js bedgraph and CRAM tracks."""

import json

from nicegui import ui

from genetics_viz.utils.data import get_data_store, get_static_prefix
from genetics_viz.utils.data_availability import check_sample_availability
from genetics_viz.utils.sharding import get_sample_url


def show_sample_dialog(sample_id: str) -> None:
    """Open a fullscreen dialog to visualize a sample's bedgraph and CRAM data."""
    # Ensure IGV.js is loaded at page level (must be before dialog opens)
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>'
    )

    store = get_data_store()
    avail = check_sample_availability(store.data_dir, sample_id)

    with (
        ui.dialog().props("maximized") as dialog,
        ui.card().classes("w-full h-full"),
    ):
        with ui.column().classes("w-full h-full p-6"):
            # Header
            with ui.row().classes("items-center justify-between w-full mb-4"):
                with ui.row().classes("items-center gap-3"):
                    ui.label(f"Sample: {sample_id}").classes(
                        "text-2xl font-bold text-blue-900"
                    )
                    # Availability badges
                    _BADGES = [
                        ("bedgraph", "CRAM", avail["cram"]),
                        ("bedgraph", "Bedgraph", avail["bedgraph"]),
                        ("vaf", "VAF", avail["vaf_bedgraph"]),
                        ("dv", "DeepVariant", avail["deepvariant"]),
                        ("svs", "SVs", avail["svs"]),
                    ]
                    for _, label, present in _BADGES:
                        color = "green" if present else "grey"
                        ui.badge(label, color=color).props("outline")

                ui.button(icon="close", on_click=lambda: dialog.close()).props(
                    "flat round"
                )

            # Locus input
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.label("Locus:").classes("font-semibold")
                locus_input = (
                    ui.input(placeholder="e.g. chr1:1000000-2000000")
                    .props("outlined dense")
                    .classes("w-96")
                )

                async def navigate_to_locus():
                    locus = locus_input.value.strip() if locus_input.value else ""
                    if not locus:
                        ui.notify("Enter a locus first", type="warning")
                        return
                    ui.run_javascript(
                        f"""
                        if (window.sampleBrowser) {{
                            window.sampleBrowser.search("{locus}");
                        }}
                        """
                    )

                ui.button("Go", on_click=navigate_to_locus).props("color=blue dense")
                locus_input.on("keydown.enter", navigate_to_locus)

            # IGV container
            (
                ui.element("div")
                .props('id="sample-igv-div"')
                .classes("w-full border border-gray-300 rounded-lg")
                .style("height: 700px")
            )

            # Build tracks
            sample_url = f"{get_static_prefix()}/{get_sample_url(store.data_dir, sample_id)}/sequences"
            tracks = []

            if avail["bedgraph"]:
                tracks.append(
                    {
                        "name": f"{sample_id} Coverage",
                        "type": "wig",
                        "format": "bedgraph",
                        "url": f"{sample_url}/{sample_id}.by1000.bedgraph.gz",
                        "indexURL": f"{sample_url}/{sample_id}.by1000.bedgraph.gz.tbi",
                        "height": 150,
                        "color": "rgb(0, 0, 150)",
                        "autoscale": True,
                    }
                )

            if avail["vaf_bedgraph"]:
                tracks.append(
                    {
                        "name": f"{sample_id} VAF",
                        "type": "wig",
                        "format": "bedgraph",
                        "url": f"{sample_url}/{sample_id}.vaf.bedgraph.gz",
                        "indexURL": f"{sample_url}/{sample_id}.vaf.bedgraph.gz.tbi",
                        "height": 100,
                        "color": "rgb(150, 0, 150)",
                        "autoscale": False,
                        "min": 0,
                        "max": 1,
                        "graphType": "points",
                    }
                )

            if avail["cram"]:
                tracks.append(
                    {
                        "name": f"{sample_id} Alignment",
                        "type": "alignment",
                        "format": "cram",
                        "url": f"{sample_url}/{sample_id}.GRCh38_GIABv3.cram",
                        "indexURL": f"{sample_url}/{sample_id}.GRCh38_GIABv3.cram.crai",
                        "height": 300,
                        "displayMode": "SQUISHED",
                    }
                )

            if not tracks:
                ui.label("No visualization data available for this sample.").classes(
                    "text-lg text-gray-500 italic"
                )
            else:
                igv_config = {
                    "genome": "hg38",
                    "tracks": tracks,
                }

                ui.run_javascript(
                    f"""
                    (function waitForIgv() {{
                        var igvDiv = document.getElementById("sample-igv-div");
                        if (igvDiv && typeof igv !== 'undefined') {{
                            igv.createBrowser(igvDiv, {json.dumps(igv_config)})
                                .then(function(browser) {{
                                    window.sampleBrowser = browser;
                                }})
                                .catch(function(error) {{
                                    console.error("Error creating IGV browser:", error);
                                }});
                        }} else {{
                            setTimeout(waitForIgv, 200);
                        }}
                    }})();
                    """
                )

    dialog.open()
