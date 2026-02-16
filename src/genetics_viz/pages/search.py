"""Search page for cohort-wide variant search."""

import asyncio
import re
from typing import Any, Dict, List

import polars as pl
from nicegui import app as nicegui_app
from nicegui import context, ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.view_presets import VIEW_PRESETS, select_preset_for_config
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    get_column_group,
    get_column_sorting,
    get_display_label,
    get_dropped_columns,
    get_schema_overrides,
    reorder_columns_by_group,
)
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.score_colors import get_score_color
from genetics_viz.utils.clinvar import (
    CLINVAR_COLORS,
    format_clinvar_display,
    get_clinvar_color,
)
from genetics_viz.utils.cytobands import (
    CHROM_ORDER,
    CHROM_SIZES_MB,
    CYTOBANDS,
    GIESTAIN_COLORS,
    VALIDATION_COLORS,
    norm_chrom,
)
from genetics_viz.utils.vep import (
    VEP_CONSEQUENCES,
    format_consequence_display,
    get_consequence_color,
    get_consequence_impact,
    get_highest_consequence_term,
)





def parse_locus_query(query: str) -> Dict[str, Any]:
    """Parse locus query into filter parameters.

    Returns dict with 'type' and relevant filter parameters.
    """
    query = query.strip()

    # chr1:10000:A:GC - exact variant
    variant_pattern = r"^(chr)?(\w+):(\d+):([ACGT]+):([ACGT]+)$"
    match = re.match(variant_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "exact_variant",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "pos": int(match.group(3)),
            "ref": match.group(4).upper(),
            "alt": match.group(5).upper(),
        }

    # chr1:10000-10100 - range
    range_pattern = r"^(chr)?(\w+):(\d+)-(\d+)$"
    match = re.match(range_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "range",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "start": int(match.group(3)),
            "end": int(match.group(4)),
        }

    # chr1:10000 - exact position
    pos_pattern = r"^(chr)?(\w+):(\d+)$"
    match = re.match(pos_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "exact_position",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "pos": int(match.group(3)),
        }

    # ENSG00000164099 - gene ID
    if re.match(r"^ENSG\d+$", query, re.IGNORECASE):
        return {
            "type": "gene_id",
            "gene_id": query.upper(),
        }

    # SHANK* - wildcard gene
    if "*" in query:
        return {
            "type": "gene_wildcard",
            "pattern": query.replace("*", "").upper(),
        }

    # SHANK3 - exact gene
    return {
        "type": "gene_name",
        "gene_name": query.upper(),
    }


def filter_dataframe(df: pl.DataFrame, query_params: Dict[str, Any]) -> pl.DataFrame:
    """Filter dataframe based on parsed query parameters."""
    query_type = query_params["type"]

    if query_type == "exact_variant":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") == query_params["pos"])
            & (pl.col("REF") == query_params["ref"])
            & (pl.col("ALT") == query_params["alt"])
        )

    elif query_type == "range":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") >= query_params["start"])
            & (pl.col("POS") <= query_params["end"])
        )

    elif query_type == "exact_position":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") == query_params["pos"])
        )

    elif query_type == "gene_id":
        # VEP_Gene can contain multiple genes separated by &
        gene_id = query_params["gene_id"]
        return df.filter(pl.col("VEP_Gene").str.to_uppercase().str.contains(gene_id))

    elif query_type == "gene_wildcard":
        # VEP_SYMBOL can contain multiple symbols
        pattern = query_params["pattern"]
        return df.filter(pl.col("VEP_SYMBOL").str.to_uppercase().str.contains(pattern))

    elif query_type == "gene_name":
        # Exact match in VEP_SYMBOL (case-insensitive, as part of the field)
        gene_name = query_params["gene_name"]
        return df.filter(
            pl.col("VEP_SYMBOL").str.to_uppercase().str.contains(gene_name)
        )

    return df


# Values that represent "unknown/missing" in pedigree files
_PED_MISSING = {"", "0", "-9"}


def _pedigree_data_from_cohort(cohort_name: str) -> Dict[str, Dict[str, str]]:
    """Build pedigree lookup from the already-parsed Cohort object.

    Uses the DataStore's Cohort (which handles all header formats robustly)
    instead of re-parsing the pedigree file.

    Returns dict mapping sample_id -> {FID, Father, Mother, Sex, Phenotype}.
    """
    store = get_data_store()
    cohort = store.get_cohort(cohort_name)
    if cohort is None:
        return {}

    pedigree_data: Dict[str, Dict[str, str]] = {}
    for family in cohort.families.values():
        for sample in family.samples:
            pedigree_data[sample.sample_id] = {
                "FID": sample.family_id,
                "Father": sample.father_id or "",
                "Mother": sample.mother_id or "",
                "Sex": sample.sex or "",
                "Phenotype": sample.phenotype or "",
            }
    return pedigree_data


@ui.page("/search/{cohort_name}")
def search_cohort_page(cohort_name: str) -> None:
    """Search page for cohort-wide variant search."""
    create_header(cohort_name)

    # Add IGV.js library
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

        # Scan for wombat files
        wombat_dir = store.data_dir / "cohorts" / cohort_name / "wombat"

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            ui.label(f"🔍 Search: {cohort_name}").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            if not wombat_dir.exists():
                ui.label(f"No wombat directory found at: {wombat_dir}").classes(
                    "text-red-500 text-lg"
                )
                return

            # Scan for wombat files matching pattern
            pattern = re.compile(
                rf"{re.escape(cohort_name)}\.rare\.([^.]+)\.(.+?)\.results\.tsv$"
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
                            "display_name": wombat_config,
                        }
                    )

            if not wombat_files:
                ui.label(
                    f"No wombat files found matching pattern in: {wombat_dir}"
                ).classes("text-gray-500 text-lg italic")
                return

            # Load pedigree data from the already-parsed Cohort object
            pedigree_data = _pedigree_data_from_cohort(cohort_name)
            sample_to_family = {
                sid: ped["FID"] for sid, ped in pedigree_data.items()
            }

            # Derive unique sex and phenotype values for individual filters
            # Exclude missing-value sentinels ("", "0", "-9") from the option lists
            available_sex_values = sorted(
                {v.get("Sex", "") for v in pedigree_data.values()}
                - _PED_MISSING
            )
            available_phenotype_values = sorted(
                {v.get("Phenotype", "") for v in pedigree_data.values()}
                - _PED_MISSING
            )

            # Load genesets from params/genesets
            genesets_dir = store.data_dir / "params" / "genesets"
            available_genesets = {}
            if genesets_dir.exists():
                for geneset_file in genesets_dir.glob("*.tsv"):
                    geneset_name = geneset_file.stem
                    genes = set()
                    with open(geneset_file, "r") as f:
                        # Skip header line
                        next(f, None)
                        for line in f:
                            gene = line.strip()
                            if gene:
                                genes.add(gene.upper())
                    if genes:
                        available_genesets[geneset_name] = genes

            # State for selected genesets and impacts
            selected_genesets: Dict[str, List[str]] = {"value": []}
            selected_impacts_search: Dict[str, List[str]] = {"value": []}

            # Validation filter state (all statuses selected by default)
            selected_validations: Dict[str, List[str]] = {
                "value": ["present", "absent", "uncertain", "conflicting", "TODO"]
            }

            # Exclude filters state
            filter_exclude_lcr: Dict[str, bool] = {"value": True}
            filter_exclude_gnomad: Dict[str, bool] = {"value": True}
            filter_exclude_gnomad_wgs: Dict[str, bool] = {"value": False}

            # Individual filters state
            filter_sex: Dict[str, List[str]] = {"value": []}
            filter_phenotype: Dict[str, List[str]] = {"value": []}
            filter_has_parents: Dict[str, bool] = {"value": False}

            # Button references for visual indicators
            geneset_button_ref: Dict[str, Any] = {"button": None}
            impact_button_ref: Dict[str, Any] = {"button": None}

            # Search panel
            with ui.card().classes("w-full p-2 mb-2").props("flat bordered"):
                # Header row: title + search button (always visible)
                with ui.row().classes("items-center justify-between w-full"):
                    ui.label("Search Parameters").classes("text-lg font-semibold")
                    search_button = ui.button(
                        "Search",
                        icon="search",
                    ).props("color=blue dense")

                # Tabs
                with ui.tabs().classes("w-full").props("dense") as search_tabs:
                    variants_tab = ui.tab("Variants", icon="biotech")
                    individuals_tab = ui.tab("Individuals", icon="people")

                with ui.tab_panels(search_tabs, value=variants_tab).classes("w-full"):
                  with ui.tab_panel(variants_tab).classes("p-0 pt-1"):
                    with ui.row().classes("items-center gap-2 w-full flex-wrap"):
                        # Source dropdown
                        source_select = (
                            ui.select(
                                options=[wf["display_name"] for wf in wombat_files],
                                label="Source",
                                value=wombat_files[0]["display_name"]
                                if wombat_files
                                else None,
                            )
                            .props("outlined dense")
                            .classes("w-64")
                        )

                        # Locus input (Enter key handler set later after function definition)
                        locus_input = (
                            ui.input(
                                label="Locus (optional)",
                                placeholder="chr1:10000-10100, SHANK3, ENSG...",
                                on_change=lambda: None,  # Placeholder to enable events
                            )
                            .props("outlined dense")
                            .classes("flex-grow")
                        )

                        # Geneset filter menu
                        if available_genesets:
                            geneset_btn = (
                                ui.button("Genesets", icon="list")
                                .props(
                                    ("outline"
                                    if not selected_genesets["value"]
                                    else "unelevated color=green")
                                    + " dense"
                                )
                            )
                            geneset_button_ref["button"] = geneset_btn
                            with geneset_btn:
                                with ui.menu():
                                    ui.label("Select Genesets:").classes(
                                        "px-4 py-2 font-semibold text-sm"
                                    )
                                    ui.separator()

                                    with ui.column().classes("p-2"):
                                        geneset_checkboxes: Dict[str, Any] = {}

                                        with ui.row().classes("gap-2 mb-2"):

                                            def select_all_genesets():
                                                selected_genesets["value"] = list(
                                                    available_genesets.keys()
                                                )
                                                for cb in geneset_checkboxes.values():
                                                    cb.value = True
                                                if geneset_button_ref["button"]:
                                                    geneset_button_ref["button"].props(
                                                        remove="outline",
                                                        add="unelevated color=green",
                                                    )
                                                    geneset_button_ref["button"].update()

                                            def select_no_genesets():
                                                selected_genesets["value"] = []
                                                for cb in geneset_checkboxes.values():
                                                    cb.value = False
                                                if geneset_button_ref["button"]:
                                                    geneset_button_ref["button"].props(
                                                        remove="unelevated color=green",
                                                        add="outline",
                                                    )
                                                    geneset_button_ref["button"].update()

                                            ui.button(
                                                "All", on_click=select_all_genesets
                                            ).props("size=sm flat dense").classes("text-xs")
                                            ui.button(
                                                "None", on_click=select_no_genesets
                                            ).props("size=sm flat dense").classes("text-xs")

                                        ui.separator()

                                        for geneset_name in sorted(
                                            available_genesets.keys()
                                        ):

                                            def make_geneset_handler(gs_name):
                                                def handler(e):
                                                    if e.value:
                                                        if (
                                                            gs_name
                                                            not in selected_genesets[
                                                                "value"
                                                            ]
                                                        ):
                                                            selected_genesets[
                                                                "value"
                                                            ].append(gs_name)
                                                    else:
                                                        if (
                                                            gs_name
                                                            in selected_genesets["value"]
                                                        ):
                                                            selected_genesets[
                                                                "value"
                                                            ].remove(gs_name)
                                                    # Update button visual state
                                                    if geneset_button_ref["button"]:
                                                        if selected_genesets["value"]:
                                                            geneset_button_ref[
                                                                "button"
                                                            ].props(
                                                                remove="outline",
                                                                add="unelevated color=green",
                                                            )
                                                        else:
                                                            geneset_button_ref[
                                                                "button"
                                                            ].props(
                                                                remove="unelevated color=green",
                                                                add="outline",
                                                            )
                                                        geneset_button_ref[
                                                            "button"
                                                        ].update()

                                                return handler

                                            geneset_checkboxes[geneset_name] = ui.checkbox(
                                                f"{geneset_name} ({len(available_genesets[geneset_name])} genes)",
                                                value=False,
                                                on_change=make_geneset_handler(
                                                    geneset_name
                                                ),
                                            ).classes("text-sm")

                        # Impact filter menu
                        impact_btn = (
                            ui.button("Impacts", icon="filter_list")
                            .props("outline dense")
                        )
                        impact_button_ref["button"] = impact_btn
                        with impact_btn:
                            with ui.menu():
                                ui.label("Select Impact Types:").classes(
                                    "px-4 py-2 font-semibold text-sm"
                                )
                                ui.separator()

                                with ui.column().classes("p-2"):
                                    impact_checkboxes_search: Dict[str, Any] = {}

                                    with ui.row().classes("gap-2 mb-2 flex-wrap"):

                                        def select_all_impacts_search():
                                            selected_impacts_search["value"] = list(
                                                VEP_CONSEQUENCES.keys()
                                            )
                                            for cb in impact_checkboxes_search.values():
                                                cb.value = True
                                            if impact_button_ref["button"]:
                                                impact_button_ref["button"].props(
                                                    remove="unelevated color=orange",
                                                    add="outline",
                                                )
                                                impact_button_ref["button"].update()

                                        def select_none_impacts_search():
                                            selected_impacts_search["value"] = []
                                            for cb in impact_checkboxes_search.values():
                                                cb.value = False
                                            if impact_button_ref["button"]:
                                                impact_button_ref["button"].props(
                                                    remove="outline",
                                                    add="unelevated color=orange",
                                                )
                                                impact_button_ref["button"].update()

                                        def select_by_impact_level(level: str):
                                            selected = [
                                                cons
                                                for cons, (
                                                    imp,
                                                    _,
                                                ) in VEP_CONSEQUENCES.items()
                                                if imp == level
                                            ]
                                            selected_impacts_search["value"] = selected
                                            for (
                                                impact,
                                                cb,
                                            ) in impact_checkboxes_search.items():
                                                cb.value = impact in selected
                                            if impact_button_ref["button"]:
                                                impact_button_ref["button"].props(
                                                    remove="outline",
                                                    add="unelevated color=orange",
                                                )
                                                impact_button_ref["button"].update()

                                        ui.button(
                                            "All", on_click=select_all_impacts_search
                                        ).props("size=sm flat dense").classes("text-xs")
                                        ui.button(
                                            "None", on_click=select_none_impacts_search
                                        ).props("size=sm flat dense").classes("text-xs")
                                        ui.button(
                                            "HIGH",
                                            on_click=lambda: select_by_impact_level("HIGH"),
                                        ).props("size=sm flat dense color=red").classes(
                                            "text-xs"
                                        )
                                        ui.button(
                                            "MODERATE",
                                            on_click=lambda: select_by_impact_level(
                                                "MODERATE"
                                            ),
                                        ).props("size=sm flat dense color=orange").classes(
                                            "text-xs"
                                        )
                                        ui.button(
                                            "LOW",
                                            on_click=lambda: select_by_impact_level("LOW"),
                                        ).props(
                                            "size=sm flat dense color=yellow-8"
                                        ).classes("text-xs")
                                        ui.button(
                                            "MODIFIER",
                                            on_click=lambda: select_by_impact_level(
                                                "MODIFIER"
                                            ),
                                        ).props("size=sm flat dense color=grey").classes(
                                            "text-xs"
                                        )

                                    ui.separator()

                                    # Pre-populate with all VEP consequences
                                    with ui.column().classes("gap-1"):

                                        def make_impact_handler_search(cons):
                                            def handler(e):
                                                if e.value:
                                                    if (
                                                        cons
                                                        not in selected_impacts_search[
                                                            "value"
                                                        ]
                                                    ):
                                                        selected_impacts_search[
                                                            "value"
                                                        ].append(cons)
                                                else:
                                                    if (
                                                        cons
                                                        in selected_impacts_search["value"]
                                                    ):
                                                        selected_impacts_search[
                                                            "value"
                                                        ].remove(cons)
                                                # Update button visual state
                                                if impact_button_ref["button"]:
                                                    if len(
                                                        selected_impacts_search["value"]
                                                    ) == len(VEP_CONSEQUENCES):
                                                        impact_button_ref["button"].props(
                                                            remove="unelevated color=orange",
                                                            add="outline",
                                                        )
                                                    else:
                                                        impact_button_ref["button"].props(
                                                            remove="outline",
                                                            add="unelevated color=orange",
                                                        )
                                                    impact_button_ref["button"].update()

                                            return handler

                                        # Group by impact level for better organization
                                        for impact_level in [
                                            "HIGH",
                                            "MODERATE",
                                            "LOW",
                                            "MODIFIER",
                                        ]:
                                            consequences = [
                                                cons
                                                for cons, (
                                                    imp,
                                                    _,
                                                ) in VEP_CONSEQUENCES.items()
                                                if imp == impact_level
                                            ]
                                            if consequences:
                                                ui.label(f"{impact_level}:").classes(
                                                    "text-xs font-bold text-gray-600 mt-2"
                                                )
                                                for cons in sorted(consequences):
                                                    impact_checkboxes_search[cons] = (
                                                        ui.checkbox(
                                                            format_consequence_display(
                                                                cons
                                                            ),
                                                            value=True,
                                                            on_change=make_impact_handler_search(
                                                                cons
                                                            ),
                                                        ).classes("text-sm")
                                                    )

                                        # Initialize with all selected
                                        selected_impacts_search["value"] = list(
                                            VEP_CONSEQUENCES.keys()
                                        )

                        # Validation filter
                        create_validation_filter_menu(
                            all_statuses=["present", "absent", "uncertain", "conflicting", "TODO"],
                            filter_state=selected_validations,
                            on_change=lambda: None,  # No action needed during search parameter setup
                            label="Validation",
                            button_classes="",
                            button_size="",
                            button_props="dense",
                        )

                    # Exclude filters row
                    with ui.row().classes("items-center gap-4 w-full flex-wrap mt-1"):
                        ui.checkbox(
                            "Exclude LCR",
                            value=filter_exclude_lcr["value"],
                            on_change=lambda e: filter_exclude_lcr.update(
                                {"value": e.value}
                            ),
                        ).props("dense")
                        ui.checkbox(
                            "Exclude gnomAD filtered",
                            value=filter_exclude_gnomad["value"],
                            on_change=lambda e: filter_exclude_gnomad.update(
                                {"value": e.value}
                            ),
                        ).props("dense")
                        ui.checkbox(
                            "Exclude gnomAD WGS",
                            value=filter_exclude_gnomad_wgs["value"],
                            on_change=lambda e: filter_exclude_gnomad_wgs.update(
                                {"value": e.value}
                            ),
                        ).props("dense")

                    # Help text
                    with ui.expansion("Query Examples", icon="help").classes("mt-1").props("dense"):
                        ui.markdown("""
- **chr1:10000-10100** - All variants in range [10000, 10100] on chr1
- **chr1:10000** - Exact position 10000 on chr1
- **chr1:10000:A:GC** - Exact variant with REF=A and ALT=GC
- **SHANK3** - All variants in SHANK3 gene
- **SHANK*** - All variants in genes starting with SHANK
- **ENSG00000164099** - All variants in gene with this Ensembl ID
                        """)

                  with ui.tab_panel(individuals_tab).classes("p-0 pt-1"):
                    with ui.row().classes("items-center gap-2 flex-wrap"):
                        # Sex filter (always shown)
                        ui.select(
                            options=available_sex_values,
                            label="Sex",
                            value=filter_sex["value"],
                            multiple=True,
                            on_change=lambda e: filter_sex.update({"value": e.value or []}),
                        ).props("outlined dense use-chips").classes("w-48")

                        # Phenotype filter (always shown)
                        ui.select(
                            options=available_phenotype_values,
                            label="Phenotype",
                            value=filter_phenotype["value"],
                            multiple=True,
                            on_change=lambda e: filter_phenotype.update({"value": e.value or []}),
                        ).props("outlined dense use-chips").classes("w-48")

                        # Has parents checkbox
                        ui.checkbox(
                            "Only samples with both parents",
                            value=filter_has_parents["value"],
                            on_change=lambda e: filter_has_parents.update({"value": e.value}),
                        ).props("dense")

            # Results container
            results_container = ui.column().classes("w-full")

            # Capture client context for callbacks
            page_client = context.client

            async def perform_search():
                """Execute the search and display results."""
                results_container.clear()

                # Show progress indicator while loading
                with results_container:
                    with ui.column().classes("items-center gap-4 justify-center py-8 w-full"):
                        progress = ui.circular_progress(
                            min=0, max=100, value=0, size="xl", color="blue"
                        )
                        status_label = ui.label("Starting search...").classes(
                            "text-lg text-gray-600"
                        )

                selected_source = source_select.value
                locus_query = locus_input.value

                # Allow empty locus if genesets are selected
                if (
                    not locus_query or not locus_query.strip()
                ) and not selected_genesets["value"]:
                    results_container.clear()
                    with results_container:
                        ui.label(
                            "Please enter a locus or select at least one geneset"
                        ).classes("text-orange-600")
                    return

                # Find selected file
                selected_file = None
                for wf in wombat_files:
                    if wf["display_name"] == selected_source:
                        selected_file = wf
                        break

                if not selected_file:
                    results_container.clear()
                    with results_container:
                        ui.label("Selected source not found").classes("text-red-500")
                    return

                # Update progress: validation complete
                progress.set_value(10)
                status_label.set_text("Loading data...")
                await asyncio.sleep(0)

                try:
                    # Define function to load and process dataframe (runs in background thread)
                    def load_and_group_data():
                        # Load dataframe
                        df = pl.read_csv(
                            selected_file["file_path"],
                            separator="\t",
                            infer_schema_length=10000,
                            schema_overrides=get_schema_overrides(),
                            null_values=[".", ""],
                        )
                        _drop = get_dropped_columns() & set(df.columns)
                        if _drop:
                            df = df.drop(list(_drop))

                        # Group by variant and sample, aggregating other columns
                        grouping_cols = ["#CHROM", "POS", "REF", "ALT", "sample"]

                        # Identify columns to aggregate
                        agg_cols = [col for col in df.columns if col not in grouping_cols]

                        # Create aggregation expressions
                        agg_exprs = [pl.len().alias("n_grouped")]  # Count rows grouped
                        for col in agg_cols:
                            # Aggregate as comma-separated unique values, excluding empty/null/'.'
                            agg_exprs.append(
                                pl.col(col)
                                .cast(pl.Utf8)
                                .filter(
                                    (pl.col(col).is_not_null())
                                    & (pl.col(col).cast(pl.Utf8) != "")
                                    & (pl.col(col).cast(pl.Utf8) != ".")
                                )
                                .unique()
                                .str.join(",")
                                .alias(col)
                            )

                        # Group and aggregate
                        df = df.group_by(grouping_cols, maintain_order=True).agg(agg_exprs)
                        return df

                    # Run in background thread to avoid blocking
                    df = await asyncio.to_thread(load_and_group_data)

                    # Update progress: CSV loaded and grouped
                    progress.set_value(35)
                    status_label.set_text("Filtering data...")
                    await asyncio.sleep(0)

                    # Apply locus filter if provided
                    if locus_query and locus_query.strip():
                        # Parse query
                        query_params = parse_locus_query(locus_query)
                        # Filter dataframe
                        filtered_df = filter_dataframe(df, query_params)
                    else:
                        filtered_df = df

                    # Apply geneset filter if selected
                    if selected_genesets["value"]:
                        # Combine all genes from selected genesets
                        combined_genes = set()
                        for geneset_name in selected_genesets["value"]:
                            combined_genes.update(available_genesets[geneset_name])

                        # Filter for rows where VEP_SYMBOL contains any of the genes
                        def matches_geneset(symbol_str):
                            if not symbol_str:
                                return False
                            # VEP_SYMBOL can contain multiple genes separated by &
                            symbols = [
                                s.strip().upper() for s in str(symbol_str).split("&")
                            ]
                            return any(s in combined_genes for s in symbols)

                        # Filter using polars
                        filtered_df = filtered_df.filter(
                            pl.col("VEP_SYMBOL").map_elements(
                                matches_geneset, return_dtype=pl.Boolean
                            )
                        )

                    # Update progress: filtering complete
                    progress.set_value(60)
                    status_label.set_text("Processing variants...")
                    await asyncio.sleep(0)

                    if len(filtered_df) == 0:
                        results_container.clear()
                        with results_container:
                            ui.label("No results found").classes(
                                "text-gray-500 text-lg italic"
                            )
                        return

                    # Convert to list of dicts
                    all_rows = filtered_df.to_dicts()

                    # Load validation data
                    validation_file = store.data_dir / "validations" / "snvs.tsv"
                    validation_map = load_validation_map(validation_file, None)

                    # Yield to event loop before badge processing
                    await asyncio.sleep(0)

                    # Track unknown terms for warnings
                    unknown_consequences = set()
                    unknown_clinvar_terms = set()

                    # Add Variant column and validation status
                    for row in all_rows:
                        chrom = row.get("#CHROM", "")
                        pos = row.get("POS", "")
                        ref = row.get("REF", "")
                        alt = row.get("ALT", "")
                        sample_id = row.get("sample", "")
                        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                        row["Variant"] = variant_key

                        # Add cohort name for FID link generation
                        row["_cohort_name"] = cohort_name

                        # Add FID and Phenotype from pedigree data
                        ped_info = pedigree_data.get(sample_id, {})
                        row["FID"] = ped_info.get("FID", sample_to_family.get(sample_id, ""))
                        row["Phenotype"] = ped_info.get("Phenotype", "")

                        # Add consequence badges (from aggregated comma-separated string)
                        consequence_str = row.get("VEP_Consequence", "")
                        if consequence_str:
                            # Split by both '&' and ',' to handle aggregated values
                            consequences = []
                            for part in str(consequence_str).split(","):
                                for cons in part.split("&"):
                                    cons = cons.strip()
                                    if cons:
                                        consequences.append(cons)

                            row["ConsequenceBadges"] = []
                            seen_badges = set()  # Track unique (label, color) pairs
                            for cons in consequences:
                                # Track unknown consequences
                                if cons and cons not in VEP_CONSEQUENCES:
                                    unknown_consequences.add(cons)
                                label = format_consequence_display(cons)
                                color = get_consequence_color(cons)
                                badge_key = (label, color)
                                if badge_key not in seen_badges:
                                    seen_badges.add(badge_key)
                                    row["ConsequenceBadges"].append(
                                        {
                                            "label": label,
                                            "color": color,
                                        }
                                    )
                        else:
                            row["ConsequenceBadges"] = []

                        # Add ClinVar badges (from aggregated comma-separated string)
                        clinvar_str = row.get("VEP_CLIN_SIG", "")
                        if clinvar_str:
                            # Split by both '&' and ',' to handle aggregated values
                            clinvar_sigs = []
                            for part in str(clinvar_str).split(","):
                                for sig in part.split("&"):
                                    sig = sig.strip()
                                    if sig and sig != ".":
                                        clinvar_sigs.append(sig)

                            row["ClinVarBadges"] = []
                            seen_badges = set()  # Track unique (label, color) pairs
                            for sig in clinvar_sigs:
                                # Track unknown ClinVar terms (case-insensitive check)
                                sig_lower = sig.lower()
                                is_known = any(
                                    key.lower() == sig_lower
                                    for key in CLINVAR_COLORS.keys()
                                )
                                if sig and not is_known:
                                    unknown_clinvar_terms.add(sig)
                                label = format_clinvar_display(sig)
                                color = get_clinvar_color(sig)
                                badge_key = (label, color)
                                if badge_key not in seen_badges:
                                    seen_badges.add(badge_key)
                                    row["ClinVarBadges"].append(
                                        {
                                            "label": label,
                                            "color": color,
                                        }
                                    )
                        else:
                            row["ClinVarBadges"] = []

                        # Add gene badges with color coding based on genesets
                        gene_scorer = get_gene_scorer()

                        # Process VEP_SYMBOL
                        symbol_str = row.get("VEP_SYMBOL", "")
                        if symbol_str:
                            symbols = [
                                s.strip()
                                for s in str(symbol_str).split(",")
                                if s.strip()
                            ]
                            row["GeneBadges"] = []
                            for symbol in symbols:
                                color = gene_scorer.get_gene_color(symbol)
                                tooltip = gene_scorer.get_gene_tooltip(symbol)
                                row["GeneBadges"].append(
                                    {
                                        "label": symbol,
                                        "color": color,
                                        "tooltip": tooltip,
                                    }
                                )
                        else:
                            row["GeneBadges"] = []

                        # Process VEP_Gene (ENSG IDs)
                        gene_str = row.get("VEP_Gene", "")
                        if gene_str:
                            genes = [
                                g.strip() for g in str(gene_str).split(",") if g.strip()
                            ]
                            row["VEP_Gene_badges"] = []
                            for gene in genes:
                                color = gene_scorer.get_gene_color(gene)
                                tooltip = gene_scorer.get_gene_tooltip(gene)
                                row["VEP_Gene_badges"].append(
                                    {
                                        "label": gene,
                                        "color": color,
                                        "tooltip": tooltip,
                                    }
                                )
                        else:
                            row["VEP_Gene_badges"] = []

                        add_validation_status_to_row(
                            row, validation_map, variant_key, sample_id
                        )

                        # Add continuous score badges
                        # Iterate over row columns and check if they have score configs
                        # Use list() to create a copy of items to avoid "dictionary changed size during iteration" error
                        for col_name, value_str in list(row.items()):
                            if value_str and value_str != ".":
                                try:
                                    value = float(value_str)
                                    badge_info = get_score_color(col_name, value)
                                    if badge_info:
                                        row[f"{col_name}_badge"] = {
                                            "label": f"{value:.3f}",
                                            "color": badge_info["color"],
                                            "tooltip": f"{col_name}: {value:.3f} ({badge_info['label']})"
                                        }
                                except (ValueError, TypeError):
                                    pass  # Skip invalid values or non-numeric columns

                    # Update progress: badge processing complete
                    progress.set_value(85)
                    status_label.set_text("Rendering table...")
                    await asyncio.sleep(0)

                    # Display warnings for unknown terms
                    if unknown_consequences:
                        ui.notify(
                            f"⚠️ Unknown VEP consequence terms found: {', '.join(sorted(unknown_consequences))}. "
                            "Please add to vep_consequences.yaml",
                            type="warning",
                            timeout=10000,
                            position="top",
                        )

                    if unknown_clinvar_terms:
                        ui.notify(
                            f"⚠️ Unknown ClinVar terms found: {', '.join(sorted(unknown_clinvar_terms))}. "
                            "Please add to clinvar_colors.yaml",
                            type="warning",
                            timeout=10000,
                            position="top",
                        )

                    # Get all columns
                    all_columns = list(filtered_df.columns)
                    if "Variant" not in all_columns:
                        all_columns.insert(0, "Variant")
                    # Ensure FID column is in the list
                    if "FID" not in all_columns:
                        all_columns.append("FID")
                    # Ensure Phenotype column is in the list
                    if "Phenotype" not in all_columns:
                        all_columns.append("Phenotype")
                    # Ensure Validation column is in the list
                    if "Validation" not in all_columns:
                        all_columns.append("Validation")

                    # Group same-group columns together
                    all_columns = reorder_columns_by_group(all_columns)

                    # Default visible columns
                    default_visible = [
                        "Variant",
                        "VEP_Consequence",
                        "VEP_SYMBOL",
                        "VEP_CLIN_SIG",
                        "fafmax_faf95_max_genomes",
                        "FID",
                        "Phenotype",
                        "sample",
                        "sample_gt",
                        "father_gt",
                        "mother_gt",
                        "Validation",
                    ]
                    # Auto-select preset based on wombat config name
                    wombat_config = selected_file["wombat_config"]
                    initial_preset = select_preset_for_config(wombat_config, VIEW_PRESETS)
                    selected_preset = {"name": initial_preset["name"]}

                    # Override with preset columns if available
                    preset_columns = initial_preset.get("columns", [])
                    initial_selected = [col for col in preset_columns if col in all_columns]

                    selected_cols = {
                        "value": initial_selected if initial_selected else [col for col in default_visible if col in all_columns]
                    }

                    # Table state for persistence across refreshes
                    table_state: Dict[str, Any] = {"sorting": [], "page": 0}

                    # Apply individual filters (sex, phenotype, has-parents)
                    if filter_sex["value"]:
                        all_rows = [
                            r for r in all_rows
                            if pedigree_data.get(r.get("sample", ""), {}).get("Sex", "") in filter_sex["value"]
                        ]

                    if filter_phenotype["value"]:
                        all_rows = [
                            r for r in all_rows
                            if pedigree_data.get(r.get("sample", ""), {}).get("Phenotype", "") in filter_phenotype["value"]
                        ]

                    if filter_has_parents["value"]:
                        def _has_parents(ped: Dict[str, str]) -> bool:
                            father = ped.get("Father", "")
                            mother = ped.get("Mother", "")
                            return (
                                father not in _PED_MISSING
                                and mother not in _PED_MISSING
                            )

                        all_rows = [
                            r for r in all_rows
                            if _has_parents(pedigree_data.get(r.get("sample", ""), {}))
                        ]

                    # Apply impact filter if some impacts are deselected
                    if selected_impacts_search["value"] and set(
                        selected_impacts_search["value"]
                    ) != set(VEP_CONSEQUENCES.keys()):
                        # Filter rows where at least one consequence matches selected impacts
                        filtered_rows = []
                        for row in all_rows:
                            consequence_str = row.get("VEP_Consequence", "")
                            if consequence_str:
                                # Split by both "," (aggregated) and "&" (compound VEP)
                                consequences = []
                                for part in str(consequence_str).split(","):
                                    for c in part.split("&"):
                                        c = c.strip()
                                        if c:
                                            consequences.append(c)
                                # Keep row if any of its consequences is in selected impacts
                                if any(
                                    c in selected_impacts_search["value"]
                                    for c in consequences
                                ):
                                    filtered_rows.append(row)
                        all_rows = filtered_rows

                    # Update progress: ready to display
                    progress.set_value(100)
                    status_label.set_text("Complete!")
                    await asyncio.sleep(0)

                    # Clear progress indicator and show results
                    results_container.clear()
                    with results_container:

                        @ui.refreshable
                        def render_results_table():
                            # No additional filtering - use all_rows as is
                            rows = all_rows.copy()

                            # Apply exclude filters
                            if filter_exclude_lcr["value"]:
                                rows = [
                                    r
                                    for r in rows
                                    if not (
                                        r.get("LCR")
                                        and "true" in str(r.get("LCR", "")).lower()
                                    )
                                ]

                            if filter_exclude_gnomad["value"]:
                                rows = [
                                    r for r in rows if not r.get("genomes_filters")
                                ]

                            if filter_exclude_gnomad_wgs["value"]:
                                rows = [
                                    r
                                    for r in rows
                                    if not r.get("fafmax_faf95_max_genomes")
                                ]

                            # Apply validation filter
                            if selected_validations["value"]:
                                rows = [
                                    row
                                    for row in rows
                                    if row.get("Validation", "") in selected_validations["value"]
                                    or (
                                        "TODO" in selected_validations["value"]
                                        and not row.get("Validation")
                                    )
                                ]



                            def get_columns():
                                cols: List[Dict[str, Any]] = [
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
                                    if col == "Validation":
                                        col_def["cellType"] = "validation"
                                    elif col == "VEP_Consequence":
                                        col_def["cellType"] = "badge_list"
                                        col_def["badgesField"] = "ConsequenceBadges"
                                    elif col == "VEP_CLIN_SIG":
                                        col_def["cellType"] = "badge_list"
                                        col_def["badgesField"] = "ClinVarBadges"
                                    elif col == "VEP_SYMBOL":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "GeneBadges"
                                    elif col == "VEP_Gene":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "VEP_Gene_badges"
                                    elif col == "FID":
                                        col_def["cellType"] = "link"
                                        col_def["href"] = "/cohort/{_cohort_name}/family/{FID}"
                                    else:
                                        col_def["cellType"] = "score_badge"
                                    apply_width_constraints(col_def, col)
                                    cols.append(col_def)
                                return cols

                            # Reference to the DataTable for column visibility updates
                            search_dt: Dict[str, Any] = {"ref": None}

                            def _apply_col_visibility():
                                if search_dt["ref"]:
                                    visible = ["actions"] + list(selected_cols["value"])
                                    search_dt["ref"].set_column_visibility(visible)

                            with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
                                ui.label(f"Results ({len(rows)} rows)").classes(
                                    "text-lg font-semibold text-blue-700"
                                )

                                # Preset dropdown
                                preset_select = ui.select(
                                    options={p["name"]: p["name"] for p in VIEW_PRESETS},
                                    value=selected_preset["name"],
                                    label="Preset"
                                ).classes("w-48")

                                ui.space()  # Push column selector to the right

                                # Column selector dialog
                                col_dialog, _sync_col_selector = build_column_selector(
                                    all_columns=all_columns,
                                    selected_cols=selected_cols,
                                    on_visibility_change=_apply_col_visibility,
                                    presets=VIEW_PRESETS,
                                )
                                ui.button(
                                    "Columns", icon="view_column",
                                    on_click=col_dialog.open,
                                ).props("outline color=blue size=sm")

                                # --- Stats button + dialog ---
                                def show_stats_dialog(current_rows=rows):
                                    from collections import Counter

                                    # Deduplicate by (#CHROM, POS, REF, ALT)
                                    seen: set = set()
                                    unique_variants: list = []
                                    for r in current_rows:
                                        key = (
                                            r.get("#CHROM", ""),
                                            r.get("POS", ""),
                                            r.get("REF", ""),
                                            r.get("ALT", ""),
                                        )
                                        if key not in seen:
                                            seen.add(key)
                                            unique_variants.append(r)

                                    # Classify variant type
                                    for r in unique_variants:
                                        ref = str(r.get("REF", ""))
                                        alt = str(r.get("ALT", ""))
                                        r["_is_snv"] = len(ref) == 1 and len(alt) == 1

                                    chrom_order = CHROM_ORDER
                                    chrom_sizes_mb = CHROM_SIZES_MB
                                    validation_colors = VALIDATION_COLORS

                                    # Filter state
                                    type_filter = {"snv": True, "indel": True}
                                    show_ideogram: Dict[str, bool] = {"value": False}
                                    _containers: Dict[str, Any] = {"charts": None, "ideo": None}

                                    with ui.dialog().props(
                                        "full-width"
                                    ) as stats_dialog, ui.card().classes("w-full"):
                                        with ui.column().classes("w-full p-4"):
                                            # Header
                                            with ui.row().classes(
                                                "items-center justify-between w-full mb-2"
                                            ):
                                                with ui.row().classes("items-center gap-3"):
                                                    ui.label("Variant Statistics").classes(
                                                        "text-xl font-bold text-blue-900"
                                                    )
                                                    subtitle_label = ui.label("").classes(
                                                        "text-sm text-gray-500"
                                                    )
                                                    ideogram_btn = ui.button(
                                                        "Ideogram",
                                                    ).props(
                                                        "outline color=blue size=sm dense no-caps"
                                                    )
                                                    snv_cb = ui.checkbox(
                                                        "SNVs", value=True
                                                    ).props("dense").classes("text-sm")
                                                    indel_cb = ui.checkbox(
                                                        "Indels", value=True
                                                    ).props("dense").classes("text-sm")
                                                ui.button(
                                                    icon="close",
                                                    on_click=lambda: stats_dialog.close(),
                                                ).props("flat round")

                                            @ui.refreshable
                                            def render_stats_content():
                                                # Filter variants by type
                                                filtered = [
                                                    r for r in unique_variants
                                                    if (type_filter["snv"] and r["_is_snv"])
                                                    or (type_filter["indel"] and not r["_is_snv"])
                                                ]
                                                snv_n = sum(1 for r in filtered if r["_is_snv"])
                                                indel_n = len(filtered) - snv_n
                                                subtitle_label.text = (
                                                    f"{len(filtered)} unique variants "
                                                    f"({snv_n} SNVs, {indel_n} Indels)"
                                                )

                                                # Chromosome distribution stacked by validation
                                                chrom_validation: Dict[str, Dict[str, int]] = {
                                                    c: {} for c in chrom_order
                                                }
                                                for r in filtered:
                                                    chrom = norm_chrom(r.get("#CHROM", ""))
                                                    status = r.get("Validation", "") or "TODO"
                                                    if chrom in chrom_validation:
                                                        chrom_validation[chrom][status] = (
                                                            chrom_validation[chrom].get(status, 0) + 1
                                                        )
                                                all_statuses: List[str] = []
                                                for c in chrom_order:
                                                    for s in chrom_validation[c]:
                                                        if s not in all_statuses:
                                                            all_statuses.append(s)

                                                # Consequence distribution
                                                consequence_counts = Counter(
                                                    get_highest_consequence_term(
                                                        str(r.get("VEP_Consequence", ""))
                                                    )
                                                    for r in filtered
                                                )
                                                # Validation distribution
                                                validation_counts = Counter(
                                                    r.get("Validation", "") or "TODO"
                                                    for r in filtered
                                                )

                                                # Scatter data for ideogram
                                                scatter_data: List[List[Any]] = []
                                                for r in filtered:
                                                    chrom = norm_chrom(r.get("#CHROM", ""))
                                                    pos = r.get("POS", 0)
                                                    try:
                                                        pos_mb = round(float(pos) / 1_000_000, 2)
                                                    except (ValueError, TypeError):
                                                        continue
                                                    if chrom in chrom_sizes_mb:
                                                        status = r.get("Validation", "") or "TODO"
                                                        scatter_data.append([pos_mb, chrom, status])

                                                # --- Charts container ---
                                                _containers["charts"] = ui.column().classes("w-full")
                                                _containers["charts"].set_visibility(not show_ideogram["value"])
                                                with _containers["charts"]:
                                                    ui.label("Variants per Chromosome").classes(
                                                        "text-lg font-semibold text-gray-800 mt-2"
                                                    )
                                                    stacked_series = [
                                                        {
                                                            "name": status,
                                                            "type": "bar",
                                                            "stack": "total",
                                                            "data": [
                                                                chrom_validation[c].get(status, 0)
                                                                for c in chrom_order
                                                            ],
                                                            "itemStyle": {
                                                                "color": validation_colors.get(
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
                                                                "data": chrom_order,
                                                                "name": "Chromosome",
                                                            },
                                                            "yAxis": {"type": "value", "name": "Count"},
                                                            "series": stacked_series,
                                                        }
                                                    ).classes("w-full h-64")

                                                    with ui.row().classes("w-full gap-4 flex-wrap mt-4"):
                                                        # Consequence pie chart
                                                        with ui.column().classes("flex-1 min-w-[400px]"):
                                                            ui.label(
                                                                "Consequence Distribution (highest per variant)"
                                                            ).classes("text-lg font-semibold text-gray-800")
                                                            cons_data = [
                                                                {
                                                                    "name": format_consequence_display(cons),
                                                                    "value": count,
                                                                    "itemStyle": {
                                                                        "color": VEP_CONSEQUENCES.get(
                                                                            cons, ("", "#6b7280")
                                                                        )[1]
                                                                    },
                                                                }
                                                                for cons, count in consequence_counts.most_common()
                                                            ]
                                                            ui.echart(
                                                                {
                                                                    "tooltip": {"trigger": "item"},
                                                                    "series": [{
                                                                        "type": "pie",
                                                                        "radius": "70%",
                                                                        "data": cons_data,
                                                                        "label": {"formatter": "{b}: {c} ({d}%)"},
                                                                    }],
                                                                }
                                                            ).classes("w-full h-80")

                                                        # Validation pie chart
                                                        with ui.column().classes("flex-1 min-w-[400px]"):
                                                            ui.label(
                                                                "Validation Status Distribution"
                                                            ).classes("text-lg font-semibold text-gray-800")
                                                            val_data = [
                                                                {
                                                                    "name": st,
                                                                    "value": cnt,
                                                                    "itemStyle": {
                                                                        "color": validation_colors.get(st, "#6b7280")
                                                                    },
                                                                }
                                                                for st, cnt in validation_counts.most_common()
                                                            ]
                                                            ui.echart(
                                                                {
                                                                    "tooltip": {"trigger": "item"},
                                                                    "series": [{
                                                                        "type": "pie",
                                                                        "radius": "70%",
                                                                        "data": val_data,
                                                                        "label": {"formatter": "{b}: {c} ({d}%)"},
                                                                    }],
                                                                }
                                                            ).classes("w-full h-80")

                                                # --- Ideogram container ---
                                                _containers["ideo"] = ui.column().classes("w-full")
                                                _containers["ideo"].set_visibility(show_ideogram["value"])
                                                with _containers["ideo"]:
                                                    svg_w = 1800
                                                    lbl_w = 50
                                                    plot_w = svg_w - lbl_w - 20
                                                    row_h = 16
                                                    tri_h = 6
                                                    row_gap = tri_h + 4
                                                    svg_h = len(chrom_order) * (row_h + row_gap) + 60
                                                    max_mb = max(chrom_sizes_mb.values())

                                                    svg_parts = [
                                                        f'<svg viewBox="0 0 {svg_w} {svg_h}" '
                                                        f'xmlns="http://www.w3.org/2000/svg" '
                                                        f'preserveAspectRatio="xMinYMin meet" '
                                                        f'style="font-family: sans-serif; width: 100%; height: auto;">'
                                                    ]

                                                    axis_y = len(chrom_order) * (row_h + row_gap)
                                                    for mb_val in range(0, 260, 50):
                                                        gx = lbl_w + (mb_val / max_mb) * plot_w
                                                        svg_parts.append(
                                                            f'<line x1="{gx:.1f}" y1="0" '
                                                            f'x2="{gx:.1f}" y2="{axis_y}" '
                                                            f'stroke="#e5e7eb" stroke-width="0.5" '
                                                            f'stroke-dasharray="3,3"/>'
                                                        )
                                                        svg_parts.append(
                                                            f'<text x="{gx:.1f}" y="{axis_y + 14}" '
                                                            f'text-anchor="middle" font-size="12" '
                                                            f'fill="#6b7280">{mb_val}</text>'
                                                        )
                                                    svg_parts.append(
                                                        f'<text x="{svg_w / 2}" y="{axis_y + 32}" '
                                                        f'text-anchor="middle" font-size="13" '
                                                        f'fill="#6b7280">Position (Mb)</text>'
                                                    )

                                                    for ci, chrom in enumerate(chrom_order):
                                                        bar_y = ci * (row_h + row_gap) + tri_h
                                                        bands = CYTOBANDS.get(chrom, [])
                                                        cs = chrom_sizes_mb.get(chrom, 0)
                                                        total_w = (cs / max_mb) * plot_w

                                                        svg_parts.append(
                                                            f'<text x="{lbl_w - 6}" y="{bar_y + row_h * 0.75}" '
                                                            f'text-anchor="end" font-size="12" '
                                                            f'fill="#374151">{chrom}</text>'
                                                        )
                                                        for band in bands:
                                                            bx = lbl_w + (band["start"] / max_mb) * plot_w
                                                            bw = max(
                                                                ((band["end"] - band["start"]) / max_mb) * plot_w,
                                                                0.5,
                                                            )
                                                            color = GIESTAIN_COLORS.get(band["stain"], "#e5e7eb")
                                                            svg_parts.append(
                                                                f'<rect x="{bx:.1f}" y="{bar_y}" '
                                                                f'width="{bw:.1f}" height="{row_h}" '
                                                                f'fill="{color}"/>'
                                                            )
                                                        svg_parts.append(
                                                            f'<rect x="{lbl_w}" y="{bar_y}" '
                                                            f'width="{total_w:.1f}" height="{row_h}" '
                                                            f'fill="none" stroke="#9ca3af" '
                                                            f'stroke-width="0.5" rx="3"/>'
                                                        )

                                                    for sd in scatter_data:
                                                        v_mb, v_chrom, v_status = sd
                                                        if v_chrom not in chrom_order:
                                                            continue
                                                        v_idx = chrom_order.index(v_chrom)
                                                        bar_y = v_idx * (row_h + row_gap) + tri_h
                                                        vx = lbl_w + (v_mb / max_mb) * plot_w
                                                        v_color = validation_colors.get(v_status, "#94a3b8")
                                                        tw = 5
                                                        svg_parts.append(
                                                            f'<polygon points="{vx - tw:.1f},{bar_y - tri_h} '
                                                            f'{vx + tw:.1f},{bar_y - tri_h} '
                                                            f'{vx:.1f},{bar_y}" '
                                                            f'fill="{v_color}" opacity="0.9"/>'
                                                        )
                                                        svg_parts.append(
                                                            f'<line x1="{vx:.1f}" y1="{bar_y}" '
                                                            f'x2="{vx:.1f}" y2="{bar_y + row_h}" '
                                                            f'stroke="{v_color}" stroke-width="1.5" '
                                                            f'opacity="0.85"/>'
                                                        )

                                                    legend_y = axis_y + 40
                                                    legend_x = lbl_w
                                                    for v_status in all_statuses:
                                                        v_color = validation_colors.get(v_status, "#94a3b8")
                                                        svg_parts.append(
                                                            f'<rect x="{legend_x}" y="{legend_y}" '
                                                            f'width="12" height="12" rx="2" fill="{v_color}"/>'
                                                        )
                                                        svg_parts.append(
                                                            f'<text x="{legend_x + 16}" y="{legend_y + 10}" '
                                                            f'font-size="12" fill="#374151">{v_status}</text>'
                                                        )
                                                        legend_x += len(v_status) * 8 + 32

                                                    svg_parts.append("</svg>")
                                                    ui.html(
                                                        "\n".join(svg_parts),
                                                        sanitize=False,
                                                    ).classes("w-full")

                                            render_stats_content()

                                            # Toggle ideogram/charts view
                                            def toggle_ideogram(_e=None):
                                                show_ideogram["value"] = not show_ideogram["value"]
                                                if _containers["charts"]:
                                                    _containers["charts"].set_visibility(
                                                        not show_ideogram["value"]
                                                    )
                                                if _containers["ideo"]:
                                                    _containers["ideo"].set_visibility(
                                                        show_ideogram["value"]
                                                    )
                                                if show_ideogram["value"]:
                                                    ideogram_btn.props(
                                                        remove="outline", add="unelevated"
                                                    )
                                                else:
                                                    ideogram_btn.props(
                                                        remove="unelevated", add="outline"
                                                    )
                                                ideogram_btn.update()

                                            ideogram_btn.on_click(toggle_ideogram)

                                            # SNV / Indel filter handler
                                            def on_type_filter_change(_e=None):
                                                type_filter["snv"] = snv_cb.value
                                                type_filter["indel"] = indel_cb.value
                                                render_stats_content.refresh()

                                            snv_cb.on_value_change(on_type_filter_change)
                                            indel_cb.on_value_change(on_type_filter_change)

                                    stats_dialog.open()

                                ui.button(
                                    "Stats", icon="bar_chart", on_click=show_stats_dialog
                                ).props("outline color=blue size=sm")

                            # Preset change handler
                            def on_preset_change(e):
                                """Handle preset selection change."""
                                preset_name = e.value
                                preset = next((p for p in VIEW_PRESETS if p["name"] == preset_name), None)
                                if not preset:
                                    return

                                # Filter columns to only those available in the data
                                available = [col for col in preset.get("columns", []) if col in all_columns]

                                selected_cols["value"] = available
                                selected_preset["name"] = preset_name
                                _apply_col_visibility()
                                _sync_col_selector()

                            # Connect preset change handler
                            preset_select.on_value_change(on_preset_change)

                            # Handle view variant click
                            def on_view_variant(e):
                                row_data = e.get("row", {})
                                variant_str = row_data.get("Variant", "")
                                sample_id = row_data.get("sample", "")

                                # Get family from sample
                                family_id = sample_to_family.get(sample_id)

                                if not family_id:
                                    ui.notify(
                                        f"Could not find family for sample {sample_id}",
                                        type="warning",
                                    )
                                    return

                                try:
                                    parts = variant_str.split(":")
                                    if len(parts) == 4:
                                        chrom, pos, ref, alt = parts

                                        # Create variant data dict
                                        variant_data = dict(row_data)

                                        # Callback to refresh validation status
                                        def on_save(validation_status: str):
                                            # Reload validation map
                                            validation_map_updated = (
                                                load_validation_map(
                                                    validation_file, None
                                                )
                                            )
                                            # Update validation status for all rows
                                            for row in all_rows:
                                                v_key = row.get("Variant", "")
                                                s_id = row.get("sample", "")
                                                add_validation_status_to_row(
                                                    row,
                                                    validation_map_updated,
                                                    v_key,
                                                    s_id,
                                                )
                                            # Refresh table using captured client context
                                            with page_client:
                                                ui.timer(
                                                    0.1,
                                                    render_results_table.refresh,
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
                                    ui.notify(
                                        f"Error parsing variant: {ex}", type="warning"
                                    )

                            # Restore table state (sorting / page) across refreshes
                            saved_sorting = table_state.get("sorting", [])
                            if saved_sorting:
                                col_id = saved_sorting[0]["id"]
                                desc = saved_sorting[0].get("desc", False)
                                col_def = next(
                                    (c for c in get_columns() if c.get("id") == col_id), {}
                                )
                                sort_field = col_def.get("sortField", col_id)
                                sort_type = col_def.get("sorting", "")
                                if sort_type == "genomic":
                                    from genetics_viz.utils.column_names import genomic_sort_key
                                    rows.sort(
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
                                    rows.sort(key=_num_key, reverse=desc)
                                else:
                                    rows.sort(
                                        key=lambda r: (
                                            r.get(sort_field) is None,
                                            r.get(sort_field, ""),
                                        ),
                                        reverse=desc,
                                    )

                            # Create table with all columns, initial visibility from preset
                            search_dt["ref"] = DataTable(
                                columns=get_columns(),
                                rows=rows,
                                row_key="Variant",
                                pagination={"rowsPerPage": 50},
                                visible_columns=["actions"] + list(selected_cols["value"]),
                                on_row_action=on_view_variant,
                                initial_sorting=saved_sorting,
                                initial_page=table_state.get("page", 0),
                                state_holder=table_state,
                            )

                        render_results_table()

                except Exception as e:
                    import traceback

                    results_container.clear()
                    with results_container:
                        ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
                        ui.label("Traceback:").classes("text-red-500 font-semibold")
                        ui.label(traceback.format_exc()).classes(
                            "text-red-500 text-xs font-mono whitespace-pre"
                        )

            # Set up handlers after function definition
            search_button.on_click(perform_search)
            locus_input.on("keydown.enter", perform_search)

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
