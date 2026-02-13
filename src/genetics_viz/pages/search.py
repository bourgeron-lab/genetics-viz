"""Search page for cohort-wide variant search."""

import asyncio
import csv
import re
from pathlib import Path
from typing import Any, Dict, List

import polars as pl
import yaml
from nicegui import app as nicegui_app
from nicegui import context, ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.pages.cohort.components.wombat_tab import (
    VIEW_PRESETS,
    select_preset_for_config,
)
from genetics_viz.utils.column_names import (
    get_column_group,
    get_display_label,
    reorder_columns_by_group,
)
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.score_colors import get_score_color

# Load VEP Consequence data from YAML
def _load_vep_consequences() -> Dict[str, tuple]:
    """Load VEP consequences from YAML config file."""
    config_path = Path(__file__).parent.parent / "config" / "vep_consequences.yaml"
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)
    # Convert to dict with (impact, color) tuples
    return {term: (info["impact"], info["color"]) for term, info in data.items()}


VEP_CONSEQUENCES = _load_vep_consequences()


def get_consequence_impact(consequence: str) -> str:
    """Get impact level for a consequence term."""
    return VEP_CONSEQUENCES.get(consequence, ("MODIFIER", "#636363"))[0]


def get_consequence_color(consequence: str) -> str:
    """Get color for a consequence term."""
    return VEP_CONSEQUENCES.get(consequence, ("MODIFIER", "#636363"))[1]


def format_consequence_display(consequence: str) -> str:
    """Format consequence for display: remove _variant suffix and replace _ with space."""
    display = consequence.replace("_variant", "").replace("_", " ")
    return display


# Load ClinVar significance colors from YAML
def _load_clinvar_colors() -> Dict[str, str]:
    """Load ClinVar colors from YAML config file."""
    config_path = Path(__file__).parent.parent / "config" / "clinvar_colors.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


CLINVAR_COLORS = _load_clinvar_colors()


def get_clinvar_color(significance: str) -> str:
    """Get color for a ClinVar significance term (case-insensitive)."""
    if not significance:
        return "#757575"  # Default to gray
    # Case-insensitive lookup
    sig_lower = significance.lower()
    for key, color in CLINVAR_COLORS.items():
        if key.lower() == sig_lower:
            return color
    return "#757575"  # Default to gray if not found


def format_clinvar_display(significance: str) -> str:
    """Format ClinVar significance for display: replace _ with space."""
    return significance.replace("_", " ")





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


def load_sample_to_family_map(pedigree_file: Path) -> Dict[str, str]:
    """Load mapping from sample ID to family ID from pedigree file."""
    sample_to_family = {}

    if not pedigree_file.exists():
        return sample_to_family

    with open(pedigree_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Handle case where there's no header
        if reader.fieldnames and not reader.fieldnames[0].lower().startswith("fid"):
            f.seek(0)
            lines = f.readlines()
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    fid, iid = parts[0], parts[1]
                    sample_to_family[iid] = fid
        elif reader.fieldnames:
            # Map column names (case-insensitive)
            fieldnames_lower = {fn.lower(): fn for fn in reader.fieldnames}
            fid_col = None
            iid_col = None

            for possible in ["fid", "family_id", "familyid", "family"]:
                if possible in fieldnames_lower:
                    fid_col = fieldnames_lower[possible]
                    break

            for possible in ["iid", "individual_id", "sample_id", "sample"]:
                if possible in fieldnames_lower:
                    iid_col = fieldnames_lower[possible]
                    break

            if fid_col and iid_col:
                for row in reader:
                    sample_to_family[row[iid_col]] = row[fid_col]

    return sample_to_family


def load_pedigree_data(pedigree_file: Path) -> Dict[str, Dict[str, str]]:
    """Load full pedigree data from pedigree file.

    Returns dict mapping sample ID to pedigree info (FID, Phenotype, etc.)
    """
    pedigree_data = {}

    if not pedigree_file.exists():
        return pedigree_data

    with open(pedigree_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Handle case where there's no header
        if reader.fieldnames and not reader.fieldnames[0].lower().startswith("fid"):
            f.seek(0)
            lines = f.readlines()
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    fid, iid, phenotype = parts[0], parts[1], parts[5]
                    pedigree_data[iid] = {
                        "FID": fid,
                        "Phenotype": phenotype,
                    }
        elif reader.fieldnames:
            # Map column names (case-insensitive)
            fieldnames_lower = {fn.lower(): fn for fn in reader.fieldnames}
            fid_col = None
            iid_col = None
            phenotype_col = None

            for possible in ["fid", "family_id", "familyid", "family"]:
                if possible in fieldnames_lower:
                    fid_col = fieldnames_lower[possible]
                    break

            for possible in ["iid", "individual_id", "sample_id", "sample"]:
                if possible in fieldnames_lower:
                    iid_col = fieldnames_lower[possible]
                    break

            for possible in ["phenotype", "pheno", "status", "affected"]:
                if possible in fieldnames_lower:
                    phenotype_col = fieldnames_lower[possible]
                    break

            if fid_col and iid_col:
                for row in reader:
                    pedigree_data[row[iid_col]] = {
                        "FID": row[fid_col],
                        "Phenotype": row.get(phenotype_col, "") if phenotype_col else "",
                    }

    return pedigree_data


@ui.page("/search/{cohort_name}")
def search_cohort_page(cohort_name: str) -> None:
    """Search page for cohort-wide variant search."""
    create_header()

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

            # Load pedigree data for family and phenotype info
            pedigree_file = (
                store.data_dir / "cohorts" / cohort_name / f"{cohort_name}.pedigree.tsv"
            )
            sample_to_family = load_sample_to_family_map(pedigree_file)
            pedigree_data = load_pedigree_data(pedigree_file)

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

            # Button references for visual indicators
            geneset_button_ref: Dict[str, Any] = {"button": None}
            impact_button_ref: Dict[str, Any] = {"button": None}

            # Search panel
            with ui.card().classes("w-full p-4 mb-4"):
                ui.label("Search Parameters").classes("text-xl font-semibold mb-4")

                with ui.row().classes("items-end gap-4 w-full flex-wrap"):
                    # Source dropdown
                    source_select = (
                        ui.select(
                            options=[wf["display_name"] for wf in wombat_files],
                            label="Source",
                            value=wombat_files[0]["display_name"]
                            if wombat_files
                            else None,
                        )
                        .props("outlined")
                        .classes("w-64")
                    )

                    # Locus input (Enter key handler set later after function definition)
                    locus_input = (
                        ui.input(
                            label="Locus (optional)",
                            placeholder="chr1:10000-10100, SHANK3, ENSG...",
                            on_change=lambda: None,  # Placeholder to enable events
                        )
                        .props("outlined")
                        .classes("flex-grow")
                    )

                    # Geneset filter menu
                    if available_genesets:
                        geneset_btn = (
                            ui.button("Genesets", icon="list")
                            .props(
                                "outline"
                                if not selected_genesets["value"]
                                else "unelevated color=green"
                            )
                            .classes("h-14")
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
                        .props("outline")
                        .classes("h-14")
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
                        button_size="h-14",
                    )

                    # Search button (handler set later after function definition)
                    search_button = ui.button(
                        "Search",
                        icon="search",
                    ).props("color=blue").classes("h-14")

                # Exclude filters row
                with ui.row().classes("items-center gap-4 w-full flex-wrap mt-2"):
                    ui.checkbox(
                        "Exclude LCR",
                        value=filter_exclude_lcr["value"],
                        on_change=lambda e: filter_exclude_lcr.update(
                            {"value": e.value}
                        ),
                    )
                    ui.checkbox(
                        "Exclude gnomAD filtered",
                        value=filter_exclude_gnomad["value"],
                        on_change=lambda e: filter_exclude_gnomad.update(
                            {"value": e.value}
                        ),
                    )
                    ui.checkbox(
                        "Exclude gnomAD WGS",
                        value=filter_exclude_gnomad_wgs["value"],
                        on_change=lambda e: filter_exclude_gnomad_wgs.update(
                            {"value": e.value}
                        ),
                    )

                # Help text
                with ui.expansion("Query Examples", icon="help").classes("mt-2"):
                    ui.markdown("""
- **chr1:10000-10100** - All variants in range [10000, 10100] on chr1
- **chr1:10000** - Exact position 10000 on chr1
- **chr1:10000:A:GC** - Exact variant with REF=A and ALT=GC
- **SHANK3** - All variants in SHANK3 gene
- **SHANK*** - All variants in genes starting with SHANK
- **ENSG00000164099** - All variants in gene with this Ensembl ID
                    """)

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
                            infer_schema_length=100,
                            schema_overrides={"sex": pl.Utf8},
                            null_values=[".", ""],
                        )

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

                    # Apply impact filter if some impacts are deselected
                    if selected_impacts_search["value"] and set(
                        selected_impacts_search["value"]
                    ) != set(VEP_CONSEQUENCES.keys()):
                        # Filter rows where at least one consequence matches selected impacts
                        filtered_rows = []
                        for row in all_rows:
                            consequence_str = row.get("VEP_Consequence", "")
                            if consequence_str:
                                consequences = [
                                    c.strip() for c in str(consequence_str).split("&")
                                ]
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

                                # Compact column selector button
                                with ui.button(
                                    "+ column", icon="view_column"
                                ).props("outline color=blue size=sm"):
                                    with ui.menu():
                                        ui.label("Show/Hide Columns:").classes(
                                            "px-4 py-2 font-semibold text-sm"
                                        )
                                        ui.separator()

                                        with ui.column().classes("p-2"):
                                            col_checkboxes = {}

                                            def _sync_col_checkboxes():
                                                for col, cb in col_checkboxes.items():
                                                    cb.value = col in selected_cols["value"]

                                            with ui.row().classes("gap-2 mb-2"):

                                                def col_select_all():
                                                    selected_cols["value"] = (
                                                        all_columns.copy()
                                                    )
                                                    _apply_col_visibility()
                                                    _sync_col_checkboxes()

                                                def col_select_none():
                                                    selected_cols["value"] = []
                                                    _apply_col_visibility()
                                                    _sync_col_checkboxes()

                                                ui.button(
                                                    "All", on_click=col_select_all
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )
                                                ui.button(
                                                    "None", on_click=col_select_none
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )

                                            def handle_col_change(col_name, is_checked):
                                                if (
                                                    is_checked
                                                    and col_name
                                                    not in selected_cols["value"]
                                                ):
                                                    selected_cols["value"].append(
                                                        col_name
                                                    )
                                                elif (
                                                    not is_checked
                                                    and col_name
                                                    in selected_cols["value"]
                                                ):
                                                    selected_cols["value"].remove(
                                                        col_name
                                                    )

                                                # Reorder to match all_columns order
                                                selected_cols["value"] = [
                                                    col for col in all_columns
                                                    if col in selected_cols["value"]
                                                ]

                                                _apply_col_visibility()

                                            for col in all_columns:
                                                col_checkboxes[col] = ui.checkbox(
                                                    get_display_label(col),
                                                    value=col in selected_cols["value"],
                                                    on_change=lambda e,
                                                    c=col: handle_col_change(
                                                        c, e.value
                                                    ),
                                                ).classes("text-sm")

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
                                _sync_col_checkboxes()

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

                            # Create table with all columns, initial visibility from preset
                            search_dt["ref"] = DataTable(
                                columns=get_columns(),
                                rows=rows,
                                row_key="Variant",
                                pagination={"rowsPerPage": 50},
                                visible_columns=["actions"] + list(selected_cols["value"]),
                                on_row_action=on_view_variant,
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
