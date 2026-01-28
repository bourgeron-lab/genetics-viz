"""Search page for cohort-wide variant search."""

import csv
import re
from pathlib import Path
from typing import Any, Dict, List

import polars as pl
import yaml
from nicegui import app as nicegui_app
from nicegui import context, ui

from genetics_viz.components.header import create_header
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.data import get_data_store

# Same table slot as wombat_tab
SEARCH_TABLE_SLOT = r"""
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
            <template v-else-if="col.name === 'VEP_Consequence'">
                <div style="display: flex; flex-wrap: wrap; gap: 4px;">
                    <q-badge v-for="(badge, idx) in (props.row.ConsequenceBadges || [])" :key="idx" 
                             :style="'background-color: ' + badge.color + '; color: white; font-size: 0.875em; padding: 4px 8px;'">
                        {{ badge.label }}
                    </q-badge>
                </div>
            </template>
            <template v-else-if="col.name === 'VEP_CLIN_SIG'">
                <div style="display: flex; flex-wrap: wrap; gap: 4px;">
                    <q-badge v-for="(badge, idx) in (props.row.ClinVarBadges || [])" :key="idx" 
                             :style="'background-color: ' + badge.color + '; color: white; font-size: 0.875em; padding: 4px 8px;'">
                        {{ badge.label }}
                    </q-badge>
                </div>
            </template>
            <template v-else>
                {{ col.value }}
            </template>
        </q-td>
    </q-tr>
"""


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


def get_display_label(col: str) -> str:
    """Get display label for column."""
    if col == "fafmax_faf95_max_genomes":
        return "gnomAD 4.1 WGS"
    elif col == "nhomalt_genomes":
        return "gnomAD 4.1 nhomalt WGS"
    elif col == "VEP_CLIN_SIG":
        return "ClinVar"
    elif col.startswith("VEP_"):
        return col[4:]
    else:
        return col


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

            # Load pedigree for sample->family mapping
            pedigree_file = (
                store.data_dir / "cohorts" / cohort_name / f"{cohort_name}.pedigree.tsv"
            )
            sample_to_family = load_sample_to_family_map(pedigree_file)

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

                    # Locus input
                    locus_input = (
                        ui.input(
                            label="Locus (optional)",
                            placeholder="chr1:10000-10100, SHANK3, ENSG...",
                            on_change=lambda: None,  # Placeholder to enable events
                        )
                        .props("outlined")
                        .classes("flex-grow")
                    )
                    # Add Enter key handler
                    locus_input.on("keydown.enter", lambda: perform_search())

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

                    # Search button
                    ui.button(
                        "Search",
                        icon="search",
                        on_click=lambda: perform_search(),
                    ).props("color=blue").classes("h-14")

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

            def perform_search():
                """Execute the search and display results."""
                results_container.clear()

                # Show spinner while loading
                with results_container:
                    with ui.row().classes("items-center gap-4 justify-center py-8"):
                        ui.spinner(size="lg", color="blue")
                        ui.label("Searching...").classes("text-lg text-gray-600")

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
                    with results_container:
                        ui.label("Selected source not found").classes("text-red-500")
                    return

                try:
                    # Load dataframe
                    df = pl.read_csv(
                        selected_file["file_path"],
                        separator="\t",
                        infer_schema_length=100,
                    )

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

                        # Add FID from sample mapping
                        row["FID"] = sample_to_family.get(sample_id, "")

                        # Add consequence badges
                        consequence_str = row.get("VEP_Consequence", "")
                        if consequence_str:
                            consequences = [
                                c.strip() for c in str(consequence_str).split("&")
                            ]
                            row["ConsequenceBadges"] = []
                            for cons in consequences:
                                # Track unknown consequences
                                if cons and cons not in VEP_CONSEQUENCES:
                                    unknown_consequences.add(cons)
                                row["ConsequenceBadges"].append(
                                    {
                                        "label": format_consequence_display(cons),
                                        "color": get_consequence_color(cons),
                                    }
                                )
                        else:
                            row["ConsequenceBadges"] = []

                        # Add ClinVar badges
                        clinvar_str = row.get("VEP_CLIN_SIG", "")
                        if clinvar_str:
                            clinvar_sigs = [
                                c.strip()
                                for c in str(clinvar_str).split("&")
                                if c.strip() and c.strip() != "."
                            ]
                            row["ClinVarBadges"] = []
                            for sig in clinvar_sigs:
                                # Track unknown ClinVar terms (case-insensitive check)
                                sig_lower = sig.lower()
                                is_known = any(
                                    key.lower() == sig_lower
                                    for key in CLINVAR_COLORS.keys()
                                )
                                if sig and not is_known:
                                    unknown_clinvar_terms.add(sig)
                                row["ClinVarBadges"].append(
                                    {
                                        "label": format_clinvar_display(sig),
                                        "color": get_clinvar_color(sig),
                                    }
                                )
                        else:
                            row["ClinVarBadges"] = []

                        add_validation_status_to_row(
                            row, validation_map, variant_key, sample_id
                        )

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
                    # Ensure Validation column is in the list
                    if "Validation" not in all_columns:
                        all_columns.append("Validation")

                    # Default visible columns
                    default_visible = [
                        "Variant",
                        "VEP_Consequence",
                        "VEP_SYMBOL",
                        "VEP_CLIN_SIG",
                        "fafmax_faf95_max_genomes",
                        "FID",
                        "sample",
                        "sample_gt",
                        "father_gt",
                        "mother_gt",
                        "Validation",
                    ]
                    selected_cols = {
                        "value": [col for col in default_visible if col in all_columns]
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

                    # Clear spinner and show results
                    results_container.clear()
                    with results_container:

                        @ui.refreshable
                        def render_results_table():
                            # No additional filtering - use all_rows as is
                            rows = all_rows.copy()

                            # Prepare columns
                            visible_cols = selected_cols["value"]

                            def get_columns():
                                cols: List[Dict[str, Any]] = [
                                    {"name": "actions", "label": "", "field": "actions"}
                                ]
                                cols.extend(
                                    [
                                        {
                                            "name": col,
                                            "label": get_display_label(col),
                                            "field": col,
                                            "sortable": True,
                                            "align": "left",
                                        }
                                        for col in visible_cols
                                    ]
                                )
                                return cols

                            with ui.row().classes("items-center gap-4 mt-4 mb-2"):
                                ui.label(f"Results ({len(rows)} rows)").classes(
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
                                            col_checkboxes = {}

                                            with ui.row().classes("gap-2 mb-2"):

                                                def col_select_all():
                                                    selected_cols["value"] = (
                                                        all_columns.copy()
                                                    )
                                                    for cb in col_checkboxes.values():
                                                        cb.value = True
                                                    render_results_table.refresh()

                                                def col_select_none():
                                                    selected_cols["value"] = []
                                                    for cb in col_checkboxes.values():
                                                        cb.value = False
                                                    render_results_table.refresh()

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
                                                render_results_table.refresh()

                                            for col in all_columns:
                                                col_checkboxes[col] = ui.checkbox(
                                                    get_display_label(col),
                                                    value=col in selected_cols["value"],
                                                    on_change=lambda e,
                                                    c=col: handle_col_change(
                                                        c, e.value
                                                    ),
                                                ).classes("text-sm")

                            # Create table
                            results_table = (
                                ui.table(
                                    columns=get_columns(),
                                    rows=rows,
                                    row_key="Variant",
                                    pagination={"rowsPerPage": 50},
                                )
                                .classes("w-full")
                                .props("dense flat")
                            )

                            results_table.add_slot("body", SEARCH_TABLE_SLOT)

                            # Handle view variant click
                            def on_view_variant(e):
                                row_data = e.args
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

                            results_table.on("view_variant", on_view_variant)

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

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
