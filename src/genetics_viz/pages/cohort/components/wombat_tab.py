"""Wombat tab component for family page."""

import re
from pathlib import Path
from typing import Any, Callable, Dict, List

import polars as pl
import yaml
from nicegui import ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.column_names import (
    get_column_group,
    get_display_label,
    reorder_columns_by_group,
)
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.score_colors import get_score_color


# Load VEP Consequence data from YAML
def _load_vep_consequences() -> Dict[str, tuple]:
    """Load VEP consequences from YAML config file."""
    config_path = (
        Path(__file__).parent.parent.parent.parent / "config" / "vep_consequences.yaml"
    )
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)
    # Convert to dict with (impact, color) tuples
    return {term: (info["impact"], info["color"]) for term, info in data.items()}


VEP_CONSEQUENCES = _load_vep_consequences()

# Cytoband stain colors for ideogram rendering
GIESTAIN_COLORS = {
    "gneg": "#f5f5f5",
    "gpos25": "#c8c8c8",
    "gpos50": "#969696",
    "gpos75": "#646464",
    "gpos100": "#323232",
    "acen": "#d92f27",
    "gvar": "#646464",
    "stalk": "#969696",
}


def _load_cytobands() -> Dict[str, list]:
    """Load cytoband data from TSV file, grouped by chromosome."""
    config_path = (
        Path(__file__).parent.parent.parent.parent / "config" / "cytobands_hg38.tsv"
    )
    bands: Dict[str, list] = {}
    if not config_path.exists():
        return bands
    with open(config_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0].replace("chr", "").upper()
            if chrom == "M":
                chrom = "MT"
            start_mb = int(parts[1]) / 1_000_000
            end_mb = int(parts[2]) / 1_000_000
            name = parts[3]
            stain = parts[4]
            bands.setdefault(chrom, []).append(
                {"start": start_mb, "end": end_mb, "name": name, "stain": stain}
            )
    return bands


CYTOBANDS = _load_cytobands()

# Create consequence priority map (lower number = higher priority)
VEP_CONSEQUENCE_PRIORITY = {
    term: idx for idx, term in enumerate(_load_vep_consequences().keys())
}


def get_consequence_priority(consequence: str) -> int:
    """Get priority for a consequence term (lower = higher priority)."""
    return VEP_CONSEQUENCE_PRIORITY.get(consequence, 9999)


def get_highest_priority_consequence(consequence_str: str) -> int:
    """Get the highest priority (lowest number) from a comma/ampersand-separated consequence string."""
    if not consequence_str:
        return 9999

    # Split by both '&' and ',' to handle aggregated values
    consequences = []
    for part in str(consequence_str).split(","):
        for cons in part.split("&"):
            cons = cons.strip()
            if cons:
                consequences.append(cons)

    if not consequences:
        return 9999

    # Return the lowest priority number (highest priority consequence)
    return min(get_consequence_priority(cons) for cons in consequences)


def get_highest_consequence_term(consequence_str: str) -> str:
    """Get the highest-priority consequence term from a comma/ampersand-separated string."""
    if not consequence_str:
        return "Unknown"
    consequences = []
    for part in str(consequence_str).split(","):
        for cons in part.split("&"):
            cons = cons.strip()
            if cons:
                consequences.append(cons)
    if not consequences:
        return "Unknown"
    return min(consequences, key=lambda c: VEP_CONSEQUENCE_PRIORITY.get(c, 9999))


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
    config_path = (
        Path(__file__).parent.parent.parent.parent / "config" / "clinvar_colors.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


CLINVAR_COLORS = _load_clinvar_colors()


# Load view presets from YAML
def _load_view_presets() -> List[Dict[str, Any]]:
    """Load view presets from YAML config file."""
    config_path = (
        Path(__file__).parent.parent.parent.parent / "config" / "view_presets.yaml"
    )
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)
    return data.get("presets", [])


VIEW_PRESETS = _load_view_presets()


def reload_wombat_configs():
    """Reload all wombat tab configurations from YAML files."""
    global VEP_CONSEQUENCES, CLINVAR_COLORS, VIEW_PRESETS
    VEP_CONSEQUENCES = _load_vep_consequences()
    CLINVAR_COLORS = _load_clinvar_colors()
    VIEW_PRESETS = _load_view_presets()


def select_preset_for_config(config_name: str, presets: List[Dict]) -> Dict:
    """Select the first preset whose keywords contain the config_name,
    or return the first preset if none match."""
    config_lower = config_name.lower()

    for preset in presets:
        keywords = preset.get("keywords", [])
        if any(keyword.lower() in config_lower for keyword in keywords):
            return preset

    # Return first preset as default
    return presets[0] if presets else {"name": "Default", "columns": []}


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



def render_wombat_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the Wombat tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    wombat_dir = store.data_dir / "families" / family_id / "wombat"

    if not wombat_dir.exists():
        ui.label(f"No wombat directory found at: {wombat_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Parse wombat TSV files
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
        return

    # Load genesets from params/genesets (once, shared across all configs)
    genesets_dir = store.data_dir / "params" / "genesets"
    available_genesets: Dict[str, set] = {}
    if genesets_dir.exists():
        for geneset_file in genesets_dir.glob("*.tsv"):
            geneset_name = geneset_file.stem
            genes: set = set()
            with open(geneset_file, "r") as f:
                next(f, None)  # Skip header line
                for line in f:
                    gene = line.strip()
                    if gene:
                        genes.add(gene.upper())
            if genes:
                available_genesets[geneset_name] = genes

    # Create dictionaries to store data for each wombat config
    wombat_data: Dict[str, Dict[str, Any]] = {}

    # Create subtabs for each wombat config
    with ui.tabs().classes("w-full") as wombat_subtabs:
        subtab_refs = {}
        for wf in wombat_files:
            subtab_refs[wf["wombat_config"]] = ui.tab(wf["wombat_config"])

    with ui.tab_panels(wombat_subtabs, value=list(subtab_refs.values())[0]).classes(
        "w-full"
    ):
        for wf in wombat_files:
            with ui.tab_panel(subtab_refs[wf["wombat_config"]]):
                config_name = wf["wombat_config"]

                with ui.card().classes("w-full p-4"):
                    ui.label(f"Wombat Configuration: {wf['wombat_config']}").classes(
                        "text-lg font-semibold text-blue-700 mb-2"
                    )
                    with ui.row().classes("gap-4 mb-4"):
                        ui.label("VEP Config:").classes("font-semibold")
                        ui.badge(wf["vep_config"]).props("color=indigo")
                    with ui.row().classes("gap-4"):
                        ui.label("File Path:").classes("font-semibold")
                        ui.label(str(wf["file_path"])).classes(
                            "text-sm text-gray-600 font-mono"
                        )

                # Display TSV content in a table
                try:
                    df = pl.read_csv(
                        wf["file_path"],
                        separator="\t",
                        infer_schema_length=100,
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
                            .str.concat(delimiter=",")
                            .alias(col)
                        )

                    # Group and aggregate
                    df = df.group_by(grouping_cols, maintain_order=True).agg(agg_exprs)

                    # Convert to list of dicts for NiceGUI table
                    all_rows = df.to_dicts()

                    # Store in the wombat_data dict keyed by config
                    wombat_data[config_name] = {
                        "df": df,
                        "all_rows": all_rows,
                    }

                    # Load validation data from snvs.tsv
                    validation_file = store.data_dir / "validations" / "snvs.tsv"
                    validation_map = load_validation_map(validation_file, family_id)

                    # Track unknown terms for warnings
                    unknown_consequences = set()
                    unknown_clinvar_terms = set()

                    # Add concatenated Variant column and Validation status to each row
                    for row in all_rows:
                        chrom = row.get("#CHROM", "")
                        pos = row.get("POS", "")
                        ref = row.get("REF", "")
                        alt = row.get("ALT", "")
                        sample_id = row.get("sample", "")
                        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                        row["Variant"] = variant_key

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

                        # Add consequence sort priority for custom sorting
                        row["_consequence_priority"] = get_highest_priority_consequence(
                            consequence_str
                        )

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

                    # All available columns (add Variant and Validation columns, exclude n_grouped from display)
                    all_columns = reorder_columns_by_group(
                        ["Variant"]
                        + [
                            col
                            for col in df.columns
                            if col not in ["#CHROM", "POS", "REF", "ALT", "n_grouped"]
                        ]
                        + ["Validation"]
                    )

                    wombat_data[config_name]["all_columns"] = all_columns

                    # Auto-select preset based on keywords
                    initial_preset = select_preset_for_config(config_name, VIEW_PRESETS)
                    selected_preset = {"name": initial_preset["name"]}

                    # Override initial_selected with preset columns (if available)
                    preset_columns = initial_preset.get("columns", [])
                    initial_selected = [col for col in preset_columns if col in all_columns]

                    selected_cols = {"value": initial_selected}
                    wombat_data[config_name]["selected_cols"] = selected_cols
                    wombat_data[config_name]["selected_preset"] = selected_preset

                    # Filter state for this config
                    wombat_data[config_name]["filter_exclude_lcr"] = {"value": True}
                    wombat_data[config_name]["filter_exclude_gnomad"] = {"value": True}
                    wombat_data[config_name]["filter_exclude_mnv"] = {"value": True}
                    wombat_data[config_name]["selected_genesets"] = {"value": []}
                    wombat_data[config_name]["selected_impacts"] = {
                        "value": list(VEP_CONSEQUENCES.keys())
                    }
                    wombat_data[config_name]["selected_validations"] = {
                        "value": [
                            "present",
                            "absent",
                            "uncertain",
                            "conflicting",
                            "TODO",
                        ]
                    }
                    # Mutable containers for UI element references
                    # (populated during filter panel construction, used by handlers)
                    wombat_data[config_name]["_refresh"] = {"fn": None}
                    wombat_data[config_name]["_geneset_cbs"] = {}
                    wombat_data[config_name]["_impact_cbs"] = {}

                    # Create a container for the data table
                    data_container = ui.column().classes("w-full")

                    # Capture the client context for use in callbacks
                    from nicegui import context

                    page_client = context.client

                    with data_container:

                        # Collapsible Filters panel
                        with ui.card().classes("w-full p-0 mt-2"):
                          with ui.expansion(
                              "Filters", icon="filter_list", value=False
                          ).classes("w-full").props(
                              "header-class='text-lg font-semibold text-blue-700'"
                          ):

                            with ui.row().classes(
                                "gap-6 items-start flex-wrap p-2"
                            ):
                                # Checkbox exclude filters
                                with ui.column().classes("gap-2"):
                                    ui.checkbox(
                                        "Exclude LCR",
                                        value=wombat_data[config_name][
                                            "filter_exclude_lcr"
                                        ]["value"],
                                        on_change=lambda e, cfg=config_name: (
                                            wombat_data[cfg][
                                                "filter_exclude_lcr"
                                            ].update({"value": e.value}),
                                            wombat_data[cfg]["_refresh"]["fn"](),
                                        ),
                                    )

                                    ui.checkbox(
                                        "Exclude gnomAD filtered",
                                        value=wombat_data[config_name][
                                            "filter_exclude_gnomad"
                                        ]["value"],
                                        on_change=lambda e, cfg=config_name: (
                                            wombat_data[cfg][
                                                "filter_exclude_gnomad"
                                            ].update({"value": e.value}),
                                            wombat_data[cfg]["_refresh"]["fn"](),
                                        ),
                                    )

                                    ui.checkbox(
                                        "Exclude FP due to MNV (mnv_proba > 0.5)",
                                        value=wombat_data[config_name][
                                            "filter_exclude_mnv"
                                        ]["value"],
                                        on_change=lambda e, cfg=config_name: (
                                            wombat_data[cfg][
                                                "filter_exclude_mnv"
                                            ].update({"value": e.value}),
                                            wombat_data[cfg]["_refresh"]["fn"](),
                                        ),
                                    )

                                # Dropdown menu filters
                                with ui.row().classes(
                                    "gap-2 items-center flex-wrap"
                                ):
                                    # Genesets filter
                                    if available_genesets:
                                        geneset_btn_ref: Dict[str, Any] = {
                                            "button": None
                                        }

                                        geneset_btn = ui.button(
                                            "Genesets", icon="list"
                                        ).props("outline")
                                        geneset_btn_ref["button"] = geneset_btn

                                        with geneset_btn:
                                            with ui.menu():
                                                ui.label(
                                                    "Select Genesets:"
                                                ).classes(
                                                    "px-4 py-2 font-semibold text-sm"
                                                )
                                                ui.separator()

                                                with ui.column().classes("p-2"):
                                                    with ui.row().classes(
                                                        "gap-2 mb-2"
                                                    ):

                                                        def select_all_genesets(
                                                            _e=None,
                                                            cfg=config_name,
                                                        ):
                                                            wombat_data[cfg][
                                                                "selected_genesets"
                                                            ]["value"] = list(
                                                                available_genesets.keys()
                                                            )
                                                            for (
                                                                cb
                                                            ) in wombat_data[cfg]["_geneset_cbs"].values():
                                                                cb.value = True
                                                            if geneset_btn_ref[
                                                                "button"
                                                            ]:
                                                                geneset_btn_ref[
                                                                    "button"
                                                                ].props(
                                                                    remove="outline",
                                                                    add="unelevated color=green",
                                                                )
                                                                geneset_btn_ref[
                                                                    "button"
                                                                ].update()
                                                            wombat_data[cfg]["_refresh"]["fn"]()

                                                        def select_no_genesets(
                                                            _e=None,
                                                            cfg=config_name,
                                                        ):
                                                            wombat_data[cfg][
                                                                "selected_genesets"
                                                            ]["value"] = []
                                                            for (
                                                                cb
                                                            ) in wombat_data[cfg]["_geneset_cbs"].values():
                                                                cb.value = (
                                                                    False
                                                                )
                                                            if geneset_btn_ref[
                                                                "button"
                                                            ]:
                                                                geneset_btn_ref[
                                                                    "button"
                                                                ].props(
                                                                    remove="unelevated color=green",
                                                                    add="outline",
                                                                )
                                                                geneset_btn_ref[
                                                                    "button"
                                                                ].update()
                                                            wombat_data[cfg]["_refresh"]["fn"]()

                                                        ui.button(
                                                            "All",
                                                            on_click=select_all_genesets,
                                                        ).props(
                                                            "size=sm flat dense"
                                                        ).classes("text-xs")
                                                        ui.button(
                                                            "None",
                                                            on_click=select_no_genesets,
                                                        ).props(
                                                            "size=sm flat dense"
                                                        ).classes("text-xs")

                                                    ui.separator()

                                                    for gs_name in sorted(
                                                        available_genesets.keys()
                                                    ):

                                                        def make_geneset_handler(
                                                            name,
                                                            cfg=config_name,
                                                        ):
                                                            def handler(e):
                                                                if e.value:
                                                                    if (
                                                                        name
                                                                        not in wombat_data[
                                                                            cfg
                                                                        ][
                                                                            "selected_genesets"
                                                                        ][
                                                                            "value"
                                                                        ]
                                                                    ):
                                                                        wombat_data[
                                                                            cfg
                                                                        ][
                                                                            "selected_genesets"
                                                                        ][
                                                                            "value"
                                                                        ].append(
                                                                            name
                                                                        )
                                                                else:
                                                                    if (
                                                                        name
                                                                        in wombat_data[
                                                                            cfg
                                                                        ][
                                                                            "selected_genesets"
                                                                        ][
                                                                            "value"
                                                                        ]
                                                                    ):
                                                                        wombat_data[
                                                                            cfg
                                                                        ][
                                                                            "selected_genesets"
                                                                        ][
                                                                            "value"
                                                                        ].remove(
                                                                            name
                                                                        )
                                                                if geneset_btn_ref[
                                                                    "button"
                                                                ]:
                                                                    if wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_genesets"
                                                                    ][
                                                                        "value"
                                                                    ]:
                                                                        geneset_btn_ref[
                                                                            "button"
                                                                        ].props(
                                                                            remove="outline",
                                                                            add="unelevated color=green",
                                                                        )
                                                                    else:
                                                                        geneset_btn_ref[
                                                                            "button"
                                                                        ].props(
                                                                            remove="unelevated color=green",
                                                                            add="outline",
                                                                        )
                                                                    geneset_btn_ref[
                                                                        "button"
                                                                    ].update()
                                                                wombat_data[cfg]["_refresh"]["fn"]()

                                                            return handler

                                                        wombat_data[config_name]["_geneset_cbs"][
                                                            gs_name
                                                        ] = ui.checkbox(
                                                            f"{gs_name} ({len(available_genesets[gs_name])} genes)",
                                                            value=False,
                                                            on_change=make_geneset_handler(
                                                                gs_name
                                                            ),
                                                        ).classes("text-sm")

                                    # Impacts filter
                                    impact_btn_ref: Dict[str, Any] = {
                                        "button": None
                                    }
                                    impact_btn = ui.button(
                                        "Impacts", icon="filter_list"
                                    ).props("outline")
                                    impact_btn_ref["button"] = impact_btn

                                    with impact_btn:
                                        with ui.menu():
                                            ui.label(
                                                "Select Impact Types:"
                                            ).classes(
                                                "px-4 py-2 font-semibold text-sm"
                                            )
                                            ui.separator()

                                            with ui.column().classes("p-2"):
                                                with ui.row().classes(
                                                    "gap-2 mb-2 flex-wrap"
                                                ):

                                                    def select_all_impacts(
                                                        _e=None,
                                                        cfg=config_name,
                                                    ):
                                                        wombat_data[cfg][
                                                            "selected_impacts"
                                                        ]["value"] = list(
                                                            VEP_CONSEQUENCES.keys()
                                                        )
                                                        for (
                                                            cb
                                                        ) in wombat_data[cfg]["_impact_cbs"].values():
                                                            cb.value = True
                                                        if impact_btn_ref[
                                                            "button"
                                                        ]:
                                                            impact_btn_ref[
                                                                "button"
                                                            ].props(
                                                                remove="unelevated color=orange",
                                                                add="outline",
                                                            )
                                                            impact_btn_ref[
                                                                "button"
                                                            ].update()
                                                        wombat_data[cfg]["_refresh"]["fn"]()

                                                    def select_none_impacts(
                                                        _e=None,
                                                        cfg=config_name,
                                                    ):
                                                        wombat_data[cfg][
                                                            "selected_impacts"
                                                        ]["value"] = []
                                                        for (
                                                            cb
                                                        ) in wombat_data[cfg]["_impact_cbs"].values():
                                                            cb.value = False
                                                        if impact_btn_ref[
                                                            "button"
                                                        ]:
                                                            impact_btn_ref[
                                                                "button"
                                                            ].props(
                                                                remove="outline",
                                                                add="unelevated color=orange",
                                                            )
                                                            impact_btn_ref[
                                                                "button"
                                                            ].update()
                                                        wombat_data[cfg]["_refresh"]["fn"]()

                                                    def select_by_impact_level(
                                                        level: str,
                                                        *,
                                                        cfg=config_name,
                                                    ):
                                                        selected = [
                                                            cons
                                                            for cons, (
                                                                imp,
                                                                _,
                                                            ) in VEP_CONSEQUENCES.items()
                                                            if imp == level
                                                        ]
                                                        wombat_data[cfg][
                                                            "selected_impacts"
                                                        ]["value"] = selected
                                                        for (
                                                            impact,
                                                            cb,
                                                        ) in wombat_data[cfg]["_impact_cbs"].items():
                                                            cb.value = (
                                                                impact
                                                                in selected
                                                            )
                                                        if impact_btn_ref[
                                                            "button"
                                                        ]:
                                                            impact_btn_ref[
                                                                "button"
                                                            ].props(
                                                                remove="outline",
                                                                add="unelevated color=orange",
                                                            )
                                                            impact_btn_ref[
                                                                "button"
                                                            ].update()
                                                        wombat_data[cfg]["_refresh"]["fn"]()

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
                                                        "HIGH",
                                                        on_click=lambda _e=None, fn=select_by_impact_level: fn(
                                                            "HIGH"
                                                        ),
                                                    ).props(
                                                        "size=sm flat dense color=red"
                                                    ).classes("text-xs")
                                                    ui.button(
                                                        "MODERATE",
                                                        on_click=lambda _e=None, fn=select_by_impact_level: fn(
                                                            "MODERATE"
                                                        ),
                                                    ).props(
                                                        "size=sm flat dense color=orange"
                                                    ).classes("text-xs")
                                                    ui.button(
                                                        "LOW",
                                                        on_click=lambda _e=None, fn=select_by_impact_level: fn(
                                                            "LOW"
                                                        ),
                                                    ).props(
                                                        "size=sm flat dense color=yellow-8"
                                                    ).classes("text-xs")
                                                    ui.button(
                                                        "MODIFIER",
                                                        on_click=lambda _e=None, fn=select_by_impact_level: fn(
                                                            "MODIFIER"
                                                        ),
                                                    ).props(
                                                        "size=sm flat dense color=grey"
                                                    ).classes("text-xs")

                                                ui.separator()

                                                with ui.column().classes(
                                                    "gap-1"
                                                ):

                                                    def make_impact_handler(
                                                        cons,
                                                        cfg=config_name,
                                                    ):
                                                        def handler(e):
                                                            if e.value:
                                                                if (
                                                                    cons
                                                                    not in wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_impacts"
                                                                    ][
                                                                        "value"
                                                                    ]
                                                                ):
                                                                    wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_impacts"
                                                                    ][
                                                                        "value"
                                                                    ].append(
                                                                        cons
                                                                    )
                                                            else:
                                                                if (
                                                                    cons
                                                                    in wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_impacts"
                                                                    ][
                                                                        "value"
                                                                    ]
                                                                ):
                                                                    wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_impacts"
                                                                    ][
                                                                        "value"
                                                                    ].remove(
                                                                        cons
                                                                    )
                                                            if impact_btn_ref[
                                                                "button"
                                                            ]:
                                                                if len(
                                                                    wombat_data[
                                                                        cfg
                                                                    ][
                                                                        "selected_impacts"
                                                                    ][
                                                                        "value"
                                                                    ]
                                                                ) == len(
                                                                    VEP_CONSEQUENCES
                                                                ):
                                                                    impact_btn_ref[
                                                                        "button"
                                                                    ].props(
                                                                        remove="unelevated color=orange",
                                                                        add="outline",
                                                                    )
                                                                else:
                                                                    impact_btn_ref[
                                                                        "button"
                                                                    ].props(
                                                                        remove="outline",
                                                                        add="unelevated color=orange",
                                                                    )
                                                                impact_btn_ref[
                                                                    "button"
                                                                ].update()
                                                            wombat_data[cfg]["_refresh"]["fn"]()

                                                        return handler

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
                                                            if imp
                                                            == impact_level
                                                        ]
                                                        if consequences:
                                                            ui.label(
                                                                f"{impact_level}:"
                                                            ).classes(
                                                                "text-xs font-bold text-gray-600 mt-2"
                                                            )
                                                            for cons in sorted(
                                                                consequences
                                                            ):
                                                                wombat_data[config_name]["_impact_cbs"][
                                                                    cons
                                                                ] = ui.checkbox(
                                                                    format_consequence_display(
                                                                        cons
                                                                    ),
                                                                    value=True,
                                                                    on_change=make_impact_handler(
                                                                        cons
                                                                    ),
                                                                ).classes(
                                                                    "text-sm"
                                                                )

                                    # Validation filter
                                    create_validation_filter_menu(
                                        all_statuses=[
                                            "present",
                                            "absent",
                                            "uncertain",
                                            "conflicting",
                                            "TODO",
                                        ],
                                        filter_state=wombat_data[config_name][
                                            "selected_validations"
                                        ],
                                        on_change=lambda cfg=config_name: wombat_data[cfg]["_refresh"]["fn"](),
                                        label="Validation",
                                        button_classes="",
                                    )

                        @ui.refreshable
                        def render_data_table(cfg=config_name):
                            data = wombat_data[cfg]
                            df_local = data["df"]
                            all_rows_local = data["all_rows"]
                            all_columns_local = data["all_columns"]
                            selected_cols_local = data["selected_cols"]

                            # Filter rows by selected members
                            if "sample" in df_local.columns:
                                rows = [
                                    r
                                    for r in all_rows_local
                                    if r.get("sample") in selected_members["value"]
                                ]
                            else:
                                rows = all_rows_local

                            total_before_filters = len(rows)

                            # Apply checkbox exclude filters
                            if data["filter_exclude_lcr"]["value"]:
                                rows = [
                                    r
                                    for r in rows
                                    if not (
                                        r.get("LCR")
                                        and "true"
                                        in str(r.get("LCR", "")).lower()
                                    )
                                ]

                            if data["filter_exclude_gnomad"]["value"]:
                                rows = [
                                    r
                                    for r in rows
                                    if not r.get("genomes_filters")
                                ]

                            if data["filter_exclude_mnv"]["value"]:

                                def has_high_mnv_proba(row):
                                    val = row.get("mnv_proba", "")
                                    if not val:
                                        return False
                                    for part in str(val).split(","):
                                        part = part.strip()
                                        if part:
                                            try:
                                                if float(part) > 0.5:
                                                    return True
                                            except (ValueError, TypeError):
                                                pass
                                    return False

                                rows = [
                                    r
                                    for r in rows
                                    if not has_high_mnv_proba(r)
                                ]

                            # Apply geneset filter
                            selected_gs = data["selected_genesets"]["value"]
                            if selected_gs:
                                combined_genes: set = set()
                                for gs_name in selected_gs:
                                    combined_genes.update(
                                        available_genesets.get(gs_name, set())
                                    )
                                rows = [
                                    r
                                    for r in rows
                                    if any(
                                        s.strip().upper() in combined_genes
                                        for s in str(
                                            r.get("VEP_SYMBOL", "")
                                        ).split(",")
                                        if s.strip()
                                    )
                                ]

                            # Apply impact filter (only if not all selected)
                            selected_imps = data["selected_impacts"]["value"]
                            if selected_imps and set(selected_imps) != set(
                                VEP_CONSEQUENCES.keys()
                            ):

                                def row_matches_impact(row):
                                    consequence_str = row.get(
                                        "VEP_Consequence", ""
                                    )
                                    if not consequence_str:
                                        return False
                                    for part in str(consequence_str).split(
                                        ","
                                    ):
                                        for cons in part.split("&"):
                                            cons = cons.strip()
                                            if (
                                                cons
                                                and cons in selected_imps
                                            ):
                                                return True
                                    return False

                                rows = [
                                    r for r in rows if row_matches_impact(r)
                                ]

                            # Apply validation filter
                            selected_vals = data["selected_validations"][
                                "value"
                            ]
                            if selected_vals:
                                rows = [
                                    r
                                    for r in rows
                                    if r.get("Validation", "")
                                    in selected_vals
                                    or (
                                        "TODO" in selected_vals
                                        and not r.get("Validation")
                                    )
                                ]

                            def make_columns():
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
                                for col in all_columns_local:
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
                                        col_def["sortField"] = "_consequence_priority"
                                    elif col == "VEP_CLIN_SIG":
                                        col_def["cellType"] = "badge_list"
                                        col_def["badgesField"] = "ClinVarBadges"
                                    elif col == "VEP_SYMBOL":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "GeneBadges"
                                    elif col == "VEP_Gene":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "VEP_Gene_badges"
                                    else:
                                        col_def["cellType"] = "score_badge"
                                    cols.append(col_def)
                                return cols

                            def _apply_col_visibility():
                                """Push current column selection to JS table."""
                                if data.get("_dt"):
                                    visible = ["actions"] + list(selected_cols_local["value"])
                                    data["_dt"].set_column_visibility(visible)

                            with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
                                row_label = (
                                    f"Data ({len(rows)} / {total_before_filters} rows)"
                                    if len(rows) < total_before_filters
                                    else f"Data ({len(rows)} rows)"
                                )
                                ui.label(row_label).classes(
                                    "text-lg font-semibold text-blue-700"
                                )

                                # Preset selector dropdown
                                preset_select = ui.select(
                                    options={p["name"]: p["name"] for p in VIEW_PRESETS},
                                    value=data["selected_preset"]["name"],
                                    label="Preset"
                                ).classes("w-48")

                                ui.space()  # Push remaining items to the right

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
                                            with ui.row().classes("gap-2 mb-2"):
                                                checkboxes: Dict[str, Any] = {}

                                                def _sync_checkboxes():
                                                    for col, cb in checkboxes.items():
                                                        cb.value = col in selected_cols_local["value"]

                                                def select_all():
                                                    selected_cols_local["value"] = list(
                                                        all_columns_local
                                                    )
                                                    _apply_col_visibility()
                                                    _sync_checkboxes()

                                                def select_none():
                                                    selected_cols_local["value"] = []
                                                    _apply_col_visibility()
                                                    _sync_checkboxes()

                                                ui.button(
                                                    "All", on_click=select_all
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )
                                                ui.button(
                                                    "None", on_click=select_none
                                                ).props("size=sm flat dense").classes(
                                                    "text-xs"
                                                )

                                            ui.separator()

                                            for col in all_columns_local:
                                                checkboxes[col] = ui.checkbox(
                                                    col,
                                                    value=col
                                                    in selected_cols_local["value"],
                                                    on_change=lambda e,
                                                    c=col: handle_col_change(
                                                        c, e.value
                                                    ),
                                                ).classes("text-sm")

                                # Stats button
                                def show_stats_dialog(current_rows=rows):
                                    # Deduplicate by (#CHROM, POS, REF, ALT) for unique variants
                                    seen = set()
                                    unique_variants = []
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

                                    # Chromosome distribution
                                    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
                                    chrom_counts: Dict[str, int] = {c: 0 for c in chrom_order}
                                    for r in unique_variants:
                                        chrom = str(r.get("#CHROM", "")).replace("chr", "").upper()
                                        if chrom == "M":
                                            chrom = "MT"
                                        if chrom in chrom_counts:
                                            chrom_counts[chrom] += 1

                                    # Consequence distribution (highest per variant)
                                    from collections import Counter

                                    consequence_counts = Counter(
                                        get_highest_consequence_term(
                                            str(r.get("VEP_Consequence", ""))
                                        )
                                        for r in unique_variants
                                    )

                                    # Validation distribution
                                    validation_counts = Counter(
                                        r.get("Validation", "") or "TODO"
                                        for r in unique_variants
                                    )

                                    # Chromosome sizes (hg38, in Mb) for ideogram
                                    chrom_sizes_mb = {
                                        "1": 249, "2": 242, "3": 198, "4": 190, "5": 182,
                                        "6": 171, "7": 159, "8": 145, "9": 138, "10": 134,
                                        "11": 135, "12": 133, "13": 114, "14": 107, "15": 102,
                                        "16": 90, "17": 83, "18": 80, "19": 59, "20": 64,
                                        "21": 47, "22": 51, "X": 156, "Y": 57, "MT": 1,
                                    }

                                    # Scatter data for ideogram: [position_mb, chromosome]
                                    scatter_data = []
                                    for r in unique_variants:
                                        chrom = str(r.get("#CHROM", "")).replace("chr", "").upper()
                                        if chrom == "M":
                                            chrom = "MT"
                                        pos = r.get("POS", 0)
                                        try:
                                            pos_mb = round(float(pos) / 1_000_000, 2)
                                        except (ValueError, TypeError):
                                            continue
                                        if chrom in chrom_sizes_mb:
                                            scatter_data.append([pos_mb, chrom])

                                    with ui.dialog().props(
                                        "full-width"
                                    ) as stats_dialog, ui.card().classes(
                                        "w-full"
                                    ):
                                        with ui.column().classes("w-full p-4"):
                                            # Header with toggle and close
                                            show_ideogram = {"value": False}

                                            with ui.row().classes(
                                                "items-center justify-between w-full mb-2"
                                            ):
                                                with ui.row().classes(
                                                    "items-center gap-3"
                                                ):
                                                    with ui.column().classes(
                                                        "gap-0"
                                                    ):
                                                        ui.label(
                                                            "Variant Statistics"
                                                        ).classes(
                                                            "text-xl font-bold text-blue-900"
                                                        )
                                                        ui.label(
                                                            f"Based on {len(unique_variants)} unique variants "
                                                            "(deduplicated by #CHROM / POS / REF / ALT)"
                                                        ).classes(
                                                            "text-sm text-gray-500"
                                                        )
                                                    ideogram_btn = ui.button(
                                                        "Ideogram",
                                                    ).props(
                                                        "outline color=blue size=sm dense no-caps"
                                                    )
                                                ui.button(
                                                    icon="close",
                                                    on_click=lambda: stats_dialog.close(),
                                                ).props("flat round")

                                            # Charts container (visible by default)
                                            charts_container = ui.column().classes(
                                                "w-full"
                                            )
                                            with charts_container:
                                                # Chromosome bar chart
                                                ui.label(
                                                    "Variants per Chromosome"
                                                ).classes(
                                                    "text-lg font-semibold text-gray-800 mt-2"
                                                )
                                                ui.echart(
                                                    {
                                                        "tooltip": {},
                                                        "xAxis": {
                                                            "type": "category",
                                                            "data": chrom_order,
                                                            "name": "Chromosome",
                                                        },
                                                        "yAxis": {
                                                            "type": "value",
                                                            "name": "Count",
                                                        },
                                                        "series": [
                                                            {
                                                                "type": "bar",
                                                                "data": [
                                                                    chrom_counts[c]
                                                                    for c in chrom_order
                                                                ],
                                                                "itemStyle": {
                                                                    "color": "#3b82f6"
                                                                },
                                                            }
                                                        ],
                                                    }
                                                ).classes("w-full h-64")

                                                # Pie charts side by side
                                                with ui.row().classes(
                                                    "w-full gap-4 flex-wrap mt-4"
                                                ):
                                                    # Consequence pie chart
                                                    with ui.column().classes(
                                                        "flex-1 min-w-[400px]"
                                                    ):
                                                        ui.label(
                                                            "Consequence Distribution (highest per variant)"
                                                        ).classes(
                                                            "text-lg font-semibold text-gray-800"
                                                        )
                                                        cons_data = [
                                                            {
                                                                "name": format_consequence_display(
                                                                    cons
                                                                ),
                                                                "value": count,
                                                                "itemStyle": {
                                                                    "color": VEP_CONSEQUENCES.get(
                                                                        cons,
                                                                        (
                                                                            "",
                                                                            "#6b7280",
                                                                        ),
                                                                    )[1]
                                                                },
                                                            }
                                                            for cons, count in consequence_counts.most_common()
                                                        ]
                                                        ui.echart(
                                                            {
                                                                "tooltip": {
                                                                    "trigger": "item"
                                                                },
                                                                "series": [
                                                                    {
                                                                        "type": "pie",
                                                                        "radius": "70%",
                                                                        "data": cons_data,
                                                                        "label": {
                                                                            "formatter": "{b}: {c} ({d}%)"
                                                                        },
                                                                    }
                                                                ],
                                                            }
                                                        ).classes("w-full h-80")

                                                    # Validation pie chart
                                                    with ui.column().classes(
                                                        "flex-1 min-w-[400px]"
                                                    ):
                                                        ui.label(
                                                            "Validation Status Distribution"
                                                        ).classes(
                                                            "text-lg font-semibold text-gray-800"
                                                        )
                                                        validation_colors = {
                                                            "present": "#22c55e",
                                                            "absent": "#ef4444",
                                                            "uncertain": "#f59e0b",
                                                            "conflicting": "#fbbf24",
                                                            "TODO": "#6b7280",
                                                        }
                                                        val_data = [
                                                            {
                                                                "name": status,
                                                                "value": count,
                                                                "itemStyle": {
                                                                    "color": validation_colors.get(
                                                                        status,
                                                                        "#6b7280",
                                                                    )
                                                                },
                                                            }
                                                            for status, count in validation_counts.most_common()
                                                        ]
                                                        ui.echart(
                                                            {
                                                                "tooltip": {
                                                                    "trigger": "item"
                                                                },
                                                                "series": [
                                                                    {
                                                                        "type": "pie",
                                                                        "radius": "70%",
                                                                        "data": val_data,
                                                                        "label": {
                                                                            "formatter": "{b}: {c} ({d}%)"
                                                                        },
                                                                    }
                                                                ],
                                                            }
                                                        ).classes("w-full h-80")

                                            # Ideogram container (hidden by default)
                                            ideo_container = ui.column().classes(
                                                "w-full"
                                            )
                                            ideo_container.set_visibility(False)
                                            with ideo_container:
                                                # Build SVG ideogram with cytobands
                                                svg_w = 1800
                                                lbl_w = 50
                                                plot_w = svg_w - lbl_w - 20
                                                row_h = 24
                                                row_gap = 6
                                                svg_h = (
                                                    len(chrom_order)
                                                    * (row_h + row_gap)
                                                    + 40
                                                )
                                                max_mb = max(
                                                    chrom_sizes_mb.values()
                                                )

                                                svg_parts = [
                                                    f'<svg viewBox="0 0 {svg_w} {svg_h}" '
                                                    f'xmlns="http://www.w3.org/2000/svg" '
                                                    f'preserveAspectRatio="xMinYMin meet" '
                                                    f'style="font-family: sans-serif; width: 100%; height: auto;">'
                                                ]

                                                # Grid lines and x-axis labels
                                                axis_y = len(chrom_order) * (
                                                    row_h + row_gap
                                                )
                                                for mb_val in range(0, 260, 50):
                                                    gx = (
                                                        lbl_w
                                                        + (mb_val / max_mb)
                                                        * plot_w
                                                    )
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

                                                for ci, chrom in enumerate(
                                                    chrom_order
                                                ):
                                                    cy = ci * (row_h + row_gap)
                                                    bands = CYTOBANDS.get(
                                                        chrom, []
                                                    )
                                                    cs = chrom_sizes_mb.get(
                                                        chrom, 0
                                                    )
                                                    total_w = (
                                                        cs / max_mb
                                                    ) * plot_w

                                                    # Chromosome label
                                                    svg_parts.append(
                                                        f'<text x="{lbl_w - 6}" y="{cy + row_h * 0.7}" '
                                                        f'text-anchor="end" font-size="13" '
                                                        f'fill="#374151">{chrom}</text>'
                                                    )

                                                    # Cytoband rectangles
                                                    for band in bands:
                                                        bx = (
                                                            lbl_w
                                                            + (
                                                                band["start"]
                                                                / max_mb
                                                            )
                                                            * plot_w
                                                        )
                                                        bw = max(
                                                            (
                                                                (
                                                                    band["end"]
                                                                    - band[
                                                                        "start"
                                                                    ]
                                                                )
                                                                / max_mb
                                                            )
                                                            * plot_w,
                                                            0.5,
                                                        )
                                                        color = (
                                                            GIESTAIN_COLORS.get(
                                                                band["stain"],
                                                                "#e5e7eb",
                                                            )
                                                        )
                                                        # Same height for all bands including centromeres
                                                        svg_parts.append(
                                                            f'<rect x="{bx:.1f}" y="{cy}" '
                                                            f'width="{bw:.1f}" height="{row_h}" '
                                                            f'fill="{color}"/>'
                                                        )

                                                    # Chromosome outline
                                                    svg_parts.append(
                                                        f'<rect x="{lbl_w}" y="{cy}" '
                                                        f'width="{total_w:.1f}" height="{row_h}" '
                                                        f'fill="none" stroke="#9ca3af" '
                                                        f'stroke-width="0.5" rx="4"/>'
                                                    )

                                                # Variant markers (blue triangles)
                                                for sd in scatter_data:
                                                    v_mb, v_chrom = (
                                                        sd[0],
                                                        sd[1],
                                                    )
                                                    if (
                                                        v_chrom
                                                        not in chrom_order
                                                    ):
                                                        continue
                                                    v_idx = chrom_order.index(
                                                        v_chrom
                                                    )
                                                    vy = v_idx * (
                                                        row_h + row_gap
                                                    )
                                                    vx = (
                                                        lbl_w
                                                        + (v_mb / max_mb)
                                                        * plot_w
                                                    )
                                                    svg_parts.append(
                                                        f'<line x1="{vx:.1f}" y1="{vy - 2}" '
                                                        f'x2="{vx:.1f}" y2="{vy + row_h + 2}" '
                                                        f'stroke="#2563eb" stroke-width="1.5" '
                                                        f'opacity="0.8"/>'
                                                    )

                                                svg_parts.append("</svg>")
                                                ui.html(
                                                    "\n".join(svg_parts),
                                                    sanitize=False,
                                                ).classes("w-full")

                                            # Toggle handler
                                            def toggle_ideogram(
                                                _e=None,
                                                _charts=charts_container,
                                                _ideo=ideo_container,
                                                _btn=ideogram_btn,
                                                _state=show_ideogram,
                                            ):
                                                _state["value"] = not _state[
                                                    "value"
                                                ]
                                                _charts.set_visibility(
                                                    not _state["value"]
                                                )
                                                _ideo.set_visibility(
                                                    _state["value"]
                                                )
                                                if _state["value"]:
                                                    _btn.props(
                                                        remove="outline",
                                                        add="unelevated",
                                                    )
                                                else:
                                                    _btn.props(
                                                        remove="unelevated",
                                                        add="outline",
                                                    )
                                                _btn.update()

                                            ideogram_btn.on_click(
                                                toggle_ideogram
                                            )

                                    stats_dialog.open()

                                ui.button(
                                    "Stats", icon="bar_chart", on_click=show_stats_dialog
                                ).props("outline color=blue size=sm")

                            def on_view_variant(e):
                                row_data = e.get("row", {})
                                chrom = row_data.get("#CHROM", "")
                                pos = row_data.get("POS", "")
                                ref = row_data.get("REF", "")
                                alt = row_data.get("ALT", "")
                                sample_val = row_data.get("sample", "")

                                # Callback to update the Validation column in the table
                                def on_save(validation_status: str):
                                    # Reload validation data from file
                                    validation_file = (
                                        store.data_dir / "validations" / "snvs.tsv"
                                    )
                                    validation_map = load_validation_map(
                                        validation_file, family_id
                                    )
                                    # Re-add validation status to all rows
                                    for row in all_rows:
                                        chrom = row.get("#CHROM", "")
                                        pos = row.get("POS", "")
                                        ref = row.get("REF", "")
                                        alt = row.get("ALT", "")
                                        sample_id = row.get("sample", "")
                                        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                                        add_validation_status_to_row(
                                            row,
                                            validation_map,
                                            variant_key,
                                            sample_id,
                                        )
                                    # Refresh the table using the captured client context
                                    with page_client:
                                        ui.timer(
                                            0.1,
                                            render_data_table.refresh,
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
                                    sample=sample_val,
                                    variant_data=row_data,
                                    on_save_callback=on_save,
                                )

                            data["_dt"] = DataTable(
                                columns=make_columns(),
                                rows=rows,
                                row_key="Variant",
                                pagination={"rowsPerPage": 10},
                                visible_columns=["actions"] + list(selected_cols_local["value"]),
                                on_row_action=on_view_variant,
                            )

                            def on_preset_change(e):
                                """Handle preset selection change."""
                                preset_name = e.value

                                # Find the selected preset
                                preset = next((p for p in VIEW_PRESETS if p["name"] == preset_name), None)
                                if not preset:
                                    return

                                # Filter columns to only those available in the data
                                available = [col for col in preset.get("columns", [])
                                             if col in all_columns_local]

                                selected_cols_local["value"] = available
                                data["selected_preset"]["name"] = preset_name
                                _apply_col_visibility()
                                _sync_checkboxes()

                            # Connect preset change handler
                            preset_select.on_value_change(on_preset_change)

                            def handle_col_change(col_name, is_checked):
                                if (
                                    is_checked
                                    and col_name not in selected_cols_local["value"]
                                ):
                                    selected_cols_local["value"].append(col_name)
                                elif (
                                    not is_checked
                                    and col_name in selected_cols_local["value"]
                                ):
                                    selected_cols_local["value"].remove(col_name)

                                # Reorder to match all_columns_local order
                                selected_cols_local["value"] = [
                                    col for col in all_columns_local
                                    if col in selected_cols_local["value"]
                                ]

                                _apply_col_visibility()

                        # Store refresh reference for filter callbacks
                        wombat_data[config_name]["_refresh"]["fn"] = render_data_table.refresh
                        data_table_refreshers.append(render_data_table.refresh)
                        render_data_table()

                except Exception as e:
                    ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
