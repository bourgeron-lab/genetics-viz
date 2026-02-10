"""Wombat tab component for family page."""

import re
from pathlib import Path
from typing import Any, Callable, Dict, List

import polars as pl
import yaml
from nicegui import ui

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.gene_scoring import get_gene_scorer


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


# Table slot template for wombat data with view button and validation icons
WOMBAT_TABLE_SLOT = r"""
    <q-tr :props="props">
        <q-td key="actions" :props="props">
            <div style="display: flex; align-items: center; gap: 4px;">
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
                <q-badge 
                    v-if="props.row.n_grouped && props.row.n_grouped > 1"
                    :label="props.row.n_grouped.toString()"
                    color="orange"
                    style="font-size: 11px;"
                >
                    <q-tooltip>
                        {{ props.row.n_grouped }} transcripts collapsed
                        <span v-if="props.row.VEP_SYMBOL"> for genes: {{ props.row.VEP_SYMBOL }}</span>
                    </q-tooltip>
                </q-badge>
            </div>
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
            <template v-else-if="col.name === 'VEP_SYMBOL'">
                <div style="display: flex; flex-wrap: wrap; gap: 4px;">
                    <q-badge v-for="(badge, idx) in (props.row.GeneBadges || [])" :key="idx" 
                             :label="badge.label"
                             :style="'background-color: ' + badge.color + '; color: ' + (badge.color === '#ffffff' ? 'black' : 'white') + '; font-size: 0.875em; padding: 4px 8px;'">
                        <q-tooltip>{{ badge.tooltip }}</q-tooltip>
                    </q-badge>
                </div>
            </template>
            <template v-else-if="col.name === 'VEP_Gene'">
                <div style="display: flex; flex-wrap: wrap; gap: 4px;">
                    <q-badge v-for="(badge, idx) in (props.row.VEP_Gene_badges || [])" :key="idx" 
                             :label="badge.label"
                             :style="'background-color: ' + badge.color + '; color: ' + (badge.color === '#ffffff' ? 'black' : 'white') + '; font-size: 0.875em; padding: 4px 8px;'">
                        <q-tooltip>{{ badge.tooltip }}</q-tooltip>
                    </q-badge>
                </div>
            </template>
            <template v-else>
                {{ col.value }}
            </template>
        </q-td>
    </q-tr>
"""


def get_wombat_display_label(col: str) -> str:
    """Get display label for wombat column, removing VEP_ prefix and renaming gnomAD columns."""
    if col == "fafmax_faf95_max_genomes":
        return "gnomAD 4.1 WGS"
    elif col == "nhomalt_genomes":
        return "gnomAD 4.1 nhomalt WGS"
    elif col == "VEP_CLIN_SIG":
        return "ClinVar"
    elif col.startswith("VEP_"):
        return col[4:]  # Remove VEP_ prefix
    else:
        return col


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
                    all_columns = (
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

                    # Create a container for the data table
                    data_container = ui.column().classes("w-full")

                    # Capture the client context for use in callbacks
                    from nicegui import context

                    page_client = context.client

                    with data_container:

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

                            def make_columns(visible_cols):
                                cols = [
                                    {
                                        "name": "actions",
                                        "label": "",
                                        "field": "actions",
                                        "sortable": False,
                                        "align": "center",
                                    }
                                ]
                                for col in visible_cols:
                                    col_def = {
                                        "name": col,
                                        "label": get_wombat_display_label(col),
                                        "field": col,
                                        "sortable": True,
                                        "align": "left",
                                    }
                                    # Use custom sort field for VEP_Consequence based on priority
                                    if col == "VEP_Consequence":
                                        col_def["field"] = "_consequence_priority"
                                    cols.append(col_def)
                                return cols

                            with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
                                ui.label(f"Data ({len(rows)} rows)").classes(
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

                                                def select_all():
                                                    selected_cols_local["value"] = list(
                                                        all_columns_local
                                                    )
                                                    update_table()

                                                def select_none():
                                                    selected_cols_local["value"] = []
                                                    update_table()

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

                            with ui.card().classes("w-full"):
                                data_table = (
                                    ui.table(
                                        columns=make_columns(
                                            selected_cols_local["value"]
                                        ),
                                        rows=rows,
                                        pagination={
                                            "rowsPerPage": 10,
                                            "sortBy": "VEP_Consequence",
                                            "descending": False,
                                        },
                                    )
                                    .classes("w-full")
                                    .props("dense flat")
                                )

                                data_table.add_slot("body", WOMBAT_TABLE_SLOT)

                                def on_view_variant(e):
                                    row_data = e.args
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

                                data_table.on("view_variant", on_view_variant)

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
                                update_table()

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
                                update_table()

                            def update_table():
                                visible = [
                                    c
                                    for c in all_columns_local
                                    if c in selected_cols_local["value"]
                                ]
                                data_table.columns = make_columns(visible)
                                data_table.update()

                                for col, checkbox in checkboxes.items():
                                    checkbox.value = col in selected_cols_local["value"]

                        data_table_refreshers.append(render_data_table.refresh)
                        render_data_table()

                except Exception as e:
                    ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
