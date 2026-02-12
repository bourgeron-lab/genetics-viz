"""SVs tab component for family page."""

from pathlib import Path
from typing import Any, Callable, Dict, List

import polars as pl
import yaml
from nicegui import ui

from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils.gene_scoring import get_gene_scorer


# Load WisecondorX thresholds and colors from YAML
def _load_wisecondorx_config():
    config_path = (
        Path(__file__).parent.parent.parent.parent
        / "config"
        / "wisecondorx_thresholds.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


WISECONDORX_CONFIG = _load_wisecondorx_config()


# Generate table slot template with config values
def _generate_svs_table_slot():
    """Generate SVS table slot with dynamic threshold values from config."""
    robust_loss = WISECONDORX_CONFIG["robust_loss"]
    permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
    robust_gain = WISECONDORX_CONFIG["robust_gain"]
    permissive_gain = WISECONDORX_CONFIG["permissive_gain"]

    return f"""
    <q-tr :props="props">
        <q-td key="actions" :props="props">
            <q-btn 
                flat 
                dense 
                size="sm" 
                icon="visibility" 
                color="blue"
                @click="$parent.$emit('view_sv', props.row)"
            >
                <q-tooltip>View in IGV</q-tooltip>
            </q-btn>
        </q-td>
        <q-td v-for="col in props.cols.filter(c => c.name !== 'actions')" :key="col.name" :props="props">
            <template v-if="col.name === 'Validation'">
                <span v-if="col.value === 'present'" style="display: flex; align-items: center; gap: 4px;">
                    <q-icon name="check_circle" color="green" size="sm">
                        <q-tooltip>Validated as present</q-tooltip>
                    </q-icon>
                    <span v-if="props.row.ValidationInheritance === 'de novo'" style="font-weight: bold;">dnm</span>
                    <span v-else-if="props.row.ValidationInheritance === 'homozygous'" style="font-weight: bold;">hom</span>
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
            <template v-else-if="col.name === 'call'">
                <q-badge 
                    v-if="col.value === '{robust_loss["label"]}'"
                    :label="col.value"
                    style="background-color: {robust_loss["color"]}; color: white;"
                />
                <q-badge 
                    v-else-if="col.value === '{permissive_loss["label"]}'"
                    :label="col.value"
                    style="background-color: {permissive_loss["color"]}; color: white;"
                />
                <q-badge 
                    v-else-if="col.value === '{robust_gain["label"]}'"
                    :label="col.value"
                    style="background-color: {robust_gain["color"]}; color: white;"
                />
                <q-badge 
                    v-else-if="col.value === '{permissive_gain["label"]}'"
                    :label="col.value"
                    style="background-color: {permissive_gain["color"]}; color: white;"
                />
                <span v-else class="text-grey-6">{{{{ col.value }}}}</span>
            </template>
            <template v-else-if="col.name === 'ratio'">
                <span 
                    :style="'color: ' + (parseFloat(col.value) <= {robust_loss["ratio_threshold"]} ? '{robust_loss["color"]}' : parseFloat(col.value) <= {permissive_loss["ratio_threshold"]} ? '{permissive_loss["color"]}' : parseFloat(col.value) >= {robust_gain["ratio_threshold"]} ? '{robust_gain["color"]}' : parseFloat(col.value) >= {permissive_gain["ratio_threshold"]} ? '{permissive_gain["color"]}' : 'inherit') + '; font-weight: ' + ((parseFloat(col.value) <= {robust_loss["ratio_threshold"]} || parseFloat(col.value) >= {robust_gain["ratio_threshold"]}) ? 'bold' : (parseFloat(col.value) <= {permissive_loss["ratio_threshold"]} || parseFloat(col.value) >= {permissive_gain["ratio_threshold"]}) ? '600' : 'normal')"
                >
                    {{{{ col.value }}}}
                </span>
            </template>
            <template v-else-if="col.name === 'zscore'">
                <span 
                    :style="'color: ' + (parseFloat(col.value) <= {robust_loss["zscore_threshold"]} ? '{robust_loss["color"]}' : parseFloat(col.value) <= {permissive_loss["zscore_threshold"]} ? '{permissive_loss["color"]}' : parseFloat(col.value) >= {robust_gain["zscore_threshold"]} ? '{robust_gain["color"]}' : parseFloat(col.value) >= {permissive_gain["zscore_threshold"]} ? '{permissive_gain["color"]}' : 'inherit') + '; font-weight: ' + ((parseFloat(col.value) <= {robust_loss["zscore_threshold"]} || parseFloat(col.value) >= {robust_gain["zscore_threshold"]}) ? 'bold' : (parseFloat(col.value) <= {permissive_loss["zscore_threshold"]} || parseFloat(col.value) >= {permissive_gain["zscore_threshold"]}) ? '600' : 'normal')"
                >
                    {{{{ col.value }}}}
                </span>
            </template>
            <template v-else-if="col.name === 'chr:start-end'">
                <span style="display: flex; align-items: center; gap: 4px;">
                    {{{{ col.value }}}}
                    <q-icon v-if="props.row.IsCurated" name="check_circle" color="green" size="xs">
                        <q-tooltip>Curated breakpoints from validation</q-tooltip>
                    </q-icon>
                </span>
            </template>
            <template v-else-if="col.name === 'gene'">
                <template v-if="props.row.GeneBadges && props.row.GeneBadges.length > 0">
                    <q-badge 
                        v-for="(badge, idx) in props.row.GeneBadges" 
                        :key="idx"
                        :label="badge.label" 
                        :style="'background-color: ' + badge.color + '; color: ' + (badge.color === '#ffffff' ? 'black' : 'white') + '; ' + (badge.type.includes('exonic') ? 'border: 2px solid black;' : '')"
                        class="q-mr-xs q-mb-xs"
                    >
                        <q-tooltip>{{{{ badge.tooltip }}}}</q-tooltip>
                    </q-badge>
                </template>
                <template v-else>
                    <span>-</span>
                </template>
            </template>
            <template v-else-if="['genic_symbol', 'genic_ensg', 'exonic_symbol', 'exonic_ensg', 'VEP_Gene'].includes(col.name)">
                <template v-if="props.row[col.name + '_badges'] && props.row[col.name + '_badges'].length > 0">
                    <q-badge 
                        v-for="(badge, idx) in props.row[col.name + '_badges']" 
                        :key="idx"
                        :label="badge.label" 
                        :style="'background-color: ' + badge.color + '; color: ' + (badge.color === '#ffffff' ? 'black' : 'white') + '; ' + (badge.isExonic ? 'border: 2px solid black;' : '')"
                        class="q-mr-xs q-mb-xs"
                    >
                        <q-tooltip>{{{{ badge.tooltip }}}}</q-tooltip>
                    </q-badge>
                </template>
                <template v-else-if="col.value && col.value !== '-' && col.value !== ''">
                    <span>{{{{ col.value }}}}</span>
                </template>
                <template v-else>
                    <span>-</span>
                </template>
            </template>
            <template v-else>
                {{{{ col.value }}}}
            </template>
        </q-td>
    </q-tr>
"""


SVS_TABLE_SLOT = _generate_svs_table_slot()


def render_svs_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the SVs tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    svs_dir = store.data_dir / "families" / family_id / "svs"

    if not svs_dir.exists():
        ui.label(f"No SVs directory found at: {svs_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Create subtabs for different SV callers
    with ui.tabs().classes("w-full") as svs_subtabs:
        wisecondorx_tab = ui.tab("WisecondorX")

    with ui.tab_panels(svs_subtabs, value=wisecondorx_tab).classes("w-full"):
        # WisecondorX subtab
        with ui.tab_panel(wisecondorx_tab):
            render_wisecondorx_subtab(
                store=store,
                family_id=family_id,
                svs_dir=svs_dir,
                selected_members=selected_members,
                data_table_refreshers=data_table_refreshers,
                cohort_name=cohort_name,
            )


def render_wisecondorx_subtab(
    store: Any,
    family_id: str,
    svs_dir: Path,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    cohort_name: str,
) -> None:
    """Render the WisecondorX subtab content.

    Args:
        store: DataStore instance
        family_id: Family ID
        svs_dir: Path to SVs directory
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
        cohort_name: Cohort name
    """
    wisecondorx_dir = svs_dir / "wisecondorx"
    aberrations_file = wisecondorx_dir / f"{family_id}_aberrations.annotated.bed"

    if not wisecondorx_dir.exists():
        ui.label(f"No WisecondorX directory found at: {wisecondorx_dir}").classes(
            "text-gray-500 italic"
        )
        return

    if not aberrations_file.exists():
        ui.label(f"No aberrations file found at: {aberrations_file}").classes(
            "text-gray-500 italic"
        )
        return

    with ui.card().classes("w-full p-4"):
        ui.label("WisecondorX Aberrations").classes(
            "text-lg font-semibold text-blue-700 mb-2"
        )
        with ui.row().classes("gap-4"):
            ui.label("File Path:").classes("font-semibold")
            ui.label(str(aberrations_file)).classes("text-sm text-gray-600 font-mono")

        # Gene badge legend
        with ui.row().classes("gap-4 mt-3 items-center"):
            ui.label("Gene badges:").classes("text-sm font-semibold")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #ffffff; color: black; border: 2px solid black; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(exonic - black border)").classes("text-xs text-gray-600")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #ffffff; color: black; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(genic - no border)").classes("text-xs text-gray-600")
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    '<span class="q-badge" style="background-color: #8b0000; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">Gene</span>',
                    sanitize=False,
                )
                ui.label("(color indicates geneset importance)").classes(
                    "text-xs text-gray-600"
                )

        # CNV call legend
        with ui.row().classes("gap-4 mt-2 items-center"):
            ui.label("CNV calls:").classes("text-sm font-semibold")

            robust_loss = WISECONDORX_CONFIG["robust_loss"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {robust_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_loss["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≤{robust_loss['ratio_threshold']} & Z≤{robust_loss['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {permissive_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_loss["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≤{permissive_loss['ratio_threshold']} & Z≤{permissive_loss['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            robust_gain = WISECONDORX_CONFIG["robust_gain"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {robust_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_gain["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≥{robust_gain['ratio_threshold']} & Z≥{robust_gain['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

            permissive_gain = WISECONDORX_CONFIG["permissive_gain"]
            with ui.row().classes("gap-2 items-center"):
                ui.html(
                    f'<span class="q-badge" style="background-color: {permissive_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_gain["label"]}</span>',
                    sanitize=False,
                )
                ui.label(
                    f"(log2≥{permissive_gain['ratio_threshold']} & Z≥{permissive_gain['zscore_threshold']})"
                ).classes("text-xs text-gray-600")

    # Display BED file content in a table
    try:
        # Read BED file - the format uses spaces for numeric columns and tabs for gene columns
        # We need to handle this mixed delimiter format
        with open(aberrations_file, "r") as f:
            lines = f.readlines()

        if not lines:
            ui.label("File is empty").classes("text-gray-500 italic")
            return

        # Parse header - split by any whitespace
        header = lines[0].strip().split()

        # Parse data rows
        data = []
        for line in lines[1:]:
            if not line.strip():
                continue
            # Split by whitespace (spaces and tabs)
            parts = line.strip().split()
            # Expecting at least 7 columns (chr start end ratio zscore type barcode)
            # Plus 4 gene columns which may be empty
            if len(parts) >= 7:
                # Pad with empty strings if gene columns are missing
                while len(parts) < len(header):
                    parts.append("")
                data.append(parts[: len(header)])

        # Create DataFrame
        df = pl.DataFrame(
            {
                col: [row[i] if i < len(row) else "" for row in data]
                for i, col in enumerate(header)
            }
        )

        # Convert numeric columns
        numeric_cols = ["start", "end", "ratio", "zscore"]
        for col in numeric_cols:
            if col in df.columns:
                df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False))

        # Rename barcode column to sample if it exists
        if "barcode" in df.columns:
            df = df.rename({"barcode": "sample"})

        # Create chr:start-end column if chr, start, and end columns exist
        if "chr" in df.columns and "start" in df.columns and "end" in df.columns:
            # Ensure chr column starts with "chr" prefix
            df = df.with_columns(
                pl.when(pl.col("chr").cast(pl.Utf8).str.starts_with("chr"))
                .then(pl.col("chr").cast(pl.Utf8))
                .otherwise(pl.lit("chr") + pl.col("chr").cast(pl.Utf8))
                .alias("chr_prefixed")
            )

            # Create the combined column (cast to int first to avoid .0 decimals)
            df = df.with_columns(
                (
                    pl.col("chr_prefixed")
                    + ":"
                    + pl.col("start").cast(pl.Int64).cast(pl.Utf8)
                    + "-"
                    + pl.col("end").cast(pl.Int64).cast(pl.Utf8)
                ).alias("chr:start-end")
            )

            # Drop the temporary and original columns
            df = df.drop(["chr_prefixed", "chr", "start", "end"])

        # Create gene column combining genic_symbol and exonic_symbol
        if "genic_symbol" in df.columns and "exonic_symbol" in df.columns:
            # For each row, combine genes and mark which are exonic
            # Format: "gene1:exonic,gene2:genic,gene3:exonic"
            def create_gene_list(genic, exonic):
                genic_genes = (
                    set(str(genic).split(",")) if genic and str(genic) != "" else set()
                )
                exonic_genes = (
                    set(str(exonic).split(","))
                    if exonic and str(exonic) != ""
                    else set()
                )

                # Remove empty strings
                genic_genes.discard("")
                exonic_genes.discard("")

                # Create combined list with tags
                result = []
                for gene in exonic_genes:
                    result.append(f"{gene.strip()}:exonic")
                for gene in genic_genes:
                    if gene.strip() not in exonic_genes:
                        result.append(f"{gene.strip()}:genic")

                return ",".join(result) if result else ""

            # Apply the function
            df = df.with_columns(
                pl.struct(["genic_symbol", "exonic_symbol"])
                .map_elements(
                    lambda row: create_gene_list(
                        row["genic_symbol"], row["exonic_symbol"]
                    ),
                    return_dtype=pl.Utf8,
                )
                .alias("gene")
            )

            # Reorder columns to put chr:start-end and gene first
            priority_cols = (
                ["chr:start-end", "gene"] if "chr:start-end" in df.columns else ["gene"]
            )
            other_cols = [col for col in df.columns if col not in priority_cols]
            df = df.select(priority_cols + other_cols)
        elif "chr:start-end" in df.columns:
            # Just reorder if no gene columns
            other_cols = [col for col in df.columns if col != "chr:start-end"]
            df = df.select(["chr:start-end"] + other_cols)

        # Add CNV call classification based on ratio and zscore
        if "ratio" in df.columns and "zscore" in df.columns:

            def classify_cnv(ratio, zscore):
                """Classify CNV based on ratio (log2) and zscore thresholds."""
                try:
                    r = float(ratio) if ratio and str(ratio) != "" else 0
                    z = float(zscore) if zscore and str(zscore) != "" else 0

                    robust_loss = WISECONDORX_CONFIG["robust_loss"]
                    permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
                    robust_gain = WISECONDORX_CONFIG["robust_gain"]
                    permissive_gain = WISECONDORX_CONFIG["permissive_gain"]

                    # Robust calls
                    if (
                        r <= robust_loss["ratio_threshold"]
                        and z <= robust_loss["zscore_threshold"]
                    ):
                        return robust_loss["label"]
                    elif (
                        r >= robust_gain["ratio_threshold"]
                        and z >= robust_gain["zscore_threshold"]
                    ):
                        return robust_gain["label"]
                    # Permissive calls
                    elif (
                        r <= permissive_loss["ratio_threshold"]
                        and z <= permissive_loss["zscore_threshold"]
                    ):
                        return permissive_loss["label"]
                    elif (
                        r >= permissive_gain["ratio_threshold"]
                        and z >= permissive_gain["zscore_threshold"]
                    ):
                        return permissive_gain["label"]
                    else:
                        return "Below threshold"
                except:
                    return "N/A"

            df = df.with_columns(
                pl.struct(["ratio", "zscore"])
                .map_elements(
                    lambda row: classify_cnv(row["ratio"], row["zscore"]),
                    return_dtype=pl.Utf8,
                )
                .alias("call")
            )

            # Reorder to put call after chr:start-end and gene
            if "chr:start-end" in df.columns and "gene" in df.columns:
                priority_cols = ["chr:start-end", "gene", "call"]
                other_cols = [col for col in df.columns if col not in priority_cols]
                df = df.select(priority_cols + other_cols)
            elif "chr:start-end" in df.columns:
                priority_cols = ["chr:start-end", "call"]
                other_cols = [col for col in df.columns if col not in priority_cols]
                df = df.select(priority_cols + other_cols)

        # Convert to list of dicts for NiceGUI table
        all_rows = df.to_dicts()

        # Store original chr:start-end for each row before any modifications
        for row in all_rows:
            if "chr:start-end" in row:
                row["_original_locus"] = row["chr:start-end"]

        # Function to reload and apply validation data
        def reload_validations():
            """Reload validation data and update all rows."""
            # Load validation data from svs.tsv
            validation_file = store.data_dir / "validations" / "svs.tsv"
            validation_map = load_validation_map(validation_file, family_id)

            # Add Validation status to each row
            for row in all_rows:
                # Reset to original locus and clear curated flags
                if "_original_locus" in row:
                    row["chr:start-end"] = row["_original_locus"]
                row.pop("IsCurated", None)
                row.pop("OriginalLocus", None)

                # For SVs, variant_key is the original chr:start-end format
                variant_key = row.get("_original_locus", row.get("chr:start-end", ""))

                sample_id = row.get("sample", "")

                # Construct variant key in the format stored in svs.tsv
                # Format: chr:start-end:type (e.g., chr1:1000-2000:del or chr1:1000-2000:dup)
                # First, determine type from call
                sv_call = row.get("call", "")
                if (
                    "GAIN" in str(sv_call).upper()
                    or "gain" in str(sv_call).lower()
                    or "Gain" in str(sv_call)
                ):
                    sv_type = "dup"
                elif (
                    "LOSS" in str(sv_call).upper()
                    or "loss" in str(sv_call).lower()
                    or "Loss" in str(sv_call)
                ):
                    sv_type = "del"
                else:
                    # Try to infer from ratio
                    ratio = row.get("ratio", 0)
                    try:
                        ratio_val = float(ratio) if ratio else 0
                        sv_type = "dup" if ratio_val > 0 else "del"
                    except (ValueError, TypeError):
                        sv_type = "del"  # Default to deletion

                # Construct the full variant key
                full_variant_key = f"{variant_key}:{sv_type}"

                add_validation_status_to_row(
                    row, validation_map, full_variant_key, sample_id
                )

                # Check if there are "present" validations with curated boundaries
                # If so, update the chr:start-end display to show curated values
                # Store original coordinates separately for dialog opening
                map_key = (full_variant_key, sample_id)
                if map_key in validation_map:
                    validations = validation_map[map_key]
                    # Find present validations with curated boundaries (not ignored)
                    present_with_curated = [
                        v
                        for v in validations
                        if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
                    ]
                    if present_with_curated:
                        # Sort by timestamp (most recent first)
                        present_with_curated.sort(key=lambda v: v[6], reverse=True)
                        most_recent = present_with_curated[0]
                        curated_start = most_recent[4]
                        curated_end = most_recent[5]

                        # Parse original chr:start-end
                        parts = variant_key.split(":")
                        if len(parts) == 2:
                            chrom = parts[0]
                            range_parts = parts[1].split("-")
                            if len(range_parts) == 2:
                                orig_start = range_parts[0]
                                orig_end = range_parts[1]

                                # Store original coordinates for dialog opening
                                row["OriginalLocus"] = variant_key

                                # Use curated values if provided, otherwise keep original
                                new_start = (
                                    curated_start if curated_start else orig_start
                                )
                                new_end = curated_end if curated_end else orig_end

                                # Update the display value and mark as curated
                                row["chr:start-end"] = f"{chrom}:{new_start}-{new_end}"
                                row["IsCurated"] = True

            # Add gene badge information for all rows
            gene_scorer = get_gene_scorer()
            for row in all_rows:
                # Process main gene column
                gene_str = row.get("gene", "")
                if gene_str and gene_str != "-":
                    # Parse gene string format: "SYMBOL:type,SYMBOL2:type"
                    gene_badges = []
                    for gene_part in str(gene_str).split(","):
                        if ":" in gene_part:
                            symbol = gene_part.split(":")[0].strip()
                            gene_type = gene_part.split(":")[1].strip()
                        else:
                            symbol = gene_part.strip()
                            gene_type = ""

                        if symbol:
                            score, _ = gene_scorer.get_gene_score_and_sets(symbol)
                            color = gene_scorer.get_gene_color(symbol)
                            tooltip = gene_scorer.get_gene_tooltip(symbol)
                            gene_badges.append(
                                {
                                    "label": symbol,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "type": gene_type,
                                    "score": score,
                                }
                            )

                    # Sort by score (descending)
                    gene_badges.sort(key=lambda x: x["score"], reverse=True)

                    # Limit to first 6 genes and add "+X genes" indicator if needed
                    total_genes = len(gene_badges)
                    if total_genes > 6:
                        gene_badges = gene_badges[:6]
                        # Add a "+X genes" badge
                        remaining_count = total_genes - 6
                        gene_badges.append(
                            {
                                "label": f"+{remaining_count} genes",
                                "color": "#9e9e9e",  # grey color
                                "tooltip": f"{remaining_count} more genes",
                                "type": "",
                            }
                        )

                    row["GeneBadges"] = gene_badges
                else:
                    row["GeneBadges"] = []

                # Process genic_symbol, exonic_symbol columns
                for col_name in ["genic_symbol", "exonic_symbol"]:
                    col_value = row.get(col_name, "")
                    if col_value and col_value != "-":
                        badges = []
                        symbols = [
                            s.strip() for s in str(col_value).split(",") if s.strip()
                        ]
                        for symbol in symbols:
                            color = gene_scorer.get_gene_color(symbol)
                            tooltip = gene_scorer.get_gene_tooltip(symbol)
                            badges.append(
                                {
                                    "label": symbol,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "isExonic": col_name == "exonic_symbol",
                                }
                            )
                        row[f"{col_name}_badges"] = badges
                    else:
                        row[f"{col_name}_badges"] = []

                # Process genic_ensg, exonic_ensg, VEP_Gene columns (ENSG IDs)
                for col_name in ["genic_ensg", "exonic_ensg", "VEP_Gene"]:
                    col_value = row.get(col_name, "")
                    if col_value and col_value != "-":
                        badges = []
                        ensgs = [
                            e.strip() for e in str(col_value).split(",") if e.strip()
                        ]
                        for ensg in ensgs:
                            color = gene_scorer.get_gene_color(ensg)
                            tooltip = gene_scorer.get_gene_tooltip(ensg)
                            badges.append(
                                {
                                    "label": ensg,
                                    "color": color,
                                    "tooltip": tooltip,
                                    "isExonic": col_name in ["exonic_ensg"],
                                }
                            )
                        row[f"{col_name}_badges"] = badges
                    else:
                        row[f"{col_name}_badges"] = []

        # Initial load of validations
        reload_validations()

        # Get all columns (add Validation column)
        all_columns = list(df.columns) + ["Validation"]

        # All columns visible by default except gene ID and symbol columns and type
        unchecked_columns = {
            "genic_ensg",
            "exonic_ensg",
            "genic_symbol",
            "exonic_symbol",
            "type",
        }
        selected_cols = {
            "value": [col for col in all_columns if col not in unchecked_columns]
        }

        # Define all possible call values
        all_call_values = [
            "Robust LOSS",
            "Robust GAIN",
            "Permissive LOSS",
            "Permissive Gain",
            "Below threshold",
        ]
        # Default: all selected except "Below threshold"
        selected_calls = {
            "value": [
                call for call in all_call_values if call != "Below threshold"
            ]
        }

        # Create a container for the data table
        data_container = ui.column().classes("w-full")

        # Capture the client context for use in callbacks
        from nicegui import context

        page_client = context.client

        with data_container:

            @ui.refreshable
            def render_data_table():
                # Filter rows by selected members if 'sample' column exists
                if "sample" in df.columns:
                    rows = [
                        r
                        for r in all_rows
                        if r.get("sample") in selected_members["value"]
                    ]
                else:
                    rows = all_rows

                # Filter rows by selected call values if 'call' column exists
                if "call" in df.columns:
                    rows = [
                        r
                        for r in rows
                        if r.get("call") in selected_calls["value"]
                    ]

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
                    cols.extend(
                        [
                            {
                                "name": col,
                                "label": col,
                                "field": col,
                                "sortable": True,
                                "align": "left",
                            }
                            for col in visible_cols
                        ]
                    )
                    return cols

                with ui.row().classes("items-center gap-4 mt-4 mb-2"):
                    ui.label(f"Data ({len(rows)} rows)").classes(
                        "text-lg font-semibold text-blue-700"
                    )

                    # Column selector
                    with ui.button("Select Columns", icon="view_column").props(
                        "outline color=blue"
                    ):
                        with ui.menu():
                            ui.label("Show/Hide Columns:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2"):
                                    checkboxes: Dict[str, Any] = {}

                                    def select_all():
                                        selected_cols["value"] = list(all_columns)
                                        update_table()

                                    def select_none():
                                        selected_cols["value"] = []
                                        update_table()

                                    ui.button("All", on_click=select_all).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button("None", on_click=select_none).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")

                                ui.separator()

                                for col in all_columns:
                                    checkboxes[col] = ui.checkbox(
                                        col,
                                        value=col in selected_cols["value"],
                                        on_change=lambda e, c=col: handle_col_change(
                                            c, e.value
                                        ),
                                    ).classes("text-sm")

                    # Call filter
                    with ui.button("Filter Call", icon="filter_list").props(
                        "outline color=blue"
                    ):
                        with ui.menu():
                            ui.label("Filter by Call:").classes(
                                "px-4 py-2 font-semibold text-sm"
                            )
                            ui.separator()

                            with ui.column().classes("p-2"):
                                with ui.row().classes("gap-2 mb-2"):
                                    call_checkboxes: Dict[str, Any] = {}

                                    def select_all_calls():
                                        selected_calls["value"] = list(all_call_values)
                                        update_call_filter()

                                    def select_none_calls():
                                        selected_calls["value"] = []
                                        update_call_filter()

                                    ui.button("All", on_click=select_all_calls).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")
                                    ui.button("None", on_click=select_none_calls).props(
                                        "size=sm flat dense"
                                    ).classes("text-xs")

                                ui.separator()

                                for call_value in all_call_values:
                                    call_checkboxes[call_value] = ui.checkbox(
                                        call_value,
                                        value=call_value in selected_calls["value"],
                                        on_change=lambda e, c=call_value: handle_call_change(
                                            c, e.value
                                        ),
                                    ).classes("text-sm")

                with ui.card().classes("w-full"):
                    data_table = (
                        ui.table(
                            columns=make_columns(selected_cols["value"]),
                            rows=rows,
                            pagination={"rowsPerPage": 10},
                        )
                        .classes("w-full")
                        .props("dense flat")
                    )

                    data_table.add_slot("body", SVS_TABLE_SLOT)

                    def on_view_sv(e):
                        row_data = e.args
                        locus = row_data.get("chr:start-end", "")
                        sample_id = row_data.get("sample", "")

                        if not locus or not sample_id:
                            ui.notify(
                                "Missing locus or sample information", type="warning"
                            )
                            return

                        # Parse locus (format: chr:start-end)
                        # Use original locus if available (for curated variants)
                        try:
                            locus_to_parse = row_data.get("OriginalLocus", locus)
                            parts = locus_to_parse.split(":")
                            if len(parts) == 2:
                                chrom = parts[0]
                                range_parts = parts[1].split("-")
                                if len(range_parts) == 2:
                                    start = range_parts[0]
                                    end = range_parts[1]

                                    # Define refresh callback that reloads validations
                                    def on_save():
                                        reload_validations()
                                        render_data_table.refresh()

                                    # Show SV dialog with refresh callback
                                    show_sv_dialog(
                                        cohort_name=cohort_name,
                                        family_id=family_id,
                                        chrom=chrom,
                                        start=start,
                                        end=end,
                                        sample=sample_id,
                                        sv_data=row_data,
                                        on_validation_saved=on_save,
                                    )
                                else:
                                    ui.notify(
                                        "Invalid locus format. Expected chr:start-end",
                                        type="warning",
                                    )
                            else:
                                ui.notify(
                                    "Invalid locus format. Expected chr:start-end",
                                    type="warning",
                                )
                        except Exception as ex:
                            ui.notify(f"Error parsing locus: {ex}", type="warning")

                    data_table.on("view_sv", on_view_sv)

                def handle_col_change(col_name, is_checked):
                    if is_checked and col_name not in selected_cols["value"]:
                        selected_cols["value"].append(col_name)
                    elif not is_checked and col_name in selected_cols["value"]:
                        selected_cols["value"].remove(col_name)

                    # Reorder to match all_columns order
                    selected_cols["value"] = [
                        col for col in all_columns
                        if col in selected_cols["value"]
                    ]

                    update_table()

                def update_table():
                    visible = [c for c in all_columns if c in selected_cols["value"]]
                    data_table.columns = make_columns(visible)
                    data_table.update()

                    for col, checkbox in checkboxes.items():
                        checkbox.value = col in selected_cols["value"]

                def handle_call_change(call_value, is_checked):
                    if is_checked and call_value not in selected_calls["value"]:
                        selected_calls["value"].append(call_value)
                    elif not is_checked and call_value in selected_calls["value"]:
                        selected_calls["value"].remove(call_value)
                    update_call_filter()

                def update_call_filter():
                    # Refresh the entire data table to apply the call filter
                    render_data_table.refresh()
                    # Update checkbox states
                    for call, checkbox in call_checkboxes.items():
                        checkbox.value = call in selected_calls["value"]

            data_table_refreshers.append(render_data_table.refresh)
            render_data_table()

    except Exception as e:
        ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
