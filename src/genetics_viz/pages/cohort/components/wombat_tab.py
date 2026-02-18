"""Wombat tab component for family page."""

import re
from typing import Any, Callable, Dict, List

import polars as pl
from nicegui import ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    get_column_group,
    get_column_sorting,
    get_display_label,
    get_dropped_columns,
    get_schema_overrides,
    reorder_columns_by_group,
)
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.genesets import load_genesets
from genetics_viz.utils.score_colors import get_score_color
from genetics_viz.utils.clinvar import (
    CLINVAR_COLORS,
    format_clinvar_display,
    get_clinvar_color,
)
from genetics_viz.utils.vep import (
    VEP_CONSEQUENCES,
    VEP_CONSEQUENCE_PRIORITY,
    format_consequence_display,
    get_consequence_color,
    get_highest_consequence_term,
    get_highest_priority_consequence,
)
from genetics_viz.utils.view_presets import VIEW_PRESETS, select_preset_for_config
from genetics_viz.utils.cytobands import (
    CHROM_ORDER,
    CHROM_SIZES_MB,
    CYTOBANDS,
    GIESTAIN_COLORS,
    VALIDATION_COLORS,
    norm_chrom,
)



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
    available_genesets = load_genesets(store.data_dir)

    # Create dictionaries to store data for each wombat config
    wombat_data: Dict[str, Dict[str, Any]] = {}

    # Create subtabs for each wombat config
    with ui.tabs().classes("w-full") as wombat_subtabs:
        subtab_refs = {}
        for wf in wombat_files:
            subtab_refs[wf["wombat_config"]] = ui.tab(wf["wombat_config"])

    with ui.tab_panels(wombat_subtabs, value=list(subtab_refs.values())[0]).props(
        "keep-alive"
    ).classes("w-full"):
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
                    wombat_data[config_name]["_table_state"] = {"sorting": [], "page": 0}
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
                                        "sorting": get_column_sorting(col),
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
                                    apply_width_constraints(col_def, col)
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

                                # Column selector dialog
                                col_dialog, _sync_col_selector = build_column_selector(
                                    all_columns=all_columns_local,
                                    selected_cols=selected_cols_local,
                                    on_visibility_change=_apply_col_visibility,
                                    presets=VIEW_PRESETS,
                                )
                                ui.button(
                                    "Columns", icon="view_column",
                                    on_click=col_dialog.open,
                                ).props("outline color=blue size=sm")

                                # Stats button
                                def show_stats_dialog(current_rows=rows):
                                    from collections import Counter

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

                                    # Classify variant type
                                    for r in unique_variants:
                                        ref = str(r.get("REF", ""))
                                        alt = str(r.get("ALT", ""))
                                        r["_is_snv"] = len(ref) == 1 and len(alt) == 1

                                    # Use shared constants from cytobands module
                                    chrom_order = CHROM_ORDER
                                    chrom_sizes_mb = CHROM_SIZES_MB
                                    validation_colors = VALIDATION_COLORS

                                    # Filter state (persists across refreshes)
                                    type_filter = {"snv": True, "indel": True}
                                    show_ideogram = {"value": False}
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

                            def on_view_variant(e):
                                row_data = e.get("row", {})
                                chrom = row_data.get("#CHROM", "")
                                pos = row_data.get("POS", "")
                                ref = row_data.get("REF", "")
                                alt = row_data.get("ALT", "")
                                sample_val = row_data.get("sample", "")

                                # Callback to update the Validation column in the table
                                # NOTE: use data["all_rows"] and data["_refresh"]["fn"]
                                # instead of bare `all_rows` / `render_data_table.refresh`
                                # to avoid late-binding closure over the loop variable.
                                def on_save(validation_status: str):
                                    # Reload validation data from file
                                    validation_file = (
                                        store.data_dir / "validations" / "snvs.tsv"
                                    )
                                    validation_map = load_validation_map(
                                        validation_file, family_id
                                    )
                                    # Re-add validation status to all rows
                                    for row in data["all_rows"]:
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
                                    refresh_fn = data["_refresh"]["fn"]
                                    with page_client:
                                        ui.timer(
                                            0.1,
                                            refresh_fn,
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

                            # Restore table state (sorting / page) across refreshes
                            table_state = data["_table_state"]
                            saved_sorting = table_state.get("sorting", [])

                            # Pre-sort rows using saved sorting so JS gets them in order
                            if saved_sorting:
                                col_id = saved_sorting[0]["id"]
                                desc = saved_sorting[0].get("desc", False)
                                col_def = next(
                                    (c for c in make_columns() if c.get("id") == col_id), {}
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

                            data["_dt"] = DataTable(
                                columns=make_columns(),
                                rows=rows,
                                row_key="Variant",
                                pagination={"rowsPerPage": 10},
                                visible_columns=["actions"] + list(selected_cols_local["value"]),
                                on_row_action=on_view_variant,
                                initial_sorting=saved_sorting,
                                initial_page=table_state.get("page", 0),
                                state_holder=table_state,
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
                                _sync_col_selector()

                            # Connect preset change handler
                            preset_select.on_value_change(on_preset_change)

                        # Store refresh reference for filter callbacks
                        wombat_data[config_name]["_refresh"]["fn"] = render_data_table.refresh
                        data_table_refreshers.append(render_data_table.refresh)
                        render_data_table()

                except Exception as e:
                    ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")
