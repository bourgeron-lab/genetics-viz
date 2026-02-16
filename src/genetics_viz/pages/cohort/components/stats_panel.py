"""Statistics panel for the cohort page — box/bar plots of variant counts."""

import re
import statistics
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import polars as pl
from nicegui import run, ui

from genetics_viz.models import Cohort, DataStore



def _discover_wombat_configs(
    store: DataStore, cohort: Cohort
) -> List[str]:
    """Scan the cohort wombat directory for available config names."""
    wombat_dir = store.data_dir / "cohorts" / cohort.name / "wombat"
    if not wombat_dir.exists():
        return []

    pattern = re.compile(
        rf"{re.escape(cohort.name)}\.rare\.([^.]+)\.(.+?)\.results\.tsv$"
    )
    configs: Set[str] = set()
    for tsv in wombat_dir.glob("*.tsv"):
        m = pattern.match(tsv.name)
        if m:
            configs.add(m.group(2))
    return sorted(configs)


def _load_geneset_genes(store: DataStore) -> Dict[str, Set[str]]:
    """Load geneset gene lists from params/genesets/*.tsv."""
    genesets_dir = store.data_dir / "params" / "genesets"
    result: Dict[str, Set[str]] = {}
    if not genesets_dir.exists():
        return result
    for geneset_file in genesets_dir.glob("*.tsv"):
        genes: Set[str] = set()
        with open(geneset_file) as f:
            next(f, None)  # skip header
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene.upper())
        if genes:
            result[geneset_file.stem] = genes
    return result


def _find_wombat_file(
    store: DataStore, cohort: Cohort, wombat_config: str
) -> Path | None:
    """Find the cohort-level wombat TSV for a given config."""
    wombat_dir = store.data_dir / "cohorts" / cohort.name / "wombat"
    if not wombat_dir.exists():
        return None
    pattern = re.compile(
        rf"{re.escape(cohort.name)}\.rare\.[^.]+\."
        rf"{re.escape(wombat_config)}\.results\.tsv$"
    )
    for f in wombat_dir.glob("*.tsv"):
        if pattern.match(f.name):
            return f
    return None


def _count_variants_per_sample(
    tsv_path: Path,
    geneset_genes: Set[str] | None,
) -> Dict[str, int]:
    """Count unique variants per sample in a cohort wombat file.

    If geneset_genes is provided, only count variants where VEP_SYMBOL
    overlaps the gene set. Returns {sample_id: variant_count}.
    """
    df = pl.read_csv(
        tsv_path,
        separator="\t",
        infer_schema_length=0,
        null_values=[".", ""],
    )

    # Ensure required columns exist
    required = {"#CHROM", "POS", "REF", "ALT", "sample"}
    if not required.issubset(set(df.columns)):
        return {}

    # Filter by geneset genes if provided
    if geneset_genes and "VEP_SYMBOL" in df.columns:
        df = df.filter(
            pl.col("VEP_SYMBOL")
            .cast(pl.Utf8)
            .fill_null("")
            .str.split(",")
            .list.eval(
                pl.element().str.to_uppercase().is_in(list(geneset_genes))
            )
            .list.any()
        )

    # Count unique variants per sample
    variant_counts = (
        df.select(["#CHROM", "POS", "REF", "ALT", "sample"])
        .unique()
        .group_by("sample")
        .agg(pl.len().alias("count"))
    )

    return {
        str(row["sample"]): row["count"]
        for row in variant_counts.iter_rows(named=True)
    }



def render_stats_panel(
    store: DataStore, cohort: Cohort, filtered_state: Dict[str, Any]
) -> None:
    """Render the statistics panel with wombat/geneset selectors and plots."""
    import plotly.graph_objects as go

    # Discover available wombat configs
    wombat_configs = _discover_wombat_configs(store, cohort)
    if not wombat_configs:
        ui.label("No wombat files found").classes("text-gray-500 italic")
        return

    # Load geneset gene lists and metadata
    geneset_genes = _load_geneset_genes(store)
    geneset_names = sorted(geneset_genes.keys())

    if not geneset_names:
        ui.label("No genesets found").classes("text-gray-500 italic")
        return

    # Shared state for plot figures
    plot_state: Dict[str, Any] = {"figs": []}

    @ui.refreshable
    def plot_area() -> None:
        figs = plot_state["figs"]
        if not figs:
            return
        for fig in figs:
            ui.plotly(fig).classes("w-full")

    with ui.card().classes("w-full"):
        ui.label("Statistics").classes("text-lg font-semibold text-blue-700 mb-2")

        # Wombat config selector
        config_select = ui.select(
            options=wombat_configs,
            value=wombat_configs[0],
            label="Wombat config",
        ).classes("w-full mb-2")

        # Geneset multiselect
        geneset_select = ui.select(
            options=geneset_names,
            value=[],
            label="Genesets",
            multiple=True,
        ).props("use-chips").classes("w-full mb-2")

        # Spinner container
        spinner_container = ui.column().classes("w-full items-center")

        # Plot button
        async def on_plot() -> None:
            config = config_select.value
            if not config:
                ui.notify("Select a wombat config", type="warning")
                return
            genesets = geneset_select.value or []
            if not genesets:
                ui.notify("Select at least one geneset", type="warning")
                return

            # Find the wombat file
            tsv_path = _find_wombat_file(store, cohort, config)
            if tsv_path is None:
                ui.notify(
                    f"No wombat file found for config: {config}", type="negative"
                )
                return

            # Show spinner, clear previous plots
            plot_state["figs"] = []
            plot_area.refresh()
            spinner_container.clear()
            with spinner_container:
                ui.spinner("dots", size="lg")
                ui.label("Loading variants...").classes("text-gray-500 text-sm")

            # Build sample -> phenotype map from currently filtered individuals
            individuals = filtered_state["individuals"]
            sample_pheno: Dict[str, str] = {
                ind["Sample ID"]: ind["Phenotype"]
                for ind in individuals
                if ind["Phenotype"] != "-"
            }
            phenotype_values = sorted(set(sample_pheno.values()))

            # Phenotype colors
            default_colors = [
                "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            ]
            pheno_colors = {
                p: default_colors[i % len(default_colors)]
                for i, p in enumerate(phenotype_values)
            }

            # Collect per-geneset data, split by median
            # box_data: genesets with median > 1 (show variant count distribution)
            # bar_data: genesets with median <= 1 (show % carriers)
            box_data: List[Tuple[str, Dict[str, List[int]]]] = []
            bar_data: List[Tuple[str, Dict[str, float]]] = []

            for gs_name in genesets:
                genes = geneset_genes.get(gs_name)
                counts = await run.io_bound(
                    _count_variants_per_sample, tsv_path, genes
                )

                # Compute per-phenotype lists
                pheno_counts: Dict[str, List[int]] = {}
                for pheno in phenotype_values:
                    pheno_counts[pheno] = [
                        counts.get(sid, 0)
                        for sid, p in sample_pheno.items()
                        if p == pheno
                    ]

                # Overall median across all samples
                all_counts = [
                    c for vals in pheno_counts.values() for c in vals
                ]
                median_val = statistics.median(all_counts) if all_counts else 0

                if median_val > 1:
                    box_data.append((gs_name, pheno_counts))
                else:
                    pheno_pcts: Dict[str, float] = {}
                    for pheno, vals in pheno_counts.items():
                        if vals:
                            carriers = sum(1 for c in vals if c > 0)
                            pheno_pcts[pheno] = carriers / len(vals) * 100
                        else:
                            pheno_pcts[pheno] = 0.0
                    bar_data.append((gs_name, pheno_pcts))

            figs: List[go.Figure] = []

            # Box plot for high-count genesets
            if box_data:
                fig_box = go.Figure()
                for pheno in phenotype_values:
                    x_vals: List[int] = []
                    y_vals: List[str] = []
                    for gs_name, pheno_counts in box_data:
                        for c in pheno_counts.get(pheno, []):
                            x_vals.append(c)
                            y_vals.append(gs_name)
                    if x_vals:
                        fig_box.add_trace(
                            go.Box(
                                x=x_vals,
                                y=y_vals,
                                name=pheno,
                                orientation="h",
                                marker_color=pheno_colors.get(pheno, "#999"),
                            )
                        )
                fig_box.update_layout(
                    title="Variant count per sample",
                    xaxis_title="Variant count",
                    boxmode="group",
                    showlegend=True,
                    margin=dict(l=120, r=20, t=40, b=40),
                    height=max(300, len(box_data) * 120),
                )
                figs.append(fig_box)

            # Bar plot for low-count genesets (% carriers)
            if bar_data:
                fig_bar = go.Figure()
                for pheno in phenotype_values:
                    fig_bar.add_trace(
                        go.Bar(
                            x=[pcts.get(pheno, 0) for _, pcts in bar_data],
                            y=[gs_name for gs_name, _ in bar_data],
                            name=pheno,
                            orientation="h",
                            marker_color=pheno_colors.get(pheno, "#999"),
                        )
                    )
                fig_bar.update_layout(
                    title="% carriers per geneset",
                    xaxis_title="% carriers",
                    barmode="group",
                    showlegend=True,
                    margin=dict(l=120, r=20, t=40, b=40),
                    height=max(300, len(bar_data) * 80),
                )
                figs.append(fig_bar)

            # Hide spinner, show plots
            spinner_container.clear()
            plot_state["figs"] = figs
            plot_area.refresh()

        ui.button("Plot", icon="show_chart", on_click=on_plot).props(
            "color=blue"
        ).classes("mb-4")

        plot_area()
