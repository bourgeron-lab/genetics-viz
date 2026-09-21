"""Cohort-wide Polygenic Scores tab for the cohort page.

Violin plots of the cohort's z-scores for one selected trait, grouped by
phenotype, sex or predicted ancestry region, optionally split by sex within
each group, and optionally restricted to a subset of regions. Everything
follows the filters applied to the individuals table.

Violins are Plotly rather than ECharts, which has no violin series and would
mean hand-rolling a KDE. The sibling Statistics tab is already Plotly, so the
two sit consistently side by side.
"""

import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.utils.ancestry import (
    NULL_VALUES,
    AncestryFiles,
    Variant,
    find_cohort_ancestry_files,
    get_cohort_ancestry_dir,
    get_region_color,
)
from genetics_viz.utils.pedigree_labels import (
    PHENO_ORDER,
    SEX_ORDER,
    normalize_pheno,
    normalize_sex,
    pheno_label,
    sex_label,
)
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_pgs_data", "render_cohort_pgs_tab"]

#: Columns of the z-score file that are not traits.
_NON_TRAIT_COLUMNS = {"IID", "family_id"}

_GROUP_PHENOTYPE = "Phenotype"
_GROUP_SEX = "Sex"
_GROUP_REGION = "Region"
_GROUP_OPTIONS = [_GROUP_PHENOTYPE, _GROUP_SEX, _GROUP_REGION]

#: Matches the phenotype palette of the sibling Statistics panel.
_PHENO_COLORS = {"2": "#d62728", "1": "#2ca02c", "-9": "#7f7f7f"}
_SEX_COLORS = {"1": "#1f77b4", "2": "#e377c2", "-9": "#7f7f7f"}

#: Sex colours for the split halves, used whatever the grouping axis is.
_SPLIT_SIDES = (("1", "negative"), ("2", "positive"))


def probe_cohort_pgs_data(cohort: Any) -> bool:
    """Check whether cohort-level PGS z-scores exist for this cohort."""
    files = find_cohort_ancestry_files(cohort)
    return files is not None and files.pgs_zscore_path is not None


def _to_float(value: Any) -> Optional[float]:
    """Coerce a cell to a float, returning None when it is not numeric.

    A trait that is entirely ``nan`` in the source comes back from polars as an
    all-null string column, so the values reaching here are not always numbers.
    """
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if result != result else result  # drop NaN


def _load(
    files: AncestryFiles,
) -> Tuple[Optional[Dict[str, Dict]], List[str], List[str]]:
    """Read the z-scores with the region predictions joined on.

    Returns ``(by_sample, traits, regions)`` where ``by_sample`` maps a sample
    ID to its z-score row plus a ``region`` key.
    """
    df = read_tsv_or_none(
        files.pgs_zscore_path, infer_schema_length=10000, null_values=NULL_VALUES
    )
    if df is None or len(df) == 0 or "IID" not in df.columns:
        return None, [], []

    by_sample: Dict[str, Dict] = {
        row["IID"]: row for row in df.to_dicts() if row.get("IID")
    }

    regions: List[str] = []
    if files.ancestry_path is not None:
        try:
            anc = read_tsv_or_none(files.ancestry_path, null_values=NULL_VALUES)
            if anc is not None and "region" in anc.columns:
                for row in anc.to_dicts():
                    entry = by_sample.get(row.get("IID"))
                    if entry is not None:
                        entry["region"] = row.get("region")
                regions = sorted({r for r in anc["region"].to_list() if r})
        except Exception:
            logger.warning("Failed to read %s", files.ancestry_path, exc_info=True)

    traits = sorted(c for c in df.columns if c not in _NON_TRAIT_COLUMNS)
    return by_sample, traits, regions


def _group_of(kind: str, sample: Dict, individual: Dict) -> str:
    """Return the grouping key of one sample for the chosen axis."""
    if kind == _GROUP_SEX:
        return normalize_sex(individual.get("Sex"))
    if kind == _GROUP_REGION:
        return sample.get("region") or "—"
    return normalize_pheno(individual.get("Phenotype"))


def _group_order(kind: str, present: set) -> List[str]:
    """Order the groups for display, keeping only those actually present."""
    if kind == _GROUP_SEX:
        return [g for g in SEX_ORDER if g in present]
    if kind == _GROUP_REGION:
        return sorted(present)
    return [g for g in PHENO_ORDER if g in present]


def _group_label(kind: str, group: str, n: int) -> str:
    """Axis label for a group, carrying its sample count."""
    if kind == _GROUP_SEX:
        base = sex_label(group)
    elif kind == _GROUP_REGION:
        base = group
    else:
        base = pheno_label(group)
    return f"{base}<br>n={n}"


def _group_color(kind: str, group: str) -> str:
    """Colour for a group on the chosen axis."""
    if kind == _GROUP_SEX:
        return _SEX_COLORS.get(group, "#7f7f7f")
    if kind == _GROUP_REGION:
        return get_region_color(group)
    return _PHENO_COLORS.get(group, "#7f7f7f")


def render_cohort_pgs_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide Polygenic Scores tab.

    Args:
        store: DataStore instance
        cohort: Cohort instance
        filtered_state: Dict with 'individuals' key holding the filtered rows
        refreshers: List to append the refresh callback to
    """
    import plotly.graph_objects as go

    ancestry_dir = get_cohort_ancestry_dir(cohort)
    files = find_cohort_ancestry_files(cohort)
    if files is None or files.pgs_zscore_path is None:
        ui.label(f"No cohort-level PGS z-scores in: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    state: Dict[str, Any] = {
        "files": files,
        "by_sample": None,
        "traits": [],
        "regions": [],
    }

    def load_variant(variant: Optional[Variant] = None) -> None:
        selected = find_cohort_ancestry_files(cohort, variant) or files
        state["files"] = selected
        try:
            by_sample, traits, regions = _load(selected)
        except Exception:
            logger.warning(
                "Failed to read cohort PGS for %s", cohort.name, exc_info=True
            )
            by_sample, traits, regions = None, [], []
        state["by_sample"] = by_sample
        state["traits"] = traits
        state["regions"] = regions

    load_variant()
    if state["by_sample"] is None or not state["traits"]:
        ui.label("No PGS z-scores in file").classes("text-gray-500 italic")
        return

    controls: Dict[str, Any] = {
        "trait": state["traits"][0],
        "group_by": _GROUP_PHENOTYPE,
        "split_sex": False,
        "regions": list(state["regions"]),
    }

    with ui.card().classes("w-full"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Polygenic Scores").classes("text-lg font-semibold text-blue-700")
            provenance = ui.label(state["files"].label).classes("text-xs text-gray-500")

        variants = state["files"].all_variants
        if len(variants) > 1:

            def on_variant(event: Any) -> None:
                bundle, tag = event.value.split(" / ", 1)
                load_variant((bundle, tag))
                provenance.text = state["files"].label
                controls["regions"] = list(state["regions"])
                render_plot.refresh()

            ui.select(
                options=[f"{b} / {t}" for b, t in variants],
                value=f"{state['files'].bundle_version} / {state['files'].filter_tag}",
                label="Bundle / filters",
                on_change=on_variant,
            ).props("outlined dense").classes("w-full mb-2")

        ui.select(
            options=state["traits"],
            value=controls["trait"],
            label="Trait",
            with_input=True,
            on_change=lambda e: (
                controls.update({"trait": e.value}),
                render_plot.refresh(),
            ),
        ).props("outlined dense").classes("w-full mb-2")

        with ui.row().classes("items-center gap-4 w-full"):
            ui.select(
                options=_GROUP_OPTIONS,
                value=controls["group_by"],
                label="Group by",
                on_change=lambda e: (
                    controls.update({"group_by": e.value}),
                    render_plot.refresh(),
                ),
            ).props("outlined dense").classes("flex-1")

            split_box = ui.checkbox(
                "Split by sex",
                value=controls["split_sex"],
                on_change=lambda e: (
                    controls.update({"split_sex": e.value}),
                    render_plot.refresh(),
                ),
            )

        if state["regions"]:
            ui.select(
                options=state["regions"],
                value=controls["regions"],
                label="Regions",
                multiple=True,
                on_change=lambda e: (
                    controls.update({"regions": e.value or []}),
                    render_plot.refresh(),
                ),
            ).props("use-chips outlined dense").classes("w-full mt-2")

        @ui.refreshable
        def render_plot() -> None:
            # Splitting by sex while grouping by sex would just halve each
            # group into itself.
            group_by = controls["group_by"]
            split = controls["split_sex"] and group_by != _GROUP_SEX
            split_box.set_enabled(group_by != _GROUP_SEX)

            by_sample = state["by_sample"]
            trait = controls["trait"]
            allowed_regions = set(controls["regions"])

            # Pair each filtered individual with its z-score row.
            paired: List[Tuple[Dict, Dict]] = []
            for ind in filtered_state.get("individuals", []):
                sample = by_sample.get(ind["Sample ID"])
                if sample is None:
                    continue
                if state["regions"] and (sample.get("region") or "—") not in (
                    allowed_regions
                ):
                    continue
                if _to_float(sample.get(trait)) is None:
                    continue
                paired.append((sample, ind))

            if not paired:
                ui.label(
                    "No individuals with a score for this trait match the"
                    " current filters"
                ).classes("text-gray-500 italic")
                return

            # group -> (sex -> values)
            grouped: Dict[str, Dict[str, List[float]]] = {}
            for sample, ind in paired:
                group = _group_of(group_by, sample, ind)
                sex = normalize_sex(ind.get("Sex"))
                value = _to_float(sample.get(trait))
                grouped.setdefault(group, {}).setdefault(sex, []).append(value)

            order = _group_order(group_by, set(grouped))
            counts = {g: sum(len(v) for v in grouped[g].values()) for g in order}
            labels = {g: _group_label(group_by, g, counts[g]) for g in order}

            ui.label(f"{len(paired)} individuals with a {trait} score").classes(
                "text-xs text-gray-500"
            )

            fig = go.Figure()
            if split:
                # One trace per sex, the two halves of a shared violin.
                for sex, side in _SPLIT_SIDES:
                    xs: List[str] = []
                    ys: List[float] = []
                    for group in order:
                        for value in grouped[group].get(sex, []):
                            xs.append(labels[group])
                            ys.append(value)
                    if not ys:
                        continue
                    fig.add_trace(
                        go.Violin(
                            x=xs,
                            y=ys,
                            name=sex_label(sex),
                            side=side,
                            line_color=_SEX_COLORS.get(sex, "#7f7f7f"),
                            points="all",
                            jitter=0.3,
                            pointpos=0,
                            box_visible=True,
                            meanline_visible=True,
                            scalemode="count",
                        )
                    )
                fig.update_layout(violinmode="overlay", violingap=0)
            else:
                for group in order:
                    values = [v for vals in grouped[group].values() for v in vals]
                    fig.add_trace(
                        go.Violin(
                            x=[labels[group]] * len(values),
                            y=values,
                            name=labels[group].replace("<br>", " "),
                            line_color=_group_color(group_by, group),
                            points="all",
                            jitter=0.3,
                            pointpos=0,
                            box_visible=True,
                            meanline_visible=True,
                            scalemode="count",
                        )
                    )
                fig.update_layout(violinmode="group")

            fig.update_layout(
                title=f"{trait} by {group_by.lower()}",
                yaxis_title="z-score",
                showlegend=split,
                margin=dict(l=50, r=20, t=50, b=40),
                height=420,
            )
            ui.plotly(fig).classes("w-full")

        render_plot()
        refreshers.append(render_plot.refresh)
