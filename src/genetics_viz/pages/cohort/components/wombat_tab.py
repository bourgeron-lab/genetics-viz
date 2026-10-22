"""Wombat tab component for family page.

One sub-tab per Wombat config of the family, each built when first shown:
the rows and columns come from ``wombat_rows``, the filters from
``wombat_filters``, the Stats dialog from ``wombat_stats``.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import polars as pl
from nicegui import ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.filters import notify_unknown_terms
from genetics_viz.components.page_link import PageLink, decode_sorting, encode_sorting
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.models import Family
from genetics_viz.pages.cohort.components.wombat_filters import (
    LINK_KEYS,
    filters_from_link,
    filters_to_link,
    WombatFilters,
    apply_filters,
    render_filter_panel,
)
from genetics_viz.pages.cohort.components.wombat_rows import (
    add_curation_status,
    build_rows,
    find_wombat_files,
    load_config_frame,
    presort_rows,
    table_columns,
    variant_key,
    wombat_columns,
    wombat_file_pattern,
)
from genetics_viz.pages.cohort.components.wombat_stats import show_wombat_stats_dialog
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.export import export_filename, is_export_column, rows_to_tsv
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.genesets import geneset_sources
from genetics_viz.utils.sample_counts import cohort_wombat_file
from genetics_viz.utils.sample_qc import denovo_risk, flag_title, load_family_qc
from genetics_viz.utils.sample_qc_acks import ack_summary, load_acks
from genetics_viz.utils.sharding import get_family_path
from genetics_viz.utils.view_presets import VIEW_PRESETS, select_preset_for_config
from genetics_viz.utils.wombat import aggregate_by_variant

logger = logging.getLogger(__name__)


def probe_wombat_data(data_dir: Path, family_id: str) -> bool:
    """Check whether any wombat result files exist for this family."""
    wombat_dir = get_family_path(data_dir, family_id) / "wombat"
    if not wombat_dir.exists():
        return False
    pattern = wombat_file_pattern(family_id)
    return any(pattern.match(f.name) for f in wombat_dir.glob("*.tsv"))


#: Wombat configs whose calls are de novo mutations.
_DENOVO_CONFIG = re.compile(r"de_?novo", re.IGNORECASE)


def _render_denovo_qc_warning(
    store: Any, family_id: str, on_open_qc: Optional[Callable[[], None]]
) -> None:
    """A warning above a de novo table when sample QC doubts the trio.

    A parent who is not the parent, a swapped, duplicated or contaminated
    sample: each turns inherited variants into apparent de novo calls. Only the
    family's own QC files are read (a few KB, cached).

    A flag a curator acknowledged as explained (known MZ twins, say) is left
    out. One acknowledged as a confirmed issue stays: the calls are still
    unreliable, and the warning says it was confirmed.
    """
    try:
        qc = load_family_qc(store.data_dir, family_id)
    except Exception:
        logger.warning("Failed to read sample QC of %s", family_id, exc_info=True)
        return
    if qc is None:
        return
    acks = load_acks(store.data_dir)
    risky = [
        f
        for f in denovo_risk(
            qc.flags, [r.get("sample") or "" for r in qc.sample_rows()]
        )
        if f not in acks or acks[f].status != "explained"
    ]
    if not risky:
        return
    show_detail = can_write()
    with (
        ui.row()
        .classes(
            "w-full items-start gap-2 p-3 mb-2 rounded border border-amber-400"
            " bg-amber-50 no-wrap"
        )
        .mark("denovo-qc-warning")
    ):
        ui.icon("warning", color="amber-9").classes("mt-0.5")
        with ui.column().classes("gap-0 min-w-0 flex-1"):
            ui.label(
                "Sample QC: de novo calls of this family may be unreliable."
            ).classes("font-semibold text-amber-900")
            if show_detail:
                for f in risky:
                    ack = acks.get(f)
                    ui.label(
                        f"{flag_title(f, True)} — {' – '.join(f.samples)}"
                        + (f": {f.evidence}" if f.evidence else "")
                        + (f" ({ack_summary(f, ack, True)})" if ack else "")
                    ).classes("text-sm text-amber-900")
            else:
                ui.label(
                    "Sample QC flags a pedigree or sample issue in this family"
                    " (relationship, duplicate or contamination)."
                ).classes("text-sm text-amber-900")
        if on_open_qc is not None:
            ui.button("Sample QC", on_click=on_open_qc).props(
                "flat dense no-caps color=amber-10"
            )


@dataclass
class _ConfigView:
    """One config's sub-tab: its rows and columns, and what the user chose."""

    has_sample: bool
    all_rows: List[Dict[str, Any]]
    all_columns: List[str]
    selected_cols: Dict[str, List[str]]
    selected_preset: Dict[str, str]
    filters: WombatFilters = field(default_factory=WombatFilters)
    # Sorting and page, kept across refreshes
    table_state: Dict[str, Any] = field(
        default_factory=lambda: {"sorting": [], "page": 0}
    )
    table: Optional[DataTable] = None


class _WombatLink:
    """The Wombat tab's part of the page link.

    It holds the shown config (``wcfg``) and that config's filters, sorting
    (``wsort``) and page (``wpage``). The other configs keep their own state
    while hidden, but only the shown one is in the link.
    """

    def __init__(self, link: Optional[PageLink], default_cfg: str, shown: str):
        self.link = link
        self.default_cfg = default_cfg
        self.current = shown
        self.views: Dict[str, _ConfigView] = {}

    def initial_state(self, cfg: str) -> tuple[WombatFilters, Dict[str, Any]]:
        """The filters and table state a config opens with.

        The link's apply to the config it names, once: a config opened later
        starts from the defaults.
        """
        if self.link is not None and cfg == self.current and not self.views:
            page = self.link.get("wpage") or "0"
            return filters_from_link(self.link), {
                "sorting": decode_sorting(self.link.get("wsort")),
                "page": int(page) if page.isdigit() else 0,
            }
        return WombatFilters(), {"sorting": [], "page": 0}

    def register(self, cfg: str, view: _ConfigView) -> None:
        self.views[cfg] = view

    def show(self, cfg: str) -> None:
        self.current = cfg
        self.write()

    def changed(self, cfg: str) -> None:
        if cfg == self.current:
            self.write()

    def write(self) -> None:
        if self.link is None:
            return
        self.link.set("wcfg", self.current, default=self.default_cfg)
        view = self.views.get(self.current)
        if view is None:
            self.link.drop(*LINK_KEYS, "wsort", "wpage")
        else:
            filters_to_link(view.filters, self.link)
            self.link.set("wsort", encode_sorting(view.table_state.get("sorting", [])))
            self.link.set("wpage", view.table_state.get("page", 0), default=0)
        self.link.push()


def render_wombat_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    family: Optional[Family] = None,
    on_open_qc: Optional[Callable[[], None]] = None,
    link: Optional[PageLink] = None,
) -> None:
    """Render the Wombat tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
        family: The family's pedigree, for the MOI column's sex and parental
            phenotypes
        on_open_qc: Opens the family's Sample QC tab, from the de novo warning
        link: The page link: the shown config, its filters, sorting and page
    """
    wombat_dir = get_family_path(store.data_dir, family_id) / "wombat"

    if not wombat_dir.exists():
        ui.label(f"No wombat directory found at: {wombat_dir}").classes(
            "text-gray-500 italic"
        )
        return

    wombat_files = find_wombat_files(wombat_dir, family_id)
    if not wombat_files:
        ui.label(
            f"No wombat TSV files found matching pattern in: {wombat_dir}"
        ).classes("text-gray-500 italic")
        return

    # Gene sets for the filter menu (shared across all configs): the user's
    # tiers, GeneTrek's lists and the data directory's params/genesets/*.tsv
    available_genesets = geneset_sources(store.data_dir, get_gene_scorer())

    # Pre-load DataFrames to determine which configs have data
    config_dfs: Dict[str, pl.DataFrame] = {}
    for wf in wombat_files:
        df = load_config_frame(wf["file_path"])
        if df is not None:
            config_dfs[wf["wombat_config"]] = df

    # Create subtabs — disable empty ones, add variant count badges
    with ui.tabs().classes("w-full").mark("wombat-subtabs") as wombat_subtabs:
        subtab_refs = {}
        for wf in wombat_files:
            cname = wf["wombat_config"]
            if cname in config_dfs:
                label = f"{cname} ({len(config_dfs[cname])})"
            else:
                label = cname
            tab = ui.tab(cname, label=label)
            if cname not in config_dfs:
                tab.props("disable")
            subtab_refs[cname] = tab

    # Pick first enabled sub-tab as default
    first_active_name = next(
        (
            wf["wombat_config"]
            for wf in wombat_files
            if wf["wombat_config"] in config_dfs
        ),
        wombat_files[0]["wombat_config"],
    )
    # A shared link names the config to show, when it has variants
    default_cfg = first_active_name
    linked_cfg = link.get("wcfg") if link is not None else None
    if linked_cfg in config_dfs:
        first_active_name = linked_cfg
    first_active = subtab_refs[first_active_name]
    wombat_link = _WombatLink(link, default_cfg, first_active_name)

    # Each sub-tab is built the first time it is shown, as the family page's
    # tabs are: a config nobody opens costs nothing. Keyed by config name,
    # which is the tab name the "update:model-value" event reports.
    unbuilt_panels: Dict[str, tuple] = {}

    with (
        ui.tab_panels(wombat_subtabs, value=first_active)
        .props("keep-alive")
        .classes("w-full")
    ):
        for wf in wombat_files:
            config_name = wf["wombat_config"]
            with ui.tab_panel(subtab_refs[config_name]):
                if config_name not in config_dfs:
                    ui.label("No variants found").classes("text-gray-500 italic mt-4")
                    continue

                with ui.element("div").classes("w-full") as panel_slot:
                    with ui.column().classes("w-full items-center py-16"):
                        ui.spinner(size="xl", color="blue")

            def build_panel(wf=wf, config_name=config_name) -> None:
                _render_config_panel(
                    store,
                    family_id,
                    cohort_name,
                    wf,
                    config_dfs[config_name],
                    selected_members,
                    data_table_refreshers,
                    available_genesets,
                    family,
                    on_open_qc,
                    wombat_link,
                )

            unbuilt_panels[config_name] = (panel_slot, build_panel)

    def build_subtab(config_name: str) -> None:
        entry = unbuilt_panels.pop(config_name, None)
        if entry is None:
            return
        panel_slot, build_panel = entry
        panel_slot.clear()
        with panel_slot:
            build_panel()

    def on_subtab_change(e) -> None:
        wombat_link.show(e.args)
        if e.args in unbuilt_panels:
            ui.timer(0.1, lambda: build_subtab(e.args), once=True)

    wombat_subtabs.on("update:model-value", on_subtab_change)
    build_subtab(first_active_name)


def _render_config_panel(
    store: Any,
    family_id: str,
    cohort_name: str,
    wf: Dict[str, Any],
    df: pl.DataFrame,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    available_genesets: List[Any],
    family: Optional[Family],
    on_open_qc: Optional[Callable[[], None]],
    wombat_link: _WombatLink,
) -> None:
    """One config's sub-tab: the config card, the filters and the table."""
    config_name = wf["wombat_config"]
    if _DENOVO_CONFIG.search(config_name):
        _render_denovo_qc_warning(store, family_id, on_open_qc)

    with ui.card().classes("w-full p-4"):
        ui.label(f"Wombat Configuration: {config_name}").classes(
            "text-lg font-semibold text-blue-700 mb-2"
        )
        with ui.row().classes("gap-4 mb-4"):
            ui.label("VEP Config:").classes("font-semibold")
            ui.badge(wf["vep_config"]).props("color=indigo")
        with ui.row().classes("gap-4"):
            ui.label("File Path:").classes("font-semibold")
            ui.label(str(wf["file_path"])).classes("text-sm text-gray-600 font-mono")

    # Display TSV content in a table
    try:
        # One row per variant and sample, other columns joined.
        df = aggregate_by_variant(df)
        all_rows, unknown_consequences, unknown_clinvar = build_rows(
            df, store.data_dir, family_id, family
        )
        notify_unknown_terms(unknown_consequences, unknown_clinvar)

        all_columns = wombat_columns(df.columns)
        # Auto-select preset based on keywords, its columns that the data has
        initial_preset = select_preset_for_config(config_name, VIEW_PRESETS)
        filters, table_state = wombat_link.initial_state(config_name)
        view = _ConfigView(
            has_sample="sample" in df.columns,
            all_rows=all_rows,
            all_columns=all_columns,
            selected_cols={
                "value": [
                    c for c in initial_preset.get("columns", []) if c in all_columns
                ]
            },
            selected_preset={"name": initial_preset["name"]},
            filters=filters,
            table_state=table_state,
        )
        wombat_link.register(config_name, view)

        # Capture the client context for use in callbacks
        from nicegui import context

        page_client = context.client

        with ui.column().classes("w-full"):

            def filters_changed() -> None:
                # A filter change shows the first page again
                view.table_state["page"] = 0
                render_data_table.refresh()
                wombat_link.changed(config_name)

            render_filter_panel(
                view.filters, available_genesets, on_change=filters_changed
            )

            @ui.refreshable
            def render_data_table() -> None:
                _render_table(
                    view,
                    store,
                    family_id,
                    cohort_name,
                    wf,
                    selected_members,
                    available_genesets,
                    on_curation_saved=refresh_curation,
                    on_table_change=lambda: wombat_link.changed(config_name),
                )

            def refresh_curation() -> None:
                """Reload validation + diagnostic data for all rows and refresh table."""
                add_curation_status(view.all_rows, store.data_dir, family_id)
                with page_client:
                    for rfn in data_table_refreshers:
                        rfn()

            data_table_refreshers.append(render_data_table.refresh)
            render_data_table()

    except pl.exceptions.NoDataError:
        ui.label("No data available").classes("text-gray-500 italic mt-4")
    except Exception as e:
        ui.label(f"Error reading file: {e}").classes("text-red-500 mt-4")


def _render_table(
    view: _ConfigView,
    store: Any,
    family_id: str,
    cohort_name: str,
    wf: Dict[str, Any],
    selected_members: Dict[str, List[str]],
    available_genesets: List[Any],
    on_curation_saved: Callable[[], None],
    on_table_change: Callable[[], None],
) -> None:
    """The count, preset, Columns/Stats/Download row, then the table itself."""
    config_name = wf["wombat_config"]

    # Filter rows by selected members
    if view.has_sample:
        rows = [
            r for r in view.all_rows if r.get("sample") in selected_members["value"]
        ]
    else:
        rows = view.all_rows
    total_before_filters = len(rows)
    rows = apply_filters(rows, view.filters, available_genesets)
    columns = table_columns(view.all_columns)

    def _apply_col_visibility():
        """Push current column selection to JS table."""
        if view.table is not None:
            visible = ["actions"] + list(view.selected_cols["value"])
            view.table.set_column_visibility(visible)

    with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
        row_label = (
            f"Data ({len(rows)} / {total_before_filters} rows)"
            if len(rows) < total_before_filters
            else f"Data ({len(rows)} rows)"
        )
        ui.label(row_label).classes("text-lg font-semibold text-blue-700")

        # Preset selector dropdown
        preset_select = ui.select(
            options={p["name"]: p["name"] for p in VIEW_PRESETS},
            value=view.selected_preset["name"],
            label="Preset",
        ).classes("w-48")

        ui.space()  # Push remaining items to the right

        # Column selector dialog
        col_dialog, _sync_col_selector = build_column_selector(
            all_columns=view.all_columns,
            selected_cols=view.selected_cols,
            on_visibility_change=_apply_col_visibility,
            presets=VIEW_PRESETS,
        )
        ui.button(
            "Columns",
            icon="view_column",
            on_click=col_dialog.open,
        ).props("outline color=blue size=sm")

        def show_stats_dialog(current_rows=rows):
            show_wombat_stats_dialog(
                current_rows,
                selected_members["value"],
                cohort_wombat_file(
                    store.data_dir, cohort_name, wf["vep_config"], config_name
                ),
                store.cohorts.get(cohort_name),
                wf["file_path"],
            )

        ui.button(
            "Stats",
            icon="bar_chart",
            on_click=show_stats_dialog,
        ).props("outline color=blue size=sm")

        def export() -> str:
            """The rows as filtered and sorted, columns shown."""
            if view.table is None:
                return ""
            shown = [c for c in view.selected_cols["value"] if is_export_column(c)]
            return rows_to_tsv(view.table.rows, shown)

        download_button(
            [
                Download(
                    "variants",
                    "The variants shown, with the columns shown",
                    lambda: export_filename(family_id, "wombat", config_name),
                    export,
                )
            ],
            mark="download-wombat",
        )

    def on_row_action(e):
        action = e.get("action", "")
        row_data = e.get("row", {})
        sample_val = row_data.get("sample", "")
        key = variant_key(row_data)

        if action == "open_diagnostic":
            show_diagnostic_dialog(
                family_id=family_id,
                variant_key=key,
                gene=row_data.get("VEP_SYMBOL", ""),
                impact=row_data.get("VEP_Consequence", ""),
                sample=sample_val,
                variant_type="snv",
                diagnostic_file=store.data_dir / "diagnostics" / "snvs.tsv",
                on_save_callback=lambda s: on_curation_saved(),
            )
        else:
            # Default: open variant dialog
            show_variant_dialog(
                cohort_name=cohort_name,
                family_id=family_id,
                chrom=row_data.get("#CHROM", ""),
                pos=row_data.get("POS", ""),
                ref=row_data.get("REF", ""),
                alt=row_data.get("ALT", ""),
                sample=sample_val,
                variant_data=row_data,
                on_save_callback=lambda s: on_curation_saved(),
            )

    # Restore table state (sorting / page) across refreshes, the rows
    # pre-sorted so JS gets them in order
    saved_sorting = view.table_state.get("sorting", [])
    presort_rows(rows, saved_sorting, columns)

    view.table = DataTable(
        columns=columns,
        rows=rows,
        row_key="Variant",
        pagination={"rowsPerPage": 10},
        visible_columns=["actions"] + list(view.selected_cols["value"]),
        on_row_action=on_row_action,
        initial_sorting=saved_sorting,
        initial_page=view.table_state.get("page", 0),
        state_holder=view.table_state,
        on_state_change=on_table_change,
    )

    def on_preset_change(e):
        """Handle preset selection change."""
        preset_name = e.value
        preset = next((p for p in VIEW_PRESETS if p["name"] == preset_name), None)
        if not preset:
            return
        # Filter columns to only those available in the data
        view.selected_cols["value"] = [
            col for col in preset.get("columns", []) if col in view.all_columns
        ]
        view.selected_preset["name"] = preset_name
        _apply_col_visibility()
        _sync_col_selector()

    preset_select.on_value_change(on_preset_change)
