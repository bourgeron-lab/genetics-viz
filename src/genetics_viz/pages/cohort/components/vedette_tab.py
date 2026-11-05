"""SVs tab, Vedette sub-tab: the family's Vedette calls, one row per carrier.

The family export (``svs/vedette/<FID>.vedette.tsv``, ``.qc.tsv``, ``.json``;
format in Vedette's ``docs/EXPORT.md``) is read off-thread, then filtered in
Python: members (the page's selection), presets and the Filters card (type,
kind, own call, size, cohort carrier fraction, carrier GQ, impact level,
inheritance, genesets, validation). The filters are in the page link (``v…``
keys). Members Vedette flags in its QC file are listed above the table.

A row opens the SV dialog with the Vedette source (``vedette_details.py``):
Vedette's details card, base-pair read windows when the call is precise, and
validations saved to ``validations/vedette.tsv``. Diagnostics are shared with
every SV (``diagnostics/svs.tsv``), keyed the same way (``chr:start-end:del``).
"""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.utils.diagnostic_snapshot import vedette_snapshot
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.filters import (
    create_geneset_filter_menu,
    create_validation_filter_menu,
)
from genetics_viz.components.page_link import PageLink
from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.vedette_details import vedette_source
from genetics_viz.utils.export import export_filename, is_export_column, rows_to_tsv
from genetics_viz.utils.gene_scoring import (
    build_gene_badges,
    get_gene_scorer,
    tier_legend_html,
)
from genetics_viz.utils.genesets import geneset_sources, selected_symbols
from genetics_viz.utils.genetrek import search_url
from genetics_viz.utils.validation_files import validation_path
from genetics_viz.utils.validation_status import FILTER_STATUSES, filter_by_status
from genetics_viz.utils.vedette import (
    LINK_KEYS,
    VEDETTE_CONFIG,
    ExportInfo,
    MemberQC,
    VedetteFilters,
    apply_filters,
    carrier_rows,
    default_preset,
    family_vedette_file,
    filters_from_params,
    filters_to_params,
    impact_level_codes,
    inheritance_codes,
    inheritance_info,
    kind_codes,
    kind_info,
    level_info,
    own_call_codes,
    own_call_info,
    preset_names,
    qc_flag_label,
    read_export_info,
    read_member_qc,
    read_vedette_table,
)

logger = logging.getLogger(__name__)

#: (id, header, shown by default). Columns an export lacks stay empty.
_COLUMNS: List[Tuple[str, str, bool]] = [
    ("chr:start-end", "Locus", True),
    ("svtype", "Type", True),
    ("kind", "Kind", True),
    ("size", "Size", True),
    ("sample", "Carrier", True),
    ("GT", "GT", True),
    ("CN", "CN", True),
    ("GQ", "GQ", True),
    ("inheritance", "Inheritance", True),
    ("DNQ", "DNQ", False),
    ("family", "Family (GT·CN)", True),
    ("own_call", "Own call", True),
    ("support", "CN support", False),
    ("impact_score", "Impact", True),
    ("impact_level", "Level", True),
    ("genes", "Genes", True),
    ("impact_reason", "Impact reason", True),
    ("cohort_carrier_frac", "Cohort carrier fraction", True),
    ("ref_carrier_frac", "Reference carrier fraction", False),
    ("sqmax", "SQ max", False),
    ("n_sites", "Cohort sites", False),
    ("mosaic_cn", "Mosaic CN", False),
    ("lcr_ovl", "LCR overlap", False),
    ("segdup_ovl", "Segdup overlap", False),
    ("filter", "Filter", False),
    ("rescue", "Rescue", False),
    ("qc_alerts", "QC alerts", False),
    ("Validation", "Validation", True),
    ("Diagnostic", "Diagnostic", True),
    ("site_id", "Row id", False),
    # Exports before Vedette v1.4.0 only.
    ("cohort_af", "Cohort AF", False),
    ("cohort_ncarrier", "Cohort carriers", False),
]
_NUMBER_COLUMNS = {
    "size", "CN", "GQ", "DNQ", "impact_score", "cohort_carrier_frac",
    "ref_carrier_frac", "sqmax", "n_sites", "mosaic_cn", "lcr_ovl", "segdup_ovl",
    "cohort_af", "cohort_ncarrier",
}  # fmt: skip
#: Coded columns drawn as a badge, and the row field holding its colour.
_BADGE_COLUMNS = {
    "svtype": "_type_color",
    "kind": "_kind_color",
    "inheritance": "_inheritance_color",
    "own_call": "_own_call_color",
    "support": "_support_color",
    "impact_level": "_impact_color",
}
_INITIAL_SORTING = [{"id": "impact_score", "desc": True}]


def _member_labels(store: Any, cohort_name: str, family_id: str) -> Dict[str, str]:
    """``father`` / ``mother`` / ``child`` per member, from the pedigree."""
    cohort = store.get_cohort(cohort_name) if cohort_name else None
    if cohort is None:
        return {}
    members = cohort.get_family_members(family_id)
    fathers = {m.get("Father") for m in members}
    mothers = {m.get("Mother") for m in members}
    labels: Dict[str, str] = {}
    for m in members:
        sid = m["Sample ID"]
        if sid in fathers:
            labels[sid] = "father"
        elif sid in mothers:
            labels[sid] = "mother"
        elif m.get("Father") not in (None, "", "0", "-") or m.get("Mother") not in (
            None,
            "",
            "0",
            "-",
        ):
            labels[sid] = "child"
    return labels


def _load(
    data_dir: Path, family_id: str, labels: Dict[str, str], scorer: Any
) -> Tuple[List[Dict[str, Any]], ExportInfo, List[MemberQC]]:
    """Carrier rows with their gene badges, the export's JSON and the members'
    QC (blocking).

    A family has thousands of carrier rows: built here, off-thread. ``scorer``
    is the user's, taken on the event loop.
    """
    info = read_export_info(data_dir, family_id)
    df = read_vedette_table(family_vedette_file(data_dir, family_id))
    rows = carrier_rows(df, labels=labels, members=info.members)
    for row in rows:
        # Exonic genes bordered, as for WisecondorX calls.
        row["GeneBadges"] = build_gene_badges(
            row["_genes"], scorer, bordered=row["_exonic"]
        )
        row["_genetrek_url"] = search_url([row["chr:start-end"]])
        row["_original_locus"] = row["chr:start-end"]
    return rows, info, read_member_qc(data_dir, family_id)


def _badge_html(text: str, color: str) -> str:
    return (
        f'<span class="q-badge" style="background-color: {color}; color: white; '
        f'padding: 2px 6px; border-radius: 4px; font-size: 12px;">{text}</span>'
    )


def render_vedette_subtab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    link: Optional[PageLink] = None,
) -> None:
    """The Vedette sub-tab: loads the family table off-thread, then builds."""
    path = family_vedette_file(store.data_dir, family_id)
    container = ui.column().classes("w-full")
    with container:
        with ui.column().classes("w-full items-center py-12"):
            ui.spinner(size="lg", color="blue")
            ui.label("Loading Vedette calls...").classes("text-gray-600")

    async def load() -> None:
        labels = _member_labels(store, cohort_name, family_id)
        scorer = get_gene_scorer()  # reads the session: on the event loop
        try:
            rows, info, qc = await asyncio.to_thread(
                _load, store.data_dir, family_id, labels, scorer
            )
        except Exception as e:  # a degraded state, shown, never raised
            logger.warning("Cannot read %s", path, exc_info=True)
            container.clear()
            with container:
                ui.label(f"Error reading {path}: {e}").classes("text-red-600")
            return
        container.clear()
        with container:
            _build(
                store,
                family_id,
                cohort_name,
                path,
                rows,
                info,
                qc,
                selected_members,
                data_table_refreshers,
                link,
            )

    ui.timer(0.1, load, once=True)


def _build(
    store: Any,
    family_id: str,
    cohort_name: str,
    path: Path,
    all_rows: List[Dict[str, Any]],
    info: ExportInfo,
    qc: List[MemberQC],
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    link: Optional[PageLink],
) -> None:
    from nicegui import context

    page_client = context.client
    data_dir = store.data_dir
    scorer = get_gene_scorer()

    with ui.card().classes("w-full p-4"):
        ui.label("Vedette calls").classes("text-lg font-semibold text-blue-700 mb-2")
        with ui.row().classes("gap-4"):
            ui.label("File Path:").classes("font-semibold")
            ui.label(str(path)).classes("text-sm text-gray-600 font-mono")
        _render_versions(info)
        ui.html(tier_legend_html(scorer, exonic_border=True), sanitize=False).classes(
            "mt-3"
        )
        with ui.row().classes("gap-2 mt-2 items-center flex-wrap"):
            for title, codes, describe in (
                ("Impact:", impact_level_codes(), level_info),
                ("Inheritance:", inheritance_codes(), inheritance_info),
                ("Kind:", kind_codes(), kind_info),
                ("Own call:", own_call_codes(), own_call_info),
            ):
                ui.label(title).classes("text-sm font-semibold ml-2")
                for code in codes:
                    shown_info = describe(code)
                    ui.html(
                        _badge_html(shown_info["label"], shown_info["color"]),
                        sanitize=False,
                    )
    _render_member_qc(qc)

    if not all_rows:
        ui.label("No Vedette calls in this family").classes("text-gray-500 italic")
        return

    by_id = {r["id"]: r for r in all_rows}

    validation_file = validation_path(data_dir, "Vedette")
    diagnostic_file = data_dir / "diagnostics" / "svs.tsv"

    def reload_curation() -> None:
        """Validation and diagnostic status, and the curated locus, per row."""
        vmap = load_validation_map(validation_file, family_id)
        dmap = load_diagnostic_map(diagnostic_file, family_id)
        for row in all_rows:
            key, sample = row["_variant_key"], row["sample"]
            add_validation_status_to_row(row, vmap, key, sample)
            add_diagnostic_status_to_row(row, dmap, key, sample)
            row["chr:start-end"] = row["_original_locus"]
            row["IsCurated"] = False
            row["_curated_tooltip"] = ""
            curated = [
                v
                for v in vmap.get((key, sample), [])
                if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
            ]
            if curated:
                newest = max(curated, key=lambda v: v[6])
                start = newest[4] or str(row["start"])
                end = newest[5] or str(row["end"])
                row["chr:start-end"] = f"{row['contig']}:{start}-{end}"
                row["IsCurated"] = True
                row["_curated_tooltip"] = (
                    f"Original: {row['_original_locus']}\n"
                    f"Curated: {row['chr:start-end']}"
                )

    reload_curation()

    filters = (
        filters_from_params(link.get, link.get_list)
        if link is not None
        else VedetteFilters.preset(default_preset())
    )
    genesets = geneset_sources(data_dir, scorer)
    geneset_state = {"value": list(filters.genesets)}
    validation_state = {"value": list(filters.validations or FILTER_STATUSES)}
    shown = {"value": [c for c, _, on in _COLUMNS if on]}
    table_ref: Dict[str, Any] = {"dt": None}

    def push_link() -> None:
        if link is None:
            return
        link.drop(*LINK_KEYS)
        for key, value in filters_to_params(filters).items():
            if key == "vpre":
                link.set(key, value, default=default_preset())
            elif isinstance(value, list):
                link.set(key, value, default=[])
            else:
                link.set(key, value)
        link.push()

    def changed() -> None:
        filters.genesets = list(geneset_state["value"])
        statuses = validation_state["value"]
        filters.validations = (
            [] if set(statuses) >= set(FILTER_STATUSES) else list(statuses)
        )
        push_link()
        presets_row.refresh()
        render_table.refresh()

    def apply_preset(name: str) -> None:
        preset = VedetteFilters.preset(name)
        for attr in VedetteFilters._PRESET_FIELDS:
            setattr(filters, attr, getattr(preset, attr))
        filter_controls.refresh()
        changed()

    @ui.refreshable
    def presets_row() -> None:
        with ui.row().classes("items-center gap-2 mt-4"):
            ui.label("Presets:").classes("text-sm font-semibold")
            for name in preset_names():
                cfg = VEDETTE_CONFIG["presets"][name]
                active = filters.matches_preset(name)
                ui.button(
                    cfg.get("label", name), on_click=lambda _, n=name: apply_preset(n)
                ).props(
                    "unelevated color=blue size=sm"
                    if active
                    else "outline color=blue size=sm"
                ).tooltip(cfg.get("tooltip", "")).mark(f"vedette-preset-{name}")

    presets_row()

    def set_value(attr: str, cast: Callable[[Any], Any]) -> Callable[[Any], None]:
        def handler(e) -> None:
            value = e.value
            setattr(filters, attr, None if value in (None, "") else cast(value))
            changed()

        return handler

    def toggle(attr: str, code: str) -> Callable[[Any], None]:
        def handler(e) -> None:
            current = list(getattr(filters, attr))
            if e.value and code not in current:
                current.append(code)
            elif not e.value and code in current:
                current.remove(code)
            setattr(filters, attr, current)
            changed()

        return handler

    def code_menu(
        label: str,
        attr: str,
        codes: List[str],
        describe: Callable[[str], Dict[str, str]],
        hint: str = "",
    ) -> None:
        """A button whose menu ticks codes of one coded column."""
        with ui.button(label, icon="filter_list").props("outline color=blue size=sm"):
            with ui.menu(), ui.column().classes("p-2 gap-0"):
                if hint:
                    ui.label(hint).classes("text-xs text-gray-500")
                for code in codes:
                    ui.checkbox(
                        describe(code)["label"],
                        value=code in getattr(filters, attr),
                        on_change=toggle(attr, code),
                    ).classes("text-sm").mark(f"vedette-{attr}-{code}")

    @ui.refreshable
    def filter_controls() -> None:
        with ui.row().classes("gap-4 items-center flex-wrap p-2"):
            with ui.row().classes("gap-1 items-center"):
                for t in VEDETTE_CONFIG.get("types", {}):
                    ui.checkbox(
                        t, value=t in filters.types, on_change=toggle("types", t)
                    ).mark(f"vedette-type-{t}")
            ui.number(
                "Min size (bp)", value=filters.min_size, min=0,
                on_change=set_value("min_size", int),
            ).props("dense outlined debounce=600").classes("w-32")  # fmt: skip
            ui.number(
                "Max size (bp)", value=filters.max_size, min=0,
                on_change=set_value("max_size", int),
            ).props("dense outlined debounce=600").classes("w-32")  # fmt: skip
            ui.number(
                "Max cohort carrier fraction", value=filters.max_freq, min=0, max=1,
                step=0.001, on_change=set_value("max_freq", float),
            ).props("dense outlined debounce=600").classes("w-48").mark(
                "vedette-max-freq"
            )  # fmt: skip
            ui.number(
                "Min carrier GQ", value=filters.min_gq, min=0,
                on_change=set_value("min_gq", int),
            ).props("dense outlined debounce=600").classes("w-32").mark(
                "vedette-min-gq"
            )  # fmt: skip
            code_menu("Kind", "kinds", kind_codes(), kind_info)
            code_menu(
                "Own call", "own_calls", own_call_codes(), own_call_info,
                hint="None ticked: any",
            )  # fmt: skip
            with ui.button("Impact", icon="filter_list").props(
                "outline color=blue size=sm"
            ):
                with ui.menu(), ui.column().classes("p-2 gap-0"):
                    for code in impact_level_codes():
                        ui.checkbox(
                            level_info(code)["label"],
                            value=code in filters.impact_levels,
                            on_change=toggle("impact_levels", code),
                        ).classes("text-sm")
            with ui.button("Inheritance", icon="filter_list").props(
                "outline color=blue size=sm"
            ):
                with ui.menu(), ui.column().classes("p-2 gap-0"):
                    ui.label("None ticked: any").classes("text-xs text-gray-500")
                    for code in [*inheritance_codes(), ""]:
                        ui.checkbox(
                            inheritance_info(code)["label"] if code else "not given",
                            value=code in filters.inheritance,
                            on_change=toggle("inheritance", code),
                        ).classes("text-sm").mark(f"vedette-inh-{code or 'none'}")
            create_geneset_filter_menu(
                geneset_state, genesets, on_change=changed, button_props="size=sm"
            )
            create_validation_filter_menu(
                all_statuses=list(FILTER_STATUSES),
                filter_state=validation_state,
                on_change=changed,
                label="Validation",
                button_classes="",
                button_props="size=sm",
            )

    with ui.card().classes("w-full p-0 mt-2"):
        with (
            ui.expansion("Filters", icon="filter_list", value=False)
            .classes("w-full")
            .props("header-class='text-lg font-semibold text-blue-700'")
        ):
            filter_controls()

    def make_columns() -> List[Dict[str, Any]]:
        cols: List[Dict[str, Any]] = [
            {
                "id": "actions",
                "header": "",
                "cellType": "action",
                "actionName": "view_sv",
                "actionIcon": "visibility",
                "actionColor": "#1976d2",
                "actionTooltip": "View reads and validate",
                "sortable": False,
            }
        ]
        for col, header, _ in _COLUMNS:
            if col not in shown["value"]:
                continue
            c: Dict[str, Any] = {"id": col, "header": header, "sortable": True}
            if col == "chr:start-end":
                c.update(
                    cellType="curated_locus",
                    curatedField="IsCurated",
                    tooltipField="_curated_tooltip",
                )
            elif col in _BADGE_COLUMNS:
                c.update(cellType="badge", colorField=_BADGE_COLUMNS[col])
            elif col == "family":
                c.update(cellType="badge_list", badgesField="family_badges")
            elif col == "genes":
                c.update(cellType="gene_badge", badgesField="GeneBadges")
            elif col == "Validation":
                c["cellType"] = "validation"
            elif col == "Diagnostic":
                c["cellType"] = "diagnostic"
            elif col in _NUMBER_COLUMNS:
                c["cellType"] = "number"
            cols.append(c)
        return cols

    def on_saved() -> None:
        reload_curation()
        with page_client:
            for refresh in data_table_refreshers:
                refresh()

    def on_row_action(e: Dict[str, Any]) -> None:
        sent = e.get("row", {}) or {}
        row = by_id.get(sent.get("id"), sent)
        if not row.get("sample"):
            ui.notify("Missing carrier", type="warning")
            return
        if e.get("action") == "open_diagnostic":
            cohort = store.get_cohort(cohort_name) if cohort_name else None
            family = cohort.families.get(family_id) if cohort else None
            show_diagnostic_dialog(
                details=vedette_snapshot(
                    row,
                    source=f"Vedette {info.caller}".strip(),
                    rows=all_rows,
                    phenotypes={
                        s.sample_id: s.phenotype
                        for s in (family.samples if family else [])
                    },
                ),
                family_id=family_id,
                variant_key=row["_variant_key"],
                gene=", ".join(row.get("_genes", [])),
                impact="LOSS" if row["svtype"] == "DEL" else "GAIN",
                sample=row["sample"],
                variant_type="sv",
                diagnostic_file=diagnostic_file,
                on_save_callback=lambda _s: on_saved(),
            )
            return
        show_sv_dialog(
            cohort_name=cohort_name,
            family_id=family_id,
            chrom=row["contig"],
            start=str(row["start"]),
            end=str(row["end"]),
            sample=row["sample"],
            sv_data=row,
            on_validation_saved=on_saved,
            source=vedette_source(data_dir, family_id, row, info),
        )

    def export() -> str:
        table = table_ref.get("dt")
        if table is None:
            return ""
        columns = [c for c in shown["value"] if is_export_column(c)]
        return rows_to_tsv(table.rows, columns)

    @ui.refreshable
    def render_table() -> None:
        symbols = (
            selected_symbols(genesets, filters.genesets) if filters.genesets else None
        )
        rows = apply_filters(
            all_rows,
            filters,
            members=selected_members["value"],
            gene_symbols=symbols,
        )
        rows = filter_by_status(rows, filters.validations)
        rows.sort(key=lambda r: -(r.get("impact_score") or -99.0))
        n_events = len({r["site_id"] for r in rows})
        with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
            ui.label(f"Data ({len(rows)} rows, {n_events} events)").classes(
                "text-lg font-semibold text-blue-700"
            ).mark("vedette-count")
            ui.space()
            with ui.button("Columns", icon="view_column").props(
                "outline color=blue size=sm"
            ):
                with ui.menu(), ui.column().classes("p-2 gap-0"):
                    for col, header, _ in _COLUMNS:
                        ui.checkbox(
                            header,
                            value=col in shown["value"],
                            on_change=lambda e, c=col: set_shown(c, e.value),
                        ).classes("text-sm")
            download_button(
                [
                    Download(
                        "calls",
                        "The calls shown, with the columns shown",
                        lambda: export_filename(family_id, "vedette"),
                        export,
                    )
                ],
                mark="download-vedette",
            )
        with ui.card().classes("w-full"):
            table_ref["dt"] = DataTable(
                columns=make_columns(),
                rows=rows,
                pagination={"rowsPerPage": 25},
                on_row_action=on_row_action,
                initial_sorting=_INITIAL_SORTING,
            )

    def set_shown(col: str, value: bool) -> None:
        current = set(shown["value"])
        current = current | {col} if value else current - {col}
        shown["value"] = [c for c, _, _ in _COLUMNS if c in current]
        render_table.refresh()

    data_table_refreshers.append(render_table.refresh)
    render_table()


def _render_versions(info: ExportInfo) -> None:
    """Caller, model and schema versions, and a warning on a newer schema."""
    parts = [
        ("Caller", info.caller),
        ("Model", info.model_version),
        ("Schema", "" if info.schema_version is None else str(info.schema_version)),
        ("Cohort run", info.cohort),
        ("Written", info.generated.replace("T", " ")),
    ]
    with ui.row().classes("gap-4 items-center flex-wrap"):
        for label, value in parts:
            if value:
                ui.label(f"{label}:").classes("font-semibold")
                ui.label(value).classes("text-sm text-gray-600 font-mono")
    if info.newer_schema:
        ui.label(
            f"This export uses schema {info.schema_version}, newer than the one "
            "this app was written for: some columns may be missing or misread."
        ).classes("text-sm text-orange-700").mark("vedette-schema-warning")
    elif info.schema_version is None:
        ui.label(
            "Export written before Vedette v1.4.0: one row per site, coordinates "
            "one base to the right, truncated gene lists."
        ).classes("text-sm text-gray-500 italic")


def _render_member_qc(qc: List[MemberQC]) -> None:
    """The members Vedette's QC flags, above the table."""
    flagged = [m for m in qc if m.flags]
    if not flagged:
        return
    with (
        ui.card()
        .classes("w-full p-3 mt-2 bg-amber-50 border border-amber-300")
        .mark("vedette-qc")
    ):
        ui.label("Vedette QC flags: read these members' calls with care").classes(
            "font-semibold text-amber-900"
        )
        for m in flagged:
            facts = []
            if m.carrier_excess is not None:
                facts.append(f"{m.carrier_excess:g}× the cohort median of PASS sites")
            if m.denovo_candidates is not None:
                facts.append(f"{m.denovo_candidates} de novo candidates")
            if m.wave_score is not None:
                facts.append(f"wave score {m.wave_score:g}")
            if m.karyotype:
                facts.append(f"karyotype {m.karyotype}")
            with ui.row().classes("items-center gap-2 flex-wrap"):
                ui.label(m.sample).classes("font-mono font-semibold")
                for code in m.flags:
                    ui.badge(code).props("color=amber-9").tooltip(qc_flag_label(code))
                ui.label(" · ".join(facts)).classes("text-sm text-gray-700")
