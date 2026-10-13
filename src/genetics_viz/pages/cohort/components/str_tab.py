"""STR tab for the family page: short tandem repeats, one sub-tab per caller.

Only ExpansionHunter (``EH``) exists today; another caller (STRling) adds a
sibling sub-tab next to it.

The EH sub-tab pivots the family's genotype table to one row per locus, with a
column per member, so a child's allele sits beside their parents'. The cohort's
tables, when present, add the founders' distribution and flag the loci where
much of the cohort crosses a threshold - ARX_2 and XYLT1 do, which says more
about the catalog or thresholds than about the individuals. Clicking a locus
opens that distribution and the members' REViewer images.
"""

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

import polars as pl
from nicegui import ui

from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.data import get_static_prefix
from genetics_viz.utils.diagnostic_files import diagnostic_path
from genetics_viz.utils.sharding import get_sample_path, get_sample_url
from genetics_viz.utils.str_eh import (
    FLAGGED,
    STATUSES,
    SUPPORT,
    CohortContext,
    Member,
    diagnostic_badge,
    eh_path,
    family_locus_rows,
    family_members,
    find_eh_labels,
    get_family_eh_dir,
    load_cohort_context,
    locus_distribution,
    locus_gene,
    member_details,
    member_key,
    probe_eh_data,
    read_eh_table,
    reviewer_path,
    reviewer_reviewed,
    status_info,
    thresholds_text,
)
from genetics_viz.utils.str_plots import distribution_option
from genetics_viz.utils.str_validation import (
    INHERITANCES,
    VALIDATIONS,
    read_str_validations,
    save_str_validation,
    str_validation_file,
    str_variant_key,
)

logger = logging.getLogger(__name__)

__all__ = ["probe_eh_data", "render_str_tab"]


def render_str_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the STR tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name ("" on the standalone family page)
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    eh_dir = get_family_eh_dir(store.data_dir, family_id)
    labels = find_eh_labels(eh_dir, family_id)

    n_flagged: Optional[int] = None
    if labels:
        try:
            fam = read_eh_table(eh_path(eh_dir, family_id, labels[0], "genotypes"))
            n_flagged = (
                0
                if fam is None
                else fam.filter(pl.col("status").is_in(FLAGGED))["locus"].n_unique()
            )
        except Exception:
            logger.warning(
                "Failed to read the EH table of %s", family_id, exc_info=True
            )

    with ui.tabs().classes("w-full").mark("str-tabs") as str_subtabs:
        eh_label = "EH" if not n_flagged else f"EH ({n_flagged} flagged)"
        eh_tab = ui.tab("EH", label=eh_label).tooltip("ExpansionHunter")
        if not labels:
            eh_tab.props("disable")

    with ui.tab_panels(str_subtabs, value=eh_tab).classes("w-full"):
        with ui.tab_panel(eh_tab):
            if not labels:
                ui.label(f"No ExpansionHunter results in: {eh_dir}").classes(
                    "text-gray-500 italic"
                )
            else:
                render_eh_subtab(
                    store=store,
                    family_id=family_id,
                    cohort_name=cohort_name,
                    labels=labels,
                    selected_members=selected_members,
                    data_table_refreshers=data_table_refreshers,
                )


def _legend() -> None:
    with ui.row().classes("items-center gap-2 text-xs text-gray-600 mt-2 flex-wrap"):
        for key in STATUSES:
            info = status_info(key)
            ui.label(info["label"]).classes("px-2 py-0.5 rounded text-white").style(
                f"background-color: {info['color']}"
            )
        for support in SUPPORT.values():
            if support.get("code"):
                ui.label(f"·{support['code']}: {support['help']}")
        ui.label(
            "Badges are the members' alleles in repeat units, coloured by class;"
            " hover one for its CI, reads and cohort percentile. ↑: longer than"
            " both parents' alleles."
        )


def render_eh_subtab(
    store: Any,
    family_id: str,
    cohort_name: str,
    labels: List[str],
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """The ExpansionHunter table of one family."""
    eh_dir = get_family_eh_dir(store.data_dir, family_id)
    cohort = store.cohorts.get(cohort_name) if cohort_name else None
    state: Dict[str, Any] = {
        "label": labels[0],
        "fam": None,
        "members": [],
        "context": None,
        "context_loading": cohort is not None,
        "hide_normal": False,
    }

    def load_family() -> bool:
        path = eh_path(eh_dir, family_id, state["label"], "genotypes")
        try:
            fam = read_eh_table(path)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            fam = None
        state["fam"] = fam
        state["members"] = family_members(fam) if fam is not None else []
        return fam is not None and len(fam) > 0

    async def load_context() -> None:
        if cohort is None:
            return
        state["context_loading"] = True
        label = state["label"]
        try:
            context = await asyncio.to_thread(load_cohort_context, cohort, label)
        except Exception:
            logger.warning("Failed to read the cohort EH tables", exc_info=True)
            context = None
        if label != state["label"]:
            return  # the catalog changed while reading
        state["context"] = context
        state["context_loading"] = False
        content.refresh()

    if not load_family():
        ui.label("The ExpansionHunter table holds no genotypes").classes(
            "text-gray-500 italic"
        )
        return

    with ui.row().classes("items-center gap-4 mb-2 w-full"):
        ui.label("ExpansionHunter").classes("text-lg font-semibold text-blue-700")
        if len(labels) > 1:

            def on_label(e: Any) -> None:
                state["label"] = e.value
                state["context"] = None
                load_family()
                content.refresh()
                ui.timer(0.05, load_context, once=True)

            ui.select(
                options=labels,
                value=state["label"],
                label="Catalog",
                on_change=on_label,
            ).props("outlined dense").classes("w-48")
        else:
            ui.label(f"catalog {state['label']}").classes("text-xs text-gray-500")
        ui.space()
        ui.switch(
            "Hide loci where every member is normal",
            value=False,
            on_change=lambda e: (
                state.update({"hide_normal": e.value}),
                content.refresh(),
            ),
        ).props("dense")

    @ui.refreshable
    def content() -> None:
        fam: pl.DataFrame = state["fam"]
        members: List[Member] = state["members"]
        context: Optional[CohortContext] = state["context"]
        selected = [m.sample for m in members if m.sample in selected_members["value"]]

        if cohort is None:
            ui.label("No cohort context on this page.").classes("text-xs text-gray-500")
        elif state["context_loading"]:
            # The table waits for the cohort: built now, it would be replaced
            # before the browser mounted it, and its script would time out.
            with ui.row().classes("items-center gap-2"):
                ui.spinner(size="sm")
                ui.label(
                    f"Reading cohort {cohort.name}'s ExpansionHunter tables…"
                ).classes("text-xs text-gray-500")
            return
        elif context is None or not context.available:
            ui.label(
                f"No cohort-level ExpansionHunter results for {cohort.name}"
                f" (catalog {state['label']}): showing the family alone."
            ).classes("text-orange-600 text-sm")
        else:
            n = context.n_founders
            ui.label(
                f"Cohort {context.name}"
                + (f": compared with {n} founders" if n else "")
                + " (members with no genotyped parent)."
            ).classes("text-xs text-gray-500")

        if not selected:
            ui.label("No members selected").classes("text-gray-500 italic")
            return

        validations = _validation_statuses(store.data_dir, family_id, fam)
        diagnostics = _diagnostic_statuses(store.data_dir, family_id, fam)
        rows = family_locus_rows(
            fam, members, selected, context, validations, diagnostics
        )
        if state["hide_normal"]:
            rows = [r for r in rows if r["worst_rank"] >= 1]

        index = {m.sample: i for i, m in enumerate(members)}
        num = {"sortable": True, "sorting": "numerical"}
        columns: List[Dict[str, Any]] = [
            {
                "id": "actions",
                "header": "",
                "cellType": "action",
                "actionName": "view_locus",
                "actionIcon": "insights",
                "actionColor": "#1976d2",
                "actionTooltip": "Cohort distribution and REViewer",
                "sortable": False,
            },
            {"id": "locus", "header": "Locus", "group": "Locus", "sortable": True,
             "filter": {"type": "text"}},
            {"id": "disease", "header": "Disease", "group": "Locus", "sortable": True,
             "filter": {"type": "text"}, "maxWidth": 220},
            {"id": "inheritance", "header": "Inh.", "group": "Locus", "sortable": True},
            {"id": "motif", "header": "Motif", "group": "Locus", "sortable": False,
             "maxWidth": 110},
            {"id": "thresholds", "header": "Thresholds", "group": "Locus",
             "sortable": False},
            {"id": "worst", "header": "Status", "cellType": "badge_list",
             "sortField": "worst_rank", **num},
        ]  # fmt: skip
        for sample in selected:
            key = member_key(index[sample])
            columns.append(
                {
                    "id": key,
                    "header": members[index[sample]].column,
                    "group": "Members",
                    "cellType": "badge_list",
                    "sortField": f"{key}_max",
                    **num,
                }
            )
        columns.append(
            {
                "id": "flags",
                "header": "Flags",
                "cellType": "badge_list",
                "sortable": False,
            }
        )
        if context is not None and context.loci is not None:
            columns += [
                {"id": "founders_median", "header": "Median", "group": "Cohort founders", **num},
                {"id": "founders_p999", "header": "p99.9", "group": "Cohort founders", **num},
                {"id": "founders_max", "header": "Max", "group": "Cohort founders", **num},
                {"id": "cohort_flagged", "header": "Flagged", "group": "Cohort founders", **num},
            ]  # fmt: skip

        def on_row_action(e: Dict[str, Any]) -> None:
            locus = (e.get("row") or {}).get("locus")
            if locus:
                show_eh_locus_dialog(
                    store,
                    family_id,
                    fam,
                    members,
                    selected,
                    state["context"],
                    state["label"],
                    locus,
                    on_saved=content.refresh,
                )

        DataTable(
            columns=columns,
            rows=rows,
            row_key="id",
            pagination={"rowsPerPage": 60},
            on_row_action=on_row_action,
            initial_sorting=[{"id": "worst", "desc": True}],
        )
        _legend()

    content()
    data_table_refreshers.append(content.refresh)
    if cohort is not None:
        ui.timer(0.05, load_context, once=True)


def _curated_statuses(
    cmap: Dict[Tuple[str, str], Any],
    fam: pl.DataFrame,
    add_status: Callable[[Dict[str, Any], Any, str, str], None],
    field: str,
) -> Dict[Tuple[str, str], str]:
    """``(locus, sample)`` -> the curators' aggregated status, ignored rows left out.

    ``cmap`` is a validation or diagnostic map keyed by the STR variant key;
    ``add_status`` the loader's matching row helper, which writes ``field``.
    """
    if not cmap:
        return {}
    out: Dict[Tuple[str, str], str] = {}
    cols = ["locus", "chrom", "pos", "end", "sample"]
    for r in fam.select(cols).unique().to_dicts():
        key = str_variant_key(r["chrom"], r["pos"], r["end"], r["locus"])
        if (key, r["sample"]) in cmap:
            row: Dict[str, Any] = {}
            add_status(row, cmap, key, r["sample"])
            if row[field]:
                out[(r["locus"], r["sample"])] = row[field]
    return out


def _validation_statuses(
    data_dir: Any, family_id: str, fam: pl.DataFrame
) -> Dict[Tuple[str, str], str]:
    vmap = load_validation_map(str_validation_file(data_dir), family_id)
    return _curated_statuses(vmap, fam, add_validation_status_to_row, "Validation")


def _diagnostic_statuses(
    data_dir: Any, family_id: str, fam: pl.DataFrame
) -> Dict[Tuple[str, str], str]:
    dmap = load_diagnostic_map(diagnostic_path(data_dir, "STR"), family_id)
    return _curated_statuses(dmap, fam, add_diagnostic_status_to_row, "Diagnostic")


def _render_validation(
    data_dir: Any,
    family_id: str,
    locus: str,
    info: Dict[str, Any],
    details: List[Dict[str, Any]],
    dialog: Any,
    on_saved: Optional[Callable[[], None]],
) -> None:
    """The locus's validation history, and the form for curators."""
    key = str_variant_key(info.get("chrom"), info.get("pos"), info.get("end"), locus)
    ui.label("Validation").classes("text-base font-semibold text-blue-700 mt-2")
    history = read_str_validations(data_dir, family_id, key)
    if history:
        ui.table(
            columns=[
                {"name": field, "label": label, "field": field, "align": "left"}
                for field, label in (
                    ("Sample", "Sample"),
                    ("Validation", "Validation"),
                    ("Inheritance", "Inheritance"),
                    ("CuratedRepeats", "Curated repeats"),
                    ("User", "User"),
                    ("Comment", "Comment"),
                    ("Date", "Date"),
                )
            ],  # fmt: skip
            rows=[
                {
                    **h,
                    "Date": (h.get("Timestamp") or "")[:10],
                    "Validation": h.get("Validation", "")
                    + (" (ignored)" if h.get("Ignore") == "1" else ""),
                    "id": i,
                }
                for i, h in enumerate(history)
            ],
            row_key="id",
        ).props("dense flat").classes("w-full")
    else:
        ui.label("No validation recorded at this locus.").classes(
            "text-sm text-gray-500 italic"
        )

    if not can_write():
        ui.label("Curators and administrators can record validations.").classes(
            "text-xs text-gray-500"
        )
        return

    worst = max(details, key=lambda r: status_info(r.get("status"))["rank"])
    with ui.row().classes("items-end gap-3 w-full flex-wrap mt-2"):
        member = (
            ui.select(
                [r["sample"] for r in details], value=worst["sample"], label="Member"
            )
            .props("outlined dense")
            .classes("w-40")
        )
        validation = (
            ui.select(VALIDATIONS, value="present", label="Validation")
            .props("outlined dense")
            .classes("w-36")
        )
        inheritance = (
            ui.select(INHERITANCES, value="unknown", label="Inheritance")
            .props("outlined dense")
            .classes("w-40")
        )
        repeats = (
            ui.input("Curated repeats (optional)")
            .props("outlined dense")
            .classes("w-44")
        )
        comment = ui.input("Comment").props("outlined dense").classes("flex-1 min-w-48")

        def save() -> None:
            try:
                save_str_validation(
                    data_dir,
                    family_id=family_id,
                    variant_key=key,
                    sample=member.value,
                    user=get_current_user() or "",
                    validation=validation.value,
                    inheritance=inheritance.value,
                    curated_repeats=repeats.value or "",
                    comment=(comment.value or "").strip(),
                )
            except (PermissionError, ValueError) as ex:
                ui.notify(str(ex), type="negative")
                return
            except OSError as ex:
                ui.notify(f"Could not save the validation: {ex}", type="negative")
                return
            ui.notify(
                f"STR validation saved: {member.value} {validation.value} at {locus}",
                type="positive",
            )
            dialog.close()
            if on_saved is not None:
                on_saved()

        ui.button("Save validation", icon="save", on_click=save).props(
            "color=primary"
        ).mark("str-validation-save")


def _render_diagnostics(
    data_dir: Any,
    family_id: str,
    locus: str,
    info: Dict[str, Any],
    details: List[Dict[str, Any]],
    changed: Dict[str, bool],
) -> None:
    """Each member's diagnostic at the locus, and the button to its review dialog.

    The review dialog is the SNV/SV one: history, ignore toggles, and a save
    form shown to curators only. Saving refreshes this section in place and
    marks ``changed``, so the family table is redrawn once the locus dialog
    closes rather than from under it.
    """
    key = str_variant_key(info.get("chrom"), info.get("pos"), info.get("end"), locus)
    path = diagnostic_path(data_dir, "STR")
    ui.label("Diagnostic").classes("text-base font-semibold text-blue-700 mt-2")
    # The review dialogs live here, outside the refreshable section: a refresh
    # from an ignore toggle would otherwise delete the dialog that sent it.
    anchor = ui.element("div")

    @ui.refreshable
    def section() -> None:
        dmap = load_diagnostic_map(path, family_id)
        with ui.row().classes("gap-2 flex-wrap w-full"):
            for r in details:
                sample = r["sample"]
                row: Dict[str, Any] = {}
                add_diagnostic_status_to_row(row, dmap, key, sample)
                with ui.row().classes(
                    "items-center gap-2 border rounded px-2 py-1 no-wrap"
                ):
                    ui.label(sample).classes("text-sm font-medium")
                    if row["Diagnostic"]:
                        badge = diagnostic_badge(row["Diagnostic"])
                        ui.label(row["Diagnostic"]).classes(
                            "text-xs text-white rounded px-1"
                        ).style(f"background-color: {badge['color']}")
                    else:
                        ui.label("none").classes("text-xs text-gray-400 italic")

                    def make_open(r: Dict[str, Any]) -> Callable[[], None]:
                        def handler() -> None:
                            with anchor:
                                show_diagnostic_dialog(
                                    family_id=family_id,
                                    variant_key=key,
                                    gene=locus_gene(locus),
                                    impact=status_info(r.get("status"))["key"],
                                    sample=r["sample"],
                                    variant_type="str",
                                    diagnostic_file=path,
                                    on_save_callback=on_change,
                                )

                        return handler

                    ui.button(icon="biotech", on_click=make_open(r)).props(
                        "flat dense round size=sm"
                    ).tooltip("Diagnostic review").mark(
                        f"str-diagnostic-open str-diagnostic-{sample}"
                    )

    def on_change(_action: str) -> None:
        changed["value"] = True
        section.refresh()

    section()


def _no_reviewer_reason(store: Any, row: Dict[str, Any], label: str) -> str:
    """Why a member has no REViewer image at a locus.

    The pipeline draws flagged genotypes only, so an unflagged member's
    missing image is the expected case, not an anomaly.
    """
    if not reviewer_reviewed(store.data_dir, row["sample"], label):
        return "REViewer not run yet"
    if row.get("status") not in FLAGGED:
        return "not flagged, not drawn"
    return "no REViewer image"


def show_eh_locus_dialog(
    store: Any,
    family_id: str,
    fam: pl.DataFrame,
    members: List[Member],
    selected: List[str],
    context: Optional[CohortContext],
    label: str,
    locus: str,
    on_saved: Optional[Callable[[], None]] = None,
) -> None:
    """One locus: the cohort distribution, each member's call, validations and
    diagnostics, REViewer."""
    details = [
        r for r in member_details(fam, locus, context) if r["sample"] in selected
    ]
    if not details:
        return
    info = details[0]
    # A diagnostic saved in the review dialog redraws the family table on close.
    changed = {"value": False}

    def on_close(e: Any) -> None:
        if not e.value and changed["value"] and on_saved is not None:
            changed["value"] = False
            on_saved()

    with ui.dialog().props("full-width") as dialog, ui.card().classes("w-full"):
        dialog.on_value_change(on_close)
        with ui.row().classes("w-full items-center"):
            ui.label(f"{locus}").classes("text-lg font-semibold text-blue-700")
            ui.label(
                f"{info.get('chrom')}:{info.get('pos')}-{info.get('end')}"
                f" · {info.get('motif') or ''} · {info.get('disease') or 'no disease'}"
                f" · {info.get('inheritance') or '–'} · {thresholds_text(info)}"
            ).classes("text-sm text-gray-600")
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")

        if context is not None and context.genotypes is not None:
            counts = locus_distribution(context.genotypes, locus)
            marked = [
                (
                    r["sample"],
                    [
                        int(a)
                        for a in (r.get("allele1"), r.get("allele2"))
                        if a is not None
                    ],
                    r.get("status") or "",
                )
                for r in details
            ]
            option = distribution_option(
                counts, info, marked, title=f"{locus} in cohort {context.name}"
            )
            height = 420 + 22 * len(marked)
            # shrink-0: the dialog card is a flex column and would squeeze the
            # chart below its height to fit the table under it.
            ui.echart(option).classes("w-full shrink-0").style(f"height: {height}px")
        else:
            ui.label(
                "No cohort genotype table: the founders' distribution is unavailable."
            ).classes("text-gray-500 italic text-sm")

        columns = [
            {"name": "sample", "label": "Member", "field": "sample", "align": "left"},
            {"name": "status", "label": "Status", "field": "status", "align": "left"},
            {"name": "genotype", "label": "Alleles", "field": "alleles", "align": "left"},
            {"name": "ci", "label": "CI", "field": "ci", "align": "left"},
            {"name": "support", "label": "Support", "field": "support", "align": "left"},
            {"name": "reads", "label": "Spanning / flanking / in-repeat", "field": "reads",
             "align": "left"},
            {"name": "coverage", "label": "Coverage", "field": "coverage"},
            {"name": "filter", "label": "Filter", "field": "filter"},
            {"name": "parents", "label": "Parents' max", "field": "parents"},
            {"name": "pct", "label": "Cohort pct", "field": "pct"},
        ]  # fmt: skip

        def pair(r: Dict[str, Any], suffix: str) -> str:
            values = [r.get(f"allele1{suffix}"), r.get(f"allele2{suffix}")]
            return " / ".join(str(v) for v in values if v is not None) or "–"

        def reads(r: Dict[str, Any]) -> str:
            parts = []
            for i in (1, 2):
                if r.get(f"allele{i}") is None:
                    continue
                parts.append(
                    f"{r.get(f'allele{i}_spanning') or 0}/"
                    f"{r.get(f'allele{i}_flanking') or 0}/"
                    f"{r.get(f'allele{i}_inrepeat') or 0}"
                )
            return " · ".join(parts) or "–"

        ui.table(
            columns=columns,
            rows=[
                {
                    "sample": r["sample"],
                    "status": status_info(r.get("status"))["label"],
                    "alleles": pair(r, ""),
                    "ci": pair(r, "_ci"),
                    "support": pair(r, "_support"),
                    "reads": reads(r),
                    "coverage": (
                        f"{r['locus_coverage']:.1f}x"
                        if r.get("locus_coverage") is not None
                        else "–"
                    ),
                    "filter": r.get("filter") or "–",
                    "parents": r.get("parents_max_allele") or "–",
                    "pct": "–"
                    if r.get("cohort_percentile") is None
                    else r["cohort_percentile"],
                }
                for r in details
            ],
            row_key="sample",
        ).props("dense flat").classes("w-full")

        _render_validation(
            store.data_dir, family_id, locus, info, details, dialog, on_saved
        )
        _render_diagnostics(store.data_dir, family_id, locus, info, details, changed)

        ui.label("REViewer").classes("text-base font-semibold text-blue-700 mt-2")
        prefix = get_static_prefix()
        for r in details:
            sample = r["sample"]
            path = reviewer_path(store.data_dir, sample, label, locus)
            if path is None:
                ui.label(f"{sample}: {_no_reviewer_reason(store, r, label)}").classes(
                    "text-sm text-gray-500 italic"
                )
                continue
            relative = path.relative_to(get_sample_path(store.data_dir, sample))
            url = (
                f"{prefix}/{get_sample_url(store.data_dir, sample)}"
                f"/{relative.as_posix()}"
            )
            with ui.expansion(f"{sample} — REViewer", icon="image").classes("w-full"):
                # loading=lazy: the images are up to a megabyte each, and the
                # browser would fetch them all even while collapsed. w-full: the
                # expansion's column aligns flex-start, so the wrapper would
                # shrink to the SVG's default 300px and the image with it.
                ui.html(
                    f'<img src="{url}" loading="lazy" alt="REViewer {sample} {locus}"'
                    ' style="width: 100%">',
                    sanitize=False,
                ).classes("w-full")

    dialog.open()
