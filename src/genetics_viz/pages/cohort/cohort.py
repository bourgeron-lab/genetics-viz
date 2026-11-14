"""Cohort detail page - displays pedigree members with filters."""

import asyncio
import csv
import logging
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from fastapi import Request
from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.page_link import PageLink
from genetics_viz.components.tab_badge import TabPills
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.models import is_sequenced_parent
from genetics_viz.pages.cohort.components.cohort_ancestry_tab import (
    probe_cohort_ancestry_data,
    render_cohort_ancestry_tab,
)
from genetics_viz.pages.cohort.components.cohort_pgs_tab import (
    probe_cohort_pgs_data,
    render_cohort_pgs_tab,
)
from genetics_viz.pages.cohort.components.cohort_roh_tab import (
    probe_cohort_roh_data,
    render_cohort_roh_tab,
)
from genetics_viz.pages.cohort.components.cohort_sample_qc_tab import (
    probe_cohort_qc_data,
    render_cohort_sample_qc_tab,
)
from genetics_viz.pages.cohort.components.cohort_denovo_tab import (
    probe_cohort_denovo,
    render_cohort_denovo_tab,
)
from genetics_viz.pages.cohort.components.cohort_str_tab import (
    probe_cohort_eh_data,
    render_cohort_str_tab,
)
from genetics_viz.pages.cohort.components.cohort_upd_tab import (
    UPD_ANY,
    cohort_standing,
    probe_cohort_upd_data,
    render_cohort_upd_tab,
    upd_column_fields,
    upd_filter_options,
)
from genetics_viz.pages.cohort.components.docs_panel import render_params_doc_panel
from genetics_viz.pages.cohort.components.params_panel import render_params_panel
from genetics_viz.pages.cohort.components.stats_panel import render_stats_panel
from genetics_viz.pages.cohort.components.status_panel import (
    render_status_panel,
    status_headline,
)
from genetics_viz.utils.auth import can_write, check_auth
from genetics_viz.utils.cohort_state import read_cohort_state
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.diagnostic_files import diagnostic_paths
from genetics_viz.utils.ghfc_ngs_docs import fetch_params_doc
from genetics_viz.utils.pedigree_labels import dummy_badges
from genetics_viz.utils.review_counts import (
    filter_pills,
    qc_cohort_pills,
    qc_worst_by_individual,
    denovo_cohort_pills,
    denovo_cohort_status,
    str_cohort_open_flags,
    str_cohort_pills,
    upd_cohort_pills,
    workflow_pills,
)
from genetics_viz.utils.genome_qc import cohort_flags as genome_flags
from genetics_viz.utils.sample_qc import (
    CATEGORIES,
    category_badges,
    flags_by_key,
    load_cohort_qc,
)
from genetics_viz.utils.sample_qc_acks import (
    acknowledged_titles,
    load_acks,
    muted_color,
    split_acknowledged,
)
from genetics_viz.utils.upd_calls import load_cohort_upd

logger = logging.getLogger(__name__)

#: QC filter options besides the categories.
_QC_ANY = "Any issue"
_QC_NONE = "No issue"
_QC_UNCHECKED = "Not checked"
_QC_ACKED = "Acknowledged"

# Diagnostic priority: higher value wins
_DIAG_PRIORITY = {"pathogenic": 3, "uncertain": 2, "benign": 1}


def _loading(message: str) -> None:
    """Placeholder shown by a lazy tab panel until its content is loaded."""
    with ui.column().classes("w-full items-center justify-center py-16"):
        ui.spinner(size="xl", color="blue")
        ui.label(message).classes("text-lg text-gray-600 mt-4")


def _render_workflow_tab(cohort: Any) -> None:
    """Pipeline status, then the cohort's parameters beside their reference.

    The status panel is open on arrival: it is the reason to come to this tab.
    The two panels below it each scroll internally rather than letting the page
    grow — the reference alone is 563 lines — and carry ``min-w-0`` so a wide
    table inside one cannot drag the row past the page.
    """
    with ui.column().classes("w-full min-w-0 gap-4 p-4"):
        headline = status_headline(cohort)
        with ui.expansion(value=True).classes(
            "w-full min-w-0 border border-gray-200 rounded-md"
        ) as status_expansion:
            with status_expansion.add_slot("header"):
                with ui.row().classes("items-center gap-2 no-wrap w-full min-w-0"):
                    ui.label("Pipeline status").classes(
                        "text-sm font-semibold text-gray-700 whitespace-nowrap"
                    )
                    ui.badge(
                        headline.get("label", ""), color=headline.get("color", "grey")
                    ).classes("text-[10px]")
                    ui.label(headline.get("caption", "")).classes(
                        "text-xs text-gray-500 truncate min-w-0"
                    )
            # The expansion header above already carries the title, the
            # status chip and the step summary.
            render_status_panel(cohort, title=None)

        # flex-wrap so the two panels stack rather than crushing each other on a
        # narrow window; min-w-0 so each can actually shrink to its share.
        with ui.row().classes("w-full min-w-0 gap-4 flex-wrap items-start"):
            with ui.column().classes("flex-1 min-w-0 basis-[24rem]"):
                render_params_panel(cohort)
            with ui.column().classes("flex-1 min-w-0 basis-[24rem]"):
                render_params_doc_panel()


#: Link key of each individuals-table column filter, and its empty value.
_FILTER_LINK = {
    "FID": ("fid", ""),
    "Sample ID": ("sid", ""),
    "Sex": ("sex", []),
    "Phenotype": ("pheno", []),
    "Diagnostic": ("diag", []),
    "QC": ("qc", []),
    "UPD": ("upd", []),
    "Father": ("father", False),
    "Mother": ("mother", False),
}

#: The page's tabs, in order: label, link value (``tab=``), icon, tooltip.
_TABS: List[Tuple[str, str, str, str]] = [
    ("Individuals", "individuals", "groups", ""),
    ("Statistics", "stats", "bar_chart", "Carriers per geneset"),
    ("Ancestry", "ancestry", "public", ""),
    ("Polygenic Scores", "pgs", "functions", ""),
    ("ROH", "roh", "linear_scale", "Runs of homozygosity"),
    ("UPD", "upd", "content_copy", "Uniparental disomy"),
    (
        "De novo",
        "denovo",
        "new_releases",
        "De novo SNVs and indels: per-child QC, flags, literature (denovo step)",
    ),
    ("STR", "str", "repeat", "Short tandem repeats"),
    (
        "Sample QC",
        "qc",
        "fact_check",
        "Relatedness, sex and contamination (somalier); DNA source and library",
    ),
    ("Workflow", "workflow", "account_tree", "Pipeline run status and parameters"),
]
_SLUGS = {label: slug for label, slug, _icon, _tip in _TABS}
_DEFAULT_TAB = "individuals"

#: Values a filter chip lists before ``+N``.
_CHIP_VALUES = 3

#: Keeps ``--gv-header-height`` at the fixed page header's height, which
#: grows when its menus wrap: the sticky tab row sits right under it.
_HEADER_HEIGHT_JS = """
(() => {
    const header = document.querySelector('.q-header');
    if (!header) return;
    const set = () => document.documentElement.style.setProperty(
        '--gv-header-height', header.offsetHeight + 'px');
    set();
    new ResizeObserver(set).observe(header);
})();
"""


#: Ten tabs with icons and badges fit a 1280 px window only with Quasar's
#: 16 px tab padding cut to 10 (and mixed-case labels): about 1,150 px with
#: the De novo and Workflow pills, of the 1,185 px the row gets.
_TABS_CSS = ".cohort-tabs .q-tab { padding: 0 10px; }"


def _linked_tab(link: PageLink) -> Optional[str]:
    """The tab a link opens. Older links name it ``panel`` or ``view=workflow``."""
    if link.get("tab"):
        return link.get("tab")
    if link.get("view") == "workflow":
        return "workflow"
    return link.get("panel")


def _filter_active(value: Any) -> bool:
    return bool(value.strip()) if isinstance(value, str) else bool(value)


def _filter_chips(
    filters: Dict[str, Any], columns: Sequence[Dict[str, Any]]
) -> List[Tuple[str, str]]:
    """``(column id, chip text)`` of each filter in use, in column order."""
    chips: List[Tuple[str, str]] = []
    for col in columns:
        spec = col.get("filter")
        value = filters.get(col["id"])
        if not spec or not _filter_active(value):
            continue
        if spec.get("type") == "checkbox":
            text = spec.get("label") or col["header"]
        elif isinstance(value, list):
            text = f"{col['header']}: " + ", ".join(
                str(v) for v in value[:_CHIP_VALUES]
            )
            if len(value) > _CHIP_VALUES:
                text += f" +{len(value) - _CHIP_VALUES}"
        else:
            text = f"{col['header']}: {str(value).strip()}"
        chips.append((col["id"], text))
    return chips


def _count_text(shown: int, total: int, filtered: bool) -> str:
    noun = "individual" if total == 1 else "individuals"
    if not filtered:
        return f"{total:,} {noun}"
    return f"{shown:,} of {total:,} {noun}"


def _filter_individuals(
    individuals: List[Dict[str, Any]], filters: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """The individuals the table's column filters leave."""
    filtered = individuals

    # Text filters: a case-insensitive substring
    for col in ("FID", "Sample ID"):
        text = (filters.get(col) or "").strip().lower()
        if text:
            filtered = [ind for ind in filtered if text in ind[col].lower()]

    # Multiselect filters: any of the selected values
    for col in ("Sex", "Phenotype"):
        wanted = set(filters.get(col) or [])
        if wanted:
            filtered = [ind for ind in filtered if ind[col] in wanted]

    # Diagnostic, with "NA" for no diagnostic
    diag_vals = set(filters.get("Diagnostic") or [])
    if diag_vals:
        has_na = "NA" in diag_vals
        filtered = [
            ind
            for ind in filtered
            if ind["Diagnostic"] in diag_vals - {"NA"}
            or (has_na and not ind["Diagnostic"])
        ]

    # Sample QC (categories, any, none, unchecked)
    wanted_qc = set(filters.get("QC") or [])
    if wanted_qc:
        filtered = [
            ind
            for ind in filtered
            if wanted_qc & set(ind.get("_qc_cats", []))
            or (_QC_ANY in wanted_qc and ind.get("_qc_rank", -1) > 0)
        ]

    # UPD (confidences, review, any, none)
    wanted_upd = set(filters.get("UPD") or [])
    if wanted_upd:
        filtered = [
            ind
            for ind in filtered
            if wanted_upd & set(ind.get("_upd_cats", []))
            or (UPD_ANY in wanted_upd and ind.get("_upd_rank", -1) > 0)
        ]

    # Has father / mother: a sequenced one, so not a dummy
    for col in ("Father", "Mother"):
        if filters.get(col):
            filtered = [ind for ind in filtered if is_sequenced_parent(ind[col])]
    return filtered


def _individual_columns(
    cohort_name: str,
    sex_values: List[str],
    phenotype_values: List[str],
    diag_options: List[str],
    *,
    with_qc: bool,
    with_upd: bool,
) -> List[Dict[str, Any]]:
    """The individuals table's columns; Sample QC and UPD when the cohort ran them."""
    columns: List[Dict[str, Any]] = [
        {
            "id": "FID",
            "header": "Family ID",
            "cellType": "link",
            "href": f"/cohort/{cohort_name}/family/{{FID}}",
            "sortable": True,
            "minWidth": 250,
            "filter": {"type": "text", "placeholder": "Filter..."},
        },
        {
            "id": "Sample ID",
            "header": "Sample ID",
            "sortable": True,
            "minWidth": 90,
            "filter": {"type": "text", "placeholder": "Filter..."},
        },
        {
            "id": "Sex",
            "header": "Sex",
            "sortable": True,
            "filter": {
                "type": "multiselect",
                "options": sex_values,
                "placeholder": "All",
            },
        },
        {
            "id": "Phenotype",
            "header": "Phenotype",
            "sortable": True,
            "filter": {
                "type": "multiselect",
                "options": phenotype_values,
                "placeholder": "All",
            },
        },
        {
            "id": "Father",
            "header": "Father ID",
            "cellType": "badge_list",
            "badgesField": "_father_badges",
            "sortable": True,
            "minWidth": 90,
            "filter": {"type": "checkbox", "label": "Has sequenced father"},
        },
        {
            "id": "Mother",
            "header": "Mother ID",
            "cellType": "badge_list",
            "badgesField": "_mother_badges",
            "sortable": True,
            "minWidth": 90,
            "filter": {"type": "checkbox", "label": "Has sequenced mother"},
        },
        {
            "id": "Diagnostic",
            "header": "Diagnostic",
            "sortable": True,
            "cellType": "badge",
            "colorField": "_diag_color",
            "filter": {
                "type": "multiselect",
                "options": diag_options,
                "placeholder": "All",
            },
        },
    ]
    if with_qc:
        columns.append(
            {
                "id": "QC",
                "header": "Sample QC",
                "cellType": "badge_list",
                "sortable": True,
                "sortField": "_qc_rank",
                "sorting": "numerical",
                "filter": {
                    "type": "multiselect",
                    "options": [_QC_ANY]
                    + [c.get("label", k) for k, c in CATEGORIES.items()]
                    + [_QC_NONE, _QC_ACKED, _QC_UNCHECKED],
                    "placeholder": "All",
                },
            }
        )
    if with_upd:
        columns.append(
            {
                "id": "UPD",
                "header": "UPD",
                "cellType": "badge_list",
                "sortable": True,
                "sortField": "_upd_rank",
                "sorting": "numerical",
                "filter": {
                    "type": "multiselect",
                    "options": upd_filter_options(),
                    "placeholder": "All",
                },
            }
        )
    return columns


def _filters_from_link(link: PageLink) -> Dict[str, Any]:
    """The table's column filter values a shared link sets."""
    filters: Dict[str, Any] = {}
    for col, (key, empty) in _FILTER_LINK.items():
        if key not in link.params:
            continue
        if isinstance(empty, list):
            filters[col] = link.get_list(key) or []
        elif isinstance(empty, bool):
            filters[col] = link.get_bool(key, False)
        else:
            filters[col] = link.get(key) or ""
    return filters


def _filters_to_link(filters: Dict[str, Any], link: PageLink) -> None:
    for col, (key, empty) in _FILTER_LINK.items():
        value = filters.get(col)
        link.set(key, value if value else None, default=empty)


@ui.page("/cohort/{cohort_name}")
def cohort_page(cohort_name: str, request: Request) -> None:
    """Render the cohort detail page.

    The open tabs and the individuals table's filters are in the URL query
    (``components/page_link.py``), so the address is a link to this view.
    """
    if redirect := check_auth():
        return redirect
    link = PageLink.from_request(request)
    create_header(cohort_name)

    try:
        store = get_data_store()
        cohort = store.get_cohort(cohort_name)

        if cohort is None:
            # The directory may well exist and simply have no usable pedigree,
            # in which case saying "not found" would send someone looking for
            # the wrong thing.
            stub = next(
                (s for s in store.incomplete_cohorts if s.name == cohort_name), None
            )
            with ui.column().classes("w-full max-w-6xl mx-auto p-6 gap-2"):
                if stub is None:
                    ui.label(f"Cohort not found: {cohort_name}").classes(
                        "text-xl text-red-500"
                    )
                else:
                    ui.label(cohort_name).classes("text-xl font-bold text-amber-900")
                    ui.label(
                        "Unreadable pedigree file"
                        if stub.reason == "unreadable"
                        else "Missing pedigree file"
                    ).classes("text-amber-800 font-semibold")
                    ui.label(
                        "This cohort cannot be explored until"
                        f" {stub.expected_pedigree.name} is in place."
                    ).classes("text-sm text-gray-600")
                    if stub.error:
                        ui.label(stub.error).classes(
                            "text-sm text-amber-900 italic font-mono"
                        )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        # Load diagnostics per sample (highest priority wins)
        sample_diag: Dict[str, str] = defaultdict(str)
        for _dtype, diag_file in diagnostic_paths(store.data_dir):
            if not diag_file.exists():
                continue
            with open(diag_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
                for row in reader:
                    if row.get("Ignore", "0") == "1":
                        continue
                    sid = row.get("Sample", "")
                    diag = row.get("Diagnostic", "")
                    if not sid or not diag:
                        continue
                    current = sample_diag[sid]
                    if _DIAG_PRIORITY.get(diag, 0) > _DIAG_PRIORITY.get(current, 0):
                        sample_diag[sid] = diag

        # Sample QC per (family, sample), when the cohort's step ran. The
        # tables are tens of KB and cached on their stat.
        try:
            cohort_qc = load_cohort_qc(cohort)
        except Exception:
            cohort_qc = None
        # Flags a curator acknowledged leave the badges and "Any issue"; they
        # get one grey badge, which the filter's "Acknowledged" option finds.
        qc_acks = load_acks(store.data_dir)
        # The genome QC flags (cohort level only) count with them.
        qc_active, qc_acked = (
            split_acknowledged(list(cohort_qc.flags) + genome_flags(cohort), qc_acks)
            if cohort_qc
            else ([], [])
        )
        qc_by_key = flags_by_key(qc_active)
        qc_acked_by_key = flags_by_key(qc_acked)
        qc_checked = (
            {
                ((r.get("family_id") or ""), (r.get("sample") or ""))
                for r in cohort_qc.sample_rows()
            }
            if cohort_qc
            else set()
        )
        show_qc_detail = can_write()

        def qc_fields(fid: str, sid: str) -> Dict[str, Any]:
            if cohort_qc is None:
                return {}
            flags = qc_by_key.get((fid, sid), [])
            acked = qc_acked_by_key.get((fid, sid), [])
            if (fid, sid) not in qc_checked and not flags and not acked:
                return {"QC": "", "_qc_cats": [_QC_UNCHECKED], "_qc_rank": -1}
            badges = category_badges(flags, show_qc_detail)
            if acked:
                badges.append(
                    {
                        "label": _QC_ACKED,
                        "color": muted_color(),
                        "tooltip": "\n".join(
                            [
                                "Acknowledged by a curator:",
                                *(
                                    f"• {t}"
                                    for t in acknowledged_titles(
                                        acked, qc_acks, show_qc_detail
                                    )
                                ),
                            ]
                        ),
                    }
                )
            return {
                "QC": "",
                "QC_badges": badges,
                "_qc_cats": [b["label"] for b in badges]
                + ([] if flags else [_QC_NONE]),
                "_qc_rank": max((f.rank for f in flags), default=0),
            }

        # The upd step's events per (family, sample), when the cohort ran it:
        # a few hundred KB, cached on their stat.
        try:
            cohort_upd = load_cohort_upd(cohort)
            upd_fields = upd_column_fields(
                cohort_upd, cohort_standing(store.data_dir) if cohort_upd else {}
            )
        except Exception:
            cohort_upd = None
            upd_fields = upd_column_fields(None, {})

        # Build flat list of all individuals from pedigree
        all_individuals: List[Dict[str, Any]] = []
        for family in cohort.families.values():
            for sample in family.samples:
                all_individuals.append(
                    {
                        "FID": sample.family_id,
                        "Sample ID": sample.sample_id,
                        "Sex": sample.sex or "-",
                        "Phenotype": sample.phenotype or "-",
                        "Father": sample.father_id or "-",
                        "Mother": sample.mother_id or "-",
                        # A dummy parent (never sequenced) shows as a muted badge.
                        "_father_badges": dummy_badges(sample.father_id),
                        "_mother_badges": dummy_badges(sample.mother_id),
                        "Diagnostic": sample_diag.get(sample.sample_id, ""),
                        "_diag_color": {
                            "pathogenic": "#dc2626",
                            "uncertain": "#f59e0b",
                            "benign": "#16a34a",
                            "conflicting": "#f59e0b",
                        }.get(sample_diag.get(sample.sample_id, ""), ""),
                        **qc_fields(sample.family_id, sample.sample_id),
                        **upd_fields(sample.family_id, sample.sample_id),
                    }
                )

        # Collect unique values for multiselect filters
        phenotype_values = sorted(
            {ind["Phenotype"] for ind in all_individuals if ind["Phenotype"] != "-"}
        )
        sex_values = sorted(
            {ind["Sex"] for ind in all_individuals if ind["Sex"] != "-"}
        )
        # Diagnostic filter: include actual values + "NA" for undiagnosed
        diag_values_raw = sorted(
            {ind["Diagnostic"] for ind in all_individuals if ind["Diagnostic"]}
        )
        diag_filter_options = diag_values_raw + ["NA"]

        # Which analyses the cohort has: a tab without data is disabled, and
        # Workflow is absent, not disabled, with no parameters file (there is
        # no pipeline to describe).
        has_data: Dict[str, bool] = {
            "Individuals": True,
            "Statistics": True,
            "Ancestry": probe_cohort_ancestry_data(cohort),
            "Polygenic Scores": probe_cohort_pgs_data(cohort),
            "ROH": probe_cohort_roh_data(cohort),
            "UPD": probe_cohort_upd_data(cohort),
            "STR": probe_cohort_eh_data(cohort),
            "De novo": probe_cohort_denovo(cohort),
            "Sample QC": probe_cohort_qc_data(cohort),
            "Workflow": cohort.params_file is not None,
        }

        columns = _individual_columns(
            cohort_name,
            sex_values,
            phenotype_values,
            diag_filter_options,
            with_qc=cohort_qc is not None,
            with_upd=cohort_upd is not None,
        )
        n_total = len(all_individuals)
        # Each individual's worst unacknowledged QC flag, for the tab's badge.
        qc_worst = qc_worst_by_individual(qc_active, has_data["UPD"])
        show_review = can_write()

        linked_filters = _filters_from_link(link)
        linked_slug = _linked_tab(link)
        initial = next(
            (
                label
                for label, slug, _icon, _tip in _TABS
                if slug == linked_slug and has_data[label]
            ),
            "Individuals",
        )
        if "panel" in link.params or "view" in link.params:
            # An older link: its tab is now ``tab``.
            link.drop("panel", "view")
            link.set("tab", _SLUGS[initial], default=_DEFAULT_TAB)
            link.push()

        # Shared filtered state for the tabs that follow the table filters
        filtered_state: Dict[str, Any] = {"individuals": list(all_individuals)}
        # Refresh callbacks of the tabs that follow the table filters. Only
        # tabs that have actually been loaded register here, so a filter
        # change never forces work for a tab never opened.
        refreshers: List[Callable[[], None]] = []
        state: Dict[str, Any] = {
            "filters": {},
            "tab": initial,
            "str_flags": [],
            "denovo": ({}, []),
        }
        dt_ref: Dict[str, Any] = {"dt": None}

        ui.run_javascript(_HEADER_HEIGHT_JS)
        ui.add_css(_TABS_CSS)

        with ui.column().classes("w-full px-6 pt-4 pb-6 gap-0"):
            # Cohort header
            with ui.row().classes("items-center gap-4 mb-3"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")

            def apply_filters(filters: Dict[str, Any], to_table: bool = False) -> None:
                """Show the individuals ``filters`` leave, on every tab.

                ``to_table`` when the change comes from the filter bar: the
                table's filter row is set to match.
                """
                state["filters"] = dict(filters)
                _filters_to_link(filters, link)
                link.push()
                filtered = _filter_individuals(all_individuals, filters)
                filtered_state["individuals"] = filtered
                dt = dt_ref["dt"]
                if dt is not None:
                    if to_table:
                        dt.set_filters(filters)
                    dt.update_data(filtered)
                filter_bar.refresh()
                update_pills()
                # Re-render the tabs that follow the filters
                for refresher in refreshers:
                    refresher()

            def on_filter(e: Dict[str, Any]) -> None:
                apply_filters(e.get("filters") or {})

            def make_clear(col: str) -> Callable[[Any], None]:
                def handler(e: Any) -> None:
                    if e.value:
                        return
                    filters = dict(state["filters"])
                    empty = _FILTER_LINK[col][1]
                    filters[col] = list(empty) if isinstance(empty, list) else empty
                    apply_filters(filters, to_table=True)

                return handler

            def show_tab(label: str) -> None:
                tabs.set_value(label)

            @ui.refreshable
            def filter_bar() -> None:
                """The count, a removable chip per filter, and a way to the filters."""
                chips = _filter_chips(state["filters"], columns)
                shown = len(filtered_state["individuals"])
                with (
                    ui.row()
                    .classes(
                        "w-full min-w-0 items-center gap-2 px-3 py-1"
                        " bg-blue-50 border-b border-blue-100"
                    )
                    .mark("cohort-filter-bar")
                ):
                    ui.icon("filter_alt", size="18px").classes("text-blue-700")
                    ui.label(_count_text(shown, n_total, bool(chips))).classes(
                        "text-sm font-semibold text-blue-900 whitespace-nowrap"
                    ).mark("cohort-count")
                    for col, text in chips:
                        ui.chip(
                            text,
                            removable=True,
                            color="white",
                            text_color="blue-10",
                            on_value_change=make_clear(col),
                        ).props("dense").classes("text-xs border border-blue-200").mark(
                            "cohort-filter-chip"
                        )
                    if chips:
                        ui.button(
                            "Clear all", on_click=lambda: apply_filters({}, True)
                        ).props("flat dense no-caps size=sm color=blue-8").mark(
                            "cohort-clear-filters"
                        )
                    ui.space()
                    if state["tab"] != "Individuals":
                        ui.button(
                            "Edit filters" if chips else "Filter individuals",
                            icon="tune",
                            on_click=lambda: show_tab("Individuals"),
                        ).props("flat dense no-caps size=sm color=blue-8").mark(
                            "cohort-edit-filters"
                        )
                    elif not chips:
                        ui.label(
                            "Filter with the boxes under the column headers:"
                            " every tab but Workflow follows."
                        ).classes("text-xs text-gray-500")

            # The tab row and the filter bar stay in view under the page
            # header, so a long tab never loses sight of its filter.
            with (
                ui.column()
                .classes("w-full gap-0 sticky z-10 bg-white")
                .style("top: var(--gv-header-height, 0px)")
            ):
                with (
                    ui.tabs(value=initial)
                    .props("inline-label no-caps align=left")
                    .classes("w-full border-b border-gray-200 cohort-tabs")
                    .mark("cohort-tabs") as tabs
                ):
                    tab_objects: Dict[str, Any] = {}
                    for label, _slug, icon, tip in _TABS:
                        if label == "Workflow" and not has_data["Workflow"]:
                            continue
                        tab = ui.tab(label, icon=icon)
                        if tip:
                            tab.tooltip(tip)
                        if not has_data[label]:
                            tab.props("disable")
                        tab_objects[label] = tab

                # Badges: the filter for everyone, as Sample QC and the run;
                # UPD and STR count what only a curator settles.
                pill_tabs = ["Individuals", "Sample QC", "De novo", "Workflow"]
                if show_review:
                    pill_tabs += ["UPD", "STR"]
                pills: Dict[str, TabPills] = {
                    label: TabPills(tab_objects[label], f"{_SLUGS[label]}-tab-pills")
                    for label in pill_tabs
                    if label in tab_objects and has_data[label]
                }

                with (
                    ui.element("div").classes("w-full").mark("cohort-filter-box")
                ) as filter_bar_box:
                    filter_bar()
                filter_bar_box.set_visibility(initial != "Workflow")

            def update_pills() -> None:
                rows = filtered_state["individuals"]
                keys = {(r["FID"], r["Sample ID"]) for r in rows}
                n_filters = len(_filter_chips(state["filters"], columns))
                pills["Individuals"].set(filter_pills(len(rows), n_total, n_filters))
                if "Sample QC" in pills:
                    pills["Sample QC"].set(
                        qc_cohort_pills(qc_worst, keys, has_data["UPD"])
                    )
                if "UPD" in pills:
                    pills["UPD"].set(
                        upd_cohort_pills(
                            (r.get("_upd_open") or [] for r in rows), len(keys)
                        )
                    )
                if "STR" in pills:
                    pills["STR"].set(str_cohort_pills(state["str_flags"], keys))
                if "De novo" in pills:
                    verdicts, issues = state["denovo"]
                    pills["De novo"].set(denovo_cohort_pills(verdicts, issues, keys))

            async def load_badges() -> None:
                """The STR flags and the run record, read once off the loop."""
                if "STR" in pills:
                    try:
                        state["str_flags"] = await asyncio.to_thread(
                            str_cohort_open_flags, cohort, store.data_dir
                        )
                    except Exception:
                        logger.warning("Failed to count STR flags", exc_info=True)
                if "De novo" in pills:
                    try:
                        state["denovo"] = await asyncio.to_thread(
                            denovo_cohort_status, cohort, store.data_dir, show_review
                        )
                    except Exception:
                        logger.warning("Failed to read the de novo QC", exc_info=True)
                if "Workflow" in pills:
                    try:
                        run = await asyncio.to_thread(read_cohort_state, cohort.path)
                        headline = await asyncio.to_thread(status_headline, cohort)
                        pills["Workflow"].set(
                            workflow_pills(run, headline.get("caption", ""))
                        )
                    except Exception:
                        logger.warning("Failed to read the run record", exc_info=True)
                update_pills()

            # Each tab reads its files on first view, not on first paint.
            openers: Dict[str, Callable[[], None]] = {}

            def lazy_content(
                label: str, render: Callable[..., None], missing: str, loading: str
            ) -> None:
                loaded = {"value": False}
                available = has_data[label]

                @ui.refreshable
                def content() -> None:
                    if not available:
                        ui.label(missing).classes("text-gray-500 italic")
                    elif loaded["value"]:
                        render(store, cohort, filtered_state, refreshers)
                    else:
                        _loading(loading)

                content()

                def opener() -> None:
                    if not available or loaded["value"]:
                        return

                    def load() -> None:
                        loaded["value"] = True
                        content.refresh()

                    ui.timer(0.1, load, once=True)

                openers[label] = opener

            # Reading the run record hashes two files and fetching the
            # reference is a network call, so the tab waits until it is opened.
            workflow_state = {"loaded": False}

            @ui.refreshable
            def workflow_content() -> None:
                if not workflow_state["loaded"]:
                    _loading("Loading workflow...")
                    return
                _render_workflow_tab(cohort)

            def open_workflow() -> None:
                if workflow_state["loaded"]:
                    return

                async def load() -> None:
                    # Warm both caches off-thread, so the synchronous renderers
                    # only read from cache: one hashes the pedigree and params
                    # file, the other fetches 31 KB over the network, and this
                    # event loop is shared by every session.
                    await asyncio.to_thread(read_cohort_state, cohort.path)
                    await asyncio.to_thread(fetch_params_doc)
                    workflow_state["loaded"] = True
                    workflow_content.refresh()

                ui.timer(0.1, load, once=True)

            openers["Workflow"] = open_workflow

            # Each renderer opens its own card, so the panels carry no border
            # of their own and the chrome is not doubled.
            panel_classes = "w-full px-0 pt-3 pb-0"
            with ui.tab_panels(tabs, value=initial).classes("w-full"):
                with ui.tab_panel(tab_objects["Individuals"]).classes(panel_classes):
                    dt_ref["dt"] = DataTable(
                        columns=columns,
                        rows=all_individuals,
                        row_key="Sample ID",
                        pagination={"rowsPerPage": 20},
                        on_filter=on_filter,
                        initial_filters=linked_filters,
                    )
                # A shared link's filters apply before the tabs are built, so
                # they open on the filtered individuals.
                if linked_filters:
                    apply_filters(linked_filters)

                with ui.tab_panel(tab_objects["Statistics"]).classes(panel_classes):
                    render_stats_panel(store, cohort, filtered_state)

                for label, render, what, loading in (
                    ("Ancestry", render_cohort_ancestry_tab, "ancestry", "ancestry"),
                    (
                        "Polygenic Scores",
                        render_cohort_pgs_tab,
                        "polygenic score",
                        "polygenic scores",
                    ),
                    ("ROH", render_cohort_roh_tab, "ROH", "runs of homozygosity"),
                    ("UPD", render_cohort_upd_tab, "UPD", "uniparental disomy"),
                    ("De novo", render_cohort_denovo_tab, "de novo", "de novo calls"),
                    ("STR", render_cohort_str_tab, "STR", "short tandem repeats"),
                    (
                        "Sample QC",
                        render_cohort_sample_qc_tab,
                        "sample QC",
                        "sample QC",
                    ),
                ):
                    with ui.tab_panel(tab_objects[label]).classes(panel_classes):
                        lazy_content(
                            label,
                            render,
                            f"No cohort-level {what} results",
                            f"Loading {loading}...",
                        )

                if "Workflow" in tab_objects:
                    with ui.tab_panel(tab_objects["Workflow"]).classes("w-full p-0"):
                        workflow_content()

            def on_tab_change(e: Any) -> None:
                label = (
                    e.value.props["name"] if isinstance(e.value, ui.tab) else e.value
                )
                state["tab"] = label
                link.set("tab", _SLUGS.get(label), default=_DEFAULT_TAB)
                link.push()
                filter_bar_box.set_visibility(label != "Workflow")
                filter_bar.refresh()
                opener = openers.get(label)
                if opener is not None:
                    opener()

            tabs.on_value_change(on_tab_change)
            if initial in openers:
                openers[initial]()
            update_pills()
            ui.timer(0.2, load_badges, once=True)

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
