"""Cohort detail page - displays pedigree members with filters."""

import asyncio
import csv
import logging
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from fastapi import Request
from nicegui import app, ui

from genetics_viz.components.header import create_header
from genetics_viz.components.tab_groups import GROUPED_TABS_PROPS, tab_group
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
from genetics_viz.pages.cohort.components.status_panel import (
    render_status_panel,
    status_headline,
)
from genetics_viz.utils.auth import can_write, check_auth
from genetics_viz.utils.cohort_state import read_cohort_state
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.diagnostic_files import diagnostic_paths
from genetics_viz.utils.family_rows import (
    DIAG_ORDER,
    FamilyFacts,
    compact_qc_cell,
    family_columns,
    family_row,
    filter_families,
    load_family_facts,
    save_family_order,
)
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

#: Link key of each families-table column filter, and its empty value.
_FAMILY_FILTER_LINK = {
    "FID": ("ffid", ""),
    "Structure": ("fstruct", []),
    "Trios": ("ftrio", False),
    "Affected": ("fmulti", False),
    "Ancestry": ("fanc", ""),
    "FROH": ("froh", []),
    "QC": ("fqc", []),
    "DenovoQC": ("fdn", []),
    "Review": ("frev", False),
    "Diagnostics": ("fdiag", []),
    "Conclusion": ("fconc", False),
}

#: The page's tab row by group: label, link value (``tab=``), tooltip. The
#: groups follow a curator's work: browse the cohort, run the checks, then
#: the views by variant type and the genetic background.
_TAB_GROUPS: List[Tuple[str, List[Tuple[str, str, str]]]] = [
    (
        "Browse",
        [
            ("Families", "families", "One row per family: structure, QC, review"),
            ("Individuals", "individuals", ""),
        ],
    ),
    (
        "Checks",
        [
            (
                "Sample QC",
                "qc",
                "Relatedness, sex and contamination (somalier); DNA source and library",
            ),
            (
                "De novo",
                "denovo",
                "De novo SNVs and indels: per-child QC, flags, literature"
                " (denovo step)",
            ),
        ],
    ),
    (
        "Variants",
        [
            ("STR", "str", "Short tandem repeats"),
            ("UPD", "upd", "Uniparental disomy"),
        ],
    ),
    (
        "Background",
        [
            ("Ancestry", "ancestry", ""),
            ("ROH", "roh", "Runs of homozygosity"),
            ("Polygenic Scores", "pgs", ""),
        ],
    ),
]
#: Workflow is about the run, not the cohort: a button on the title line, not
#: a tab, though its panel and its ``tab=workflow`` link work as a tab's.
_WORKFLOW = ("Workflow", "workflow", "Pipeline run status and parameters")
#: Every tab, the row's then Workflow: label, link value, tooltip.
_TABS: List[Tuple[str, str, str]] = [
    tab for _caption, tabs in _TAB_GROUPS for tab in tabs
] + [_WORKFLOW]
_SLUGS = {label: slug for label, slug, _tip in _TABS}
_DEFAULT_TAB = "families"
#: The tabs holding a filterable table.
_TABLE_TABS = ("Families", "Individuals")

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


def _count_text(
    shown: int,
    total: int,
    filtered: bool,
    one: str = "individual",
    many: str = "individuals",
) -> str:
    noun = one if total == 1 else many
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


def _filters_from_link(
    link: PageLink, spec: Dict[str, Tuple[str, Any]] = _FILTER_LINK
) -> Dict[str, Any]:
    """A table's column filter values a shared link sets (``spec``: link keys)."""
    filters: Dict[str, Any] = {}
    for col, (key, empty) in spec.items():
        if key not in link.params:
            continue
        if isinstance(empty, list):
            filters[col] = link.get_list(key) or []
        elif isinstance(empty, bool):
            filters[col] = link.get_bool(key, False)
        else:
            filters[col] = link.get(key) or ""
    return filters


def _filters_to_link(
    filters: Dict[str, Any],
    link: PageLink,
    spec: Dict[str, Tuple[str, Any]] = _FILTER_LINK,
) -> None:
    for col, (key, empty) in spec.items():
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
        # Per family, each diagnosed variant (kind, key) and its best verdict.
        family_diag: Dict[str, Dict[Tuple[str, str], str]] = defaultdict(dict)
        for dtype, diag_file in diagnostic_paths(store.data_dir):
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
                    variants = family_diag[row.get("FID") or ""]
                    key = (dtype, row.get("Variant") or "")
                    if _DIAG_PRIORITY.get(diag, 0) > _DIAG_PRIORITY.get(
                        variants.get(key, ""), 0
                    ):
                        variants[key] = diag

        def diagnostics_of(fid: str) -> Tuple[int, str]:
            """A family's diagnosed variants and the strongest verdict."""
            verdicts = list(family_diag.get(fid, {}).values())
            strongest = next((v for v in DIAG_ORDER if v in verdicts), "")
            return len(verdicts), strongest

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

        # The same flags by family (a cross-family pair counts for both).
        family_flags: Dict[str, List[Any]] = defaultdict(list)
        family_acked: Dict[str, List[Any]] = defaultdict(list)
        for flags, by_family in ((qc_active, family_flags), (qc_acked, family_acked)):
            for flag in flags:
                for fid in dict.fromkeys(k[0] for k in flag.keys):
                    by_family[fid].append(flag)
        qc_checked_families = {fid for fid, _sid in qc_checked}

        def qc_fields(fid: str, sid: str) -> Dict[str, Any]:
            if cohort_qc is None:
                return {}
            return qc_cell(
                qc_by_key.get((fid, sid), []),
                qc_acked_by_key.get((fid, sid), []),
                (fid, sid) in qc_checked,
            )

        def family_qc_fields(fid: str) -> Dict[str, Any]:
            if cohort_qc is None:
                return {}
            flags, acked = family_flags.get(fid, []), family_acked.get(fid, [])
            return compact_qc_cell(
                qc_cell(flags, acked, fid in qc_checked_families),
                flags,
                acked,
                fid,
                _QC_ACKED,
            )

        def qc_cell(
            flags: List[Any], acked: List[Any], checked: bool
        ) -> Dict[str, Any]:
            """The Sample QC cell of an individual or a family."""
            if not checked and not flags and not acked:
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
            "Families": True,
            "Individuals": True,
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

        # One row per family. The pedigree, Sample QC, UPD and diagnostics
        # columns are filled now; ancestry, FROH and de novo come from the
        # cohort tables, read off the loop once the page is drawn, and the
        # STR flags with the tab badges.
        qc_options = next(
            (c["filter"]["options"] for c in columns if c["id"] == "QC"), []
        )
        family_cols = family_columns(
            cohort_name,
            qc_options=qc_options,
            diag_options=[d for d in DIAG_ORDER if d in diag_values_raw],
            with_ancestry=has_data["Ancestry"],
            with_roh=has_data["ROH"],
            with_denovo=has_data["De novo"],
            curator=show_review,
        )
        upd_open_by_family: Dict[str, List[str]] = defaultdict(list)
        for ind in all_individuals:
            upd_open_by_family[ind["FID"]].extend(ind.get("_upd_open") or [])
        family_data: Dict[str, Any] = {"facts": FamilyFacts(), "rows": []}

        def build_family_rows() -> None:
            str_open: Dict[str, List[int]] = defaultdict(lambda: [0, 0])
            for (fid, _sid), _locus, severity in state["str_flags"]:
                str_open[fid][0] += 1
                str_open[fid][1] += severity == "high"
            family_data["rows"] = [
                family_row(
                    cohort.families[fid],
                    family_data["facts"],
                    qc=family_qc_fields(fid),
                    upd_open=upd_open_by_family.get(fid, []),
                    str_open=tuple(str_open.get(fid, (0, 0))),
                    diagnostics=diagnostics_of(fid),
                    curator=show_review,
                    show_relationship=show_review,
                )
                for fid in sorted(cohort.families)
            ]

        n_families = len(cohort.families)

        linked_filters = _filters_from_link(link)
        linked_family_filters = _filters_from_link(link, _FAMILY_FILTER_LINK)
        linked_slug = _linked_tab(link)
        initial = next(
            (
                label
                for label, slug, _tip in _TABS
                if slug == linked_slug and has_data[label]
            ),
            "Families",
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
            "family_filters": {},
            "families": [],
            # The table "Edit filters" goes back to: the last one filtered.
            "filter_tab": "Families",
            "tab": initial,
            "str_flags": [],
            "denovo": ({}, []),
        }
        build_family_rows()
        state["families"] = list(family_data["rows"])
        dt_ref: Dict[str, Any] = {"dt": None, "families": None}

        ui.run_javascript(_HEADER_HEIGHT_JS)

        with ui.column().classes("w-full px-6 pt-4 pb-6 gap-0"):
            # Cohort header
            with ui.row().classes("w-full items-center gap-4 mb-3"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")
                ui.space()
                workflow_button: Optional[ui.button] = None
                if has_data["Workflow"]:
                    workflow_button = (
                        ui.button("Workflow", on_click=lambda: show_tab("Workflow"))
                        .props("outline no-caps color=blue-8")
                        .tooltip(_WORKFLOW[2])
                        .mark("cohort-workflow")
                    )

            def recompute(
                to_individuals: bool = False, to_families: bool = False
            ) -> None:
                """Apply both tables' filters, to both tables and every tab.

                The scope is the individuals the individuals filters leave, in
                the families the families filters leave. A family is listed when
                it passes its filters and, with individuals filters on, one of
                its members passes them. ``to_*`` when the change comes from the
                filter bar: that table's filter row is set to match.
                """
                ind_filters = state["filters"]
                fam_filters = state["family_filters"]
                by_individual = _filter_individuals(all_individuals, ind_filters)
                families = filter_families(family_data["rows"], fam_filters, _QC_ANY)
                if _filter_chips(fam_filters, family_cols):
                    kept = {r["FID"] for r in families}
                    individuals = [i for i in by_individual if i["FID"] in kept]
                else:
                    individuals = by_individual
                if _filter_chips(ind_filters, columns):
                    with_member = {i["FID"] for i in by_individual}
                    families = [r for r in families if r["FID"] in with_member]
                filtered_state["individuals"] = individuals
                state["families"] = families
                for key, rows, filters, push in (
                    ("dt", individuals, ind_filters, to_individuals),
                    ("families", families, fam_filters, to_families),
                ):
                    dt = dt_ref[key]
                    if dt is not None:
                        if push:
                            dt.set_filters(filters)
                        # The sort is Python's: keep the one the header shows.
                        dt.update_data(dt.sort_rows(list(rows)))
                remember_order()
                filter_bar.refresh()
                update_pills()
                # Re-render the tabs that follow the filters
                for refresher in refreshers:
                    refresher()

            def remember_order() -> None:
                """The family page's arrows walk the families shown, as sorted.

                Forgotten while the table is neither filtered nor sorted: the
                arrows then walk every family by ID, as they always did.
                """
                dt = dt_ref["families"]
                rows = dt.rows if dt is not None else state["families"]
                narrowed = _filter_chips(state["family_filters"], family_cols) or (
                    _filter_chips(state["filters"], columns)
                )
                sorted_ = dt is not None and bool(dt.sorting)
                save_family_order(
                    app.storage.user,
                    cohort_name,
                    [r["FID"] for r in rows] if narrowed or sorted_ else None,
                )

            def apply_filters(filters: Dict[str, Any], to_table: bool = False) -> None:
                """Set the individuals table's filters (from it, or the bar)."""
                state["filters"] = dict(filters)
                _filters_to_link(filters, link)
                link.push()
                recompute(to_individuals=to_table)

            def apply_family_filters(
                filters: Dict[str, Any], to_table: bool = False
            ) -> None:
                """Set the families table's filters (from it, or the bar)."""
                state["family_filters"] = dict(filters)
                _filters_to_link(filters, link, _FAMILY_FILTER_LINK)
                link.push()
                recompute(to_families=to_table)

            def on_filter(e: Dict[str, Any]) -> None:
                state["filter_tab"] = "Individuals"
                apply_filters(e.get("filters") or {})

            def on_family_filter(e: Dict[str, Any]) -> None:
                state["filter_tab"] = "Families"
                apply_family_filters(e.get("filters") or {})

            def clear_all() -> None:
                state["filters"], state["family_filters"] = {}, {}
                _filters_to_link({}, link)
                _filters_to_link({}, link, _FAMILY_FILTER_LINK)
                link.push()
                recompute(to_individuals=True, to_families=True)

            def make_clear(col: str, family: bool = False) -> Callable[[Any], None]:
                spec = _FAMILY_FILTER_LINK if family else _FILTER_LINK

                def handler(e: Any) -> None:
                    if e.value:
                        return
                    key = "family_filters" if family else "filters"
                    filters = dict(state[key])
                    empty = spec[col][1]
                    filters[col] = list(empty) if isinstance(empty, list) else empty
                    if family:
                        apply_family_filters(filters, to_table=True)
                    else:
                        apply_filters(filters, to_table=True)

                return handler

            def show_tab(label: str) -> None:
                tabs.set_value(label)

            @ui.refreshable
            def filter_bar() -> None:
                """The count, a removable chip per filter, and a way to the filters."""
                family_chips = _filter_chips(state["family_filters"], family_cols)
                individual_chips = _filter_chips(state["filters"], columns)
                chips = [(c, t, True) for c, t in family_chips] + [
                    (c, t, False) for c, t in individual_chips
                ]
                shown = len(filtered_state["individuals"])
                count = " · ".join(
                    (
                        _count_text(
                            len(state["families"]),
                            n_families,
                            bool(chips),
                            "family",
                            "families",
                        ),
                        _count_text(shown, n_total, bool(chips)),
                    )
                )
                with (
                    ui.row()
                    .classes(
                        "w-full min-w-0 items-center gap-2 px-3 py-1"
                        " bg-blue-50 border-b border-blue-100"
                    )
                    .mark("cohort-filter-bar")
                ):
                    ui.icon("filter_alt", size="18px").classes("text-blue-700")
                    ui.label(count).classes(
                        "text-sm font-semibold text-blue-900 whitespace-nowrap"
                    ).mark("cohort-count")
                    for col, text, family in chips:
                        ui.chip(
                            text,
                            removable=True,
                            color="white",
                            text_color="blue-10",
                            on_value_change=make_clear(col, family),
                        ).props("dense").classes("text-xs border border-blue-200").mark(
                            "cohort-filter-chip"
                        )
                    if chips:
                        ui.button("Clear all", on_click=clear_all).props(
                            "flat dense no-caps size=sm color=blue-8"
                        ).mark("cohort-clear-filters")
                    ui.space()
                    if state["tab"] not in _TABLE_TABS:
                        target = state["filter_tab"]
                        ui.button(
                            "Edit filters" if chips else f"Filter {target.lower()}",
                            icon="tune",
                            on_click=lambda: show_tab(target),
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
                    .props(GROUPED_TABS_PROPS)
                    .classes("w-full border-b border-gray-200")
                    .mark("cohort-tabs") as tabs
                ):
                    tab_objects: Dict[str, Any] = {}
                    for index, (caption, group) in enumerate(_TAB_GROUPS):
                        with tab_group(caption, first=index == 0):
                            for label, _slug, tip in group:
                                tab = ui.tab(label)
                                if tip:
                                    tab.tooltip(tip)
                                if not has_data[label]:
                                    tab.props("disable")
                                tab_objects[label] = tab

                # Badges: the filter for everyone, as Sample QC and the run;
                # UPD and STR count what only a curator settles.
                pill_tabs = ["Families", "Individuals", "Sample QC", "De novo"]
                if show_review:
                    pill_tabs += ["UPD", "STR"]
                pills: Dict[str, TabPills] = {
                    label: TabPills(tab_objects[label], f"{_SLUGS[label]}-tab-pills")
                    for label in pill_tabs
                    if has_data[label]
                }
                if workflow_button is not None:
                    pills["Workflow"] = TabPills(workflow_button, "workflow-tab-pills")

                with (
                    ui.element("div").classes("w-full").mark("cohort-filter-box")
                ) as filter_bar_box:
                    filter_bar()
                filter_bar_box.set_visibility(initial != "Workflow")

            def update_pills() -> None:
                rows = filtered_state["individuals"]
                keys = {(r["FID"], r["Sample ID"]) for r in rows}
                n_filters = len(_filter_chips(state["filters"], columns)) + len(
                    _filter_chips(state["family_filters"], family_cols)
                )
                pills["Individuals"].set(filter_pills(len(rows), n_total, n_filters))
                pills["Families"].set(
                    filter_pills(
                        len(state["families"]),
                        n_families,
                        n_filters,
                        "family",
                        "families",
                    )
                )
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
                if state["str_flags"]:
                    build_family_rows()
                    recompute()
                else:
                    update_pills()

            async def load_families() -> None:
                """The Families table's ancestry, FROH and de novo columns."""
                try:
                    family_data["facts"] = await asyncio.to_thread(
                        load_family_facts, cohort, store.data_dir, show_review
                    )
                except Exception:
                    logger.warning("Failed to read the family facts", exc_info=True)
                    return
                build_family_rows()
                recompute()

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
                with ui.tab_panel(tab_objects["Families"]).classes(panel_classes):
                    dt_ref["families"] = DataTable(
                        columns=family_cols,
                        rows=family_data["rows"],
                        row_key="FID",
                        pagination={"rowsPerPage": 20},
                        on_filter=on_family_filter,
                        initial_filters=linked_family_filters,
                        on_state_change=lambda: remember_order(),
                    )
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
                # they open on the filtered families and individuals.
                if linked_filters or linked_family_filters:
                    state["filters"] = dict(linked_filters)
                    state["family_filters"] = dict(linked_family_filters)
                    if linked_family_filters:
                        state["filter_tab"] = "Families"
                    elif linked_filters:
                        state["filter_tab"] = "Individuals"
                    recompute()

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

                if workflow_button is not None:
                    with ui.tab_panel("Workflow").classes("w-full p-0"):
                        workflow_content()

            def on_tab_change(e: Any) -> None:
                label = (
                    e.value.props["name"] if isinstance(e.value, ui.tab) else e.value
                )
                state["tab"] = label
                link.set("tab", _SLUGS.get(label), default=_DEFAULT_TAB)
                link.push()
                filter_bar_box.set_visibility(label != "Workflow")
                mark_workflow(label == "Workflow")
                filter_bar.refresh()
                opener = openers.get(label)
                if opener is not None:
                    opener()

            def mark_workflow(shown: bool) -> None:
                """The Workflow button looks pressed while its panel is shown."""
                if workflow_button is None:
                    return
                if shown:
                    workflow_button.props(remove="outline", add="unelevated")
                else:
                    workflow_button.props(remove="unelevated", add="outline")

            tabs.on_value_change(on_tab_change)
            mark_workflow(initial == "Workflow")
            if initial in openers:
                openers[initial]()
            update_pills()
            ui.timer(0.2, load_badges, once=True)
            ui.timer(0.2, load_families, once=True)
            # A bare page forgets an order left by an earlier, filtered visit.
            remember_order()

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
