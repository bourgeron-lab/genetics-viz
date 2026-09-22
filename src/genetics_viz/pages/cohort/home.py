"""Home page - displays available cohorts and quick search."""

import asyncio
import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.sample_dialog import show_sample_dialog
from genetics_viz.models import Cohort, CohortStub
from genetics_viz.pages.cohort.components.status_panel import (
    render_status_summary,
    status_headline,
)
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.pedigree_labels import (
    PHENO_LABELS,
    PHENO_ORDER,
    normalize_pheno,
)
from genetics_viz.utils.data_availability import (
    check_family_availability,
    check_sample_availability,
)
from genetics_viz.utils.pedigree import load_family_pedigree
from genetics_viz.utils.sharding import get_family_path, get_sample_path

# Diagnostic priority for "highest" resolution
_DIAG_PRIORITY = {"pathogenic": 3, "uncertain": 2, "benign": 1}

# Cohort card geometry. Three separate things were eating the card's width:
#
# ``p-0``: nicegui-card carries its own 16px padding, which doubled up with the
# 16px Quasar puts on every q-card__section — 32px a side, leaving 256px of a
# 320px card for content.
#
# ``w-full`` on each section: nicegui-card is a flex column with
# ``align-items: flex-start``, so a section is sized by its content instead of
# stretching. The sections came out 274px wide in a 320px card, which is why
# the phenotype table's right edge never lined up with the card's.
#
# ``shrink-0``: flex-shrink defaults to 1, so a card could be compressed below
# its own 20rem on a narrow row. This is for a consistent card width, not an
# overflow fix -- with min-w-0 throughout, a shrunken card's contents shrink
# with it rather than spilling.
#
# ``pl-2 pr-3``: the 4px accent border on the left is inside the border box and
# eats into the content, so equal padding would leave the left inset 4px wider
# than the right.
_CARD_BASE = "w-80 shrink-0 p-0"
_SECTION_X = "w-full pl-2 pr-3 py-2"


def _compute_cohort_stats(cohort: Cohort, data_dir: Path) -> Dict[str, Any]:
    """Compute per-phenotype counts and diagnostic breakdown for a cohort.

    Returns dict with:
        per_pheno: {pheno_key: {"n": int, "pathogenic": int, "uncertain": int}}
    """
    # Build sample_id -> normalized phenotype map
    sample_pheno: Dict[str, str] = {}
    for fam in cohort.families.values():
        for s in fam.samples:
            sample_pheno[s.sample_id] = normalize_pheno(s.phenotype)

    # Load diagnostics per sample (highest priority wins)
    sample_diag: Dict[str, str] = defaultdict(str)
    for fname in ["snvs.tsv", "svs.tsv"]:
        diag_file = data_dir / "diagnostics" / fname
        if not diag_file.exists():
            continue
        with open(diag_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Ignore", "0") == "1":
                    continue
                sample = row.get("Sample", "")
                diag = row.get("Diagnostic", "")
                if sample not in sample_pheno or not diag:
                    continue
                current = sample_diag[sample]
                if _DIAG_PRIORITY.get(diag, 0) > _DIAG_PRIORITY.get(current, 0):
                    sample_diag[sample] = diag

    # Bucket by phenotype
    per_pheno: Dict[str, Dict[str, int]] = {
        key: {"n": 0, "pathogenic": 0, "uncertain": 0} for key in PHENO_ORDER
    }
    for sample_id, pheno in sample_pheno.items():
        bucket = per_pheno[pheno]
        bucket["n"] += 1
        diag = sample_diag.get(sample_id, "")
        if diag == "pathogenic":
            bucket["pathogenic"] += 1
        elif diag == "uncertain":
            bucket["uncertain"] += 1

    return {"per_pheno": per_pheno}


def _gather_headlines(cohorts: List[Cohort]) -> Dict[str, Dict[str, Any]]:
    """Read every cohort's run-status headline in one offloaded pass.

    One stat and one small JSON read per cohort, but on a network-mounted data
    directory that still adds up, and the page handler must not block the
    single-threaded event loop for it.
    """
    return {
        cohort.name: status_headline(cohort)
        for cohort in cohorts
        if cohort.params_file is not None
    }


def _render_status_expansion(cohort: Cohort, headline: Dict[str, Any]) -> None:
    """Folded pipeline status on a cohort card.

    The card navigates on click, so this swallows its own clicks — otherwise
    opening the expansion would leave the page.
    """
    loaded = {"value": False}

    # min-w-0 clamps this subtree's min-content contribution: the card is a
    # fixed 20rem, and without it the expansion's intrinsic width forced the
    # card section wider than the card itself (and clipped the phenotype
    # table's last column on the way).
    with (
        ui.element("div")
        .classes("w-full min-w-0 mt-3")
        .on("click", None, js_handler="(e) => e.stopPropagation()")
    ):
        with (
            ui.expansion()
            .classes("w-full min-w-0 border border-gray-200 rounded-md")
            .props("dense dense-toggle") as exp
        ):
            with exp.add_slot("header"):
                # Stacked, not one row: the card is 20rem wide and a caption
                # like "3/8 steps complete, 1 not measured" pushes a single
                # row -- and with it the whole expansion -- past the card edge.
                with ui.column().classes("gap-0 w-full min-w-0"):
                    with ui.row().classes("items-center gap-2 no-wrap"):
                        ui.label("Status").classes(
                            "text-xs font-semibold text-gray-700 whitespace-nowrap"
                        )
                        ui.badge(
                            headline.get("label", ""),
                            color=headline.get("color", "grey"),
                        ).classes("text-[10px]")
                    ui.label(headline.get("caption", "")).classes(
                        "text-[11px] text-gray-500 w-full truncate"
                    )

            # The per-step rows re-read the record and may hash the pedigree,
            # so they wait until someone actually opens this.
            @ui.refreshable
            def body() -> None:
                if loaded["value"]:
                    render_status_summary(cohort)
                else:
                    with ui.row().classes("items-center gap-2 py-2"):
                        ui.spinner(size="sm", color="blue")
                        ui.label("Reading run record…").classes("text-xs text-gray-500")

            body()

        def on_open(e: Any) -> None:
            if not e.value or loaded["value"]:
                return

            def load() -> None:
                loaded["value"] = True
                body.refresh()

            ui.timer(0.1, load, once=True)

        exp.on_value_change(on_open)


def _render_incomplete_card(stub: CohortStub) -> None:
    """A cohort directory with no usable pedigree: shown, but not explorable."""
    unreadable = stub.reason == "unreadable"
    with ui.card().classes(
        f"{_CARD_BASE} border-l-4 border-amber-500 bg-amber-50 cursor-default"
    ):
        with ui.card_section().classes(_SECTION_X):
            with ui.row().classes("items-center gap-2 no-wrap"):
                ui.icon("warning", color="amber-8").classes("text-xl")
                ui.label(stub.name).classes(
                    "text-xl font-bold text-amber-900 truncate min-w-0"
                )

        with ui.card_section().classes(f"{_SECTION_X} pt-0"):
            # A pedigree that is absent and one that will not parse need
            # different fixes, so they do not share a message.
            ui.label(
                "Unreadable pedigree file" if unreadable else "Missing pedigree file"
            ).classes("text-sm font-semibold text-amber-800")
            ui.label(
                "This cohort cannot be explored until its pedigree is in place."
            ).classes("text-xs text-gray-600")
            if stub.error:
                ui.label(stub.error).classes("text-xs text-amber-900 italic mt-1")

        with ui.card_section().classes(f"{_SECTION_X} bg-amber-100"):
            ui.label(f"⚠ {stub.expected_pedigree.name}").classes(
                "text-xs text-amber-900 font-mono"
            )


@ui.page("/")
async def home_page() -> None:
    """Render the home/welcome page."""
    if redirect := check_auth():
        return redirect
    create_header()

    # Pre-load IGV.js for sample visualization dialog
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>'
    )

    with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
        ui.label("Welcome to Genetics-Viz").classes(
            "text-3xl font-bold mb-2 text-blue-900"
        )
        ui.label("Select a cohort to explore").classes("text-lg text-gray-600 mb-6")

        try:
            store = get_data_store()

            # Quick Search section
            with ui.card().classes("w-full mb-6 bg-gray-50"):
                with ui.column().classes("p-4 gap-3"):
                    with ui.row().classes("items-center gap-2"):
                        ui.icon("search", color="blue").classes("text-xl")
                        ui.label("Quick Search").classes(
                            "text-lg font-semibold text-blue-800"
                        )

                    search_result = ui.column().classes("w-full")

                    with ui.row().classes("items-end gap-4 w-full"):
                        sample_input = (
                            ui.input(
                                label="Sample Barcode",
                                placeholder="e.g. B00EYRL",
                            )
                            .props("outlined dense")
                            .classes("w-64")
                        )

                        async def search_sample():
                            barcode = (
                                sample_input.value.strip() if sample_input.value else ""
                            )
                            if not barcode:
                                ui.notify("Enter a barcode", type="warning")
                                return
                            search_result.clear()
                            with search_result:
                                with ui.row().classes("items-center gap-2"):
                                    ui.spinner(size="sm")
                                    ui.label("Searching...").classes("text-gray-500")

                            sample_path = get_sample_path(store.data_dir, barcode)
                            exists = await asyncio.to_thread(sample_path.is_dir)

                            search_result.clear()
                            with search_result:
                                if exists:
                                    avail = await asyncio.to_thread(
                                        check_sample_availability,
                                        store.data_dir,
                                        barcode,
                                    )
                                    with ui.card().classes(
                                        "w-full border-l-4 border-green-500"
                                    ):
                                        with ui.row().classes("items-center gap-3 p-3"):
                                            ui.icon(
                                                "check_circle", color="green"
                                            ).classes("text-xl")
                                            ui.label(
                                                f"Sample found: {barcode}"
                                            ).classes("font-semibold text-green-800")
                                            for key, label in [
                                                ("cram", "CRAM"),
                                                ("bedgraph", "Bedgraph"),
                                                ("vaf_bedgraph", "VAF"),
                                                ("deepvariant", "DeepVariant"),
                                            ]:
                                                color = (
                                                    "green"
                                                    if avail.get(key)
                                                    else "grey"
                                                )
                                                ui.badge(label, color=color).props(
                                                    "outline"
                                                    if not avail.get(key)
                                                    else ""
                                                ).classes("text-xs")
                                            ui.space()
                                            if avail.get("bedgraph") or avail.get(
                                                "cram"
                                            ):
                                                ui.button(
                                                    "Visualize",
                                                    icon="visibility",
                                                    on_click=lambda: show_sample_dialog(
                                                        barcode
                                                    ),
                                                ).props("color=blue dense")
                                else:
                                    with ui.card().classes(
                                        "w-full border-l-4 border-amber-500"
                                    ):
                                        with ui.row().classes("items-center gap-2 p-3"):
                                            ui.icon("warning", color="amber").classes(
                                                "text-xl"
                                            )
                                            ui.label(
                                                f"Sample not found: {barcode}"
                                            ).classes("text-amber-800")

                        ui.button(
                            "Search Sample",
                            on_click=search_sample,
                            icon="person_search",
                        ).props("color=blue dense")
                        sample_input.on("keydown.enter", search_sample)

                    with ui.row().classes("items-end gap-4 w-full"):
                        family_input = (
                            ui.input(
                                label="Family ID",
                                placeholder="e.g. C0733-011-068",
                            )
                            .props("outlined dense")
                            .classes("w-64")
                        )

                        async def search_family():
                            fid = (
                                family_input.value.strip() if family_input.value else ""
                            )
                            if not fid:
                                ui.notify("Enter a family ID", type="warning")
                                return
                            search_result.clear()
                            with search_result:
                                with ui.row().classes("items-center gap-2"):
                                    ui.spinner(size="sm")
                                    ui.label("Searching...").classes("text-gray-500")

                            family_path = get_family_path(store.data_dir, fid)
                            exists = await asyncio.to_thread(family_path.is_dir)

                            search_result.clear()
                            with search_result:
                                if exists:
                                    ped_file = family_path / f"{fid}.pedigree.tsv"
                                    member_count = 0
                                    if ped_file.exists():
                                        members = await asyncio.to_thread(
                                            load_family_pedigree, ped_file
                                        )
                                        member_count = len(members)

                                    avail = await asyncio.to_thread(
                                        check_family_availability,
                                        store.data_dir,
                                        fid,
                                        [],
                                    )
                                    with ui.card().classes(
                                        "w-full border-l-4 border-green-500"
                                    ):
                                        with ui.row().classes("items-center gap-3 p-3"):
                                            ui.icon(
                                                "check_circle", color="green"
                                            ).classes("text-xl")
                                            ui.label(f"Family found: {fid}").classes(
                                                "font-semibold text-green-800"
                                            )
                                            if member_count:
                                                ui.badge(
                                                    f"{member_count} members"
                                                ).props("color=blue")
                                            for key, label in [
                                                ("pedigree", "Pedigree"),
                                                ("vcfs", "VCFs"),
                                                ("wombat", "Wombat"),
                                                ("wisecondorx", "WisecondorX"),
                                            ]:
                                                color = (
                                                    "green"
                                                    if avail.get(key)
                                                    else "grey"
                                                )
                                                ui.badge(label, color=color).props(
                                                    "outline"
                                                    if not avail.get(key)
                                                    else ""
                                                ).classes("text-xs")
                                            ui.space()
                                            ui.button(
                                                "Open Family Page",
                                                icon="open_in_new",
                                                on_click=lambda: ui.navigate.to(
                                                    f"/family/{fid}"
                                                ),
                                            ).props("color=blue dense")
                                else:
                                    with ui.card().classes(
                                        "w-full border-l-4 border-amber-500"
                                    ):
                                        with ui.row().classes("items-center gap-2 p-3"):
                                            ui.icon("warning", color="amber").classes(
                                                "text-xl"
                                            )
                                            ui.label(
                                                f"Family not found: {fid}"
                                            ).classes("text-amber-800")

                        ui.button(
                            "Search Family",
                            on_click=search_family,
                            icon="family_restroom",
                        ).props("color=blue dense")
                        family_input.on("keydown.enter", search_family)

            if not store.cohorts and not store.incomplete_cohorts:
                with ui.card().classes("w-full p-6 bg-yellow-50"):
                    ui.label("⚠️ No cohorts found").classes(
                        "text-xl font-semibold text-yellow-800"
                    )
                    ui.label(
                        f"No cohort directories were found in: {store.cohorts_dir}"
                    ).classes("text-gray-600")
                    ui.label(
                        "Each cohort is a directory under cohorts/, holding a"
                        " <name>.pedigree.tsv file."
                    ).classes("text-gray-500 text-sm")
                return

            cohorts = sorted(store.cohorts.values(), key=lambda c: c.name)

            # One offloaded pass rather than one blocking read per card.
            headlines = await asyncio.to_thread(_gather_headlines, cohorts)

            # Display cohorts as cards
            ui.label("Available Cohorts").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.row().classes("w-full flex-wrap gap-4"):
                for cohort in cohorts:
                    cohort_name = cohort.name  # Capture name for lambda
                    stats = _compute_cohort_stats(cohort, store.data_dir)
                    with (
                        ui.card()
                        .classes(
                            "cursor-pointer hover:shadow-lg transition-shadow"
                            f" {_CARD_BASE} border-l-4 border-blue-500"
                        )
                        .on(
                            "click",
                            lambda _, n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                        )
                    ):
                        with ui.card_section().classes(_SECTION_X):
                            ui.label(cohort.name).classes(
                                "text-xl font-bold text-blue-700"
                            )

                        with ui.card_section().classes(f"{_SECTION_X} pt-0"):
                            # Top-line summary: Families + Samples
                            with ui.row().classes("gap-6 mb-3"):
                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_families)).classes(
                                        "text-2xl font-bold text-blue-600"
                                    )
                                    ui.label("Families").classes(
                                        "text-xs text-gray-500"
                                    )

                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_samples)).classes(
                                        "text-2xl font-bold text-green-600"
                                    )
                                    ui.label("Samples").classes("text-xs text-gray-500")

                            # Per-phenotype breakdown table
                            per_pheno = stats["per_pheno"]
                            non_empty = [
                                p for p in PHENO_ORDER if per_pheno[p]["n"] > 0
                            ]
                            if non_empty:
                                with ui.element("div").classes(
                                    "border border-gray-300 rounded-md"
                                    " overflow-hidden w-full"
                                ):
                                    with ui.grid(columns=4).classes(
                                        "gap-x-3 gap-y-1 text-xs w-full px-2 py-1"
                                    ):
                                        # Header row (separator below via border-b)
                                        _hdr = (
                                            "text-gray-500 font-semibold"
                                            " border-b border-gray-200 pb-1"
                                            " whitespace-nowrap"
                                        )
                                        ui.label("Pheno").classes(_hdr)
                                        ui.label("N").classes(_hdr + " text-right")
                                        ui.label("Pat/Unc").classes(
                                            _hdr + " text-right"
                                        )
                                        ui.label("%").classes(_hdr + " text-right")
                                        # Data rows
                                        for pheno_key in non_empty:
                                            bucket = per_pheno[pheno_key]
                                            n = bucket["n"]
                                            pat = bucket["pathogenic"]
                                            unc = bucket["uncertain"]
                                            diag_pct = (
                                                (pat + unc) / n * 100 if n > 0 else 0
                                            )

                                            ui.label(PHENO_LABELS[pheno_key]).classes(
                                                "text-gray-700 whitespace-nowrap"
                                            )
                                            ui.label(str(n)).classes(
                                                "text-gray-700 text-right"
                                            )
                                            with ui.row().classes(
                                                "gap-0.5 justify-end items-baseline"
                                            ):
                                                ui.label(str(pat)).classes(
                                                    "text-red-600 font-semibold"
                                                )
                                                ui.label("/").classes("text-gray-400")
                                                ui.label(str(unc)).classes(
                                                    "text-amber-500 font-semibold"
                                                )
                                            ui.label(f"{diag_pct:.0f}%").classes(
                                                "text-gray-700 text-right"
                                            )

                            # Pipeline status, folded. Only cohorts with a
                            # workflow parameters file have a pipeline to
                            # report on.
                            headline = headlines.get(cohort_name)
                            if headline is not None:
                                _render_status_expansion(cohort, headline)
                            else:
                                ui.label(
                                    f"No workflow parameters file"
                                    f" ({cohort_name}.params.yml)"
                                ).classes("text-xs text-gray-400 italic mt-3")

                        with ui.card_section().classes(f"{_SECTION_X} bg-gray-50"):
                            ui.label(f"📄 {cohort.pedigree_file.name}").classes(
                                "text-xs text-gray-400"
                            )

            if store.incomplete_cohorts:
                ui.label("Incomplete Cohorts").classes(
                    "text-2xl font-semibold mt-8 mb-1 text-amber-800"
                )
                ui.label(
                    "Directories under cohorts/ with no readable pedigree file."
                ).classes("text-sm text-gray-600 mb-4")

                with ui.row().classes("w-full flex-wrap gap-4"):
                    for stub in sorted(store.incomplete_cohorts, key=lambda s: s.name):
                        _render_incomplete_card(stub)

        except RuntimeError as e:
            ui.label(f"Error: {e}").classes("text-red-500")
