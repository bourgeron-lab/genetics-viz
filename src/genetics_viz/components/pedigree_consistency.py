"""Pedigree disagreements across cohorts: the family banner, the lists, the explain dialog.

Every role sees that a family is also listed elsewhere (another cohort, another
ID), with links. Curators also see what differs, and can mark a finding as
explained (``utils/pedigree_consistency``). The scan reads every cohort
pedigree, so it runs off the event loop; it is cached on the files' stats.
"""

import asyncio
import logging
from typing import Any, Callable, Iterable, List, Mapping, Optional

from nicegui import ui

from genetics_viz.utils import pedigree_consistency as pc
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.sample_qc import severity_info

logger = logging.getLogger(__name__)


def family_href(cohort: str, family_id: str) -> str:
    return f"/cohort/{cohort}/family/{family_id}"


def _place_links(places: Iterable[pc.Place]) -> None:
    for cohort, fid in places:
        ui.link(f"{fid} in {cohort}", family_href(cohort, fid)).classes(
            "text-xs text-blue-700 whitespace-nowrap"
        )


def _color(finding: pc.Finding) -> str:
    return (
        severity_info(finding.severity)["color"]
        if finding.severity != "info"
        else "#9e9e9e"
    )


def render_finding(
    finding: pc.Finding,
    ack: Optional[pc.Ack],
    *,
    curator: bool,
    here: Optional[str] = None,
    cohort: Optional[str] = None,
    on_explain: Optional[Callable[[pc.Finding, Optional[pc.Ack]], None]] = None,
) -> None:
    """One finding: its label, the other places (linked), and for curators the
    details and the explanation. ``here`` is the family the page shows."""
    info = finding.info
    places = finding.others(here, cohort) if here else list(finding.places)
    if here and not (
        finding.code.startswith("SAME_FID") or finding.code == "FID_COLLISION"
    ):
        # An alias: the other IDs; this ID in other cohorts has its own line.
        places = [(c, f) for c, f in places if f != here] or places
    with ui.column().classes("w-full gap-0").mark(f"consistency-{finding.code}"):
        with ui.row().classes("items-center gap-2 wrap"):
            ui.badge(info.get("label", finding.code), color=_color(finding)).classes(
                "text-xs"
            ).tooltip(info.get("text", ""))
            if ack is not None:
                ui.badge("Explained", color="grey").classes("text-xs").tooltip(
                    f"{ack.user}, {ack.date}" + (f": {ack.note}" if curator else "")
                )
            ui.label("Also:" if here else "Families:").classes("text-xs text-gray-600")
            _place_links(places)
            if curator and on_explain is not None and finding.acknowledgeable:
                ui.button(
                    "Withdraw" if ack else "Explain",
                    on_click=lambda f=finding, a=ack: on_explain(f, a),
                ).props("flat dense no-caps size=sm color=blue-8").mark(
                    "consistency-explain"
                )
        if curator:
            for line in finding.details:
                ui.label(line).classes("text-xs text-gray-700 font-mono pl-2")
            if ack is not None and ack.note:
                ui.label(f"Explained: {ack.note}").classes(
                    "text-xs text-gray-500 italic pl-2"
                )


def explain_dialog(
    anchor: ui.element,
    data_dir: Any,
    finding: pc.Finding,
    ack: Optional[pc.Ack],
    on_saved: Callable[[], None],
) -> None:
    """Mark a finding explained (with a note), or withdraw the explanation."""
    if not can_write():
        return
    anchor.clear()
    with anchor:
        dialog = ui.dialog()
    with dialog, ui.card().classes("w-[32rem] max-w-full").mark("consistency-dialog"):
        ui.label(finding.info.get("label", finding.code)).classes(
            "text-lg font-semibold text-blue-900"
        )
        ui.label(finding.info.get("text", "")).classes("text-sm text-gray-700")
        for line in finding.details:
            ui.label(line).classes("text-xs font-mono text-gray-700")
        note = None
        if ack is None:
            note = (
                ui.input(
                    label="Why it is expected",
                    placeholder="e.g. the ADHD cohort holds only the trio",
                )
                .props("outlined dense")
                .classes("w-full")
                .mark("consistency-note")
            )
        else:
            ui.label(f"Explained by {ack.user}, {ack.date}: {ack.note}").classes(
                "text-sm text-gray-600"
            )

        def save() -> None:
            try:
                pc.save_ack(
                    data_dir,
                    finding,
                    note.value if note is not None else "",
                    user=get_current_user() or "",
                    status=pc.EXPLAINED if ack is None else pc.WITHDRAWN,
                )
            except (PermissionError, ValueError, OSError) as e:
                ui.notify(str(e), type="negative")
                return
            dialog.close()
            on_saved()

        with ui.row().classes("w-full justify-end gap-2"):
            ui.button("Cancel", on_click=dialog.close).props("flat no-caps")
            ui.button(
                "Mark explained" if ack is None else "Withdraw the explanation",
                on_click=save,
            ).props("unelevated no-caps color=primary").mark("consistency-save")
    dialog.open()


class FindingsPanel:
    """Findings loaded off-thread into a refreshable box, with the explain dialog.

    ``select`` picks the findings to show from the data directory's; an empty
    selection draws nothing (or ``empty``, when given).
    """

    def __init__(
        self,
        store: Any,
        select: Callable[[pc.Consistency], List[pc.Finding]],
        render: Callable[
            [List[pc.Finding], Mapping[str, pc.Ack], "FindingsPanel"], None
        ],
    ) -> None:
        self.store = store
        self.select = select
        self.render = render
        self.curator = can_write()
        self.findings: Optional[List[pc.Finding]] = None
        self.acks: Mapping[str, pc.Ack] = {}
        self._anchor = ui.element()
        with ui.element("div").classes("w-full") as self._box:

            @ui.refreshable
            def draw() -> None:
                if self.findings:
                    self.render(self.findings, self.acks, self)

            draw()
        self._draw = draw
        ui.timer(0.1, self.update, once=True)

    async def update(self) -> None:
        try:
            consistency = await asyncio.to_thread(pc.scan, self.store)
            self.acks = await asyncio.to_thread(pc.load_acks, self.store.data_dir)
            self.findings = self.select(consistency)
        except Exception:
            logger.warning("Failed to compare the cohort pedigrees", exc_info=True)
            self.findings = []
        self._draw.refresh()

    def explain(self, finding: pc.Finding, ack: Optional[pc.Ack]) -> None:
        def saved() -> None:
            with self._box:
                ui.timer(0.05, self.update, once=True)

        explain_dialog(self._anchor, self.store.data_dir, finding, ack, saved)


def render_family_banner(store: Any, family_id: str, cohort: Optional[str]) -> None:
    """Under the family title: where else the family is listed, and what differs."""

    def render(
        findings: List[pc.Finding], acks: Mapping[str, pc.Ack], panel: FindingsPanel
    ) -> None:
        open_ = pc.open_findings(findings, acks)
        worst = max((f for f in open_), key=lambda f: f.rank, default=None)
        tone = (
            "border-red-300 bg-red-50"
            if worst is not None and worst.severity == "high"
            else "border-amber-300 bg-amber-50"
            if worst is not None
            else "border-gray-200 bg-gray-50"
        )
        with (
            ui.column()
            .classes(f"w-full gap-1 border rounded-md p-2 mb-3 {tone}")
            .mark("family-consistency")
        ):
            if panel.curator and open_:
                ui.label(
                    "This family's pedigree disagrees with another cohort's or"
                    " another family ID's. Make the cohorts agree in the pedigree"
                    " editor (Sample QC tab), or mark the finding explained."
                ).classes("text-xs text-gray-700")
            for f in findings:
                render_finding(
                    f,
                    acks.get(f.key),
                    curator=panel.curator,
                    here=family_id,
                    cohort=cohort,
                    on_explain=panel.explain,
                )

    FindingsPanel(store, lambda c: c.for_family(family_id), render)


def _list(
    findings: List[pc.Finding], acks: Mapping[str, pc.Ack], panel: FindingsPanel
) -> None:
    open_ = pc.open_findings(findings, acks)
    with ui.column().classes("w-full gap-2"):
        ui.label(
            f"{len(open_)} to review, {len(findings) - len(open_)} explained or"
            " informative"
        ).classes("text-sm text-gray-700").mark("consistency-count")
        for f in findings:
            with ui.card().classes("w-full p-2"):
                render_finding(
                    f, acks.get(f.key), curator=panel.curator, on_explain=panel.explain
                )


def render_cohort_findings(store: Any, cohort_name: str) -> None:
    """Curators, on the Workflow tab: the findings involving the cohort."""

    def render(
        findings: List[pc.Finding], acks: Mapping[str, pc.Ack], panel: FindingsPanel
    ) -> None:
        open_ = pc.open_findings(findings, acks)
        with (
            ui.column()
            .classes("w-full gap-2 border border-amber-300 rounded-md p-3 bg-amber-50")
            .mark("workflow-consistency")
        ):
            ui.label(
                f"Pedigrees that disagree with another cohort's ({len(open_)} to"
                " review). ghfc-ngs rebuilds a shared family from the rows of"
                " the cohort it runs: make them agree before a run."
            ).classes("text-sm font-semibold text-amber-900")
            _list(findings, acks, panel)

    FindingsPanel(
        store,
        lambda c: [f for f in c.findings if cohort_name in f.cohorts],
        render,
    )


def render_all_findings(store: Any) -> None:
    """The data directory's findings, for the consistency page."""
    FindingsPanel(store, lambda c: list(c.findings), _list)
