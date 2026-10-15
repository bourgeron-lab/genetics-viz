"""Profile card: how the current user's gene badges are coloured.

Every GeneTrek geneset is assigned to a tier (1-4), "Show" (listed in the badge
tooltip, no colour) or "Hidden". A gene takes the colour of the best tier its
genesets reach (see ``utils/gene_scoring``). Edits apply to a live preview and
are saved server-side on demand (``utils/user_preferences``).
"""

from __future__ import annotations

import asyncio
from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.config_model import get_genetrek_url
from genetics_viz.utils.auth import is_admin
from genetics_viz.utils.gene_scoring import (
    GeneScorer,
    build_gene_badges,
    gene_badges_html,
    get_gene_scorer,
)
from genetics_viz.utils.genetrek import (
    GENETREK_CONFIG,
    GenesetInfo,
    get_snapshot,
    refresh_snapshot,
)
from genetics_viz.utils.user_preferences import (
    N_TIERS,
    SHOW,
    GenePrefs,
    get_default_gene_prefs,
    has_own_gene_prefs,
    reset_gene_prefs,
    save_default_gene_prefs,
    save_gene_prefs,
)

HIDDEN = "hidden"
ASSIGNMENT_OPTIONS: Dict[Any, str] = {
    **{t: f"Tier {t}" for t in range(1, N_TIERS + 1)},
    SHOW: "Show",
    HIDDEN: "Hidden",
}
_SOURCE_TEXT = {
    "user": "your own settings",
    "admin default": "the default set by an administrator",
    "packaged default": "the default shipped with the app",
}


def _draft_from(scorer: GeneScorer) -> Dict[str, Any]:
    """Editable state: one assignment per catalog set, keyed by its key."""
    sets = {
        gs.key: scorer.assignments.get(gs.name, HIDDEN)
        for gs in scorer.snapshot.catalog
    }
    stale = {k: scorer.prefs.sets[k] for k in scorer.stale_keys}
    return {"sets": sets, "colors": list(scorer.prefs.tier_colors), "stale": stale}


def _prefs_from(draft: Dict[str, Any], source: str = "user") -> GenePrefs:
    sets = {k: v for k, v in draft["sets"].items() if v != HIDDEN}
    # Keys GeneTrek no longer has are kept: the list may come back.
    sets.update(draft["stale"])
    return GenePrefs(tuple(draft["colors"]), sets, source)


def render_gene_prefs_card(username: str) -> None:
    """The "Gene badges (GeneTrek)" card of the profile page."""
    with ui.card().classes("w-full p-6 mt-6"):
        with ui.row().classes("w-full items-center gap-2"):
            ui.label("Gene badges").classes("text-xl font-semibold")
            ui.link("GeneTrek", get_genetrek_url(), new_tab=True).classes("text-sm")
        _card_body(username)


@ui.refreshable
def _card_body(username: str) -> None:
    snapshot = get_snapshot()
    scorer = get_gene_scorer(username)

    with ui.row().classes("w-full items-center gap-2"):
        ui.label(snapshot.status_text()).classes("text-sm text-gray-600")
        if is_admin():

            async def refresh_now() -> None:
                ui.notify("Checking GeneTrek…")
                await asyncio.to_thread(refresh_snapshot, force=True)
                _card_body.refresh()

            ui.button("Refresh now", icon="refresh", on_click=refresh_now).props(
                "flat dense size=sm"
            )

    if not snapshot.catalog:
        ui.label("The genesets are listed once GeneTrek data has been loaded.").classes(
            "text-sm text-orange-700 mt-2"
        )
        return

    ui.label(
        "Each geneset can colour badges (Tier 1 is the strongest colour), only be"
        " listed in the badge tooltip (Show), or be left out (Hidden). A gene takes"
        " the colour of the best tier among its genesets."
    ).classes("text-sm text-gray-700 mt-2")
    ui.label(f"Currently using {_SOURCE_TEXT[scorer.prefs.source]}.").classes(
        "text-xs text-gray-500"
    )

    draft = _draft_from(scorer)

    # -- tier colours and preview ------------------------------------------
    with ui.row().classes("w-full items-end gap-3 mt-3"):
        for i in range(N_TIERS):

            def make_color_handler(index: int):
                def handler(e) -> None:
                    if e.value:
                        draft["colors"][index] = e.value
                        preview.refresh()

                return handler

            ui.color_input(
                f"Tier {i + 1}",
                value=draft["colors"][i],
                on_change=make_color_handler(i),
            ).props("dense outlined").classes("w-32")

    @ui.refreshable
    def preview() -> None:
        draft_scorer = GeneScorer(snapshot, _prefs_from(draft, "draft"))
        genes = list(GENETREK_CONFIG.get("preview_genes") or [])
        with ui.row().classes("items-center gap-2 mt-2"):
            ui.label("Preview:").classes("text-sm font-semibold")
            ui.html(
                gene_badges_html(
                    build_gene_badges(genes, draft_scorer, cap=None, sort=False)
                ),
                sanitize=False,
            )

    preview()

    # -- geneset assignments -----------------------------------------------
    show_secondary = {"value": False}

    @ui.refreshable
    def geneset_list() -> None:
        by_category: Dict[str, List[GenesetInfo]] = {}
        for gs in snapshot.catalog:
            assigned = draft["sets"].get(gs.key, HIDDEN) != HIDDEN
            if gs.is_primary or assigned or show_secondary["value"]:
                by_category.setdefault(gs.category or "Other", []).append(gs)
        for category, sets in by_category.items():
            n_used = sum(draft["sets"].get(gs.key, HIDDEN) != HIDDEN for gs in sets)
            with (
                ui.expansion(
                    f"{category} ({n_used}/{len(sets)} used)", value=n_used > 0
                )
                .classes("w-full")
                .props("dense")
            ):
                for gs in sets:
                    _geneset_row(gs, draft, preview)

    with ui.row().classes("w-full items-center mt-4"):
        ui.label("Genesets").classes("text-lg font-semibold")
        ui.space()

        def toggle_secondary(e) -> None:
            show_secondary["value"] = e.value
            geneset_list.refresh()

        ui.switch(
            "Show older versions and sub-lists", on_change=toggle_secondary
        ).props("dense")

    geneset_list()

    if draft["stale"]:
        ui.label(
            "No longer in GeneTrek, kept in your settings: "
            + ", ".join(sorted(draft["stale"]))
        ).classes("text-xs text-orange-700 mt-2")

    # -- actions -----------------------------------------------------------
    def save() -> None:
        save_gene_prefs(username, _prefs_from(draft))
        ui.notify("Gene badge settings saved", type="positive")
        _card_body.refresh()

    def reset() -> None:
        reset_gene_prefs(username)
        ui.notify("Back to the default settings", type="info")
        _card_body.refresh()

    def make_default() -> None:
        save_default_gene_prefs(_prefs_from(draft, "admin default"))
        ui.notify("Saved as the default for every user", type="positive")
        _card_body.refresh()

    with ui.row().classes("gap-2 mt-4"):
        ui.button("Save", icon="save", on_click=save).props("color=blue unelevated")
        if has_own_gene_prefs(username):
            ui.button("Reset to default", icon="restart_alt", on_click=reset).props(
                "outline"
            )
        if is_admin():
            ui.button(
                "Use as default for everyone", icon="groups", on_click=make_default
            ).props("outline color=red").tooltip(
                f"Currently: {get_default_gene_prefs().source}"
            )


def _geneset_row(gs: GenesetInfo, draft: Dict[str, Any], preview) -> None:
    def on_change(e) -> None:
        draft["sets"][gs.key] = e.value
        preview.refresh()

    with ui.row().classes("w-full items-center gap-3 no-wrap py-1"):
        ui.select(
            ASSIGNMENT_OPTIONS,
            value=draft["sets"].get(gs.key, HIDDEN),
            on_change=on_change,
        ).props("dense outlined options-dense").classes("w-28 shrink-0").mark(
            f"geneset-{gs.key}"
        )
        with ui.column().classes("gap-0 min-w-0"):
            with ui.row().classes("items-center gap-1"):
                ui.label(gs.label).classes("text-sm font-medium")
                if gs.is_rollup:
                    ui.badge("roll-up", color="grey-5").props("outline").tooltip(
                        "Union of: " + ", ".join(gs.derived_from)
                    )
                if not gs.is_primary:
                    ui.badge("secondary", color="grey-5").props("outline")
            meta = f"{gs.key} · {gs.n_genes:,} genes"
            if gs.subcategory:
                meta = f"{gs.subcategory} · {meta}"
            ui.label(meta).classes("text-xs text-gray-500 truncate w-full")
        if gs.description:
            ui.icon("info", size="xs").classes("text-gray-400 shrink-0").tooltip(
                gs.description
            )
