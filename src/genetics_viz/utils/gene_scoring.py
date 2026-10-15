"""Gene badges: colour, tooltip and GeneTrek link, from the user's tiers.

Each user assigns GeneTrek genesets to one of four tiers, or marks them "show"
(listed in the tooltip, no colour), in their profile
(``utils/user_preferences.GenePrefs``). A gene takes the colour of the best
tier any of its genesets is in. There is no arithmetic: a sum would
double-count GeneTrek's coexisting versions of one list (``hc_ndd_v1`` …
``hc_ndd_v7``) and its roll-ups next to their parts.

A scorer belongs to one user. ``get_gene_scorer()`` reads the user from
``app.storage.user``, so call it on the event loop (in a page handler or a
timer callback), never inside ``asyncio.to_thread``; pass the scorer into any
threaded code instead.
"""

from __future__ import annotations

import html
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Tuple

from genetics_viz.utils.genetrek import (
    GENETREK_CONFIG,
    GenesetInfo,
    GenetrekSnapshot,
    gene_url,
    get_snapshot,
    search_url,
)
from genetics_viz.utils.user_preferences import (
    N_TIERS,
    SHOW,
    Assignment,
    GenePrefs,
    get_gene_prefs,
)

NEUTRAL_COLOR = str(GENETREK_CONFIG.get("neutral_color") or "#eeeeee")
OVERFLOW_COLOR = "#9e9e9e"
#: Badges shown per cell before the rest collapse into "+N genes".
BADGE_CAP = 6


def text_color_for(background: str) -> str:
    """Black or white, whichever reads on ``background`` (a ``#rrggbb``)."""
    try:
        r, g, b = (int(background.lstrip("#")[i : i + 2], 16) / 255 for i in (0, 2, 4))
    except (ValueError, IndexError):
        return "black"

    def lin(c: float) -> float:
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4

    luminance = 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
    # Contrast with white beats contrast with black below about 0.18.
    return "black" if luminance > 0.18 else "white"


class GeneScorer:
    """One user's view of GeneTrek: tiers, colours, tooltips, links."""

    def __init__(self, snapshot: GenetrekSnapshot, prefs: GenePrefs):
        self.snapshot = snapshot
        self.prefs = prefs
        #: Catalog name -> assignment, preference keys resolved via aliases.
        self.assignments: Dict[str, Assignment] = {}
        #: Preference keys GeneTrek no longer has.
        self.stale_keys: List[str] = []
        for key, assignment in prefs.sets.items():
            gs = snapshot.geneset(key)
            if gs is None:
                self.stale_keys.append(key)
            else:
                self.assignments[gs.name] = assignment
        self._tooltips: Dict[str, str] = {}

    # -- lookups -----------------------------------------------------------

    def resolve(self, identifier: Optional[str]) -> Optional[str]:
        return self.snapshot.resolve(identifier)

    def sets_of(self, identifier: str) -> FrozenSet[str]:
        symbol = self.resolve(identifier)
        return self.snapshot.members.get(symbol, frozenset()) if symbol else frozenset()

    def tier(self, identifier: str) -> Optional[int]:
        tiers = [
            a
            for name in self.sets_of(identifier)
            if isinstance(a := self.assignments.get(name), int)
        ]
        return min(tiers) if tiers else None

    def rank(self, identifier: str) -> int:
        """Sort key: higher is better; 0 for a gene in no tiered set."""
        tier = self.tier(identifier)
        return 0 if tier is None else N_TIERS + 1 - tier

    def tier_color(self, tier: int) -> str:
        return self.prefs.tier_colors[tier - 1]

    def color(self, identifier: str) -> str:
        tier = self.tier(identifier)
        return NEUTRAL_COLOR if tier is None else self.tier_color(tier)

    def text_color(self, identifier: str) -> str:
        return text_color_for(self.color(identifier))

    def url(self, identifier: str) -> str:
        return gene_url(self.resolve(identifier) or identifier)

    def tier_symbols(self, max_tier: int) -> FrozenSet[str]:
        """Symbols in any set assigned to tiers 1..``max_tier``."""
        out: set = set()
        for name, assignment in self.assignments.items():
            if isinstance(assignment, int) and assignment <= max_tier:
                out.update(self.snapshot.set_members.get(name, ()))
        return frozenset(out)

    def sets_by_tier(self) -> Dict[int, List[GenesetInfo]]:
        """The user's tiered sets, for legends."""
        out: Dict[int, List[GenesetInfo]] = {t: [] for t in range(1, N_TIERS + 1)}
        for gs in self.snapshot.catalog:
            assignment = self.assignments.get(gs.name)
            if isinstance(assignment, int):
                out[assignment].append(gs)
        return out

    # -- tooltip -----------------------------------------------------------

    def tooltip(self, identifier: str) -> str:
        cached = self._tooltips.get(identifier)
        if cached is None:
            cached = self._tooltips[identifier] = self._build_tooltip(identifier)
        return cached

    def _build_tooltip(self, identifier: str) -> str:
        snap = self.snapshot
        if not snap.ok:
            return f"{identifier}\n{snap.status_text()}"
        symbol = self.resolve(identifier)
        if symbol is None:
            return f"{identifier}: not in GeneTrek"

        info = snap.info.get(symbol)
        head = symbol if symbol == identifier else f"{identifier} → {symbol}"
        lines = [f"{head} — {info.name}" if info and info.name else head]

        by_tier: Dict[int, List[str]] = {}
        shown: List[str] = []
        others = 0
        for gs in snap.catalog:  # catalog order is GeneTrek's display order
            if gs.name not in snap.members.get(symbol, ()):
                continue
            assignment = self.assignments.get(gs.name)
            if isinstance(assignment, int):
                by_tier.setdefault(assignment, []).append(gs.label)
            elif assignment == SHOW:
                shown.append(gs.label)
            else:
                others += 1
        for tier in sorted(by_tier):
            lines.append(f"Tier {tier}: {', '.join(by_tier[tier])}")
        if shown:
            lines.append(f"Also in: {', '.join(shown)}")
        if not by_tier and not shown:
            lines.append("In none of your genesets")
        if others:
            plural = "s" if others > 1 else ""
            lines.append(f"+ {others} other GeneTrek set{plural}, click to open")
        if info and info.values:
            lines.append(" · ".join(f"{label} {value}" for label, value in info.values))
        if info and info.inheritance:
            modes = sorted({mode for _, mode in info.inheritance})
            lines.append(f"Inheritance: {' / '.join(modes)}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Per-user scorer cache
# ---------------------------------------------------------------------------

#: (username, prefs revision, snapshot identity) -> scorer. The snapshot is
#: held by the scorer, so its id() cannot be reused while the entry lives.
_scorers: Dict[Tuple[Any, ...], GeneScorer] = {}
_MAX_SCORERS = 64


def clear_scorer_cache() -> None:
    _scorers.clear()


def _current_username() -> Optional[str]:
    from genetics_viz.utils.auth import get_current_user

    try:
        return get_current_user()
    except RuntimeError:  # no request context: tests, startup
        return None


def get_gene_scorer(username: Optional[str] = None) -> GeneScorer:
    """The scorer for ``username``, by default the logged-in user."""
    if username is None:
        username = _current_username()
    prefs = get_gene_prefs(username)
    snapshot = get_snapshot()
    key = (username, prefs.revision, id(snapshot))
    scorer = _scorers.get(key)
    if scorer is None or scorer.snapshot is not snapshot:
        if len(_scorers) >= _MAX_SCORERS:
            _scorers.clear()
        scorer = _scorers[key] = GeneScorer(snapshot, prefs)
    return scorer


# ---------------------------------------------------------------------------
# Badges
# ---------------------------------------------------------------------------


def split_genes(value: Any, separators: str = ",&") -> List[str]:
    """Distinct non-empty gene identifiers from a ``,``/``&`` joined cell."""
    if value is None:
        return []
    text = str(value)
    for sep in separators[1:]:
        text = text.replace(sep, separators[0])
    out: List[str] = []
    for part in text.split(separators[0]):
        part = part.strip()
        if part and part not in out:
            out.append(part)
    return out


def gene_badge(identifier: str, scorer: GeneScorer) -> Dict[str, Any]:
    color = scorer.color(identifier)
    return {
        "label": identifier,
        "color": color,
        "textColor": text_color_for(color),
        "tooltip": scorer.tooltip(identifier),
        "href": scorer.url(identifier),
        "score": scorer.rank(identifier),
    }


def build_gene_badges(
    identifiers: Iterable[str],
    scorer: GeneScorer,
    *,
    cap: Optional[int] = BADGE_CAP,
    sort: bool = True,
    bordered: Iterable[str] = (),
    border_color: str = "black",
) -> List[Dict[str, Any]]:
    """Badges for a cell: best tier first, capped, the rest as "+N genes".

    ``bordered`` genes (an SV's exonic genes) get a border. The "+N genes"
    badge links to GeneTrek's table filtered to every gene of the cell.
    """
    idents = split_genes(",".join(i for i in identifiers if i))
    bordered = set(bordered)
    badges = []
    for ident in idents:
        badge = gene_badge(ident, scorer)
        if ident in bordered:
            badge["borderColor"] = border_color
        badges.append(badge)
    if sort:
        badges.sort(key=lambda b: b["score"], reverse=True)  # stable
    if cap is not None and len(badges) > cap:
        rest = badges[cap:]
        badges = badges[:cap]
        names = [b["label"] for b in rest]
        badges.append(
            {
                "label": f"+{len(rest)} genes",
                "color": OVERFLOW_COLOR,
                "textColor": "white",
                "tooltip": (
                    f"{len(rest)} more genes: {', '.join(names)}\n"
                    "Click to open all the genes in GeneTrek"
                ),
                "href": search_url(scorer.resolve(i) or i for i in idents),
                "score": -1,
                "overflow": True,
            }
        )
    return badges


def badge_labels(badges: List[Dict[str, Any]]) -> List[str]:
    """The gene names of a badge list, without the "+N genes" badge."""
    return [b["label"] for b in badges if not b.get("overflow")]


def gene_badge_html(badge: Dict[str, Any], *, extra_style: str = "") -> str:
    """One badge as escaped HTML for ``ui.html``: a link when it has one."""
    style = (
        f"background-color: {badge['color']}; color: {badge.get('textColor', 'black')};"
        " padding: 2px 6px; border-radius: 4px; font-size: 12px;"
        " text-decoration: none; display: inline-block;"
    )
    if badge.get("borderColor"):
        style += f" border: 2px solid {badge['borderColor']};"
    style += extra_style
    title = html.escape(badge.get("tooltip") or "", quote=True)
    label = html.escape(badge["label"])
    if badge.get("href"):
        href = html.escape(badge["href"], quote=True)
        return (
            f'<a class="q-badge gene-badge" href="{href}" target="_blank"'
            f' rel="noopener noreferrer" style="{style}" title="{title}">{label}</a>'
        )
    return f'<span class="q-badge gene-badge" style="{style}" title="{title}">{label}</span>'


def gene_badges_html(badges: List[Dict[str, Any]]) -> str:
    inner = "".join(gene_badge_html(b) for b in badges)
    return f'<div style="display: flex; flex-wrap: wrap; gap: 4px;">{inner}</div>'


def tier_legend_html(scorer: GeneScorer, *, exonic_border: bool = False) -> str:
    """A legend of the user's tiers (and the exonic border, for SVs)."""
    parts = []
    for tier, sets in scorer.sets_by_tier().items():
        if not sets:
            continue
        badge = {
            "label": f"Tier {tier}",
            "color": scorer.tier_color(tier),
            "textColor": text_color_for(scorer.tier_color(tier)),
            "tooltip": ", ".join(gs.label for gs in sets),
        }
        parts.append(gene_badge_html(badge))
    neutral = {
        "label": "none",
        "color": NEUTRAL_COLOR,
        "textColor": "black",
        "tooltip": "In no tiered geneset. Set the tiers in your profile.",
    }
    parts.append(gene_badge_html(neutral))
    if exonic_border:
        exonic = dict(
            neutral, label="exonic", borderColor="black", tooltip="Exonic overlap"
        )
        parts.append(gene_badge_html(exonic))
    inner = "".join(parts)
    return (
        '<div style="display: flex; flex-wrap: wrap; gap: 4px; align-items: center;'
        ' font-size: 12px;"><a href="/profile" style="color: inherit;"'
        f' title="Set your tiers in your profile">Gene tiers:</a>{inner}</div>'
    )
