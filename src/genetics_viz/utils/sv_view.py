"""Which read windows an SV is drawn in: one around it, or one per breakpoint.

A short SV fits in one window with its flanks, and both breakpoints are seen in
one view. A long one does not fit igv.js's alignment visibility window, so each
breakpoint gets its own. The threshold is that visibility window itself
(``read_view`` in ``config/cnv.yaml``), which is also set on the igv.js CRAM
tracks, so the switch and igv.js agree.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from genetics_viz.utils.cnv_coverage import CNV_CONFIG

READ_VIEW_CONFIG: Dict[str, Any] = CNV_CONFIG.get("read_view", {})

MODES = ("auto", "single", "split")


@dataclass(frozen=True)
class Window:
    chrom: str
    start: int
    end: int

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def width(self) -> int:
        return self.end - self.start


def visibility_window_bp(cfg: Optional[Dict[str, Any]] = None) -> int:
    return int((cfg or READ_VIEW_CONFIG).get("visibility_window_bp", 30000))


def auto_mode(
    sv_type: str, start: int, end: int, cfg: Optional[Dict[str, Any]] = None
) -> str:
    """``"single"`` or ``"split"`` for an SV, by type and size."""
    cfg = cfg or READ_VIEW_CONFIG
    kind = (sv_type or "").lower()
    if kind == "ins":
        return "single"
    if kind in ("bnd", "tra"):
        return "split"
    span = max(1, int(end) - int(start))
    expansion = float(cfg.get("expansion", 1.6))
    return "single" if span * expansion <= visibility_window_bp(cfg) else "split"


def read_view_windows(
    sv_type: str,
    chrom: str,
    start: int,
    end: int,
    *,
    mode: str = "auto",
    mate: Optional[Tuple[str, int]] = None,
    precise: bool = True,
    cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[str, List[Window]]:
    """The mode actually used and its windows, left to right.

    ``mate`` is the second breakpoint of a translocation, on another
    chromosome; without it the second window is at ``end`` on ``chrom``.
    ``precise`` is False for a caller whose breakpoints are approximate
    (WisecondorX bins): each breakpoint window is then wider.
    """
    cfg = cfg or READ_VIEW_CONFIG
    start, end = int(start), int(end)
    used = auto_mode(sv_type, start, end, cfg) if mode == "auto" else mode
    flank = int(
        cfg.get("flank_bp", 1500)
        if precise
        else cfg.get("flank_bp_imprecise", cfg.get("flank_bp", 1500))
    )

    if used == "single":
        expansion = float(cfg.get("expansion", 1.6))
        span = max(1, end - start)
        width = max(int(span * expansion), 2 * flank)
        centre = (start + end) // 2
        lo = max(0, centre - width // 2)
        return used, [Window(chrom, lo, lo + width)]

    second_chrom, second_pos = mate if mate else (chrom, end)
    return used, [
        Window(chrom, max(0, start - flank), start + flank),
        Window(second_chrom, max(0, second_pos - flank), second_pos + flank),
    ]
