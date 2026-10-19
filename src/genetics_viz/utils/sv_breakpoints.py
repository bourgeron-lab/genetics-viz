"""Breakpoints suggested by the reads and the coverage around an SV.

The SV's start and end are judged separately. Evidence left of the SV's middle
speaks for its start, and evidence right of it for its end. For each side:

- **Exact position** from soft-clip clusters and split-read junctions, scored
  as clips + 2 × split reads. Evidence whose orientation does not fit the SV
  type counts for ``mismatch_weight`` of that. A deletion's start has reads
  clipped after it and its end reads clipped before it; a tandem duplication
  is the reverse. Inversions, insertions and translocations accept both.
- **Interval** from discordant pairs of the type's class. A ``+`` read lies
  before its breakpoint, within the largest expected insert of it; a ``-``
  read lies after it. The pairs' constraints intersect.
- **Copy-number check**, for copy-number-changing types only (``del``,
  ``dup``): the coverage must step the right way at the exact position.
  Inversions and other balanced events are judged on reads alone.
- **Coverage step**, for those same types: the strongest step in the
  coverage, as an approximate position when the reads give none (repeat-
  mediated events, whose reads have MAPQ 0).
- **Confidence**: strong, moderate or weak, from the above.

The coverage is the 25 bp read depth (``2 × depth / autosomal median``) where
it exists, else the 250 bp CN. Bins the control panel masks are left out:
repeats lose depth in everyone, and would otherwise read as steps.

Pure: works on the read evidence, CN tracks and depth already loaded.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from genetics_viz.utils.cnv_coverage import CNTrack
from genetics_viz.utils.sv_reads import READ_CONFIG, WindowReads

#: Pair classes that support each SV type's breakpoints.
_PAIR_CLASSES = {
    "del": {"deletion"},
    "dup": {"duplication"},
    "inv": {"inversion_pp", "inversion_mm"},
}
#: Clip side expected at (start, end), per type; absent = either is fine.
_EXPECTED_SIDES = {"del": ("right", "left"), "dup": ("left", "right")}
#: Split-read strand class expected, per type.
_EXPECTED_SPLIT = {"del": "split", "dup": "split", "inv": "split_inv"}
#: CN change from left to right of the (start, end), per type.
_CN_DIRECTION = {"del": (-1, 1), "dup": (1, -1)}


def changes_copy_number(sv_type: str) -> bool:
    return (sv_type or "").lower() in _CN_DIRECTION


@dataclass(frozen=True)
class CoverageStep:
    pos: int
    delta: float  # CN right minus left
    resolution: int  # bp


@dataclass
class SideSuggestion:
    side: str  # "start" | "end"
    pos: Optional[int] = None
    clips: int = 0
    splits: int = 0
    mismatched: int = 0
    score: float = 0.0
    confidence: Optional[str] = None  # "strong" | "moderate" | "weak"
    #: "supported" | "contradicted" | "unknown" | "n/a" (balanced type)
    coverage: str = "n/a"
    coverage_delta: Optional[float] = None
    interval: Optional[Tuple[int, int]] = None
    interval_pairs: int = 0
    outside_interval: bool = False
    step: Optional[CoverageStep] = None

    def position(self) -> Tuple[Optional[int], bool]:
        """The position to curate with, and whether it is approximate."""
        if self.pos is not None:
            return self.pos, False
        if self.step is not None:
            return self.step.pos, True
        return None, False


# ---------------------------------------------------------------------------
# Exact position from clips and split reads
# ---------------------------------------------------------------------------


def _exact(
    evidence: Sequence[Tuple[int, str, float]], tolerance: int
) -> Optional[Tuple[int, int, int, int, float]]:
    """Best ``(pos, clips, splits, mismatched, score)`` among weighted items.

    Items are ``(pos, "clip" | "split" | "clip?" | "split?", weight)``, the
    ``?`` kinds being orientation mismatches. Items within ``tolerance`` bp
    count together, at the position with the most weight.
    """
    if not evidence:
        return None
    items = sorted(evidence)
    groups: List[List[Tuple[int, str, float]]] = [[items[0]]]
    for item in items[1:]:
        if item[0] - groups[-1][-1][0] <= tolerance:
            groups[-1].append(item)
        else:
            groups.append([item])
    best = None
    for group in groups:
        score = sum(w for _, _, w in group)
        by_pos: Dict[int, float] = {}
        for p, _, w in group:
            by_pos[p] = by_pos.get(p, 0.0) + w
        pos = max(by_pos, key=by_pos.get)
        clips = sum(1 for _, k, _ in group if k == "clip")
        splits = sum(1 for _, k, _ in group if k == "split")
        mismatched = sum(1 for _, k, _ in group if k.endswith("?"))
        if best is None or score > best[4]:
            best = (pos, clips, splits, mismatched, score)
    return best


def _read_evidence(
    reads: Sequence[WindowReads],
    sv_type: str,
    middle: float,
    cfg: Dict[str, Any],
) -> Tuple[List[Tuple[int, str, float]], List[Tuple[int, str, float]]]:
    mismatch = float(cfg.get("mismatch_weight", 0.25))
    expected_sides = _EXPECTED_SIDES.get(sv_type)
    expected_split = _EXPECTED_SPLIT.get(sv_type)
    sides: Tuple[List, List] = ([], [])

    def add(pos: int, side: str, kind: str, weight: float, cls_ok: bool) -> None:
        index = 0 if pos < middle else 1
        ok = cls_ok and (expected_sides is None or side == expected_sides[index])
        sides[index].append(
            (pos, kind if ok else kind + "?", weight * (1 if ok else mismatch))
        )

    for cell in reads:
        for clip in cell.clips:
            if clip.low:
                continue
            for _ in range(clip.count):
                add(clip.pos, clip.side, "clip", 1.0, True)
        for split in cell.splits:
            if split.low:
                continue
            cls_ok = expected_split is None or split.cls == expected_split
            add(split.ja, split.ja_side, "split", 2.0, cls_ok)
            add(split.jb, split.jb_side, "split", 2.0, cls_ok)
    return sides


# ---------------------------------------------------------------------------
# Interval from discordant pairs
# ---------------------------------------------------------------------------


def _pair_constraints(
    reads: Sequence[WindowReads], sv_type: str, middle: float
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """``(lo, hi)`` bounds each supporting read puts on its side's breakpoint."""
    classes = _PAIR_CLASSES.get(sv_type)
    sides: Tuple[List, List] = ([], [])
    for cell in reads:
        reach = cell.insert_hi or cell.insert_max
        if not reach:
            continue
        for p in cell.pairs:
            if p.low or (classes is not None and p.cls not in classes):
                continue
            strands = {
                "deletion": "+-",
                "duplication": "-+",
                "inversion_pp": "++",
                "inversion_mm": "--",
            }[p.cls]
            for (r1, r2), strand in (
                ((p.x1, p.left_end), strands[0]),
                ((p.right_start, p.x2), strands[1]),
            ):
                index = 0 if r1 < middle else 1
                if strand == "+":
                    sides[index].append((r2, int(r1 + reach)))
                else:
                    sides[index].append((int(r2 - reach), r1))
    return sides


def _intersect(bounds: Sequence[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
    if not bounds:
        return None
    lo = max(b[0] for b in bounds)
    hi = min(b[1] for b in bounds)
    if lo <= hi:
        return lo, hi
    # Pairs disagree (a second event, mismapped reads): the middle half of the
    # lower and upper bounds, so one stray pair does not empty the interval.
    los = sorted(b[0] for b in bounds)
    his = sorted(b[1] for b in bounds)
    lo, hi = los[len(los) // 4], his[(3 * len(his)) // 4]
    return (lo, hi) if lo <= hi else None


# ---------------------------------------------------------------------------
# Coverage: profile, check and step
# ---------------------------------------------------------------------------


def coverage_profile(
    depth: Optional[Sequence[Sequence[float]]], cn: Optional[CNTrack]
) -> List[Tuple[float, float]]:
    """``(pos, CN-like value)``, finest available, without panel-masked bins."""
    masked_bins: List[Tuple[float, float]] = []
    if cn is not None and cn.positions and cn.bin_bp:
        half = cn.bin_bp / 2
        masked_bins = [
            (p - half, p + half) for p, m in zip(cn.positions, cn.masked) if m
        ]

    def masked(pos: float) -> bool:
        return any(lo <= pos < hi for lo, hi in masked_bins)

    if depth:
        return [(p, v) for p, v in depth if v is not None and not masked(p)]
    if cn is None:
        return []
    return [
        (p, v)
        for p, v, m in zip(cn.positions, cn.values, cn.masked)
        if v is not None and not m
    ]


def _mean(values: Sequence[float]) -> Optional[float]:
    return sum(values) / len(values) if values else None


def _check(
    profile: Sequence[Tuple[float, float]],
    pos: int,
    direction: int,
    cfg: Dict[str, Any],
) -> Tuple[str, Optional[float]]:
    width = int(cfg.get("cn_check_bp", 1000))
    left = [v for p, v in profile if pos - width <= p < pos]
    right = [v for p, v in profile if pos < p <= pos + width]
    if len(left) < 3 or len(right) < 3:
        return "unknown", None
    delta = _mean(right) - _mean(left)
    min_step = float(cfg.get("min_cn_step", 0.5))
    supported = delta * direction >= min_step
    return ("supported" if supported else "contradicted"), round(delta, 2)


def _step(
    profile: Sequence[Tuple[float, float]], direction: int, cfg: Dict[str, Any]
) -> Optional[CoverageStep]:
    """The strongest step the right way, as a two-segment mean fit."""
    min_seg = int(cfg.get("min_step_bins", 6))
    min_step = float(cfg.get("min_cn_step", 0.5))
    n = len(profile)
    if n < 2 * min_seg:
        return None
    values = [v for _, v in profile]
    prefix = [0.0]
    for v in values:
        prefix.append(prefix[-1] + v)
    best: Optional[Tuple[float, int, float]] = None
    for k in range(min_seg, n - min_seg + 1):
        left = prefix[k] / k
        right = (prefix[n] - prefix[k]) / (n - k)
        delta = right - left
        if delta * direction < min_step:
            continue
        stat = abs(delta) * math.sqrt(k * (n - k) / n)
        if best is None or stat > best[0]:
            best = (stat, k, delta)
    if best is None:
        return None
    _, k, delta = best
    pos = int(round((profile[k - 1][0] + profile[k][0]) / 2))
    resolution = int(round(abs(profile[k][0] - profile[k - 1][0]))) or 1
    return CoverageStep(pos, round(delta, 2), resolution)


# ---------------------------------------------------------------------------
# Together
# ---------------------------------------------------------------------------


def _confidence(s: SideSuggestion, cn_type: bool) -> str:
    if s.coverage == "contradicted" or s.outside_interval:
        return "weak"
    good_splits = s.splits
    good_clips = s.clips
    if (
        good_splits >= 1
        and good_clips >= 2
        and (not cn_type or s.coverage == "supported")
    ):
        return "strong"
    if good_clips >= 5 or good_splits >= 2:
        return "moderate"
    return "weak"


def suggest_breakpoints(
    reads: Sequence[Optional[WindowReads]],
    cns: Sequence[Optional[CNTrack]],
    depths: Sequence[Optional[Sequence[Sequence[float]]]],
    sv_type: str,
    start: int,
    end: int,
    cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[SideSuggestion, SideSuggestion]:
    """Start and end suggestions for one sample, from its windows' evidence."""
    cfg = cfg or READ_CONFIG
    sv_type = (sv_type or "").lower()
    middle = (int(start) + int(end)) / 2
    present = [r for r in reads if r is not None]
    cn_type = changes_copy_number(sv_type)

    evidence = _read_evidence(present, sv_type, middle, cfg)
    constraints = _pair_constraints(present, sv_type, middle)
    profile = sorted(
        point for cn, depth in zip(cns, depths) for point in coverage_profile(depth, cn)
    )
    tolerance = int(cfg.get("breakpoint_tolerance_bp", 3))
    min_score = float(cfg.get("min_breakpoint_score", 3))
    min_pairs = int(cfg.get("min_interval_pairs", 2))

    out = []
    for index, side in enumerate(("start", "end")):
        s = SideSuggestion(side)
        best = _exact(evidence[index], tolerance)
        if best is not None and best[4] >= min_score:
            s.pos, s.clips, s.splits, s.mismatched, s.score = best
            s.score = round(s.score, 2)
        # A pair spanning both breakpoint windows is seen from each: once.
        bounds = sorted(set(constraints[index]))
        if len(bounds) >= min_pairs:
            s.interval = _intersect(bounds)
            s.interval_pairs = len(bounds)
        side_profile = [pt for pt in profile if (pt[0] < middle) == (index == 0)]
        if cn_type:
            direction = _CN_DIRECTION[sv_type][index]
            if s.pos is not None:
                s.coverage, s.coverage_delta = _check(
                    side_profile, s.pos, direction, cfg
                )
            else:
                s.coverage = "unknown"
                s.step = _step(side_profile, direction, cfg)
        if s.pos is not None:
            if s.interval is not None:
                slack = tolerance + 1
                s.outside_interval = not (
                    s.interval[0] - slack <= s.pos <= s.interval[1] + slack
                )
            s.confidence = _confidence(s, cn_type)
        out.append(s)
    return out[0], out[1]
