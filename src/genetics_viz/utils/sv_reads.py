"""Read evidence around an SV, samplot style: discordant pairs, split reads, clips.

Each discordant read pair is drawn from its leftmost to its rightmost aligned
base, at a height equal to its insert size. Its class comes from the
orientation of its two reads (the left read's strand, then the right one's):

- ``+-`` is the normal orientation. It is a *deletion* signature when the
  insert is larger than the window's own proper pairs allow.
- ``-+`` (pair pointing outward) is a *duplication* signature.
- ``++`` and ``--`` are the two *inversion* signatures.
- A mate on another chromosome is *interchromosomal*.

Concordant pairs are not kept. Their insert sizes only give the expected
band (``insert_lo``-``insert_hi``) that the panel shades.

A split read (the ``SA`` tag of a primary alignment) keeps its two aligned
parts. Its height is the gap between them. Soft clips are counted per clipped
position; only positions where several reads are clipped are kept, since that
is what marks a breakpoint. The two together suggest exact breakpoints
(``utils/sv_breakpoints.py``).

Reads under ``min_mapq`` are left out unless asked for; they are then kept
and flagged ``low``, to be drawn faded. Repeat-flanked events (22q11.2) have
most of their evidence there. Read depth at ``coverage_bin_bp`` is computed
from the reads that pass ``min_mapq`` only, as mosdepth's ``-Q`` does.

Blocking (CRAM decoding over the network mount): call it off-thread. The
CRAMs need the exact reference they were aligned to (``reference_fasta`` in
the app config). Without it nothing is opened, since htslib would otherwise
download the reference sequences by MD5.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genetics_viz.utils.cnv_coverage import CNV_CONFIG
from genetics_viz.utils.sv_view import Window

READ_CONFIG: Dict[str, Any] = CNV_CONFIG.get("read_evidence", {})

PAIR_CLASSES = ("deletion", "duplication", "inversion_pp", "inversion_mm")
SPLIT_CLASSES = ("split", "split_inv")

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_SOFT_CLIP = 4


@dataclass(frozen=True)
class Pair:
    """A discordant pair: its span, its two read bodies and its insert."""

    x1: int
    x2: int
    y: int
    cls: str
    left_end: int  # end of the left read's body
    right_start: int  # start of the right read's body
    low: bool = False  # a read or its mate under min_mapq


@dataclass(frozen=True)
class SplitRead:
    """The two aligned parts of a split read, left to right, and their junction.

    A part joins the other at its clipped end: ``ja`` and ``jb`` are those
    ends, with ``ja_side`` / ``jb_side`` named like a :class:`SoftClip`'s
    side (``"right"``: clipped after the junction, ``"left"``: before it).
    A deletion joins ``a2`` to ``b1``, a tandem duplication ``b2`` to ``a1``.
    """

    a1: int
    a2: int
    b1: int
    b2: int
    y: int
    cls: str
    low: bool = False  # either part under min_mapq
    ja: int = 0
    ja_side: str = "right"
    jb: int = 0
    jb_side: str = "left"


@dataclass(frozen=True)
class SoftClip:
    pos: int
    #: ``"left"``: the bases before ``pos`` are clipped; ``"right"``: after.
    side: str
    count: int
    low: bool = False  # clipped reads under min_mapq, counted apart


@dataclass(frozen=True)
class Interchrom:
    x: int
    mate: str


@dataclass
class WindowReads:
    window: Window
    pairs: List[Pair] = field(default_factory=list)
    splits: List[SplitRead] = field(default_factory=list)
    clips: List[SoftClip] = field(default_factory=list)
    interchrom: List[Interchrom] = field(default_factory=list)
    insert_median: Optional[float] = None
    #: The expected insert range of concordant pairs (``insert_band`` percentiles).
    insert_lo: Optional[float] = None
    insert_hi: Optional[float] = None
    insert_max: float = 0.0
    n_concordant: int = 0
    #: ``[bin centre, mean depth]`` at ``coverage_bin_bp``, MAPQ-filtered.
    depth: List[List[float]] = field(default_factory=list)


def _mad(values: Sequence[float], centre: float) -> float:
    return median(abs(v - centre) for v in values) if values else 0.0


def _orientation(read) -> str:
    """``+-`` style orientation of a same-chromosome pair, left read first."""
    read_is_left = read.reference_start <= read.next_reference_start
    left_rev = read.is_reverse if read_is_left else read.mate_is_reverse
    right_rev = read.mate_is_reverse if read_is_left else read.is_reverse
    return ("-" if left_rev else "+") + ("-" if right_rev else "+")


def _ref_length(cigar: str) -> int:
    """Reference bases a CIGAR string covers."""
    return sum(int(n) for n, op in _CIGAR_RE.findall(cigar) if op in "MDN=X")


def _clipped_at_end(ops: Sequence[Tuple[str, int]]) -> bool:
    """Whether an alignment's longer clip is after it (else before it)."""
    lead = ops[0][1] if ops and ops[0][0] in "SH" else 0
    trail = ops[-1][1] if ops and ops[-1][0] in "SH" else 0
    return trail >= lead


def _junction(start: int, end: int, at_end: bool) -> Tuple[int, str]:
    return (end, "right") if at_end else (start, "left")


def _first_sa(read) -> Optional[Tuple[str, int, str, int, int, str]]:
    """``(chrom, 0-based pos, strand, ref length, MAPQ, CIGAR)`` of the first SA part."""
    try:
        sa = read.get_tag("SA")
    except KeyError:
        return None
    first = str(sa).split(";")[0].split(",")
    if len(first) < 5:
        return None
    try:
        return (
            first[0],
            int(first[1]) - 1,
            first[2],
            _ref_length(first[3]),
            int(first[4]),
            first[3],
        )
    except ValueError:
        return None


def _mate_mapq(read) -> Optional[int]:
    try:
        return int(read.get_tag("MQ"))
    except KeyError:
        return None


def collect_window(
    alignment_file,
    window: Window,
    *,
    cfg: Optional[Dict[str, Any]] = None,
    include_low_mapq: bool = False,
) -> WindowReads:
    """Classify the discordant pairs, split reads and soft clips of one window."""
    cfg = cfg or READ_CONFIG
    min_mapq = int(cfg.get("min_mapq", 1))
    min_clip = int(cfg.get("min_clip_bp", 20))
    bin_bp = max(1, int(cfg.get("coverage_bin_bp", 25)))
    n_bins = max(1, -(-window.width // bin_bp))
    depth = [0.0] * n_bins
    out = WindowReads(window=window)

    same_chrom: Dict[str, Tuple[int, int, int, str, int, int, bool]] = {}
    proper_inserts: List[int] = []
    clips: Counter = Counter()
    seen_split: set[str] = set()
    seen_inter: set[str] = set()

    try:
        reads = alignment_file.fetch(window.chrom, window.start, window.end)
    except ValueError:
        return out  # contig not in this file
    for read in reads:
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.is_duplicate
            or read.is_qcfail
        ):
            continue
        low = read.mapping_quality < min_mapq
        if low and not include_low_mapq:
            continue
        r_start = read.reference_start
        r_end = read.reference_end or r_start
        if not low:
            lo = max(r_start, window.start) - window.start
            hi = min(r_end, window.end) - window.start
            for b in range(lo // bin_bp, min(n_bins, -(-hi // bin_bp))):
                b_lo = window.start + b * bin_bp
                overlap = min(r_end, b_lo + bin_bp) - max(r_start, b_lo)
                depth[b] += overlap / bin_bp
        cigar = read.cigartuples or []
        if cigar and cigar[0][0] == _SOFT_CLIP and cigar[0][1] >= min_clip:
            if window.start <= r_start < window.end:
                clips[(r_start, "left", low)] += 1
        if cigar and cigar[-1][0] == _SOFT_CLIP and cigar[-1][1] >= min_clip:
            if window.start <= r_end < window.end:
                clips[(r_end, "right", low)] += 1

        name = read.query_name
        sa = _first_sa(read)
        if sa is not None and name not in seen_split:
            seen_split.add(name)
            sa_chrom, sa_pos, sa_strand, sa_len, sa_mapq, sa_cigar = sa
            if sa_chrom == read.reference_name:
                own_ops = [("MIDNSHP=XB"[op], n) for op, n in cigar]
                sa_ops = [(op, int(n)) for n, op in _CIGAR_RE.findall(sa_cigar)]
                own = (r_start, r_end, _clipped_at_end(own_ops))
                other = (sa_pos, sa_pos + max(1, sa_len), _clipped_at_end(sa_ops))
                (a1, a2, a_end), (b1, b2, b_end) = sorted((own, other))
                ja, ja_side = _junction(a1, a2, a_end)
                jb, jb_side = _junction(b1, b2, b_end)
                own_strand = "-" if read.is_reverse else "+"
                out.splits.append(
                    SplitRead(
                        a1,
                        a2,
                        b1,
                        b2,
                        max(1, abs(jb - ja)),
                        "split" if sa_strand == own_strand else "split_inv",
                        low or sa_mapq < min_mapq,
                        ja,
                        ja_side,
                        jb,
                        jb_side,
                    )
                )
            else:
                out.interchrom.append(
                    Interchrom(r_start, f"{sa_chrom}:{sa_pos + 1} (SA)")
                )

        if not read.is_paired or read.mate_is_unmapped:
            continue
        if read.next_reference_name != read.reference_name:
            if name not in seen_inter:
                seen_inter.add(name)
                out.interchrom.append(
                    Interchrom(
                        r_start,
                        f"{read.next_reference_name}:{read.next_reference_start + 1}",
                    )
                )
            continue
        if name in same_chrom:
            continue
        mate_start = read.next_reference_start
        read_len = read.infer_read_length() or (r_end - r_start)
        if r_start <= mate_start:
            x1, left_end = r_start, r_end
            right_start, x2 = mate_start, mate_start + read_len
        else:
            x1, left_end = mate_start, mate_start + read_len
            right_start, x2 = r_start, r_end
        x2 = max(x2, left_end)
        insert = abs(read.template_length) or (x2 - x1)
        orient = _orientation(read)
        mate_mapq = _mate_mapq(read)
        pair_low = low or (mate_mapq is not None and mate_mapq < min_mapq)
        same_chrom[name] = (x1, x2, insert, orient, left_end, right_start, pair_low)
        if read.is_proper_pair and orient == "+-" and insert <= 5000 and not low:
            proper_inserts.append(insert)

    if len(proper_inserts) >= 20:
        centre = float(median(proper_inserts))
        spread = _mad(proper_inserts, centre)
        lo_pct, hi_pct = cfg.get("insert_band", [2.5, 97.5])
        out.insert_median = centre
        out.insert_lo, out.insert_hi = (
            float(v) for v in np.percentile(proper_inserts, [lo_pct, hi_pct])
        )
        out.insert_max = centre + float(cfg.get("discordant_insert_mads", 6)) * max(
            spread, 10.0
        )
    else:
        out.insert_max = float(cfg.get("default_insert_max", 1000))

    for x1, x2, insert, orient, left_end, right_start, pair_low in same_chrom.values():
        if orient == "+-":
            if insert <= out.insert_max:
                out.n_concordant += 1
                continue
            cls = "deletion"
        elif orient == "-+":
            cls = "duplication"
        elif orient == "++":
            cls = "inversion_pp"
        else:
            cls = "inversion_mm"
        out.pairs.append(
            Pair(x1, x2, max(1, insert), cls, left_end, right_start, pair_low)
        )

    min_reads = int(cfg.get("min_clip_reads", 2))
    out.clips = sorted(
        (
            SoftClip(pos, side, n, low)
            for (pos, side, low), n in clips.items()
            if n >= min_reads
        ),
        key=lambda c: c.pos,
    )
    out.depth = [
        [window.start + (i + 0.5) * bin_bp, round(d, 2)] for i, d in enumerate(depth)
    ]
    return out


def load_read_evidence(
    crams: Sequence[Optional[Path]],
    windows: Sequence[Window],
    reference: Path,
    cfg: Optional[Dict[str, Any]] = None,
    include_low_mapq: bool = False,
) -> List[Optional[List[WindowReads]]]:
    """Per sample (``None`` without a CRAM), per window, the read evidence."""
    import pysam

    results: List[Optional[List[WindowReads]]] = []
    for cram in crams:
        if cram is None:
            results.append(None)
            continue
        with pysam.AlignmentFile(
            str(cram), "rc", reference_filename=str(reference)
        ) as af:
            results.append(
                [
                    collect_window(af, w, cfg=cfg, include_low_mapq=include_low_mapq)
                    for w in windows
                ]
            )
    return results
