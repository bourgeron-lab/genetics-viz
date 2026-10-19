"""Read-pair and split-read evidence around an SV, samplot style.

Each read pair becomes a horizontal segment from its leftmost to its rightmost
aligned base, at a height equal to its insert size, and is classed by the
orientation of its two reads (the left read's strand, then the right one's):

- ``+-`` is the normal orientation, a *deletion* signature when the insert is
  larger than the window's own proper pairs allow;
- ``-+`` (pair pointing outward) is a *duplication* signature;
- ``++`` and ``--`` are the two *inversion* signatures;
- a mate on another chromosome is *interchromosomal*.

Split reads (the ``SA`` tag of a primary alignment) join their two aligned
parts. Concordant pairs are downsampled; everything else is kept whole.

Blocking (CRAM decoding over the network mount): call it off-thread. The
CRAMs need the exact reference they were aligned to (``reference_fasta`` in
the app config); without it nothing is opened, since htslib would otherwise try
to download the reference sequences by MD5.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any, Dict, List, Optional, Sequence

from genetics_viz.utils.cnv_coverage import CNV_CONFIG
from genetics_viz.utils.sv_view import Window

READ_CONFIG: Dict[str, Any] = CNV_CONFIG.get("read_evidence", {})

#: Classes in drawing order (concordant underneath).
PAIR_CLASSES = (
    "concordant",
    "deletion",
    "duplication",
    "inversion_pp",
    "inversion_mm",
)
SPLIT_CLASSES = ("split", "split_inv")


@dataclass(frozen=True)
class Segment:
    x1: int
    x2: int
    y: int
    cls: str


@dataclass(frozen=True)
class Interchrom:
    x: int
    mate: str


@dataclass
class WindowReads:
    window: Window
    pairs: List[Segment] = field(default_factory=list)
    splits: List[Segment] = field(default_factory=list)
    interchrom: List[Interchrom] = field(default_factory=list)
    #: ``(bin start, mean depth)`` at ``coverage_bin_bp``.
    coverage: List[List[float]] = field(default_factory=list)
    insert_median: Optional[float] = None
    insert_max: float = 0.0
    n_concordant: int = 0


def _mad(values: Sequence[float], centre: float) -> float:
    return median(abs(v - centre) for v in values) if values else 0.0


def _orientation(read) -> str:
    """``+-`` style orientation of a same-chromosome pair, left read first."""
    read_is_left = read.reference_start <= read.next_reference_start
    left_rev = read.is_reverse if read_is_left else read.mate_is_reverse
    right_rev = read.mate_is_reverse if read_is_left else read.is_reverse
    return ("-" if left_rev else "+") + ("-" if right_rev else "+")


def _first_sa(read) -> Optional[tuple[str, int, str]]:
    """``(chrom, 0-based pos, strand)`` of the first supplementary part."""
    try:
        sa = read.get_tag("SA")
    except KeyError:
        return None
    first = str(sa).split(";")[0].split(",")
    if len(first) < 3:
        return None
    try:
        return first[0], int(first[1]) - 1, first[2]
    except ValueError:
        return None


def collect_window(
    alignment_file,
    window: Window,
    *,
    cfg: Optional[Dict[str, Any]] = None,
    seed: int = 0,
) -> WindowReads:
    """Classify the read pairs and split reads of one window."""
    cfg = cfg or READ_CONFIG
    min_mapq = int(cfg.get("min_mapq", 1))
    bin_bp = max(1, int(cfg.get("coverage_bin_bp", 25)))
    n_bins = max(1, -(-window.width // bin_bp))
    depth = [0.0] * n_bins
    out = WindowReads(window=window)

    same_chrom: Dict[str, tuple[int, int, int, str]] = {}  # name -> x1, x2, ins, orient
    proper_inserts: List[int] = []
    seen_split: set[str] = set()
    seen_inter: set[str] = set()

    chrom = window.chrom
    try:
        reads = alignment_file.fetch(chrom, window.start, window.end)
    except ValueError:
        return out  # contig not in this file
    for read in reads:
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.is_duplicate
            or read.is_qcfail
            or read.mapping_quality < min_mapq
        ):
            continue
        r_start, r_end = (
            read.reference_start,
            read.reference_end or read.reference_start,
        )
        lo = max(r_start, window.start) - window.start
        hi = min(r_end, window.end) - window.start
        for b in range(lo // bin_bp, min(n_bins, -(-hi // bin_bp))):
            b_lo, b_hi = window.start + b * bin_bp, window.start + (b + 1) * bin_bp
            depth[b] += (min(r_end, b_hi) - max(r_start, b_lo)) / bin_bp

        name = read.query_name
        sa = _first_sa(read)
        if sa is not None and name not in seen_split:
            seen_split.add(name)
            sa_chrom, sa_pos, sa_strand = sa
            if sa_chrom == read.reference_name:
                own_strand = "-" if read.is_reverse else "+"
                x1, x2 = sorted((r_start, sa_pos))
                x2 = max(x2, r_end)
                out.splits.append(
                    Segment(
                        x1,
                        x2,
                        max(1, x2 - x1),
                        "split" if sa_strand == own_strand else "split_inv",
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
        x1 = min(r_start, mate_start)
        x2 = max(r_end, mate_start + read_len)
        insert = abs(read.template_length) or (x2 - x1)
        orient = _orientation(read)
        same_chrom[name] = (x1, x2, insert, orient)
        if read.is_proper_pair and orient == "+-" and insert <= 5000:
            proper_inserts.append(insert)

    if len(proper_inserts) >= 20:
        centre = float(median(proper_inserts))
        spread = _mad(proper_inserts, centre)
        out.insert_median = centre
        out.insert_max = centre + float(cfg.get("discordant_insert_mads", 6)) * max(
            spread, 10.0
        )
    else:
        out.insert_max = float(cfg.get("default_insert_max", 1000))

    concordant: List[Segment] = []
    for x1, x2, insert, orient in same_chrom.values():
        if orient == "+-":
            cls = "deletion" if insert > out.insert_max else "concordant"
        elif orient == "-+":
            cls = "duplication"
        elif orient == "++":
            cls = "inversion_pp"
        else:
            cls = "inversion_mm"
        seg = Segment(x1, x2, max(1, insert), cls)
        (concordant if cls == "concordant" else out.pairs).append(seg)

    out.n_concordant = len(concordant)
    cap = int(cfg.get("max_concordant_per_window", 1500))
    if len(concordant) > cap:
        concordant = random.Random(seed).sample(concordant, cap)
    out.pairs = sorted(concordant, key=lambda s: s.x1) + out.pairs
    out.coverage = [
        [window.start + (i + 0.5) * bin_bp, round(d, 2)] for i, d in enumerate(depth)
    ]
    return out


def load_read_evidence(
    crams: Sequence[Optional[Path]],
    windows: Sequence[Window],
    reference: Path,
    cfg: Optional[Dict[str, Any]] = None,
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
            results.append([collect_window(af, w, cfg=cfg) for w in windows])
    return results
