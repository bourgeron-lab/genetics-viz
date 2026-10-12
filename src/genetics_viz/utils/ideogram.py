"""SVG building blocks for chromosome ideograms (hg38).

Shared by the search statistics ideogram and the ROH genome view. Both draw
their own overlays; what they have in common is the chromosome itself - its
Giemsa-stained bands and outline - and the Mb scale under it. Functions here
return lists of SVG element strings for the caller to join into one document.
"""

from typing import List

from .cytobands import CHROM_SIZES_MB, CYTOBANDS, GIESTAIN_COLORS


def chrom_length_mb(chrom: str) -> float:
    """Exact chromosome length in Mb, from the last cytoband.

    ``CHROM_SIZES_MB`` is rounded to the nearest Mb - chr2 is 242.19 Mb but
    listed as 242 - so a position near the end of a chromosome would overflow
    it. Falls back to the rounded size when the band table lacks the chromosome.
    """
    bands = CYTOBANDS.get(chrom)
    if bands:
        return max(b["end"] for b in bands)
    return float(CHROM_SIZES_MB.get(chrom, 0))


def chromosome_svg(
    chrom: str, x0: float, y: float, px_per_mb: float, height: float
) -> List[str]:
    """The banded body of one chromosome, with a rounded outline.

    Args:
        chrom: Normalised chromosome name (``"1"`` ... ``"22"``, ``"X"``).
        x0: Left edge, in SVG units.
        y: Top edge, in SVG units.
        px_per_mb: Horizontal scale.
        height: Bar height, in SVG units.
    """
    parts: List[str] = []
    for band in CYTOBANDS.get(chrom, []):
        bx = x0 + band["start"] * px_per_mb
        bw = max((band["end"] - band["start"]) * px_per_mb, 0.5)
        color = GIESTAIN_COLORS.get(band["stain"], "#e5e7eb")
        parts.append(
            f'<rect x="{bx:.1f}" y="{y}" width="{bw:.1f}" height="{height}" '
            f'fill="{color}"/>'
        )
    total_w = CHROM_SIZES_MB.get(chrom, 0) * px_per_mb
    parts.append(
        f'<rect x="{x0}" y="{y}" width="{total_w:.1f}" height="{height}" '
        f'fill="none" stroke="#9ca3af" stroke-width="0.5" rx="3"/>'
    )
    return parts


def mb_axis_svg(
    x0: float,
    px_per_mb: float,
    top: float,
    bottom: float,
    max_mb: float,
    step: int = 50,
) -> List[str]:
    """Dashed Mb gridlines from ``top`` to ``bottom``, labelled below ``bottom``."""
    parts: List[str] = []
    for mb_val in range(0, int(max_mb) + step, step):
        if mb_val > max_mb + step / 2:
            break
        gx = x0 + mb_val * px_per_mb
        parts.append(
            f'<line x1="{gx:.1f}" y1="{top}" x2="{gx:.1f}" y2="{bottom}" '
            f'stroke="#e5e7eb" stroke-width="0.5" stroke-dasharray="3,3"/>'
        )
        parts.append(
            f'<text x="{gx:.1f}" y="{bottom + 14}" text-anchor="middle" '
            f'font-size="12" fill="#6b7280">{mb_val}</text>'
        )
    return parts
