"""Cytoband data for ideogram rendering (hg38)."""

from pathlib import Path
from typing import Dict

# Giemsa stain colors for cytoband rendering
GIESTAIN_COLORS: Dict[str, str] = {
    "gneg": "#f5f5f5",
    "gpos25": "#c8c8c8",
    "gpos50": "#969696",
    "gpos75": "#646464",
    "gpos100": "#323232",
    "acen": "#d92f27",
    "gvar": "#646464",
    "stalk": "#969696",
}

# Canonical chromosome order and approximate sizes (Mb, hg38)
CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

CHROM_SIZES_MB: Dict[str, int] = {
    "1": 249,
    "2": 242,
    "3": 198,
    "4": 190,
    "5": 182,
    "6": 171,
    "7": 159,
    "8": 145,
    "9": 138,
    "10": 134,
    "11": 135,
    "12": 133,
    "13": 114,
    "14": 107,
    "15": 102,
    "16": 90,
    "17": 83,
    "18": 80,
    "19": 59,
    "20": 64,
    "21": 47,
    "22": 51,
    "X": 156,
    "Y": 57,
    "MT": 1,
}

# Validation status colours (avoids grays used by cytoband stains)
VALIDATION_COLORS: Dict[str, str] = {
    "present": "#22c55e",
    "absent": "#ef4444",
    "uncertain": "#f59e0b",
    "conflicting": "#fb923c",
    "TODO": "#8b5cf6",
}


def _load_cytobands() -> Dict[str, list]:
    """Load cytoband data from the bundled TSV file, grouped by chromosome."""
    config_path = Path(__file__).parent.parent / "config" / "cytobands_hg38.tsv"
    bands: Dict[str, list] = {}
    if not config_path.exists():
        return bands
    with open(config_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0].replace("chr", "").upper()
            if chrom == "M":
                chrom = "MT"
            start_mb = int(parts[1]) / 1_000_000
            end_mb = int(parts[2]) / 1_000_000
            name = parts[3]
            stain = parts[4]
            bands.setdefault(chrom, []).append(
                {"start": start_mb, "end": end_mb, "name": name, "stain": stain}
            )
    return bands


CYTOBANDS = _load_cytobands()


def norm_chrom(raw: str) -> str:
    """Normalise a chromosome string: strip 'chr' prefix, uppercase, M→MT."""
    c = str(raw).replace("chr", "").upper()
    return "MT" if c == "M" else c
