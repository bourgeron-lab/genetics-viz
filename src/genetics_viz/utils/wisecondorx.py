"""WisecondorX BED parsing, CNV classification, and color utilities."""

from pathlib import Path
from typing import Any, Dict, List, Optional

import polars as pl
import yaml


def _load_wisecondorx_config() -> Dict[str, Any]:
    """Load WisecondorX thresholds and colors from YAML config."""
    config_path = (
        Path(__file__).parent.parent / "config" / "wisecondorx_thresholds.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


WISECONDORX_CONFIG: Dict[str, Any] = _load_wisecondorx_config()


def classify_cnv(ratio: Any, zscore: Any) -> str:
    """Classify CNV based on ratio (log2) and zscore thresholds.

    Returns the label from the config (e.g. "Robust LOSS", "Permissive GAIN")
    or "Below threshold" / "N/A".
    """
    try:
        r = float(ratio) if ratio and str(ratio) != "" else 0
        z = float(zscore) if zscore and str(zscore) != "" else 0

        # Robust calls checked first
        rl = WISECONDORX_CONFIG["robust_loss"]
        if r <= rl["ratio_threshold"] and z <= rl["zscore_threshold"]:
            return rl["label"]
        rg = WISECONDORX_CONFIG["robust_gain"]
        if r >= rg["ratio_threshold"] and z >= rg["zscore_threshold"]:
            return rg["label"]
        # Permissive calls (fallback)
        pl_ = WISECONDORX_CONFIG["permissive_loss"]
        if r <= pl_["ratio_threshold"] and z <= pl_["zscore_threshold"]:
            return pl_["label"]
        pg = WISECONDORX_CONFIG["permissive_gain"]
        if r >= pg["ratio_threshold"] and z >= pg["zscore_threshold"]:
            return pg["label"]
        return "Below threshold"
    except Exception:
        return "N/A"


def build_color_thresholds(metric: str) -> List[Dict[str, Any]]:
    """Build color_scale threshold list for ratio or zscore columns.

    Args:
        metric: Either "ratio" or "zscore".
    """
    robust_loss = WISECONDORX_CONFIG["robust_loss"]
    permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
    robust_gain = WISECONDORX_CONFIG["robust_gain"]
    permissive_gain = WISECONDORX_CONFIG["permissive_gain"]

    key = f"{metric}_threshold"
    return [
        {
            "op": "<=",
            "value": robust_loss[key],
            "color": robust_loss["color"],
            "weight": "bold",
        },
        {
            "op": "<=",
            "value": permissive_loss[key],
            "color": permissive_loss["color"],
            "weight": "600",
        },
        {
            "op": ">=",
            "value": robust_gain[key],
            "color": robust_gain["color"],
            "weight": "bold",
        },
        {
            "op": ">=",
            "value": permissive_gain[key],
            "color": permissive_gain["color"],
            "weight": "600",
        },
    ]


def build_call_colors() -> Dict[str, str]:
    """Build call label -> color mapping from WisecondorX config."""
    return {
        WISECONDORX_CONFIG["robust_loss"]["label"]: WISECONDORX_CONFIG["robust_loss"][
            "color"
        ],
        WISECONDORX_CONFIG["permissive_loss"]["label"]: WISECONDORX_CONFIG[
            "permissive_loss"
        ]["color"],
        WISECONDORX_CONFIG["robust_gain"]["label"]: WISECONDORX_CONFIG["robust_gain"][
            "color"
        ],
        WISECONDORX_CONFIG["permissive_gain"]["label"]: WISECONDORX_CONFIG[
            "permissive_gain"
        ]["color"],
    }


def parse_wisecondorx_bed(file_path: Path) -> Optional[pl.DataFrame]:
    """Parse a WisecondorX BED file into a Polars DataFrame.

    Keeps chr, start, end as separate columns for overlap filtering,
    and also creates the chr:start-end combined column and the
    wisecondorX CNV call classification column.

    Used by the search page for cohort-wide WCX queries.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()

    if not lines:
        return None

    # Parse tab-separated header
    header = lines[0].strip().split("\t")

    # Parse data rows
    data = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 7:
            while len(parts) < len(header):
                parts.append("")
            data.append(parts[: len(header)])

    if not data:
        return None

    df = pl.DataFrame(
        {
            col: [row[i] if i < len(row) else "" for row in data]
            for i, col in enumerate(header)
        }
    )

    # Convert numeric columns
    for col in ["start", "end", "ratio", "zscore"]:
        if col in df.columns:
            df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False))

    # Rename barcode column to sample if it exists
    if "barcode" in df.columns:
        df = df.rename({"barcode": "sample"})

    # Normalize chr prefix
    if "chr" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("chr").cast(pl.Utf8).str.starts_with("chr"))
            .then(pl.col("chr").cast(pl.Utf8))
            .otherwise(pl.lit("chr") + pl.col("chr").cast(pl.Utf8))
            .alias("chr")
        )

    # Create chr:start-end combined column (keep chr, start, end for filtering)
    if all(c in df.columns for c in ["chr", "start", "end"]):
        df = df.with_columns(
            (
                pl.col("chr")
                + ":"
                + pl.col("start").cast(pl.Int64).cast(pl.Utf8)
                + "-"
                + pl.col("end").cast(pl.Int64).cast(pl.Utf8)
            ).alias("chr:start-end")
        )

    # Add wisecondorX CNV call classification
    if "ratio" in df.columns and "zscore" in df.columns:
        df = df.with_columns(
            pl.struct(["ratio", "zscore"])
            .map_elements(
                lambda row: classify_cnv(row["ratio"], row["zscore"]),
                return_dtype=pl.Utf8,
            )
            .alias("wisecondorX")
        )

    return df


def parse_wisecondorx_bed_for_display(file_path: Path) -> Optional[pl.DataFrame]:
    """Parse a WisecondorX BED file for display in the SVs tab.

    Unlike parse_wisecondorx_bed(), this drops chr/start/end after creating
    the combined column, adds svlen and gene columns, and uses "call" as
    the classification column name.

    Used by the svs_tab for per-family display.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()

    if not lines:
        return None

    # Parse header - split by any whitespace
    header = lines[0].strip().split()

    # Parse data rows
    data = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) >= 7:
            while len(parts) < len(header):
                parts.append("")
            data.append(parts[: len(header)])

    if not data:
        return None

    df = pl.DataFrame(
        {
            col: [row[i] if i < len(row) else "" for row in data]
            for i, col in enumerate(header)
        }
    )

    # Convert numeric columns
    for col in ["start", "end", "ratio", "zscore"]:
        if col in df.columns:
            df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False))

    # Rename barcode column to sample if it exists
    if "barcode" in df.columns:
        df = df.rename({"barcode": "sample"})

    # Create chr:start-end column and svlen, then drop coordinate columns
    if "chr" in df.columns and "start" in df.columns and "end" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("chr").cast(pl.Utf8).str.starts_with("chr"))
            .then(pl.col("chr").cast(pl.Utf8))
            .otherwise(pl.lit("chr") + pl.col("chr").cast(pl.Utf8))
            .alias("chr_prefixed")
        )

        df = df.with_columns(
            (
                pl.col("chr_prefixed")
                + ":"
                + pl.col("start").cast(pl.Int64).cast(pl.Utf8)
                + "-"
                + pl.col("end").cast(pl.Int64).cast(pl.Utf8)
            ).alias("chr:start-end")
        )

        # Compute svlen before dropping
        df = df.with_columns(
            (pl.col("end").cast(pl.Int64) - pl.col("start").cast(pl.Int64)).alias(
                "svlen"
            )
        )

        df = df.drop(["chr_prefixed", "chr", "start", "end"])

    # Create gene column combining genic_symbol and exonic_symbol
    if "genic_symbol" in df.columns and "exonic_symbol" in df.columns:

        def create_gene_list(genic: Any, exonic: Any) -> str:
            genic_genes = (
                set(str(genic).split(",")) if genic and str(genic) != "" else set()
            )
            exonic_genes = (
                set(str(exonic).split(",")) if exonic and str(exonic) != "" else set()
            )
            genic_genes.discard("")
            exonic_genes.discard("")

            result = []
            for gene in exonic_genes:
                result.append(f"{gene.strip()}:exonic")
            for gene in genic_genes:
                if gene.strip() not in exonic_genes:
                    result.append(f"{gene.strip()}:genic")
            return ",".join(result) if result else ""

        df = df.with_columns(
            pl.struct(["genic_symbol", "exonic_symbol"])
            .map_elements(
                lambda row: create_gene_list(row["genic_symbol"], row["exonic_symbol"]),
                return_dtype=pl.Utf8,
            )
            .alias("gene")
        )

        priority_cols = (
            ["chr:start-end", "gene"] if "chr:start-end" in df.columns else ["gene"]
        )
        other_cols = [col for col in df.columns if col not in priority_cols]
        df = df.select(priority_cols + other_cols)
    elif "chr:start-end" in df.columns:
        other_cols = [col for col in df.columns if col != "chr:start-end"]
        df = df.select(["chr:start-end"] + other_cols)

    # Add CNV call classification
    if "ratio" in df.columns and "zscore" in df.columns:
        df = df.with_columns(
            pl.struct(["ratio", "zscore"])
            .map_elements(
                lambda row: classify_cnv(row["ratio"], row["zscore"]),
                return_dtype=pl.Utf8,
            )
            .alias("call")
        )

        # Reorder to put call after chr:start-end and gene
        if "chr:start-end" in df.columns and "gene" in df.columns:
            priority_cols = ["chr:start-end", "gene", "call"]
            other_cols = [col for col in df.columns if col not in priority_cols]
            df = df.select(priority_cols + other_cols)
        elif "chr:start-end" in df.columns:
            priority_cols = ["chr:start-end", "call"]
            other_cols = [col for col in df.columns if col not in priority_cols]
            df = df.select(priority_cols + other_cols)

    return df
