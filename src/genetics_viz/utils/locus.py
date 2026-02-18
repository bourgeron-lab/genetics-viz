"""Locus query parsing and DataFrame filtering utilities."""

import re
from typing import Any, Dict

import polars as pl


def parse_locus_query(query: str) -> Dict[str, Any]:
    """Parse locus query into filter parameters.

    Returns dict with 'type' and relevant filter parameters.
    """
    query = query.strip()

    # chr1:10000:A:GC - exact variant
    variant_pattern = r"^(chr)?(\w+):(\d+):([ACGT]+):([ACGT]+)$"
    match = re.match(variant_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "exact_variant",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "pos": int(match.group(3)),
            "ref": match.group(4).upper(),
            "alt": match.group(5).upper(),
        }

    # chr1:10000-10100 - range
    range_pattern = r"^(chr)?(\w+):(\d+)-(\d+)$"
    match = re.match(range_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "range",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "start": int(match.group(3)),
            "end": int(match.group(4)),
        }

    # chr1:10000 - exact position
    pos_pattern = r"^(chr)?(\w+):(\d+)$"
    match = re.match(pos_pattern, query, re.IGNORECASE)
    if match:
        return {
            "type": "exact_position",
            "chrom": f"chr{match.group(2)}"
            if not match.group(1)
            else match.group(1) + match.group(2),
            "pos": int(match.group(3)),
        }

    # ENSG00000164099 - gene ID
    if re.match(r"^ENSG\d+$", query, re.IGNORECASE):
        return {
            "type": "gene_id",
            "gene_id": query.upper(),
        }

    # SHANK* - wildcard gene
    if "*" in query:
        return {
            "type": "gene_wildcard",
            "pattern": query.replace("*", "").upper(),
        }

    # SHANK3 - exact gene
    return {
        "type": "gene_name",
        "gene_name": query.upper(),
    }


def filter_dataframe(df: pl.DataFrame, query_params: Dict[str, Any]) -> pl.DataFrame:
    """Filter dataframe based on parsed query parameters."""
    query_type = query_params["type"]

    if query_type == "exact_variant":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") == query_params["pos"])
            & (pl.col("REF") == query_params["ref"])
            & (pl.col("ALT") == query_params["alt"])
        )

    elif query_type == "range":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") >= query_params["start"])
            & (pl.col("POS") <= query_params["end"])
        )

    elif query_type == "exact_position":
        return df.filter(
            (pl.col("#CHROM") == query_params["chrom"])
            & (pl.col("POS") == query_params["pos"])
        )

    elif query_type == "gene_id":
        # VEP_Gene can contain multiple genes separated by &
        gene_id = query_params["gene_id"]
        return df.filter(pl.col("VEP_Gene").str.to_uppercase().str.contains(gene_id))

    elif query_type == "gene_wildcard":
        # VEP_SYMBOL can contain multiple symbols
        pattern = query_params["pattern"]
        return df.filter(pl.col("VEP_SYMBOL").str.to_uppercase().str.contains(pattern))

    elif query_type == "gene_name":
        # Exact match in VEP_SYMBOL (case-insensitive, as part of the field)
        gene_name = query_params["gene_name"]
        return df.filter(
            pl.col("VEP_SYMBOL").str.to_uppercase().str.contains(gene_name)
        )

    return df


def filter_bed_dataframe(
    df: pl.DataFrame, query_params: Dict[str, Any], exonic: bool = False
) -> pl.DataFrame:
    """Filter WisecondorX BED dataframe based on parsed query parameters.

    Uses chr/start/end for coordinate queries (overlap logic) and
    genic_symbol/genic_ensg (or exonic_symbol/exonic_ensg when exonic=True)
    for gene queries.
    """
    query_type = query_params["type"]
    symbol_col = "exonic_symbol" if exonic else "genic_symbol"
    ensg_col = "exonic_ensg" if exonic else "genic_ensg"

    if query_type == "exact_variant":
        # No REF/ALT in BED files — fall back to position overlap
        chrom = query_params["chrom"]
        pos = float(query_params["pos"])
        return df.filter(
            (pl.col("chr") == chrom)
            & (pl.col("start") <= pos)
            & (pl.col("end") >= pos)
        )

    elif query_type == "range":
        chrom = query_params["chrom"]
        q_start = float(query_params["start"])
        q_end = float(query_params["end"])
        return df.filter(
            (pl.col("chr") == chrom)
            & (pl.col("start") <= q_end)
            & (pl.col("end") >= q_start)
        )

    elif query_type == "exact_position":
        chrom = query_params["chrom"]
        pos = float(query_params["pos"])
        return df.filter(
            (pl.col("chr") == chrom)
            & (pl.col("start") <= pos)
            & (pl.col("end") >= pos)
        )

    elif query_type == "gene_id":
        gene_id = query_params["gene_id"]
        if ensg_col in df.columns:
            return df.filter(
                pl.col(ensg_col)
                .cast(pl.Utf8)
                .str.to_uppercase()
                .str.contains(gene_id)
            )
        return df.head(0)

    elif query_type == "gene_wildcard":
        pattern = query_params["pattern"]
        if symbol_col in df.columns:
            return df.filter(
                pl.col(symbol_col)
                .cast(pl.Utf8)
                .str.to_uppercase()
                .str.contains(pattern)
            )
        return df.head(0)

    elif query_type == "gene_name":
        gene_name = query_params["gene_name"]
        if symbol_col in df.columns:
            return df.filter(
                pl.col(symbol_col)
                .cast(pl.Utf8)
                .str.to_uppercase()
                .str.contains(gene_name)
            )
        return df.head(0)

    return df
