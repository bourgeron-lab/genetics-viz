"""Wombat result tables: one row per variant and sample.

Wombat writes one row per variant, sample and transcript/annotation, so the
same variant appears several times with differing annotation columns. Both the
family Wombat tab and the cohort search collapse those rows; the logic lives
here so the two stay identical.
"""

import polars as pl

#: The columns identifying one variant in one sample.
GROUPING_COLUMNS = ["#CHROM", "POS", "REF", "ALT", "sample"]


def aggregate_by_variant(df: pl.DataFrame) -> pl.DataFrame:
    """Collapse rows to one per variant and sample.

    Every other column becomes the comma-joined distinct values it took, nulls,
    empty strings and ``.`` left out, in the order they first appear. The order
    is kept deliberately: plain ``unique()`` returns values in an arbitrary
    order, so the joined string changed between runs and polars versions.
    ``n_grouped`` counts the rows merged.
    """
    exprs = [pl.len().alias("n_grouped")]
    for col in df.columns:
        if col in GROUPING_COLUMNS:
            continue
        value = pl.col(col).cast(pl.Utf8)
        exprs.append(
            value.filter(value.is_not_null() & (value != "") & (value != "."))
            .unique(maintain_order=True)
            .str.join(",")
            .alias(col)
        )
    return df.group_by(GROUPING_COLUMNS, maintain_order=True).agg(exprs)
