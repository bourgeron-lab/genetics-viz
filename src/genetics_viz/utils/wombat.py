"""Wombat result tables: one row per variant and sample.

Wombat writes one row per variant, sample and transcript/annotation, so the
same variant appears several times with differing annotation columns. Both the
family Wombat tab and the cohort search collapse those rows; the logic lives
here so the two stay identical.
"""

import polars as pl

#: The columns identifying one variant in one sample.
GROUPING_COLUMNS = ["#CHROM", "POS", "REF", "ALT", "sample"]

#: The HGVS of one transcript, picked before the rows are merged (see
#: :func:`aggregate_by_variant`): its ID, its MANE Select RefSeq ID, c. and p.
PICKED_HGVS_COLUMNS = {
    "_hgvs_transcript": "VEP_Feature",
    "_hgvs_mane": "VEP_MANE_SELECT",
    "_hgvsc": "VEP_HGVSc",
    "_hgvsp": "VEP_HGVSp",
}

_EMPTY = ("", ".", "-")


def _present(col: str) -> pl.Expr:
    value = pl.col(col).cast(pl.Utf8)
    return value.is_not_null() & ~value.is_in(_EMPTY)


def _pick_rank(columns: list[str]) -> pl.Expr:
    """0 for the MANE Select transcript, 1 for the canonical one, else null."""
    rank = pl.lit(None, dtype=pl.Int8)
    if "VEP_CANONICAL" in columns:
        rank = (
            pl.when(pl.col("VEP_CANONICAL").cast(pl.Utf8) == "YES")
            .then(pl.lit(1, dtype=pl.Int8))
            .otherwise(rank)
        )
    if "VEP_MANE_SELECT" in columns:
        rank = (
            pl.when(_present("VEP_MANE_SELECT"))
            .then(pl.lit(0, dtype=pl.Int8))
            .otherwise(rank)
        )
    return rank


def aggregate_by_variant(df: pl.DataFrame) -> pl.DataFrame:
    """Collapse rows to one per variant and sample.

    Every other column becomes the comma-joined distinct values it took, nulls,
    empty strings and ``.`` left out, in the order they first appear. The order
    is kept deliberately: plain ``unique()`` returns values in an arbitrary
    order, so the joined string changed between runs and polars versions.
    ``n_grouped`` counts the rows merged.

    Joining each column on its own loses which value came from which
    transcript, so the :data:`PICKED_HGVS_COLUMNS` are taken first, from one
    row: the MANE Select transcript's, else the canonical one's. They are
    absent when the table has no ``VEP_HGVSc``.
    """
    pick = "VEP_HGVSc" in df.columns
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
    if not pick:
        return df.group_by(GROUPING_COLUMNS, maintain_order=True).agg(exprs)
    df = df.with_columns(
        _pick_rank(df.columns).alias("_hgvs_rank"),
        *(
            pl.lit(None, dtype=pl.Utf8).alias(source)
            for source in PICKED_HGVS_COLUMNS.values()
            if source not in df.columns
        ),
    )
    for name, source in PICKED_HGVS_COLUMNS.items():
        value = pl.col(source).cast(pl.Utf8)
        exprs.append(value.sort_by("_hgvs_rank", nulls_last=True).first().alias(name))
    exprs.append(pl.col("_hgvs_rank").min())
    out = df.group_by(GROUPING_COLUMNS, maintain_order=True).agg(exprs)
    picked = pl.col("_hgvs_rank").is_not_null()
    return out.with_columns(
        pl.when(picked & ~pl.col(name).is_in(_EMPTY)).then(pl.col(name)).alias(name)
        for name in PICKED_HGVS_COLUMNS
    ).drop("_hgvs_rank")
