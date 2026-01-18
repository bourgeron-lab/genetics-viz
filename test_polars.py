from pathlib import Path

import polars as pl

pedigree_file = Path("example_data/cohorts/test_cohort/test_cohort.pedigree.tsv")
df = pl.read_csv(pedigree_file, separator="\t", infer_schema_length=0)
print("Columns:", df.columns)
print("Sample rows:")
for row in df.head(4).iter_rows(named=True):
    print(row)
