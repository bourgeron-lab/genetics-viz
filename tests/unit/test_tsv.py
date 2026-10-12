"""TSV reading helper (utils/tsv.py)."""

import gzip

import polars as pl
import pytest

from genetics_viz.utils.tsv import read_tsv_or_none


def numeric_then_text(path, n=150, sep="\n"):
    rows = [f"S{i}\t{1000 + i}" for i in range(n)] + ["Sx\tAU-ANG-FER-MPX3"]
    path.write_text(sep.join(["IID\tfamily_id", *rows]) + sep)
    return path


def test_zero_byte_file_is_none(tmp_path):
    (tmp_path / "e.tsv").write_text("")
    assert read_tsv_or_none(tmp_path / "e.tsv") is None


def test_header_only_is_an_empty_frame(tmp_path):
    (tmp_path / "h.tsv").write_text("a\tb\n")
    df = read_tsv_or_none(tmp_path / "h.tsv")
    assert df is not None and len(df) == 0 and df.columns == ["a", "b"]


def test_missing_file_still_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_tsv_or_none(tmp_path / "absent.tsv")


def test_inference_alone_breaks_on_late_text_ids(tmp_path):
    """The trap string_columns exists for: polars infers from the first rows."""
    with pytest.raises(pl.exceptions.ComputeError):
        read_tsv_or_none(numeric_then_text(tmp_path / "f.tsv"))


def test_string_columns_keep_ids_as_text(tmp_path):
    df = read_tsv_or_none(
        numeric_then_text(tmp_path / "f.tsv"), string_columns=["family_id", "absent"]
    )
    assert df.schema["family_id"] == pl.String
    assert df["family_id"][-1] == "AU-ANG-FER-MPX3"


def test_string_columns_with_gzip(tmp_path):
    path = tmp_path / "f.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("IID\tx\n001\t1\n")
    df = read_tsv_or_none(path, string_columns=["IID"])
    assert df["IID"].to_list() == ["001"]


def test_string_columns_with_crlf_line_endings(tmp_path):
    df = read_tsv_or_none(
        numeric_then_text(tmp_path / "f.tsv", sep="\r\n"), string_columns=["family_id"]
    )
    assert df.schema["family_id"] == pl.String


def test_callers_schema_overrides_win(tmp_path):
    (tmp_path / "f.tsv").write_text("IID\tn\n1\t2\n")
    df = read_tsv_or_none(
        tmp_path / "f.tsv", string_columns=["IID"], schema_overrides={"IID": pl.Int64}
    )
    assert df.schema["IID"] == pl.Int64
