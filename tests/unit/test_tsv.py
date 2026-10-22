"""TSV reading helper (utils/tsv.py)."""

import gzip
import threading

import polars as pl
import pytest

from genetics_viz.utils.tsv import (
    append_tsv_row,
    append_tsv_rows,
    read_tsv_or_none,
    rewrite_tsv_rows,
    set_ignore,
)


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


def test_clean_tsv_field():
    from genetics_viz.utils.tsv import clean_tsv_field

    assert clean_tsv_field("a\tb\nc\r\nd\re") == "a b c d e"


HEADER = ["FID", "Variant", "Sample", "Comment", "Ignore", "Timestamp"]


def test_concurrent_first_saves_write_one_header(tmp_path):
    def race(path):
        start = threading.Barrier(8)

        def save(i):
            start.wait()
            append_tsv_row(path, HEADER, ["F1", "v", f"S{i}", "", "0", "t"])

        threads = [threading.Thread(target=save, args=(i,)) for i in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        return path.read_text().splitlines()

    # One race rarely interleaves; twenty catch an unflushed write every time.
    for n in range(20):
        lines = race(tmp_path / "validations" / f"new{n}.tsv")
        assert lines.count("\t".join(HEADER)) == 1 and len(lines) == 9


def test_append_without_rows_only_creates_the_header(tmp_path):
    path = tmp_path / "notes.tsv"
    append_tsv_rows(path, HEADER, [])
    append_tsv_rows(path, HEADER, [])
    assert path.read_text() == "\t".join(HEADER) + "\n"


def test_append_rows_writes_every_row(tmp_path):
    path = tmp_path / "v.tsv"
    append_tsv_rows(
        path,
        HEADER,
        [["F1", "v", "S1", "", "0", "t"], ["F1", "v", "S2", "x", "0", "t"]],
    )
    assert len(path.read_text().splitlines()) == 3


def test_rewrite_keeps_quotes_literal(tmp_path):
    path = tmp_path / "v.tsv"
    append_tsv_row(path, HEADER, ["F1", "v", "S1", '"a" b', "0", "t1"])
    append_tsv_row(path, HEADER, ["F1", "v", "S2", 'c "d"', "0", "t2"])
    assert set_ignore(path, "F1", "v", "S2", "t2", "1")
    lines = path.read_text().splitlines()
    assert lines[1].split("\t")[3] == '"a" b'
    assert lines[2].split("\t")[3:5] == ['c "d"', "1"]


def test_rewrite_leaves_the_file_alone_without_a_change(tmp_path):
    path = tmp_path / "v.tsv"
    assert not rewrite_tsv_rows(path, lambda rows: True)  # missing
    append_tsv_row(path, HEADER, ["F1", "v", "S1", "", "0", "t1"])
    before = path.read_text()
    assert not set_ignore(path, "F1", "v", "S1", "other-ts", "1")
    assert path.read_text() == before


def test_rewrite_can_drop_rows_and_shrink_the_file(tmp_path):
    path = tmp_path / "v.tsv"
    for i in range(3):
        append_tsv_row(
            path, HEADER, ["F1", "v", f"S{i}", "long comment " * 5, "0", "t"]
        )

    def keep_first(rows):
        del rows[1:]
        return True

    assert rewrite_tsv_rows(path, keep_first)
    assert len(path.read_text().splitlines()) == 2
