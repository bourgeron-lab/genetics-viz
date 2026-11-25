"""The Families tab's rows, columns and filter (utils/family_rows.py)."""

import datagen
from genetics_viz.models import Cohort
from genetics_viz.utils import family_rows as fr


def _demo(data_dir):
    return Cohort.from_directory(data_dir / "cohorts" / "DEMO")


def test_structure_of_each_demo_family(data_dir):
    cohort = _demo(data_dir)
    got = {fid: fr.family_structure(f) for fid, f in cohort.families.items()}
    assert got["1001"] == ("Trio", "Trio", 1)
    assert got["1002"] == ("Singleton", "Singleton", 0)
    assert got["AU-X-1"] == ("Duo", "Duo", 0)


def test_dummy_parents_are_not_a_trio(data_dir):
    datagen.add_dummy_case(data_dir)
    cohort = _demo(data_dir)
    # MZ twins S7/S8 hang from dummy parents: no sequenced parent.
    kind, label, trios = fr.family_structure(cohort.families["2001"])
    assert (kind, trios) == ("No parents", 0) and label == "No parents (2)"


def test_row_reads_the_cohort_tables(data_dir):
    datagen.add_denovo_calls(data_dir)
    cohort = _demo(data_dir)
    facts = fr.load_family_facts(cohort, data_dir, curator=True)
    assert {"ancestry", "roh", "denovo"} <= facts.found
    row = fr.family_row(
        cohort.families["1001"],
        facts,
        upd_open=["high", "medium"],
        str_open=(1, 0),
        diagnostics=(2, "pathogenic"),
        curator=True,
    )
    assert row["Ancestry"] == "EUR"
    assert row["Affected"] == "1 / 3" and row["_affected"] == 1
    # S3's verdict in the DEMO tables is low.
    assert row["DenovoQC_badges"][0]["label"] == fr.denovo.qc_info("low")["label"]
    assert row["Diagnostics_badges"][0]["label"] == "P · 2"
    # SHANK3 is strict and not validated; PLAIN is flagged, so not strict.
    labels = [b["label"] for b in row["Review_badges"]]
    assert labels == ["De novo 1", "UPD 2", "STR 1"] and row["_review_n"] == 4
    assert row["Conclusion"] == ""


def test_readers_get_no_review_or_conclusion(data_dir):
    cohort = _demo(data_dir)
    facts = fr.load_family_facts(cohort, data_dir, curator=False)
    row = fr.family_row(cohort.families["1002"], facts, upd_open=["high"])
    assert "Review" not in row and "Conclusion" not in row
    assert row["_diag_cats"] == [fr.NO_DIAGNOSTIC]
    columns = fr.family_columns("DEMO", with_roh=True)
    assert {"Review", "Conclusion"}.isdisjoint(c["id"] for c in columns)


def test_dn_open_counts_skip_validated_and_low_impact():
    rows = [
        {"family_id": "F", "child": "C", "chrom": "chr1", "pos": 1, "ref": "A",
         "alt": "T", "strict": "yes", "consequences": "stop_gained"},
        {"family_id": "F", "child": "C", "chrom": "chr1", "pos": 2, "ref": "A",
         "alt": "T", "strict": "yes", "consequences": "synonymous_variant"},
        {"family_id": "F", "child": "C", "chrom": "chr1", "pos": 3, "ref": "A",
         "alt": "T", "strict": "yes", "consequences": "missense_variant"},
        {"family_id": "F", "child": "C", "chrom": "chr1", "pos": 4, "ref": "A",
         "alt": "T", "strict": "no", "consequences": "stop_gained"},
    ]  # fmt: skip
    assert fr.dn_open_counts(rows, set()) == {"F": (2, 1)}
    assert fr.dn_open_counts(rows, {("chr1:1:A:T", "C")}) == {"F": (1, 0)}


def test_filter_families():
    rows = [
        {"FID": "A1", "_structure": "Trio", "Trios": 1, "_affected": 1,
         "_diag_cats": ["pathogenic"], "_review_n": 0, "Conclusion": "2026-11-01",
         "_qc_cats": ["Sex"], "_qc_rank": 3},
        {"FID": "B2", "_structure": "Extended", "Trios": 2, "_affected": 2,
         "_diag_cats": [fr.NO_DIAGNOSTIC], "_review_n": 3, "Conclusion": "",
         "_qc_cats": ["No issue"], "_qc_rank": 0},
    ]  # fmt: skip
    fids = lambda f: [r["FID"] for r in fr.filter_families(rows, f, "Any issue")]  # noqa: E731
    assert fids({}) == ["A1", "B2"]
    assert fids({"FID": "b"}) == ["B2"]
    assert fids({"Structure": ["Trio"]}) == ["A1"]
    assert fids({"Affected": True}) == ["B2"]
    assert fids({"Diagnostics": [fr.NO_DIAGNOSTIC]}) == ["B2"]
    assert fids({"Review": True, "Conclusion": True}) == ["B2"]
    assert fids({"QC": ["Any issue"]}) == ["A1"]


def test_compact_qc_cell_is_one_badge_plus_cell_line():
    from types import SimpleNamespace as NS

    flags = [
        NS(rank=3, severity="high", code="NOT_PARENT", keys=[("F", "S1"), ("F", "S3")]),
        NS(rank=2, severity="medium", code="CELL_LINE", keys=[("F", "S3")]),
    ]
    cell = {
        "QC": "",
        "QC_badges": [
            {"label": "Relationship in family", "color": "red", "tooltip": "S1, S3"},
            {"label": "DNA source", "color": "orange", "tooltip": "S3"},
        ],
        "_qc_cats": ["Relationship in family", "DNA source"],
        "_qc_rank": 3,
    }
    out = fr.compact_qc_cell(cell, flags, [], "F", "Acknowledged")
    labels = [b["label"] for b in out["QC_badges"]]
    assert labels == ["High · 2 members", "Cell line"]
    assert "Relationship in family: S1, S3" in out["QC_badges"][0]["tooltip"]
    assert out["_qc_cats"] == cell["_qc_cats"]  # the filter still sees them

    # Acknowledged only: the grey badge; the cell line stays.
    acked = fr.compact_qc_cell(
        {"QC": "", "QC_badges": [{"label": "Acknowledged", "color": "grey"}]},
        [],
        [flags[1]],
        "F",
        "Acknowledged",
    )
    assert [b["label"] for b in acked["QC_badges"]] == ["Acknowledged", "Cell line"]


def test_family_order_is_remembered_per_cohort():
    storage: dict = {}
    fr.save_family_order(storage, "DEMO", ["B", "A"])
    fr.save_family_order(storage, "OTHER", ["X"])
    assert fr.family_order(storage, "DEMO", "A", ["A", "B", "C"]) == (["B", "A"], True)
    # A family outside the list walks every family by ID.
    assert fr.family_order(storage, "DEMO", "C", ["C", "B", "A"]) == (
        ["A", "B", "C"],
        False,
    )
    fr.save_family_order(storage, "DEMO", None)
    assert fr.family_order(storage, "DEMO", "A", ["B", "A"]) == (["A", "B"], False)
    assert storage[fr.FAMILY_ORDER_KEY] == {"OTHER": ["X"]}
