"""Dummy pedigree members (``DUMMY_`` IDs): kept in the pedigree, never looked up.

A dummy was never sequenced; it links real members, as the shared parents of
MZ twins or the grandparents joining an aunt to a cousin. The data are
``datagen.DUMMY_PEDIGREE``: families 2001 (twins) and 2002 (three
generations, one dummy named only in a father column).
"""

from pathlib import Path

import datagen
import pytest

from genetics_viz.components.validation_loader import parents_to_validate
from genetics_viz.models import (
    Cohort,
    DataStore,
    Family,
    Sample,
    is_dummy_id,
    is_sequenced_parent,
)
from genetics_viz.utils import roh
from genetics_viz.utils.family_report import build_family_report
from genetics_viz.utils.inheritance import family_pedigree
from genetics_viz.utils.pedigree import load_family_object, load_family_pedigree
from genetics_viz.utils.pedigree_labels import DUMMY_TOOLTIP, dummy_badges
from genetics_viz.utils.report_html import render_report_html
from genetics_viz.utils.sample_counts import FOUNDERS, compare_to_cohort
from genetics_viz.utils.search_rows import (
    IndividualFilters,
    filter_individuals,
    pedigree_data_from_cohort,
)
from genetics_viz.utils.sharding import get_family_path

REAL = {"2001": ["S7", "S8"], "2002": ["S9", "S10"]}


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = datagen.build_data_dir(tmp_path)
    datagen.add_dummy_case(data)
    return data


@pytest.fixture
def cohort(data_dir) -> Cohort:
    store = DataStore(data_dir)
    store.load()
    return store.cohorts["DEMO"]


# --- Detection -------------------------------------------------------------------


@pytest.mark.parametrize(
    "sample_id, dummy",
    [
        ("DUMMY_B00IDL0_F", True),
        ("DUMMY_B00INOK_MF", True),
        ("DUMMY_", True),
        ("B00IDL0", False),
        ("dummy_B00IDL0_F", False),  # the prefix is case-sensitive
        ("XDUMMY_B00IDL0_F", False),
        ("DUMMYB00IDL0", False),
        ("", False),
        (None, False),
    ],
)
def test_a_dummy_is_its_prefix(sample_id, dummy):
    assert is_dummy_id(sample_id) is dummy


@pytest.mark.parametrize(
    "parent, sequenced",
    [("S1", True), ("DUMMY_S7_F", False), (None, False), ("0", False), ("-", False)],
)
def test_a_sequenced_parent_is_neither_missing_nor_dummy(parent, sequenced):
    assert is_sequenced_parent(parent) is sequenced


# --- The cohort model ---------------------------------------------------------------


def test_dummies_are_members_but_not_samples(cohort):
    twins = cohort.families["2001"]
    assert [s.sample_id for s in twins.samples] == ["S7", "S8"]
    assert [s.sample_id for s in twins.dummies] == ["DUMMY_S7_F", "DUMMY_S7_M"]
    assert [s.sample_id for s in twins.members] == [
        "S7", "S8", "DUMMY_S7_F", "DUMMY_S7_M",
    ]  # fmt: skip
    assert all(s.is_dummy for s in twins.dummies)
    assert not any(s.is_dummy for s in twins.samples)
    assert twins.num_samples == 2
    assert twins.get_sample("DUMMY_S7_F") is None
    # DEMO's six samples, then S7-S10: no dummy counted.
    assert cohort.num_samples == 10
    assert sum(f.num_samples for f in cohort.families.values()) == 10


def test_three_generations_and_a_dummy_with_no_row(cohort):
    family = cohort.families["2002"]
    assert [s.sample_id for s in family.samples] == ["S9", "S10"]
    dummies = {s.sample_id: s for s in family.dummies}
    assert list(dummies) == [
        "DUMMY_S9_M", "DUMMY_S9_MF", "DUMMY_S9_MM", "DUMMY_S9_F",
    ]  # fmt: skip
    # A dummy has dummy parents of its own...
    mother = dummies["DUMMY_S9_M"]
    assert (mother.father_id, mother.mother_id) == ("DUMMY_S9_MF", "DUMMY_S9_MM")
    # ...and the one named only in S9's father column is a member too: male,
    # no parents, no phenotype.
    father = dummies["DUMMY_S9_F"]
    assert (father.sex, father.father_id, father.phenotype) == ("1", None, None)
    assert father.family_id == "2002"


def test_member_rows_leave_dummies_out_unless_asked(cohort):
    members = cohort.get_family_members("2001")
    assert [m["Sample ID"] for m in members] == ["S7", "S8"]
    assert members[0]["Father"] == "DUMMY_S7_F" and members[0]["Dummy"] is False
    full = cohort.get_family_members("2001", include_dummies=True)
    assert [(m["Sample ID"], m["Dummy"]) for m in full] == [
        ("S7", False), ("S8", False), ("DUMMY_S7_F", True), ("DUMMY_S7_M", True),
    ]  # fmt: skip


def test_twins_with_dummy_parents_are_founders():
    """No sequenced parent: grouped with the founders, as for de novo counts."""
    twin = Sample("S7", "2001", father_id="DUMMY_S7_F", mother_id="DUMMY_S7_M")
    assert twin.is_founder
    assert not Sample("S3", "1001", father_id="S1", mother_id="DUMMY_S3_M").is_founder
    family = Family("2001")
    for s in (twin, Sample("DUMMY_S7_F", "2001"), Sample("DUMMY_S7_M", "2001")):
        family.add_member(s)
    assert (family.num_samples, family.num_founders) == (1, 1)


@pytest.mark.parametrize("header", ["", "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"])
def test_per_family_pedigree_splits_dummies_alike(data_dir, tmp_path, header):
    text = (get_family_path(data_dir, "2002") / "2002.pedigree.tsv").read_text()
    path = tmp_path / "2002.pedigree.tsv"
    path.write_text(header + text)
    assert [m["Sample ID"] for m in load_family_pedigree(path)] == ["S9", "S10"]
    full = load_family_pedigree(path, include_dummies=True)
    assert [m["Sample ID"] for m in full if m["Dummy"]] == [
        "DUMMY_S9_M", "DUMMY_S9_MF", "DUMMY_S9_MM", "DUMMY_S9_F",
    ]  # fmt: skip
    family = load_family_object(path, "2002")
    assert (family.num_samples, len(family.dummies)) == (2, 4)


# --- Pedigree logic resolves through dummies ---------------------------------------


def test_twins_with_dummy_parents_are_declared_siblings(data_dir):
    froh = roh.read_roh_table(
        next((get_family_path(data_dir, "2001") / "roh").glob("*.froh.tsv"))
    )
    parents = roh.parents_from_froh(froh)
    assert parents["S7"] == ("DUMMY_S7_F", "DUMMY_S7_M")
    assert roh.declared_relation("S7", "S8", parents).label == "Full sibs"


def test_dummies_are_neither_affected_nor_unaffected(cohort):
    """A dummy's phenotype is never read: its children's parents are unknown."""
    lookup = family_pedigree(cohort.families["2001"].samples)
    assert lookup("S7") == {"sex": "1", "father_pheno": None, "mother_pheno": None}
    assert lookup("DUMMY_S7_F") == {}


# --- Lookups and counts skip dummies -----------------------------------------------


def test_search_knows_no_dummy_and_wants_sequenced_parents(cohort):
    pedigree = pedigree_data_from_cohort(cohort)
    assert not any(is_dummy_id(s) for s in pedigree)
    assert pedigree["S7"]["Father"] == "DUMMY_S7_F"
    filters = IndividualFilters()
    filters.has_parents["value"] = True
    rows = [{"sample": s} for s in ("S3", "S7", "S8", "S9", "S10")]
    assert filter_individuals(rows, pedigree, filters) == [{"sample": "S3"}]


def test_cohort_comparison_counts_no_dummy(cohort):
    counts = {"S7": 3, "DUMMY_S7_F": 99}
    comparisons, missing = compare_to_cohort(counts, cohort, ["S7", "S8"])
    assert missing == []
    assert {c.sample_id: c.group.name for c in comparisons} == {
        "S7": FOUNDERS,
        "S8": FOUNDERS,
    }
    # DEMO's four founders and S7-S10, whose parents are dummies; no dummy.
    assert comparisons[0].group.n == 8


def test_validating_a_childs_call_records_no_dummy_parent():
    parents = {"father": "DUMMY_S7_F", "mother": "DUMMY_S7_M"}
    for inheritance in ("maternal", "paternal", "either"):
        assert parents_to_validate(inheritance, parents, set()) == []
    half = {"father": "S1", "mother": "DUMMY_S3_M"}
    assert parents_to_validate("either", half, set()) == ["S1"]
    # Nor the "-" placeholder the member dicts give a missing parent.
    assert parents_to_validate("either", {"father": "-", "mother": "-"}, set()) == []


def test_table_badge_marks_a_dummy_id_only():
    (badge,) = dummy_badges("DUMMY_S7_F")
    assert badge["label"] == "DUMMY_S7_F" and badge["tooltip"] == DUMMY_TOOLTIP
    assert dummy_badges("S1") == [] and dummy_badges(None) == []


# --- Family report -------------------------------------------------------------------


def test_report_lists_dummies_muted_and_reads_real_members_only(
    cohort, data_dir, app_config
):
    report = build_family_report(
        data_dir, "2002", cohort.families["2002"], cohort, "curator"
    )
    pedigree = report.section("pedigree")
    assert [r["iid"] for r in pedigree.rows] == [
        "S9", "S10", "DUMMY_S9_M", "DUMMY_S9_MF", "DUMMY_S9_MM", "DUMMY_S9_F",
    ]  # fmt: skip
    s9 = pedigree.rows[0]
    assert s9["father"] == "DUMMY_S9_F" and s9["_father_note"]
    assert not any(r.get("_dummy") for r in pedigree.rows[:2])
    for row in pedigree.rows[2:]:
        assert row["_dummy"] and row["_muted"] and row["_iid_note"]
        # No status, never "Unknown": a dummy was not assessed.
        assert row["phenotype"] == "—" and row["karyotype"] == "—"
    assert pedigree.footnotes
    # Every result section reads the sequenced members only.
    assert {r["iid"] for r in report.section("ancestry").rows} == {"S9", "S10"}
    html = render_report_html(report)
    assert html.count('class="muted dummy ') == 4
    assert "dummy, not sequenced" in html


def test_report_of_a_family_without_dummies_is_unchanged(cohort, data_dir, app_config):
    report = build_family_report(
        data_dir, "1001", cohort.families["1001"], cohort, "curator"
    )
    pedigree = report.section("pedigree")
    assert [r["iid"] for r in pedigree.rows] == ["S1", "S2", "S3"]
    assert pedigree.footnotes == []
    html = render_report_html(report)
    assert 'class="muted dummy' not in html and "not sequenced" not in html
