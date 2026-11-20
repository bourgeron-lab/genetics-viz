"""The cohort's pedigree edit history (utils/pedigree_history.py)."""

import re
from datetime import datetime
from pathlib import Path

import datagen
import pytest

from genetics_viz.utils import pedigree_edit as pe
from genetics_viz.utils import pedigree_history as ph
from genetics_viz.utils.sharding import get_family_path


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    return datagen.build_data_dir(tmp_path)


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: True)


def _source(data_dir: Path, fid: str) -> pe.FamilySource:
    return pe.load_family(
        data_dir, fid, [data_dir / "cohorts" / "DEMO" / "DEMO.pedigree.tsv"]
    )


def _edits(data_dir: Path) -> None:
    """1001: a dummy father, then S3 unaffected; AU-X-1: S6 affected."""
    src = _source(data_dir, "1001")
    draft = src.draft()
    draft.set_parent("S3", "father", None)
    draft.add_dummy_parents(["S3"], roles=["father"])
    pe.save_pedigree_edit(
        src,
        draft,
        user="alice",
        comment="S1 is not S3's father",
        flags=["NOT_PARENT"],
        now=datetime(2026, 11, 10, 9, 0),
    )
    src = _source(data_dir, "1001")
    draft = src.draft()
    draft.set_phenotype("S3", "1")
    pe.save_pedigree_edit(
        src, draft, user="bob", comment="Recoded", now=datetime(2026, 11, 12, 9, 0)
    )
    src = _source(data_dir, "AU-X-1")
    draft = src.draft()
    draft.set_phenotype("S6", "2")
    pe.save_pedigree_edit(
        src, draft, user="alice", comment="Diagnosed", now=datetime(2026, 11, 11, 9, 0)
    )


def _classes(svg: str) -> dict:
    return dict(
        (m.group(2), m.group(1))
        for m in re.finditer(r'<g class="([^"]*)"[^>]*data-id="([^"]+)"', svg)
    )


def test_families_are_grouped_last_edited_first(data_dir, curator):
    _edits(data_dir)
    history = ph.build_history(data_dir, "DEMO", user="carol")
    assert [f.family_id for f in history.families] == ["1001", "AU-X-1"]
    assert history.n_edits == 3 and history.users == ["alice", "bob"]
    first = history.families[0]
    assert [e.comment for e in first.edits] == ["S1 is not S3's father", "Recoded"]
    # The family log's record: it keeps the flags the index does not.
    assert first.edits[0].flags == ("NOT_PARENT",)
    assert history.document_title.startswith("DEMO_pedigree_history_")


def test_the_trees_span_the_first_edit_to_the_last(data_dir, curator):
    _edits(data_dir)
    family = ph.build_history(data_dir, "DEMO").families[0]
    before, after = _classes(family.before.svg), _classes(family.after.svg)
    assert "DUMMY_S3_F" not in before and "changed" not in before["S3"]
    # S3 changed twice (father, phenotype): outlined once, against the first rows.
    assert "changed" in after["S3"] and "changed" in after["DUMMY_S3_F"]
    assert "changed" not in after["S1"]
    assert family.before.width.endswith("mm")


def test_only_the_cohorts_edits(data_dir, curator):
    _edits(data_dir)
    assert ph.cohort_edits(data_dir, "OTHER") == []
    assert ph.history_html(data_dir, "OTHER") == ""


def test_missing_kept_rows_are_said_not_drawn(data_dir, curator):
    _edits(data_dir)
    for kept in (get_family_path(data_dir, "AU-X-1") / "pedigree_edits").glob(
        "*.before.tsv"
    ):
        kept.unlink()
    family = ph.build_history(data_dir, "DEMO").families[1]
    assert not family.before.svg and "missing" in family.before.missing
    assert family.after.svg


def test_the_document(data_dir, curator):
    _edits(data_dir)
    statuses = {"1001": pe.EditStatus(pe.read_edit_index(data_dir)[1], True)}
    html = ph.history_html(data_dir, "DEMO", user="carol", statuses=statuses)
    assert html.startswith("<!DOCTYPE html>")
    assert "Cohort DEMO" in html and "by carol" in html
    assert "S3: father S1 -&gt; DUMMY_S3_F" in html
    assert "2026-11-10 09:00:00" in html and "bob" in html
    assert "Sample QC flags addressed: NOT_PARENT" in html
    assert "Waiting for a run" in html
    assert html.count("<svg") == 4  # two families, before and after
    # Each drawing's ids are its own.
    assert 'id="h0-b-' in html and 'id="h1-a-' in html
