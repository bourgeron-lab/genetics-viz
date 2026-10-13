"""Diagnostic kinds (utils/diagnostic_files.py) and the STR kind's readers."""

import pytest

from genetics_viz.components import diagnostic_loader as dl
from genetics_viz.utils import str_eh
from genetics_viz.utils.diagnostic_files import (
    DIAGNOSTIC_TYPES,
    diagnostic_path,
    diagnostic_paths,
)

KEY = "chr4:3074876-3074933:HTT"


@pytest.fixture(autouse=True)
def curator(monkeypatch):
    monkeypatch.setattr(dl, "can_write", lambda: True)


def _save(data_dir, dtype, variant, sample="S3", verdict="pathogenic", ts="t1"):
    dl.save_diagnostic_entry(
        diagnostic_path(data_dir, dtype), "1001", variant, "HTT", "pathogenic",
        sample, "curator", ts, "", verdict,
    )  # fmt: skip


def test_kinds_include_str(tmp_path):
    assert DIAGNOSTIC_TYPES == ["SNV", "SV", "STR"]
    assert diagnostic_path(tmp_path, "STR") == tmp_path / "diagnostics" / "strs.tsv"
    assert [t for t, _ in diagnostic_paths(tmp_path)] == DIAGNOSTIC_TYPES


def test_family_diagnostics_read_every_kind(tmp_path):
    _save(tmp_path, "SNV", "chr1:1:A:T")
    _save(tmp_path, "STR", KEY)
    entries = dl.load_family_diagnostics(tmp_path, "1001", ["S3"])
    assert {(e["_source"], e["Variant"]) for e in entries} == {
        ("snv", "chr1:1:A:T"),
        ("str", KEY),
    }


def test_diagnostic_badge_in_the_eh_table(data_dir):
    from genetics_viz.pages.cohort.components import str_tab

    fam = str_eh.read_eh_table(
        str_eh.eh_path(str_eh.get_family_eh_dir(data_dir, "1001"), "1001",
                       "jmsa50", "genotypes")
    )  # fmt: skip
    _save(data_dir, "STR", KEY)
    statuses = str_tab._diagnostic_statuses(data_dir, "1001", fam)
    assert statuses == {("HTT", "S3"): "pathogenic"}

    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(
        fam, members, [m.sample for m in members], diagnostics=statuses
    )
    htt = next(r for r in rows if r["locus"] == "HTT")
    s3 = next(i for i, m in enumerate(members) if m.sample == "S3")
    labels = [b["label"] for b in htt[f"{str_eh.member_key(s3)}_badges"]]
    assert "Dx P" in labels


def test_locus_gene():
    assert [str_eh.locus_gene(x) for x in ("HTT", "FXN_A", "NUTM2B-AS1")] == [
        "HTT",
        "FXN",
        "NUTM2B-AS1",
    ]
