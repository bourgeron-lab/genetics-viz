"""Sample QC against a pedigree draft: the pipeline's verdict, and the fixes it supports.

``pedigree_qc`` mirrors ghfc-ngs's ``sampleqc_classify``; the first tests hold
it to the flags the pipeline wrote for ``datagen``'s families.
"""

from pathlib import Path

import datagen
import pytest

from genetics_viz.utils import pedigree_edit as pe
from genetics_viz.utils import pedigree_qc as pq
from genetics_viz.utils.sample_qc import FlagRow, PairRow, load_family_qc, pair_rows

#: Observed kinship and IBS0 ratio of each class, as KING would give them.
OBSERVED = {
    "duplicate": (0.49, 0.0),
    "parent-offspring": (0.25, 0.0),
    "full-sibs": (0.25, 0.2),
    "second-degree": (0.12, 0.6),
    "third-degree": (0.06, 0.8),
    "unrelated": (0.0, 1.0),
}


def pair(a: str, b: str, observed: str) -> PairRow:
    phi, ratio = OBSERVED[observed]
    return PairRow("F", a, "F", b, phi, ratio, 16000, None, None, observed, None, None)


def draft_of(*rows) -> pe.PedigreeDraft:
    """Rows as ``(iid, father, mother, sex)``; phenotype unaffected."""
    return pe.PedigreeDraft("F", [pe.Member(*r, "1") for r in rows])


def codes(draft, pairs, **kw):
    return {(c.a, c.b): c.code for c in pq.check_pairs(draft, pairs, **kw)}


def titles(suggestions):
    return [s.title for s in suggestions]


# --- The pipeline's verdict ------------------------------------------------------------


@pytest.mark.parametrize("fid", ["1001", "AU-X-1"])
def test_the_checks_reproduce_the_pipelines_pair_flags(tmp_path: Path, fid):
    data = datagen.build_data_dir(tmp_path)
    source = pe.load_family(data, fid, [data / "cohorts/DEMO/DEMO.pedigree.tsv"])
    qc = load_family_qc(data, fid)
    written = {
        frozenset((f.sample_a, f.sample_b)): f.code for f in qc.flags if f.sample_b
    }
    for check in pq.check_pairs(source.original, pair_rows(qc.pairs)):
        assert check.code == written.get(frozenset((check.a, check.b))), check


def test_expected_kinship_follows_the_recursive_definition():
    draft = draft_of(
        ("GF", "0", "0", "1"), ("GM", "0", "0", "2"),
        ("A", "GF", "GM", "1"), ("B", "GF", "GM", "2"),
        ("C", "A", "0", "1"), ("D", "0", "B", "2"),  # first cousins
        ("E", "C", "D", "1"),  # their child: inbred
        ("H", "GF", "0", "2"),  # half sib of A and B
    )  # fmt: skip
    k = pq.Kinship(draft)
    assert k.kinship("A", "B") == 0.25 and k.kinship("A", "H") == 0.125
    assert k.kinship("C", "D") == 0.0625
    assert k.kinship("E", "E") == 0.5 * (1 + 0.0625)
    assert k.relationship("A", "B") == "full-sibs"
    assert k.relationship("A", "H") == "half-sibs"
    assert k.relationship("C", "D") == "third-degree"


def test_a_cyclic_draft_has_no_fit():
    draft = draft_of(("A", "B", "0", "1"), ("B", "A", "0", "1"))
    assert pq.check_pairs(draft, [pair("A", "B", "parent-offspring")]) is None
    assert pq.suggest(draft, [pair("A", "B", "parent-offspring")]) == []


def test_a_upd_pair_is_no_parent_contradiction():
    draft = draft_of(("D", "0", "0", "1"), ("M", "0", "0", "2"), ("K", "D", "M", "1"))
    pairs = [pair("D", "K", "unrelated"), pair("M", "K", "parent-offspring")]
    assert codes(draft, pairs)[("D", "K")] == "NOT_PARENT"
    assert codes(draft, pairs, upd_pairs=[("D", "K")])[("D", "K")] is None


# --- Suggestions ---------------------------------------------------------------------


def test_a_parent_who_is_not_one_gets_a_dummy():
    draft = draft_of(("D", "0", "0", "1"), ("M", "0", "0", "2"), ("K", "D", "M", "1"))
    pairs = [
        pair("D", "M", "unrelated"),
        pair("D", "K", "unrelated"),
        pair("M", "K", "parent-offspring"),
    ]
    (s,) = pq.suggest(draft, pairs)
    assert s.title == "Replace K's father D with a dummy"
    assert (s.before.contradictions, s.after.contradictions) == (1, 0)
    assert s.flags == ("NOT_PARENT",) and "swap" in s.caution
    trial = draft.copy()
    s.apply(trial)
    assert trial.get("K").father == "DUMMY_K_F"


def test_a_swap_is_offered_when_it_explains_every_pair():
    # The barcode labelled D is the child, the one labelled K the father.
    draft = draft_of(("D", "0", "0", "1"), ("M", "0", "0", "2"), ("K", "D", "M", "1"))
    pairs = [
        pair("D", "M", "parent-offspring"),
        pair("D", "K", "parent-offspring"),
        pair("K", "M", "unrelated"),
    ]
    best = [s for s in pq.suggest(draft, pairs) if s.key.startswith("swap")]
    assert titles(best)[0] == "Swap D and K (sample swap)"
    assert best[0].after == pq.Fit(0, 0)


def test_children_sharing_a_mother_with_father_written_0():
    """The case the editor is for: full or half sibs is undecided."""
    draft = draft_of(
        ("M", "0", "0", "2"),
        ("K1", "0", "M", "1"), ("K2", "0", "M", "2"), ("K3", "0", "M", "2"),
    )  # fmt: skip
    pairs = [
        pair("M", "K1", "parent-offspring"),
        pair("M", "K2", "parent-offspring"),
        pair("M", "K3", "parent-offspring"),
        pair("K1", "K2", "full-sibs"),
        pair("K1", "K3", "second-degree"),
        pair("K2", "K3", "second-degree"),
    ]
    # The pipeline flags none of it.
    assert not any(pq.check_pairs(draft, pairs)[i].code for i in range(6))
    found = titles(pq.suggest(draft, pairs))
    assert "Give K1, K2 shared dummy father (full sibs)" in found
    assert "Give K1 and K3 a dummy father each (half sibs)" in found
    # Applied: K1 and K2 share DUMMY_K1_F; K3 then gets its own.
    shared = next(s for s in pq.suggest(draft, pairs) if s.key.startswith("shared"))
    shared.apply(draft)
    assert {draft.get(k).father for k in ("K1", "K2")} == {"DUMMY_K1_F"}
    nxt = titles(pq.suggest(draft, pairs))
    assert "Give K1, K2 shared dummy father (full sibs)" not in nxt
    assert "Give K3 a dummy father of their own (half sib of K1)" in nxt


def test_siblings_with_no_parents_get_shared_dummies():
    draft = draft_of(("A", "0", "0", "1"), ("B", "0", "0", "2"))
    pairs = [pair("A", "B", "full-sibs")]
    assert codes(draft, pairs)[("A", "B")] == "RELATED_BEYOND_PEDIGREE"
    (s,) = [s for s in pq.suggest(draft, pairs) if s.key.startswith("shared")]
    assert s.title == "Give A, B shared dummy parents (full sibs)"
    assert (s.before.unexplained, s.after.unexplained) == (1, 0)


def test_half_sibs_declared_full_get_a_dummy_of_their_own():
    draft = draft_of(
        ("D", "0", "0", "1"), ("M", "0", "0", "2"),
        ("A", "D", "M", "1"), ("B", "D", "M", "2"),
    )  # fmt: skip
    pairs = [
        pair("D", "A", "parent-offspring"),
        pair("M", "A", "parent-offspring"),
        pair("M", "B", "parent-offspring"),
        pair("D", "B", "unrelated"),
        pair("A", "B", "second-degree"),
        pair("D", "M", "unrelated"),
    ]
    assert codes(draft, pairs)[("A", "B")] == "HALF_SIB_NOT_FULL"
    found = pq.suggest(draft, pairs)
    # Only B's father fits D–B and A–B at once.
    assert titles(found)[0] in (
        "Replace B's father D with a dummy",
        "Give B a dummy father of their own (half sibs)",
    )
    assert found[0].after.contradictions == 0


def test_sex_from_the_karyotype_with_a_caution_when_relationships_fail():
    draft = draft_of(("D", "0", "0", "2"), ("M", "0", "0", "2"), ("K", "D", "M", "1"))
    samples = [
        {"sample": "D", "karyotype": "XY", "sex_evidence": "coverage"},
        {"sample": "M", "karyotype": "XX"},
        {"sample": "K", "karyotype": "XXY"},  # atypical: nothing suggested
    ]
    pairs = [pair("D", "K", "unrelated"), pair("M", "K", "parent-offspring")]
    sex = [s for s in pq.suggest(draft, pairs, samples) if s.kind == "sex"]
    assert titles(sex) == ["Record D as male"]
    assert "sample swap" in sex[0].caution and sex[0].flags == ("SEX_MISMATCH",)


def test_qc_context_reads_upd_flags():
    flag = FlagRow("family", "F", "D", "K", "UPD_CANDIDATE", "medium", "chr2")
    ctx = pq.QcContext(pairs=[pair("D", "K", "unrelated")], flags=[flag])
    draft = draft_of(("D", "0", "0", "1"), ("K", "D", "0", "1"))
    assert ctx.upd_pairs == [("D", "K")]
    assert ctx.fit(draft) == pq.Fit(0, 0)
    assert pq.QcContext().fit(draft) is None
