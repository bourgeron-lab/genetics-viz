"""The cohort SNVs/Indels tab's rows, filter and recurrent genes."""

from genetics_viz.utils import cohort_snvs as cs
from genetics_viz.utils.denovo import Curation


def _dnm(fid, child, pos, genes, consequence, strict="yes", flags="."):
    return {
        "family_id": fid, "child": child, "chrom": "chr1", "pos": str(pos),
        "ref": "A", "alt": "T", "var_type": "SNV", "genes": genes,
        "consequences": consequence, "strict": strict, "flags": flags,
    }  # fmt: skip


ROWS = [
    _dnm("F1", "C1", 10, "SHANK3", "stop_gained"),
    _dnm("F2", "C2", 20, "SHANK3,ARID1B", "missense_variant"),
    _dnm("F3", "C3", 30, "ARID1B", "synonymous_variant"),
    _dnm("F3", "C3", 40, "ARID1B", "frameshift_variant", strict="no",
         flags="PARENT_ALT_POSSIBLE"),
]  # fmt: skip


def test_split_genes_and_impact():
    assert cs.split_genes("SHANK3,ARID1B&SHANK3;.") == ["SHANK3", "ARID1B"]
    assert cs.impact_of("synonymous_variant&stop_gained") == "HIGH"
    assert cs.impact_of("") == ""


def test_candidate_rows_carry_curation_and_follow_the_children():
    curation = Curation(
        validations={("F1", "chr1:10:A:T", "C1"): ("present", "de novo")},
        diagnostics={("F1", "chr1:10:A:T", "C1"): "pathogenic"},
    )
    rows = cs.candidate_rows(ROWS, {("F1", "C1"), ("F2", "C2")}, curation)
    assert [r["Variant"] for r in rows] == ["chr1:10:A:T", "chr1:20:A:T"]
    first = rows[0]
    assert (first["impact"], first["Validation"], first["Diagnostic"]) == (
        "HIGH",
        "present",
        "pathogenic",
    )
    assert first["Validation_badges"][0]["label"] == "present · de novo"
    assert rows[1]["Validation"] == "TODO" and rows[1]["_genes"] == ["SHANK3", "ARID1B"]


def test_filter_defaults_to_strict_high_and_moderate():
    rows = cs.candidate_rows(ROWS, None)
    assert [r["_pos"] for r in cs.filter_candidates(rows)] == [10, 20]
    loose = cs.filter_candidates(rows, strict_only=False, impacts=())
    assert len(loose) == 4
    assert [r["_pos"] for r in cs.filter_candidates(rows, statuses=["present"])] == []


def test_recurrent_genes_need_two_families():
    rows = cs.candidate_rows(ROWS, None)
    genes = cs.recurrent_genes(cs.filter_candidates(rows))
    assert [(g["gene"], g["n_families"]) for g in genes] == [("SHANK3", 2)]
    assert genes[0]["impact"] == "HIGH" and genes[0]["families"] == "F1, F2"
    # Without the filters, ARID1B is hit in F2 and F3 (three candidates).
    every = cs.recurrent_genes(rows)
    assert [(g["gene"], g["n_candidates"]) for g in every] == [
        ("ARID1B", 3),
        ("SHANK3", 2),
    ]
