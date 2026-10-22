"""Wombat row badges and filters (utils/variant_rows.py)."""

from genetics_viz.utils import variant_rows as vr
from genetics_viz.utils.vep import VEP_CONSEQUENCES


def test_split_terms_handles_aggregates_and_transcripts():
    assert vr.split_terms("a&b, c,,.") == ["a", "b", "c"]
    assert vr.split_terms(None) == vr.split_terms("") == []


def test_vep_badges_dedupe_and_collect_unknown_terms():
    row = {
        "VEP_Consequence": "missense_variant&missense_variant,made_up_term",
        "VEP_CLIN_SIG": "pathogenic,Pathogenic,.,weird_sig",
    }
    cons, clin = set(), set()
    vr.add_vep_badges(row, cons, clin)
    assert [b["label"] for b in row["ConsequenceBadges"]][0] == "missense"
    assert len(row["ConsequenceBadges"]) == 2
    assert cons == {"made_up_term"} and clin == {"weird_sig"}
    assert len(row["ClinVarBadges"]) == 3  # case variants keep their label


def test_score_badges_only_for_scored_numeric_columns():
    row = {"mpc": "2.0", "MPC_other": ".", "sample": "S1", "POS": "12"}
    vr.add_score_badges(row)
    assert row["mpc_badge"]["label"] == "2.000"
    assert "MPC_other_badge" not in row and "POS_badge" not in row


def test_exclusion_predicates():
    assert vr.in_lcr({"LCR": "True"}) and not vr.in_lcr({"LCR": "false"})
    assert vr.gnomad_filtered({"genomes_filters": "AC0"})
    assert not vr.gnomad_filtered({"genomes_filters": None})
    assert vr.in_gnomad_wgs({"fafmax_faf95_max_genomes": 0.01})
    assert vr.likely_mnv({"mnv_proba": "0.1, 0.9"})
    assert not vr.likely_mnv({"mnv_proba": "0.1,x"}) and not vr.likely_mnv({})


def test_consequence_selection():
    row = {"VEP_Consequence": "missense_variant&splice_region_variant"}
    assert vr.consequences_selected(row, ["missense_variant"])
    assert not vr.consequences_selected(row, ["stop_gained"])
    # No selection, or every term, filters nothing - a row without terms too.
    assert vr.consequences_selected({}, [])
    assert vr.consequences_selected({}, list(VEP_CONSEQUENCES))
    assert not vr.consequences_selected({}, ["stop_gained"])
