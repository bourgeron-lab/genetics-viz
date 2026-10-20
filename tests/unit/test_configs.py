"""The YAML/TSV configuration shipped in src/genetics_viz/config.

These files are loaded once at import and edited by hand; a typo shows up as a
wrong colour or a KeyError deep in a page. Checking them here catches it on the
pull request instead.
"""

import re
from pathlib import Path

import pytest
import yaml

import genetics_viz

CONFIG = Path(genetics_viz.__file__).parent / "config"
HEX = re.compile(r"^#[0-9a-fA-F]{6}$")


def load(name):
    with open(CONFIG / name) as f:
        return yaml.safe_load(f)


def test_vep_terms_have_a_valid_impact_and_colour():
    for term, spec in load("vep_consequences.yaml").items():
        assert spec["impact"] in {"HIGH", "MODERATE", "LOW", "MODIFIER"}, term
        assert HEX.match(spec["color"]), term


def test_clinvar_and_region_colours_are_hex():
    assert all(HEX.match(c) for c in load("clinvar_colors.yaml").values())
    regions = load("ancestry.yaml")["region_colors"]
    assert "unknown" in regions and all(HEX.match(c) for c in regions.values())


def test_roh_grades_and_kinship_degrees_are_ordered():
    cfg = load("roh.yaml")
    grades = cfg["froh_grades"]
    categories = [g["category"] for g in grades]
    assert len(set(categories)) == len(categories) and "none" in categories
    thresholds = [g["threshold"] for g in grades if "threshold" in g]
    assert thresholds == sorted(thresholds, reverse=True)
    mins = [d["min"] for d in cfg["kinship_degrees"]]
    assert mins == sorted(mins, reverse=True)
    assert set(cfg["patterns"]) == {
        "typical",
        "population",
        "recent",
        "recent+population",
    }
    assert set(cfg["classes"]) == {"short", "intermediate", "long"}


def test_str_statuses_cover_the_pipeline_vocabulary():
    cfg = load("str.yaml")
    statuses = cfg["statuses"]
    # eh_family_table's classes and statuses, plus "." (no thresholds row)
    assert set(statuses) == {
        "pathogenic", "carrier", "intermediate", "above_normal", "normal", "unclassified"
    }  # fmt: skip
    assert all(HEX.match(v["color"]) for v in statuses.values())
    ranks = [v["rank"] for v in statuses.values()]
    assert len(set(ranks)) == len(ranks)
    assert set(cfg["support"]) == {"SPANNING", "FLANKING", "INREPEAT"}
    assert 0 < cfg["common_in_cohort_fraction"] < 1
    assert {"{sample}", "{label}", "{locus}"} <= set(
        re.findall(r"{\w+}", cfg["reviewer_filename"])
    )


def test_sample_qc_vocabulary_covers_the_pipeline():
    cfg = load("sample_qc.yaml")
    # sampleqc_classify's flag codes (KARYOTYPE_<k> goes by prefix).
    assert set(cfg["flags"]) == {
        "SEX_MISMATCH", "SEX_AMBIGUOUS", "SEX_UNKNOWN_IN_PED", "X_HET_INCONSISTENT",
        "LOW_INFORMATIVE_SITES", "POSSIBLE_CONTAMINATION",
        "SAMPLE_UNRELATED_TO_FAMILY", "DUPLICATE_OR_MZ", "NOT_PARENT",
        "HALF_SIB_NOT_FULL", "FULL_SIB_NOT_HALF", "UNEXPECTED_RELATED_IN_FAMILY",
        "EXPECTED_RELATED_UNRELATED", "RELATIONSHIP_MISMATCH",
        "MENDELIAN_INCONSISTENT", "CRYPTIC_RELATED", "SAMPLE_SWAP_CANDIDATE",
    }  # fmt: skip
    categories = set(cfg["categories"])
    assert {f["category"] for f in cfg["flags"].values()} <= categories
    assert cfg["karyotype_flag"]["category"] in categories
    assert set(cfg["denovo_risk"]) <= set(cfg["flags"])
    ranks = [s["rank"] for s in cfg["severities"].values()]
    assert set(cfg["severities"]) == {"high", "medium", "info"}
    assert ranks == sorted(ranks, reverse=True)
    assert set(cfg["relationships"]) == {
        "duplicate", "parent-offspring", "full-sibs", "first-degree", "half-sibs",
        "second-degree", "third-degree", "unrelated", "low_confidence",
    }  # fmt: skip
    assert set(cfg["pair_status"]) == {"ok", "mismatch", "low_confidence"}
    colors = [c["color"] for c in cfg["categories"].values()]
    colors += [v["color"] for v in cfg["relationships"].values()]
    assert all(HEX.match(c) for c in colors)


def test_wisecondorx_thresholds_bracket_zero():
    cfg = load("wisecondorx_thresholds.yaml")
    assert (
        cfg["robust_loss"]["ratio_threshold"]
        <= cfg["permissive_loss"]["ratio_threshold"]
        < 0
    )
    assert (
        0
        < cfg["permissive_gain"]["ratio_threshold"]
        <= cfg["robust_gain"]["ratio_threshold"]
    )


def test_continuous_score_points_ascend():
    for name, spec in load("continuous_scores.yaml")["scores"].items():
        thresholds = [p["threshold"] for p in spec["points"]]
        assert thresholds == sorted(thresholds), name
        assert all(HEX.match(p["color"]) and p.get("label") for p in spec["points"]), (
            name
        )


def test_column_types_are_known():
    from genetics_viz.utils.column_names import _TYPE_MAP

    for col, spec in load("column_names.yaml").items():
        assert spec.get("type", "string") in _TYPE_MAP, col


def test_genetrek_config_is_well_formed():
    cfg = load("genetrek.yaml")
    assert len(cfg["tier_colors"]) == 4
    assert all(HEX.match(c) for c in cfg["tier_colors"] + [cfg["neutral_color"]])
    for key, value in cfg["default_sets"].items():
        assert value == "show" or value in (1, 2, 3, 4), key
    assert float(cfg["refresh_hours"]) > 0


def test_pipeline_steps_are_complete():
    cfg = load("pipeline_steps.yaml")
    for step, spec in cfg["steps"].items():
        assert spec.get("label") and spec.get("unit"), step
    for status, spec in cfg["run_status"].items():
        assert {"label", "color", "icon"} <= set(spec), status


def test_view_presets_have_columns():
    presets = load("view_presets.yaml")["presets"]
    assert presets and all(p.get("name") and p.get("columns") for p in presets)


@pytest.mark.parametrize("name", ["validation_badges.yaml", "diagnostic_badges.yaml"])
def test_badge_configs_load(name):
    assert load(name)


def test_cytobands_cover_every_chromosome():
    from genetics_viz.utils.cytobands import CHROM_ORDER, CYTOBANDS

    assert set(CHROM_ORDER) - {"MT"} <= set(CYTOBANDS)
