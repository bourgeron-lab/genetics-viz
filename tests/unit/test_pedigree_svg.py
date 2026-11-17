"""Pedigree drawings through pypedigree: structure only, never pixels."""

import re
import warnings
import xml.etree.ElementTree as ET

from genetics_viz.utils.pedigree_svg import (
    PedigreePerson as P,
    pedigree_svg,
    pedigree_text,
)

SVG = "{http://www.w3.org/2000/svg}"

TRIO = [P("S1", sex="1", phenotype="1"), P("S2", sex="2", phenotype="1"),
        P("S3", "S1", "S2", "1", "2")]  # fmt: skip


def _nodes(svg: str) -> dict:
    """``data-id`` -> the node's attributes."""
    root = ET.fromstring(svg.split("?>", 1)[-1])
    return {
        g.get("data-id"): g.attrib
        for g in root.iter(f"{SVG}g")
        if "node" in (g.get("class") or "").split()
    }


def _titles(svg: str) -> dict:
    root = ET.fromstring(svg.split("?>", 1)[-1])
    return {
        g.get("data-id"): g.find(f"{SVG}title").text
        for g in root.iter(f"{SVG}g")
        if g.get("data-id")
    }


def test_trio_is_a_static_svg_under_its_prefix():
    drawing = pedigree_svg("1001", TRIO, id_prefix="pedigree-")
    assert drawing is not None
    svg = drawing.svg
    assert "<script" not in svg
    assert 'class="pypedigree"' in svg and 'width="100%"' in svg
    ids = re.findall(r'\sid="([^"]+)"', svg)
    assert ids and all(i.startswith("pedigree-") for i in ids)
    nodes = _nodes(svg)
    assert set(nodes) == {"S1", "S2", "S3"}
    assert nodes["S3"]["data-parents"] == "S1,S2"
    assert drawing.width > 0 and drawing.height > 0
    assert not drawing.has_unsampled and not drawing.has_notes


def test_one_parent_gets_a_dashed_partner_and_stays_connected():
    # pypedigree's builder left this child unconnected; the parser does not.
    drawing = pedigree_svg("AU-X-1", [P("S5", sex="1"), P("S6", father="S5", sex="2")])
    parents = _nodes(drawing.svg)["S6"]["data-parents"].split(",")
    assert parents[0] == "S5" and len(parents) == 2
    assert drawing.has_unsampled


def test_dummies_and_absent_parents_are_drawn_not_sampled_without_warning():
    people = [
        P("S9", "DUMMY_S9_F", "DUMMY_S9_M", "2", "2"),
        P("DUMMY_S9_M", sex="2", phenotype="2", dummy=True),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        drawing = pedigree_svg("2002", people)
    titles = _titles(drawing.svg)
    assert set(titles) == {"S9", "DUMMY_S9_M", "DUMMY_S9_F"}
    assert "Not sampled" in titles["DUMMY_S9_F"]  # named, but no row
    assert "Not sampled" in titles["DUMMY_S9_M"]
    # A dummy has no status, whatever its row says.
    assert "Phenotype: Unknown" in titles["DUMMY_S9_M"]
    assert "Not sampled" not in titles["S9"]
    assert drawing.has_unsampled


def test_notes_and_ids_are_escaped():
    drawing = pedigree_svg(
        "F", [P("<b>&x", sex="1", note='XXY<i>"'), P("kid", father="<b>&x")]
    )
    _nodes(drawing.svg)  # parses
    assert "<b>" not in drawing.svg and "&lt;b&gt;&amp;x" in drawing.svg
    assert drawing.has_notes


def test_whitespace_in_a_cell_cannot_shift_the_columns():
    text = pedigree_text("F A", [P("a\tb", note="two words")])
    header, row = text.splitlines()
    assert len(row.split("\t")) == len(header.split("\t"))
    assert row.startswith("F_A\ta_b\t") and row.endswith("two_words")


def test_long_ids_widen_the_layout():
    short = pedigree_svg("F", [P("A", sex="1"), P("B", sex="2"), P("C", "A", "B")])
    long = pedigree_svg(
        "F",
        [P("DUMMY_B00IDL0_MF", sex="1", dummy=True),
         P("DUMMY_B00IDL0_MM", sex="2", dummy=True),
         P("C", "DUMMY_B00IDL0_MF", "DUMMY_B00IDL0_MM")],
    )  # fmt: skip

    def span(drawing):
        xs = [float(n["data-x"]) for n in _nodes(drawing.svg).values()]
        return max(xs) - min(xs)

    # The couple sits further apart than a 16-character label is wide (about
    # 2.1 units at 10 pt), so the two labels clear each other.
    assert span(short) < 2.1 < span(long)


def test_nobody_gives_no_drawing():
    assert pedigree_svg("F", []) is None
