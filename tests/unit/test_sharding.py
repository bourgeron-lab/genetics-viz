"""Two-level sharded entity directories (utils/sharding.py)."""

from genetics_viz.utils.sharding import (
    compute_shard_prefix,
    get_entity_url_segment,
    get_family_path,
    get_sample_path,
)


def test_shard_prefix_uses_last_two_characters_uppercased():
    assert compute_shard_prefix("C000EZJ") == ("J", "Z")
    assert compute_shard_prefix("ab-c.d_e") == ("E", "D")
    assert compute_shard_prefix("AU-X-1") == ("1", "X")


def test_sharded_path_wins_when_present(tmp_path):
    (tmp_path / "families" / "1" / "0" / "1001").mkdir(parents=True)
    assert (
        get_family_path(tmp_path, "1001") == tmp_path / "families" / "1" / "0" / "1001"
    )
    assert get_entity_url_segment(tmp_path, "families", "1001") == "families/1/0/1001"


def test_flat_layout_is_the_fallback(tmp_path):
    assert get_sample_path(tmp_path, "C000EZJ") == tmp_path / "samples" / "C000EZJ"
    assert get_entity_url_segment(tmp_path, "samples", "C000EZJ") == "samples/C000EZJ"
