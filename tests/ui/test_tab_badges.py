"""The family page's review badges on its tabs, by role."""

from ui_helpers import badge_text, login


async def _open(user, username, path="/cohort/DEMO/family/1001"):
    await login(user, username)
    await user.open(path)
    await user.should_see("1001", retries=10)


async def test_curator_sees_what_waits_for_a_review(user, data_dir):
    await _open(user, "curator")
    # 1001: NOT_PARENT (high) and S3's contamination (medium) not acknowledged.
    assert await badge_text(user, "qc-tab-count") == "2"
    # The NOT_PARENT flag holds the de novo table back: no count, a warning.
    assert await badge_text(user, "wombat-tab-count") == "!"
    # Two genotypes flagged by the pipeline, none validated.
    assert await badge_text(user, "str-tab-count") == "2"


async def test_reader_sees_the_qc_badge_only(user, data_dir):
    await _open(user, "reader")
    assert await badge_text(user, "qc-tab-count") == "2"
    await user.should_not_see(marker="wombat-tab-count")
    await user.should_not_see(marker="str-tab-count")


async def test_denovo_calls_are_counted_once_the_table_is_read(user, data_dir):
    with open(data_dir / "validations" / "sampleqc.tsv", "a") as f:
        f.write(
            "\t".join(["NOT_PARENT", "1001", "S1", "S3", "explained", "Adopted",
                       "curator", "2026-11-10T10:00:00"]) + "\n"
        )  # fmt: skip
    await _open(user, "curator")
    # S3's SHANK3 missense de novo, MODERATE, not yet validated.
    assert await badge_text(user, "wombat-tab-count") == "1"
    assert await badge_text(user, "qc-tab-count") == "1"


async def test_standalone_family_page_has_the_badges(user, data_dir):
    await _open(user, "curator", "/family/1001")
    assert await badge_text(user, "qc-tab-count") == "2"
    assert await badge_text(user, "str-tab-count") == "2"
