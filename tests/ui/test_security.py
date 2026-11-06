"""Login-gated data files by directory and role, password migration, backoff."""

import sys

import pytest
import yaml
from nicegui import ui
from ui_helpers import capture_tables, login

import datagen
from genetics_viz.utils import login_throttle

#: A file under the data directory, served at /data-0/... .
PEDIGREE = "/data-0/cohorts/DEMO/DEMO.pedigree.tsv"


def _sample_file(data_dir, sample, *parts):
    """``/data-0/samples/<shards>/<sample>/<parts>``."""
    rel = datagen.sample_dir(data_dir, sample).relative_to(data_dir)
    return "/data-0/" + "/".join([rel.as_posix(), *parts])


async def test_anonymous_data_request_is_refused(user):
    response = await user.http_client.get(PEDIGREE)
    assert response.status_code == 401
    assert "S1" not in response.text


async def test_curator_data_request_supports_ranges(user):
    await login(user, "curator")
    whole = await user.http_client.get(PEDIGREE)
    assert whole.status_code == 200
    assert whole.headers["cache-control"].startswith("private")
    part = await user.http_client.get(PEDIGREE, headers={"Range": "bytes=0-9"})
    assert part.status_code == 206 and part.content == whole.content[:10]


@pytest.mark.parametrize("role", ["reader", "student"])
async def test_readers_get_only_what_the_pages_load(user, data_dir, role):
    """Tables, curation and QC files are refused; IGV tracks and images are not."""
    await login(user, role)
    for refused in (
        PEDIGREE,
        "/data-0/cohorts/DEMO/qc/DEMO.somalier_hg38.pairs.tsv",
        "/data-0/validations/snvs.tsv",
        _sample_file(data_dir, "S3", "sequences", "S3.GRCh38_GIABv3.cram.json"),
    ):
        assert (await user.http_client.get(refused)).status_code == 403, refused
    cram = _sample_file(data_dir, "S3", "sequences", "S3.GRCh38_GIABv3.cram")
    part = await user.http_client.get(cram, headers={"Range": "bytes=0-3"})
    assert part.status_code == 206 and part.content == b"CRAM"
    image = _sample_file(
        data_dir, "S3", "str", "expansionhunter", "reviewer",
        f"S3.{datagen.EH_LABEL}.HTT.svg",
    )  # fmt: skip
    assert (await user.http_client.get(image)).status_code == 200


async def test_path_outside_the_data_directory_is_not_served(user):
    # %2E%2E reaches the route undecoded by the client; ``..`` would not.
    escape = "/data-0/%2E%2E/genetics_viz.yaml"
    assert (await user.http_client.get(escape)).status_code == 401
    await login(user, "curator")
    response = await user.http_client.get(escape)
    assert response.status_code == 404 and "storage_secret" not in response.text


async def test_removed_user_loses_the_session(user):
    await login(user, "reader")
    config = sys.modules["genetics_viz.config_model"].get_config()
    config.user_list = [u for u in config.user_list if u.username != "reader"]
    assert (await user.http_client.get(PEDIGREE)).status_code == 401
    await user.open("/")
    await user.should_see("Sign in to continue")


async def test_role_change_applies_at_the_next_request(user):
    await login(user, "reader")
    assert (await user.http_client.get(PEDIGREE)).status_code == 403
    config = sys.modules["genetics_viz.config_model"].get_config()
    next(u for u in config.user_list if u.username == "reader").role = "curator"
    assert (await user.http_client.get(PEDIGREE)).status_code == 200


async def test_student_gets_no_downloads(user, monkeypatch):
    await login(user, "student")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    await user.should_not_see(marker="download-wombat")
    assert tables, "the table itself is still shown"
    await user.open("/cohort/DEMO")
    await user.should_see("DEMO", retries=30)
    await user.should_not_see(marker="download-ancestry")


async def test_reader_cannot_ignore_a_validation(user, monkeypatch, data_dir):
    """#54: the switch is shown, disabled."""
    from genetics_viz.components.validation_loader import SNV_VALIDATION_HEADER
    from genetics_viz.utils.tsv import append_tsv_row

    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    row = next(r for r in tables[0]["rows"] if r["VEP_SYMBOL"] == "SHANK3")
    append_tsv_row(
        data_dir / "validations" / "snvs.tsv",
        SNV_VALIDATION_HEADER,
        ["1001", row["Variant"], row["sample"], "curator", "", "present", "", "0",
         "2026-01-01T00:00:00"],
    )  # fmt: skip
    with user.client:
        tables[0]["on_row_action"]({"action": "view_variant", "row": row})
    await user.should_see(marker="validation-ignore", retries=30)
    switch = user.find(marker="validation-ignore").elements.pop()
    assert switch.props.get("disable") is True


@pytest.mark.parametrize("role, shown", [("reader", False), ("curator", True)])
async def test_trash_button_is_for_curators(user, data_dir, role, shown):
    """#55: moving a file to validate to the trash is for curators."""
    path = data_dir / "to_validate" / "mine.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("FID\tVariant\tSample\n1001\tchr1:1:A:T\tS1\n")
    await login(user, role)
    await user.open("/validation/file/mine")
    await user.should_see("Validating: mine", retries=30)
    if shown:
        await user.should_see(marker="validation-file-trash")
    else:
        await user.should_not_see(marker="validation-file-trash")


async def test_legacy_hash_is_upgraded_at_login(user, config_file):
    def stored(name):
        users = yaml.safe_load(config_file.read_text())["user_list"]
        return next(u["password"] for u in users if u["username"] == name)

    assert not stored("reader").startswith("scrypt$")
    await login(user, "reader")
    assert stored("reader").startswith("scrypt$")
    # The other accounts are untouched until they log in.
    assert not stored("curator").startswith("scrypt$")


async def _fail(user, username):
    user.find("Username").clear().type(username)
    user.find("Password").type("not-the-password")
    user.find(kind=ui.button, content="Sign in").click()


async def test_repeated_failures_back_off(user):
    await user.open("/login")
    for _ in range(login_throttle.FREE_FAILURES_USER + 1):
        await _fail(user, "reader")
    await _fail(user, "reader")
    await user.should_see("Too many failed attempts")
    # Even the right password is refused during the wait.
    user.find("Password").type("reader-pw")
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("Too many failed attempts")
    await user.should_not_see("Welcome to Genetics-Viz")
