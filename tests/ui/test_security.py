"""Login-gated data files, password migration and login backoff."""

import yaml
from nicegui import ui
from ui_helpers import login

from genetics_viz.utils import login_throttle

#: A file under the data directory, served at /data-0/... .
PEDIGREE = "/data-0/cohorts/DEMO/DEMO.pedigree.tsv"


async def test_anonymous_data_request_is_refused(user):
    response = await user.http_client.get(PEDIGREE)
    assert response.status_code == 401
    assert "S1" not in response.text


async def test_logged_in_data_request_supports_ranges(user):
    await login(user, "reader")
    whole = await user.http_client.get(PEDIGREE)
    assert whole.status_code == 200
    assert whole.headers["cache-control"].startswith("private")
    part = await user.http_client.get(PEDIGREE, headers={"Range": "bytes=0-9"})
    assert part.status_code == 206 and part.content == whole.content[:10]


async def test_path_outside_the_data_directory_is_not_served(user):
    # %2E%2E reaches the route undecoded by the client; ``..`` would not.
    escape = "/data-0/%2E%2E/genetics_viz.yaml"
    assert (await user.http_client.get(escape)).status_code == 401
    await login(user, "reader")
    response = await user.http_client.get(escape)
    assert response.status_code == 404 and "storage_secret" not in response.text


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
