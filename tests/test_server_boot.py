"""Start the real server as production does, and check it answers.

The user-simulation tests run the app in-process; this one runs the console
entry point in a subprocess, so it also catches packaging and startup problems
(a missing config file in the wheel, an import that only fails at boot). The CI
``package`` job runs it against the built wheel.
"""

import os
import socket
import subprocess
import sys
import time

import httpx
import pytest

pytestmark = pytest.mark.slow


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_server_boots_and_serves(config_file, tmp_path):
    port = free_port()
    # NiceGUI treats a process with PYTEST_CURRENT_TEST set as its own screen
    # test and looks for a port it never gets; the child is a plain server.
    env = {k: v for k, v in os.environ.items() if k != "PYTEST_CURRENT_TEST"}
    env["NICEGUI_STORAGE_PATH"] = str(tmp_path / ".nicegui")
    proc = subprocess.Popen(
        [sys.executable, "-m", "genetics_viz", str(config_file), "--port", str(port)],
        env=env,
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    base = f"http://127.0.0.1:{port}"
    try:
        deadline = time.monotonic() + 60
        while True:
            if proc.poll() is not None:
                pytest.fail(f"server exited early:\n{proc.stdout.read()}")
            try:
                health = httpx.get(f"{base}/health", timeout=2)
                if health.status_code == 200:
                    break
            except httpx.TransportError:
                pass
            if time.monotonic() > deadline:
                pytest.fail("server never became healthy")
            time.sleep(0.5)

        assert health.json()["status"] == "ok"
        assert "Sign in to continue" in httpx.get(f"{base}/login", timeout=10).text
        root = httpx.get(f"{base}/", timeout=10, follow_redirects=False)
        assert root.status_code in (302, 303, 307) and root.headers[
            "location"
        ].endswith("/login")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()
