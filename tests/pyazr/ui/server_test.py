#!/usr/bin/env python3
"""The web UI's server as a real process on a free port above 20000, used over HTTP with httpx, on a copy of T2:

  1. started with ``python -m pyazr.ui``, it prints its URL with the token, serves the session of the user guide's
     page "The web UI's API" (the code block is read from docs/source/user_guide/web_ui_api.rst and run), and is not
     reachable on a non-loopback address of the machine;
  2. killed (SIGKILL) while a job computes, its session worker and job worker end with it (checked by PID); started
     again on the same project, it lists the job as "interrupted" with its checkpoint;
  3. stopped (SIGTERM), it ends with every worker (checked by PID);
  4. asked to bind to a non-loopback address, it refuses and binds nothing.

Every process started here is stopped and checked gone by PID.
Run: python3 tests/pyazr/ui/server_test.py
"""
import os
import re
import signal
import socket
import subprocess
import sys
import time

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import AUTH, ROOT, TOKEN, copy_project, gone, wait_until  # noqa: E402

GUIDE = os.path.join(ROOT, "docs", "source", "user_guide", "web_ui_api.rst")


def start_server(path, port, *extra):
    from pyazr.ui.server import free_port
    port = port or free_port()
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([ROOT] + [p for p in os.environ.get("PYTHONPATH", "").split(
        os.pathsep) if p]))
    log = open(os.path.join(os.path.dirname(path), f"server-{port}.log"), "ab")
    proc = subprocess.Popen([sys.executable, "-m", "pyazr.ui", path, "--port", str(port), "--token", TOKEN, *extra],
                            stdout=subprocess.PIPE, stderr=log, env=env, text=True)
    log.close()
    return proc, port


def connect(proc, port):
    import httpx
    url = proc.stdout.readline().strip()
    assert url == f"http://127.0.0.1:{port}/?token={TOKEN}", url
    client = httpx.Client(base_url=f"http://127.0.0.1:{port}", headers=AUTH, timeout=600)

    def up():
        try:
            return client.get("/api/describe").status_code == 200
        except httpx.TransportError:
            return False

    assert wait_until(up, 60), "the server did not answer"
    return client, url


def guide_session():
    """The python code block of the user guide page, as text."""
    with open(GUIDE) as f:
        text = f.read()
    block = text.split(".. code-block:: python", 1)[1].split("\n\n", 1)[1]
    lines = []
    for line in block.splitlines():
        if line and not line.startswith("   "):
            break
        lines.append(line[3:])
    return "\n".join(lines)


def other_address():
    """A non-loopback IPv4 address of this machine, or None (no traffic is sent)."""
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("192.0.2.1", 9))
            address = s.getsockname()[0]
    except OSError:
        return None
    return None if address.startswith("127.") else address


def test_server_process_workers_die_with_it(tmp_path):
    path = copy_project("T2", tmp_path)
    started = []
    proc, port = start_server(path, 0)
    started.append(proc.pid)
    try:
        api, url = connect(proc, port)
        assert api.post("/api/evaluate").status_code == 200
        # the user guide's session, as written
        code = guide_session()
        assert code.count("\n") >= 9
        exec(compile(code, GUIDE, "exec"), {"url": f"http://127.0.0.1:{port}", "token": TOKEN})
        assert api.get("/api/project").json()["unsaved"] is False
        address = other_address()
        if address:
            with pytest.raises(OSError):
                socket.create_connection((address, port), timeout=2).close()
        job = api.post("/api/jobs", json={"kind": "fit", "method": "scipy_lbfgsb", "x0": "text"}).json()
        assert wait_until(lambda: (api.get(f"/api/jobs/{job['id']}").json().get("progress") or {}).get("calls"), 120)
        processes = api.get("/api/server").json()
        roles = {p["role"].split()[0]: p for p in processes}
        assert roles["server"]["pid"] == proc.pid and {"session", "job"} <= set(roles)
        print("\n  memory (MB, resident / peak): " + ", ".join(
            f"{p['role'].split()[0]} {p['rss_kb'] / 1024:.0f} / {p['peak_kb'] / 1024:.0f}" for p in processes))
        workers = [p["pid"] for p in processes if p["role"] != "server"]
        started += workers
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()
        t0 = time.perf_counter()
        assert wait_until(lambda: all(gone(p) for p in workers), 10), f"workers outlived the server: {workers}"
        print(f"  workers gone {time.perf_counter() - t0:.2f} s after the server was killed")

        proc, port = start_server(path, 0)
        started.append(proc.pid)
        api, _ = connect(proc, port)
        listed = {j["id"]: j for j in api.get("/api/jobs").json()}
        assert listed[job["id"]]["state"] == "interrupted"
        assert api.get(f"/api/jobs/{job['id']}").json()["best"]["x"]
        assert api.post("/api/evaluate").status_code == 200
        workers = [p["pid"] for p in api.get("/api/server").json() if p["role"] != "server"]
        started += workers
        proc.send_signal(signal.SIGTERM)
        assert proc.wait(30) is not None
        assert wait_until(lambda: all(gone(p) for p in workers), 10)
    finally:
        for pid in started:
            if not gone(pid):
                os.kill(pid, signal.SIGKILL)
        assert wait_until(lambda: all(gone(p) for p in started), 10), started


def test_non_loopback_host_refused(tmp_path):
    path = copy_project("T1", tmp_path)
    proc, port = start_server(path, 0, "--host", "0.0.0.0")
    try:
        out, _ = proc.communicate(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill()
    assert proc.returncode != 0 and not re.search(r"http://", out or "")
    with pytest.raises(OSError):
        socket.create_connection(("127.0.0.1", port), timeout=2).close()
    assert gone(proc.pid)


if __name__ == "__main__":
    uitest.run(__file__)
