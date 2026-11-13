"""Shared helpers of the web UI's tests (tests/pyazr/ui/*_test.py): copies of the test projects, an app with the
test client, the thread count, process checks, and the runner that skips (exit 77) with a reason when a package or
the engine is missing.

The thread count is the environment's (``OMP_NUM_THREADS``, ctest sets it); the tests never set it.
"""
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)

THREADS = int(os.environ.get("OMP_NUM_THREADS", "1"))
TOKEN = "test-token-of-the-ui-tests"
AUTH = {"Authorization": f"Bearer {TOKEN}"}
PROJECTS = {name: os.path.join(ROOT, "tests", "upgrade", name) for name in ("T1", "T2", "T3")}


def missing_requirements(*packages):
    """Why the tests cannot run here (None when they can): a package or the engine module missing."""
    import importlib.util
    for name in ("pytest", "fastapi", "httpx", "uvicorn", "numpy", "scipy") + packages:
        if importlib.util.find_spec(name) is None:
            return f"the package {name} is not installed (pip install pyazr[web])"
    try:
        from pyazr.described import _engine
        _engine()
    except Exception as e:  # noqa: BLE001
        return f"the engine module _azure2 is not importable ({e})"
    return None


def run(test_file):
    """Run ``test_file`` with pytest; exit 77 (ctest: skipped) when it cannot run here."""
    reason = missing_requirements()
    if reason:
        print(f"skip: {reason}")
        sys.exit(77)
    import pytest
    code = pytest.main([test_file, "-q", "-p", "no:cacheprovider", "-W", "ignore::DeprecationWarning"])
    print("ok" if code == 0 else f"FAILED (pytest exit {code})")
    sys.exit(int(code))


def copy_project(name, dest):
    """A copy of the test project ``name`` (T1, T2, T3) in ``dest``; returns the path of its .azr."""
    target = os.path.join(str(dest), name)
    shutil.copytree(PROJECTS[name], target, ignore=shutil.ignore_patterns("output", ".*.ui*", "checks"))
    for directory in ("output", "checks"):  # the program refuses to run without them
        os.makedirs(os.path.join(target, directory), exist_ok=True)
    return os.path.join(target, f"{name}.azr")


def make_client(path, threads=THREADS, **kwargs):
    """The test client of an app serving ``path`` (use it as a context manager: the workers stop at its end)."""
    from fastapi.testclient import TestClient
    from pyazr.ui.server import create_app
    app = create_app(path, threads=threads, token=TOKEN, **kwargs)
    return TestClient(app, base_url="http://127.0.0.1", headers=AUTH)


def gone(pid):
    """True when the process ``pid`` no longer runs (or is a zombie no one waited for)."""
    try:
        with open(f"/proc/{pid}/stat") as f:
            return f.read().split(")")[-1].split()[0] == "Z"
    except FileNotFoundError:
        return True
    except OSError:
        try:
            os.kill(pid, 0)
            return False
        except ProcessLookupError:
            return True


def wait_until(predicate, timeout=120.0, step=0.05):
    """Poll ``predicate`` until it is true or ``timeout`` seconds pass; returns its last value."""
    end = time.monotonic() + timeout
    while True:
        value = predicate()
        if value or time.monotonic() > end:
            return value
        time.sleep(step)


def wait_job(client, job_id, timeout=600.0):
    """The job's record once its state is final."""
    from pyazr.ui.jobs import FINAL
    rec = {}

    def final():
        rec.update(client.get(f"/api/jobs/{job_id}").json())
        return rec["state"] in FINAL

    assert wait_until(final, timeout), f"job {job_id} did not end: {rec}"
    return dict(rec)


def run_program(azr):
    """The command-line program on ``azr`` (menu 1: calculate with data, then quit) at the environment's threads;
    returns the text of its chiSquared.out."""
    binary = os.environ.get("AZURE2_BIN")
    proc = subprocess.run([binary, "--no-readline", os.path.basename(azr)], input="1\n\n\n7\n", text=True,
                          cwd=os.path.dirname(azr), capture_output=True, timeout=600)
    path = os.path.join(os.path.dirname(azr), "output", "chiSquared.out")
    assert os.path.exists(path), f"the program wrote no chiSquared.out (exit {proc.returncode}): {proc.stdout[-2000:]}"
    with open(path) as f:
        return f.read()
