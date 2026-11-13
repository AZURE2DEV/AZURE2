"""The browser smoke test of the page (docs/upgrade/web_ui.md section 9), with the headless Chrome of the machine: a
real server (``python -m pyazr.ui``) on copies of T1 and T2, the page's self-test (``#selftest``: every tab, a level,
a data set's form and its plot drawn by plotly.js and as a table, a record kind without layout hints, a new record's
form, a width typed into its box and taken back with the Undo button) read from the document Chrome dumps, and one
screenshot per tab.  Fails on any script error, console error or
failed request the self-test recorded.

Environment:
- ``CHROME``: the browser (default: google-chrome, chromium or chromium-browser on PATH); without one, prints
  ``skip: ...`` and exits 77.  Each run of it gets a fresh profile directory under the test's temporary directory.
- ``AZURE2_UI_SCREENSHOTS``: where the screenshots go (default: the temporary directory, removed afterwards);
  named ``<project>_<theme>_<width>x<height>_<view>.png``.
- ``AZURE2_UI_SHOTS=all``: both themes at 1400x900 and 1000x700 (default: light at 1400x900).
- ``AZURE2_UI_SMOKE_PROJECTS``: more projects (``.azr`` paths separated by ``:``; each one's directory is copied).

Screenshots use Chrome's old headless mode, whose viewport is the window size (the new mode's is about 90 pixels
shorter); the self-test uses the new mode.
"""
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from html import unescape

import pytest

import uitest
from uitest import TOKEN, copy_project, gone

CHROME = os.environ.get("CHROME") or next((shutil.which(n) for n in ("google-chrome", "chromium", "chromium-browser")
                                           if shutil.which(n)), None)
VIEWS = ("project", "nucleus", "nucleus/level/{level}", "data", "data/{dataset}", "fit", "results", "physics", "thm")


def extra_projects():
    return [p for p in os.environ.get("AZURE2_UI_SMOKE_PROJECTS", "").split(":") if p]


def copy_any(azr, dest):
    """A copy of the directory of ``azr`` (without the UI's and the engine's outputs) in ``dest``; its .azr path."""
    source = os.path.dirname(os.path.abspath(azr))
    target = os.path.join(str(dest), os.path.basename(source))
    shutil.copytree(source, target, ignore=shutil.ignore_patterns("output", ".*.ui*", "checks"))
    for directory in ("output", "checks"):
        os.makedirs(os.path.join(target, directory), exist_ok=True)
    return os.path.join(target, os.path.basename(azr))


class Served:
    """``python -m pyazr.ui`` on a free port; stopped (SIGINT, then SIGKILL) and its workers checked gone."""

    def __init__(self, azr, log):
        from pyazr.ui.server import free_port
        self.port = free_port()
        self.log = open(log, "w")
        self.proc = subprocess.Popen([sys.executable, "-m", "pyazr.ui", azr, "--port", str(self.port), "--token",
                                      TOKEN], cwd=uitest.ROOT, stdout=subprocess.PIPE, stderr=self.log, text=True,
                                     env=dict(os.environ, PYTHONPATH=uitest.ROOT))
        self.url = self.proc.stdout.readline().strip()
        assert self.url.startswith("http://127.0.0.1:"), f"the server printed {self.url!r}"
        assert uitest.wait_until(self._answers, 60), "the server does not answer"

    def _answers(self):
        import httpx
        try:
            return httpx.get(f"http://127.0.0.1:{self.port}/api/project", headers=uitest.AUTH).status_code == 200
        except httpx.HTTPError:
            return False

    def pids(self):
        import httpx
        rows = httpx.get(f"http://127.0.0.1:{self.port}/api/server", headers=uitest.AUTH).json()
        return [r["pid"] for r in rows]

    def stop(self):
        pids = self.pids()
        self.proc.send_signal(signal.SIGINT)
        try:
            self.proc.wait(30)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            self.proc.wait(10)
        self.log.close()
        assert uitest.wait_until(lambda: all(gone(p) for p in pids), 20), f"processes left: {pids}"


def chrome(tmp_path, url, *args, mode="new", timeout=150):
    """Run the browser once on ``url`` with a fresh profile; its standard output."""
    profile = tmp_path / f"chrome-{time.monotonic_ns()}"
    profile.mkdir()
    cmd = [CHROME, f"--headless={mode}", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
           "--hide-scrollbars", f"--user-data-dir={profile}", *args, url]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return proc.stdout


def selftest_result(dom):
    m = re.search(r'<pre id="selftest-result"[^>]*>(.*?)</pre>', dom, re.S)
    assert m, f"the self-test did not finish within the time budget: {dom[-2000:]}"
    return json.loads(unescape(m.group(1)))


def projects():
    return ["T1", "T2"] + extra_projects()


@pytest.mark.parametrize("name", projects())
def test_selftest_and_screenshots(name, tmp_path):
    azr = copy_project(name, tmp_path) if name in uitest.PROJECTS else copy_any(name, tmp_path)
    label = os.path.splitext(os.path.basename(azr))[0]
    shots = os.environ.get("AZURE2_UI_SCREENSHOTS") or str(tmp_path / "shots")
    os.makedirs(shots, exist_ok=True)
    server = Served(azr, tmp_path / "server.log")
    try:
        t0 = time.monotonic()
        dom = chrome(tmp_path, f"{server.url}#selftest?edit=1", "--window-size=1400,900",
                     "--virtual-time-budget=120000", "--dump-dom")
        result = selftest_result(dom)
        print(f"{label}: self-test in {time.monotonic() - t0:.1f} s (Chrome included); "
              f"{len(result['steps'])} views, chi2 {result['chi2']}")
        assert result["errors"] == [], result["errors"]
        assert result["failed"] == [], result["failed"]
        drawn = [s["plot"] for s in result["steps"] if s["plot"]]
        if result["chi2"] is not None:
            assert "plotly" in drawn and "table" in drawn, result["steps"]
            assert result["plotly"], "plotly.js did not load"
        for step in result["steps"]:
            assert not step["errorBoxes"], step
        edit = result["edit"]  # a width typed in and Enter, then the Undo button: the engine's chi2 moves and returns
        if "skipped" not in edit:
            assert edit["edited"] and edit["undone"], edit
            assert edit["chi2_after"] != edit["chi2_before"] and edit["chi2_back"] == edit["chi2_before"], edit
            assert edit["rev_after"] == edit["rev_before"] + 1 and edit["rev_back"] == edit["rev_after"] + 1, edit
        take_screenshots(server, tmp_path, label, shots, result)
    finally:
        server.stop()


def take_screenshots(server, tmp_path, label, shots, result):
    """One screenshot per view (the routes of VIEWS filled from the self-test's own routes)."""
    routes = [s["route"] for s in result["steps"]]
    fill = {"level": next((r.split("/", 3)[3] for r in routes if r.startswith("#/nucleus/level/")), None),
            "dataset": next((r.split("/", 2)[2] for r in routes if r.startswith("#/data/") and "?" not in r
                             and r != "#/data/+new"), None)}
    sets = [("light", 1400, 900)]
    if os.environ.get("AZURE2_UI_SHOTS") == "all":
        sets = [(t, w, h) for t in ("light", "dark") for (w, h) in ((1400, 900), (1000, 700))]
    for theme, width, height in sets:
        for view in VIEWS:
            keys = re.findall(r"\{(\w+)\}", view)
            if any(fill[k] is None for k in keys):
                continue
            route = view.format(**fill)
            name = re.sub(r"[^A-Za-z0-9]+", "-", view).strip("-")
            path = os.path.join(shots, f"{label}_{theme}_{width}x{height}_{name}.png")
            chrome(tmp_path, f"{server.url}#/{route}?theme={theme}", f"--window-size={width},{height}",
                   "--virtual-time-budget=30000", f"--screenshot={path}", mode="old")
            assert os.path.getsize(path) > 10000, f"{path}: no page drawn"


if __name__ == "__main__":
    if not CHROME:
        print("skip: no headless Chrome or Chromium found (set CHROME=/path/to/chrome)")
        sys.exit(77)
    uitest.run(__file__)
