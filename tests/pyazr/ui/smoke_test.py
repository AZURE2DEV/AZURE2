"""The browser smoke test of the page (docs/upgrade/web_ui.md section 9), with the headless Chrome of the machine: a
real server (``python -m pyazr.ui``) on copies of T1 and T2, the page's self-test (``#selftest``: every tab, a level,
a data set's form and its plot drawn by plotly.js and as a table, a record kind without layout hints, a new record's
form, a width typed into its box and taken back with the Undo button; in the Fit tab a recipe fit started with Run, its
progress seen, Stop, Accept into project and Undo) read from the document Chrome dumps, and one screenshot per tab.  Fails on any script error, console error or
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
from urllib.parse import quote

import pytest

import uitest
from uitest import TOKEN, copy_project, gone

CHROME = os.environ.get("CHROME") or next((shutil.which(n) for n in ("google-chrome", "chromium", "chromium-browser")
                                           if shutil.which(n)), None)
# Chrome's own work in a fresh profile kept out of the test: extensions, component updates, background requests.
QUIET_PROFILE = ("--disable-extensions", "--disable-component-update", "--disable-background-networking")
VIEWS = ("project", "nucleus", "nucleus/level/{label}", "nucleus/{record}", "data", "data/{dataset}", "fit", "results",
         "physics", ("physics-masked", "physics/coulomb_waves?from=-2&to=20&points=111&spacing=linear&v.l=0 1 2 3"
                     "&fields=F,G"), ("physics-stopping", "physics/stopping_power?logy=1"), "thm",
         ("thm-vertex", "thm?fn=thm_vertex"),
         # the first view's own tabs (a fresh profile opens in the first view), and two views in the expert's
         ("first-pairs", "nucleus/pair"), ("first-segments", "data?section=data"),
         ("first-effects", "data?section=targets"), ("first-calculated", "results?calculated=1"),
         ("expert-project", "project?depth=expert"), ("expert-dataset", "data/{dataset}?depth=expert"),
         # the two utilities: the rate tab (the expert's depth has its link) and a data set's covariance builder
         ("expert-rate", "rate?depth=expert"), ("covariance", "covariance/{dataset}"))


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
    # Each run has a fresh profile, in which Chrome installs and loads its component extensions a few seconds after
    # start, while the self-test runs; some of the page's requests not yet sent were then cancelled in the renderer
    # (net::ERR_ABORTED, "Failed to fetch"; neither Chrome's network log nor the server saw them): 1 in 10 runs on
    # T2, 4 in 10 before the Results tab; 0 in 25 with --disable-extensions, 0 in 25 without component updates and
    # background requests, 6 in 25 with Safe Browsing off (tests/pyazr/ui/README.md).  A browser in use has loaded
    # its extensions before the page opens, so the test's Chrome loads none.
    cmd = [CHROME, f"--headless={mode}", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
           "--hide-scrollbars", *QUIET_PROFILE, f"--user-data-dir={profile}", *args, url]
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
        dom = chrome(tmp_path, f"{server.url}#selftest?edit=1&fit=1", "--window-size=1400,900",
                     "--virtual-time-budget=120000", "--dump-dom")
        result = selftest_result(dom)
        print(f"{label}: self-test in {time.monotonic() - t0:.1f} s (Chrome included); "
              f"{len(result['steps'])} views, chi2 {result['chi2']}, "
              f"evaluation refused: {result['evaluation_refused']}")
        assert result["errors"] == [], result["errors"]
        assert result["failed"] == [], result["failed"]
        if label in ("T1", "T2"):
            assert result["chi2"] is not None, result["evaluation_refused"]
        drawn = [s["plot"] for s in result["steps"] if s["plot"]]
        if result["chi2"] is not None:
            assert "plotly" in drawn and "table" in drawn, result["steps"]
            assert result["plotly"], "plotly.js did not load"
        for step in result["steps"]:
            # a THM diagnostic the experiment does not have (no line shape, no distortion) is the engine's refusal
            refused = step["route"].startswith("#/thm?fn=") and all("(422)" in b for b in step["errorBoxes"])
            assert not step["errorBoxes"] or refused, step
        # the THM tab: listed exactly where the project has THM content; its diagnostics drew their answers
        # the first view: each tab of the old program's path marked as it opens, a data set's form with its rarely
        # set fields folded, Calculate on the Plot tab
        first = [s for s in result["steps"] if s["depth"] == "first"]
        titles = [s["tab"] for s in first if not s["route"].startswith("#/data/") and s["route"] != "Calculate"]
        assert titles == ["Project", "Particle Pairs", "Levels and Channels", "Segments", "Experimental Effects", "Fit",
                          "Plot"], titles
        assert all(s["depth"] == "expert" for s in result["steps"] if s not in first), "a view not in its depth"
        if result["chi2"] is not None:
            opened = [s for s in first if s["route"].startswith("#/data/")]
            assert opened and opened[0]["forms"] == 1 and opened[0]["folded"] == 1 and opened[0]["tab"] == "Segments", \
                opened
            assert [s["tab"] for s in first if s["route"] == "Calculate"] == ["Plot", "Plot"], first
        thm = [s for s in result["steps"] if s["route"].startswith("#/thm")]
        listed = thm[0]["tabListed"]
        assert listed == (label != "T1") or label not in ("T1", "T2"), thm[0]
        if listed:
            assert len(thm) > 1 and thm[0]["tables"] >= 1, thm
            vertex = next(s for s in thm if "fn=thm_vertex" in s["route"])
            assert vertex["physics"] and vertex["physics"]["computed"] > 0 and vertex["plot"], vertex
        # the Physics tab: every view drew its plot from the route's answer with the status line and the pyazr call;
        # the Coulomb function asked on rho from -2 has elements outside its domain, counted and shown as gaps
        physics = {s["route"]: s["physics"] for s in result["steps"] if s["route"].startswith("#/physics")}
        assert len(physics) == 5 and all(v and v["computed"] > 0 and v["python"].startswith("import numpy")
                                         for v in physics.values()), physics
        masked = next(v for r, v in physics.items() if "from=-2" in r)
        assert masked["elements"] == 444 and masked["computed"] == 4 * 100, masked  # rho = -2, -1.8, ..., 0 not
        assert "rho outside its domain" in masked["text"], masked
        # a stopping power carries SRIM's citation on the page and in its Python text (the PI's condition)
        stopping = next(v for r, v in physics.items() if "stopping_power" in r)
        assert "SRIM" in stopping["citation"] and "# cite: " + stopping["citation"] in stopping["python"], stopping
        edit = result["edit"]  # a width typed in and Enter, then the Undo button: the engine's chi2 moves and returns
        if "skipped" not in edit:
            assert edit["edited"] and edit["undone"], edit
            assert edit["chi2_after"] != edit["chi2_before"] and edit["chi2_back"] == edit["chi2_before"], edit
            assert edit["rev_after"] == edit["rev_before"] + 1 and edit["rev_back"] == edit["rev_after"] + 1, edit
            # Calculate after the edit: compared with the first view's Calculate on the text before it
            assert edit["calculated"] and edit["calc_chi2"] == [edit["chi2_before"], edit["chi2_after"]], edit
            assert edit["calc_total"].startswith("Δχ² ") and edit["calc_first"] and edit["calc_tab"] == "Plot", edit
        pair = result["pair"]  # a pair added with the Nucleus tab's form, removed with its button, both undone
        print(f"{label}: pair step {pair}")
        if "skipped" not in pair:
            assert pair["added"] and pair["removed"] and pair["undone"], pair
            assert (pair["after_add"], pair["after_remove"], pair["after_undo"]) == (
                pair["before"] + 1, pair["before"], pair["before"]), pair
            assert pair["chi2_added"] == pair["chi2_before"] == pair["chi2_back"], pair  # no level names the pair
        fit = result["fit"]  # Run, progress, Stop (or the end), Accept, Undo: the header's chi2 moves and returns
        if "skipped" not in fit:
            assert fit["progress_seen"] and fit["plot"] and fit["ended"] and fit["accepted"] and fit["undone"], fit
            assert fit["state"] in ("stopped", "converged", "ended"), fit
            assert fit["chi2_header_after"] == fit["chi2_written"] and fit["chi2_back"] == fit["chi2_before"], fit
            assert fit["rev_accepted"] == fit["rev_before"] + 1 and fit["rev_back"] == fit["rev_accepted"] + 1, fit
            print(f"{label}: fit {fit['method']} {fit['state']} after {fit['calls_seen']} calls seen; chi2 "
                  f"{fit['chi2_before']} -> fit {fit['chi2_fit']}, as written {fit['chi2_written']} -> undo "
                  f"{fit['chi2_back']}")
        take_screenshots(server, tmp_path, label, shots, result)
    finally:
        server.stop()


def take_screenshots(server, tmp_path, label, shots, result):
    """One screenshot per view (the routes of VIEWS filled from the self-test's own routes)."""
    routes = [s["route"] for s in result["steps"]]
    fill = {"label": next((r.split("/", 3)[3] for r in routes if r.startswith("#/nucleus/level/")), None),
            "dataset": next((r.split("/", 2)[2] for r in routes if r.startswith("#/data/") and "?" not in r
                             and r != "#/data/+new"), None),
            "record": next((r.split("/", 2)[2] for r in routes if r.startswith("#/nucleus/")
                            and not r.startswith("#/nucleus/level/")), None)}
    sets = [("light", 1400, 900)]
    if os.environ.get("AZURE2_UI_SHOTS") == "all":
        sets = [(t, w, h) for t in ("light", "dark") for (w, h) in ((1400, 900), (1000, 700))]
    for theme, width, height in sets:
        for view in VIEWS:
            name, view = view if isinstance(view, tuple) else (None, view)
            keys = re.findall(r"\{(\w+)\}", view)
            if any(fill[k] is None for k in keys):
                continue
            route = view.format(**fill)
            name = name or re.sub(r"[^A-Za-z0-9]+", "-", view.replace("{", "").replace("}", "")).strip("-")
            path = os.path.join(shots, f"{label}_{theme}_{width}x{height}_{name}.png")
            url = f"{server.url}#/{quote(route, safe='/?=&,.-_')}{'&' if '?' in route else '?'}theme={theme}"
            chrome(tmp_path, url, f"--window-size={width},{height}",
                   "--virtual-time-budget=30000", f"--screenshot={path}", mode="old")
            assert os.path.getsize(path) > 10000, f"{path}: no page drawn"


def opened_tab(tmp_path, url, done, timeout=60):
    """Open ``url`` in a headless Chrome that stays up, and poll its tab list (the debugging port's /json) until
    ``done(tab)`` holds for the page's tab or ``timeout`` seconds pass; the tab's last (url, title).  For what
    --dump-dom cannot show: where the browser ends after a click and the navigations that follow it."""
    import httpx
    from pyazr.ui.server import free_port
    profile = tmp_path / f"chrome-{time.monotonic_ns()}"
    profile.mkdir()
    port = free_port()
    proc = subprocess.Popen([CHROME, "--headless=new", "--disable-gpu", "--no-first-run",
                             "--no-default-browser-check", *QUIET_PROFILE, f"--user-data-dir={profile}",
                             f"--remote-debugging-port={port}", url], stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    last = {}

    def look():
        try:
            tabs = httpx.get(f"http://127.0.0.1:{port}/json", timeout=5).json()
        except (httpx.HTTPError, ValueError):
            return False
        pages = [t for t in tabs if t.get("type") == "page"]
        last.update(pages[0] if pages else {})
        return bool(pages) and done(pages[0])

    try:
        uitest.wait_until(look, timeout, step=0.2)
    finally:
        proc.terminate()
        try:
            proc.wait(20)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(10)
    return last.get("url"), last.get("title")


def test_printed_url_followed_from_another_site_opens_the_page(tmp_path):
    """The printed URL as a link on a page of another site, clicked there (how the URL reaches a user: a chat, a
    mail, a notebook in the browser): the project's page opens.  With a SameSite=Strict cookie the browser kept the
    cookie back on the redirect that follows such a click and the server answered 401 to its own link.  The link's
    page is served on ``localhost`` and the project on ``127.0.0.1``: two sites for the browser.  Then the address
    without the token in a browser that has no cookie: the page that asks for the token, not JSON."""
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    server = Served(copy_project("T1", tmp_path), tmp_path / "server.log")
    link = (f'<!doctype html><title>elsewhere</title><a id="a" href="{server.url}">the printed URL</a>'
            '<script>setTimeout(() => document.getElementById("a").click(), 200)</script>').encode()

    class Elsewhere(BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802  -- the base class's name
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(link)))
            self.end_headers()
            self.wfile.write(link)

        def log_message(self, *args):
            pass

    other = ThreadingHTTPServer(("127.0.0.1", 0), Elsewhere)
    threading.Thread(target=other.serve_forever, daemon=True).start()
    try:
        page = f"http://127.0.0.1:{server.port}/"
        url, title = opened_tab(tmp_path, f"http://localhost:{other.server_address[1]}/",
                                lambda tab: tab["url"] == page and tab["title"].endswith("— AZURE2"))
        assert (url, title) == (page, "T1.azr — AZURE2"), f"the click ended at {url!r} with the title {title!r}"
        url, title = opened_tab(tmp_path, page, lambda tab: "token" in tab["title"])
        assert (url, title) == (page, "AZURE2 — token needed"), f"without the token: {url!r}, title {title!r}"
    finally:
        other.shutdown()
        other.server_close()
        server.stop()


if __name__ == "__main__":
    if not CHROME:
        print("skip: no headless Chrome or Chromium found (set CHROME=/path/to/chrome)")
        sys.exit(77)
    uitest.run(__file__)
