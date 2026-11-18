"""The page's navigation walked as a person walks it (docs/upgrade/web_navigation.md), in a headless Chrome driven
through its DevTools protocol: the printed URL opens the first tab of the work with the overview of a first visit;
the main row holds the old setup program's tabs in its order and words, and each draws its table; a row chosen in
Levels and Channels puts its label in the address and its channels beside the list; Reload keeps the tab and the row,
Back returns to the place before; Calculate shows the Plot tab and the status line the engine's χ²; the validator's
mark opens its list; a fit started on Calculate and Fit is shown, on another tab, in the status line with its Stop,
which stops it.  Each step is a click on what a person sees (a link, a button), never a route typed in.

Written against the navigation of 2026-11-19; on the page before it (a Project tab first, no Calculate and Fit tab,
no status-line marks) it fails at its first assertion on the tab row.

Environment as ``smoke_test.py`` (``CHROME``; skipped without it); needs the package ``websockets`` (skipped without
it).  Label ``browser``.
"""
import asyncio
import base64
import json
import os
import subprocess
import sys
import time

try:
    import pytest
    import websockets
except ImportError as e:  # ctest counts exit 77 as skipped
    print(f"skip: {e.name} is not installed")
    sys.exit(77)

import uitest
from smoke_test import CHROME, QUIET_PROFILE, Served
from uitest import copy_project

MAIN_ROW = ["Particle Pairs", "Levels and Channels", "Segments", "Experimental Effects", "Calculate and Fit", "Plot"]
SHOTS = os.environ.get("AZURE2_UI_SCREENSHOTS")


class Page:
    """One tab of a headless Chrome on the DevTools protocol: evaluate, click, wait, reload, screenshot."""

    def __init__(self, ws):
        self.ws, self.n = ws, 0

    async def send(self, method, **params):
        self.n += 1
        mine = self.n
        await self.ws.send(json.dumps({"id": mine, "method": method, "params": params}))
        while True:
            msg = json.loads(await self.ws.recv())
            if msg.get("id") == mine:
                assert "error" not in msg, msg
                return msg.get("result")

    async def js(self, expr):
        r = await self.send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
        assert not r.get("exceptionDetails"), r
        return r["result"].get("value")

    async def until(self, expr, timeout=60, what=""):
        t0 = time.monotonic()
        while time.monotonic() - t0 < timeout:
            if await self.js(f"(() => {{ try {{ return Boolean({expr}); }} catch {{ return false; }} }})()"):
                return
            await asyncio.sleep(0.2)
        raise AssertionError(f"not within {timeout} s: {what or expr}; at {await self.js('location.hash')}, status "
                             f"{await self.js('document.querySelector(\".status\")?.textContent')!r}")

    async def click(self, selector, text=None):
        """Click the element a person would: the first visible match of ``selector`` (whose text is ``text``)."""
        found = await self.js(f"""(() => {{
          const all = [...document.querySelectorAll({json.dumps(selector)})].filter((e) => e.offsetParent !== null
            && ({json.dumps(text)} === null || e.textContent.trim() === {json.dumps(text)}));
          if (!all.length) return false;
          all[0].click();
          return true; }})()""")
        assert found, f"nothing to click: {selector} {text or ''}"

    async def shot(self, name):
        if SHOTS:
            r = await self.send("Page.captureScreenshot", format="png")
            with open(os.path.join(SHOTS, f"navigation_{name}.png"), "wb") as f:
                f.write(base64.b64decode(r["data"]))


def tab_titles():
    """The main row's titles as shown (the key number and the group mark left out)."""
    return ("[...document.querySelectorAll('nav.tabs:not(.subtabs-row) a.tab')].map((a) => [...a.childNodes]"
            ".filter((n) => n.nodeType === 3).map((n) => n.textContent).join('').trim())")


ACTIVE = ("[...document.querySelector('nav.tabs:not(.subtabs-row) a.tab.active')?.childNodes || []]"
          ".filter((n) => n.nodeType === 3).map((n) => n.textContent).join('').trim()")


async def walk(port, url):
    import httpx
    for _ in range(100):
        try:
            tabs = httpx.get(f"http://127.0.0.1:{port}/json", timeout=5).json()
            target = next(t for t in tabs if t.get("type") == "page")
            break
        except (httpx.HTTPError, ValueError, StopIteration):
            await asyncio.sleep(0.2)
    async with websockets.connect(target["webSocketDebuggerUrl"], max_size=50_000_000) as ws:
        p = Page(ws)
        await p.send("Page.enable")
        await p.send("Page.navigate", url=url)
        # the printed URL: the first tab of the work, its address, the overview of a first visit
        await p.until("document.querySelector('main .tab-panel table')", what="the first tab's table")
        assert await p.js(tab_titles()) == MAIN_ROW + ["Rate and ENDF", "Physics and THM", "Covariances", "Project"]
        assert await p.js("location.hash") == "#/nucleus/pair"
        await p.until("document.getElementById('overview')", 10, "the overview of a first visit")
        await p.shot("overview")
        await p.click("#overview button", "Close")
        assert not await p.js("document.getElementById('overview')")
        # each main tab of the old program's row: one click, its name marked, its table (or its empty list) drawn
        for title in MAIN_ROW:
            await p.js(f"""[...document.querySelectorAll('nav.tabs:not(.subtabs-row) a.tab')].find((a) =>
              a.dataset.main === {json.dumps(title)}).click()""")
            await p.until(f"{ACTIVE} === {json.dumps(title)} && !document.querySelector('p.tab-loading')"
                          " && document.querySelector('main .tab-panel table, main .tab-panel .plot-area, "
                          "main .tab-panel p.empty')", 90, f"the tab {title} and its table")
            await p.shot(title.replace(" ", "_"))
        # a level chosen with a click: its label in the address, its channels beside the list
        await p.js("document.querySelector('a.tab[data-main=\"Levels and Channels\"]').click()")
        await p.until("document.querySelector('[aria-label=\"levels: arrow keys move, Enter opens\"] tr')")
        before = await p.js("location.hash")
        await p.js("document.querySelectorAll('[aria-label=\"levels: arrow keys move, Enter opens\"] tr')[2].click()")
        await p.until("location.hash.startsWith('#/nucleus/level/') && document.querySelector("
                      "'.master-detail .detail [aria-label^=\"channels\"] tr')", what="the level's channels beside")
        chosen = await p.js("location.hash")
        label = await p.js("document.querySelector('.master-detail .detail h2').textContent")
        assert chosen != before and label.startswith("Level "), (chosen, label)
        await p.shot("level_chosen")
        # Reload keeps the tab and the row; Back returns to the place before
        await p.send("Page.reload")
        await p.until(f"location.hash === {json.dumps(chosen)} && {ACTIVE} === 'Levels and Channels' && "
                      "document.querySelector('.master-detail .detail h2')?.textContent === " + json.dumps(label),
                      what="the same tab and level after Reload")
        assert not await p.js("document.getElementById('overview')"), "the overview came back after it was closed"
        await p.js("history.back()")
        await p.until(f"location.hash === {json.dumps(before)}", what="Back to the place before")
        # Calculate: the Plot tab, the engine's χ² in the status line; the validator's mark and its list
        await p.click("button.calculate")
        await p.until(f"{ACTIVE} === 'Plot' && /χ² \\d/.test(document.querySelector('.status').textContent)",
                      90, "Calculate's Plot tab and the χ² in the status line")
        assert (await p.js("document.querySelector('button.validation').textContent")) == "✓ valid"
        await p.click("button[data-action='validate']")
        await p.until("!document.querySelector('.diag-panel').hidden", what="the validator's list")
        await p.shot("calculated")
        await p.click(".diag-panel button", "Close")
        # a fit: started on Calculate and Fit, followed on another tab by the status line, stopped there
        await p.js("document.querySelector('a.tab[data-main=\"Calculate and Fit\"]').click()")
        await p.until("document.querySelector('button[data-action=\"run\"]')")
        await p.click("button[data-action='run']")
        await p.until("location.hash.startsWith('#/fit/job/')", what="the started job's panel")
        await p.js("document.querySelector('a.tab[data-main=\"Segments\"]').click()")
        await p.until(f"{ACTIVE} === 'Segments' && document.querySelector('.job-mark:not([hidden]) a')"
                      "?.textContent.startsWith('job ')", 30, "the fit in the status line on another tab")
        if await p.js("Boolean(document.querySelector('.job-mark.live'))"):
            await p.shot("fit_running")
            await p.click("button[data-action='status-stop']")
        await p.until("document.querySelector('.job-mark.done, .job-mark.failed')", 120, "the fit ended")
        mark = await p.js("document.querySelector('.job-mark a').textContent")
        assert mark.startswith("job ") and ("stopped" in mark or "converged" in mark or "ended" in mark), mark
        await p.click(".job-mark a")
        await p.until(f"{ACTIVE} === 'Calculate and Fit' && location.hash.startsWith('#/fit/job/')",
                      what="the job's panel from the status line")
        return mark


def test_navigation_walk(tmp_path):
    from pyazr.ui.server import free_port
    if not CHROME:
        pytest.skip("no Chrome or Chromium")
    server = Served(copy_project("T1", tmp_path), tmp_path / "server.log")
    port = free_port()
    profile = tmp_path / "chrome"
    profile.mkdir()
    chrome = subprocess.Popen([CHROME, "--headless=new", "--disable-gpu", "--no-first-run",
                               "--no-default-browser-check", *QUIET_PROFILE, "--window-size=1400,900",
                               f"--user-data-dir={profile}", f"--remote-debugging-port={port}", "about:blank"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        mark = asyncio.run(asyncio.wait_for(walk(port, server.url), 400))
        print(f"navigation walked; the fit's mark: {mark}")
    finally:
        chrome.terminate()
        try:
            chrome.wait(20)
        except subprocess.TimeoutExpired:
            chrome.kill()
        server.stop()


if __name__ == "__main__":
    uitest.run(__file__)
