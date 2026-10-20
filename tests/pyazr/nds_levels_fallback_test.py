#!/usr/bin/env python3
"""fetch_levels falls back from IAEA LiveChart to the NNDC ENSDF API, and both give the same records.

From some networks every IAEA request is answered by a Cloudflare challenge (HTTP 403, or a 200 HTML
page); on the Notre Dame CRC cluster on 2026-10-15 LiveChart answered some requests and refused the
next ones.  ``fetch_levels`` now tries LiveChart and, when it cannot be reached or does not answer with
its CSV, reads the same ENSDF adopted levels from https://www.nndc.bnl.gov/ensdf-api/ (JSON).

Offline: the web is replaced by canned answers -- a LiveChart CSV and an ENSDF-API page for three 14N
levels, saved 2026-10-15 (tests/pyazr/data/ensdf_api_14N_levels.json; the LiveChart rows carry the energy,
uncertainty, J^pi and half-life text LiveChart returned for the same levels that day).  half_life_sec is
compared to 6 significant digits: for a level known by its width the two services differ at 1e-7.  Checks: the two sources give identical Level
records (LiveChart's J^pi list and half-life spellings included); a 403, a URLError and an HTML page
each trigger the fallback with a warning; source="livechart" re-raises; source="nndc" never asks
LiveChart; the API's pages are concatenated; a bad source is refused.

Run from anywhere:  python3 tests/pyazr/nds_levels_fallback_test.py
"""
import importlib.util
import io
import json
import os
import sys
import urllib.error
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
# load pyazr/nds.py by path: it needs numpy only, not the compiled engine
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone", os.path.join(ROOT, "pyazr", "nds.py"))
N = importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)

API_PAGE = json.load(open(os.path.join(HERE, "data", "ensdf_api_14N_levels.json")))
LIVECHART_CSV = (
    "z,n,symbol,idx,energy,unc_e,jp,half_life,half_life_sec\n"
    "7,7,N,0,0.0,0,1+,STABLE,\n"
    "7,7,N,1,2312.798,0.011,0+,68,6.8e-14\n"
    "7,7,N,24,10101.0,15,2+ 1+,12,3.80198e-20\n")
CLOUDFLARE_HTML = "<!DOCTYPE html><html><head><title>Just a moment...</title></head><body></body></html>"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


class Web:
    """Stands in for nds._http_get: LiveChart answers per `livechart`, the ENSDF API from pages."""

    def __init__(self, livechart, pages=None):
        self.livechart, self.pages, self.calls = livechart, pages or [API_PAGE], []

    def __call__(self, url, params=None, timeout=60):
        self.calls.append(url)
        if url == N.LIVECHART_URL:
            if isinstance(self.livechart, Exception):
                raise self.livechart
            return self.livechart
        if url == f"{N.ENSDF_API_URL}/levels":
            assert params["document-type"] == "adopted" and params["nuclide"]
            return json.dumps(self.pages[int(params["offset"]) // int(params["limit"])])
        raise AssertionError(f"unexpected URL {url}")


def fetch(web, source="auto", nuclide="14n"):
    N._http_get = web
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        levels = N.fetch_levels(nuclide, source=source)
    return levels, [str(x.message) for x in w]


def key(lv):
    return (lv.z, lv.n, lv.symbol, lv.index, lv.energy, lv.energy_unc, lv.jp, lv.half_life,
            float(f"{lv.half_life_sec:.6g}"))


def main():
    real = N._http_get
    try:
        lc, warn = fetch(Web(LIVECHART_CSV))
        check("LiveChart answering: used, no warning", [l.source for l in lc] == ["livechart"] * 3 and not warn,
              f"{[l.source for l in lc]} {warn}")

        web = Web(LIVECHART_CSV)
        nn, warn = fetch(web, source="nndc")
        check("source='nndc': LiveChart never asked", N.LIVECHART_URL not in web.calls and
              [l.source for l in nn] == ["nndc"] * 3 and not warn, f"{web.calls} {warn}")
        check("same records from both services", [key(a) for a in lc] == [key(b) for b in nn],
              f"\n      {[key(a) for a in lc]}\n      {[key(b) for b in nn]}")
        check("NNDC J^pi list in LiveChart's spelling ('2+,1+' -> '2+ 1+')", nn[2].jp == "2+ 1+", nn[2].jp)
        check("NNDC half_life as LiveChart's bare number ('68', width '12', 'STABLE')",
              [l.half_life for l in nn] == ["STABLE", "68", "12"], [l.half_life for l in nn])
        check("width_kev from the ENSDF total width (12 keV), None otherwise",
              [l.width_kev for l in nn] == [None, None, 12.0], [l.width_kev for l in nn])
        check("spin/parity parse as for LiveChart", (nn[1].spin, nn[1].parity) == (0, 1))

        for what, answer in (("HTTP 403 (Cloudflare challenge)",
                              urllib.error.HTTPError(N.LIVECHART_URL, 403, "Forbidden", {}, io.BytesIO(b""))),
                             ("URLError (no route)", urllib.error.URLError("unreachable")),
                             ("HTML page instead of the CSV", CLOUDFLARE_HTML)):
            web = Web(answer)
            got, warn = fetch(web)
            check(f"{what}: falls back to the ENSDF API with a warning",
                  [key(l) for l in got] == [key(l) for l in nn] and len(warn) == 1 and "NNDC" in warn[0]
                  and web.calls[0] == N.LIVECHART_URL, f"{warn} {web.calls}")

        try:
            fetch(Web(urllib.error.HTTPError(N.LIVECHART_URL, 403, "Forbidden", {}, io.BytesIO(b""))),
                  source="livechart")
            check("source='livechart' re-raises", False, "no exception")
        except urllib.error.HTTPError:
            check("source='livechart' re-raises", True)

        r = API_PAGE["results"]          # two pages of two, the first one out of energy order
        page1 = {"meta": {"total": 3, "limit": 2, "offset": 0, "next": "..."}, "results": [r[2], r[1]]}
        page2 = {"meta": {"total": 3, "limit": 2, "offset": 2, "previous": "..."}, "results": [r[0]]}
        web = Web(LIVECHART_CSV, pages=[page1, page2])
        N._http_get = web
        paged = N._levels_nndc("14n", page=2)
        check("API pages concatenated and sorted by energy", [key(l) for l in paged] == [key(l) for l in nn]
              and len(web.calls) == 2, f"{len(web.calls)} calls")

        for bad in (("14n", "iaea"), ("N-14", "auto")):
            try:
                N.fetch_levels(bad[0], source=bad[1])
                check(f"refused: fetch_levels{bad}", False, "no exception")
            except ValueError:
                check(f"refused: fetch_levels{bad}", True)
    finally:
        N._http_get = real

    print("PASS" if not failures else f"FAILED: {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
