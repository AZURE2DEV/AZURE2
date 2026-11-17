#!/usr/bin/env python3
"""EXFOR Web-API requests carry the user-agent the IAEA Nuclear Data Section gave AZURE2.

The API sits behind a firewall that refuses many clients (Cloudflare 403 or a challenge page, since
2026-10-22).  On 2026-11-09 IAEA asked AZURE2 to send ``Exfor-Azure2/1.0``.  The string is for EXFOR:
LiveChart, the NNDC ENSDF API and Crossref keep pyazr's own.

Offline: urlopen is replaced by a recorder, so this checks what would be sent, not what the server
answers.  Checks: x4list (search_exfor) and x4get (the web route of fetch_exfor) send the IAEA string;
the other three services do not.  (From dev, 9b6ecb17; dev's check of the GUI client is not here: this
build has no graphical utility.)

Run from anywhere:  python3 tests/pyazr/nds_exfor_user_agent_test.py
"""
import importlib.util
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
# load pyazr/sources/nds.py by path: it needs numpy only, not the compiled engine
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone",
                                              os.path.join(ROOT, "pyazr", "sources", "nds.py"))
N = importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)

EXPECTED = "Exfor-Azure2/1.0"
sent = []


class Reply(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def fake_urlopen(request, timeout=None):
    sent.append((request.full_url, request.get_header("User-agent")))
    return Reply(b'{"x4Datasets": []}')


N.urllib.request.urlopen = fake_urlopen
failures = []


def check(label, url_part, want_iaea):
    hits = [ua for url, ua in sent if url_part in url]
    if not hits:
        failures.append(f"{label}: no request to {url_part}")
    elif (all(ua == EXPECTED for ua in hits)) != want_iaea:
        failures.append(f"{label}: user-agent {hits[0]!r}")
    else:
        print(f"ok  {label}: {hits[0]}")


N.search_exfor(target="C-13", reaction="a,n", quantity="SIG")
check("x4list", "/exfor/x4list", True)

try:
    N._fetch_exfor_web("C3062002")
except ValueError:
    pass  # the canned reply is not a CSV; only the request matters here
check("x4get", "/exfor/x4get", True)

N._http_get(N.LIVECHART_URL, {"fields": "levels", "nuclides": "14n"})
check("LiveChart", N.LIVECHART_URL, False)
N._http_get(N.ENSDF_API_URL + "/levels", {"nuclide": "14N"})
check("NNDC ENSDF API", N.ENSDF_API_URL, False)
N._http_get(N.CROSSREF_URL + "/10.1103/x")
check("Crossref", N.CROSSREF_URL, False)

for f in failures:
    print("FAIL", f)
sys.exit(1 if failures else 0)
