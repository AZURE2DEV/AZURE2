#!/usr/bin/env python3
"""EXFOR Web-API requests carry the user-agent the IAEA Nuclear Data Section gave AZURE2.

The API sits behind a firewall that refuses many clients (Cloudflare 403 or a challenge page, since
2026-10-22).  On 2026-11-09 IAEA asked AZURE2 to send ``Exfor-Azure2/1.0``.  The string is for EXFOR:
LiveChart, the NNDC ENSDF API and Crossref keep pyazr's own.

Offline: urlopen is replaced by a recorder, so this checks what would be sent, not what the server
answers.  Checks: x4list (search_exfor) and x4get (the web route of fetch_exfor) send the IAEA string;
the other three services do not; the GUI client sends the same string on both of its requests.

Run from anywhere:  python3 tests/pyazr/nds_exfor_user_agent_test.py
"""
import importlib.util
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
# load pyazr/nds.py by path: it needs numpy only, not the compiled engine
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone", os.path.join(ROOT, "pyazr", "nds.py"))
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

# The GUI client: Qt headers are not available here, so read its source.  Every request must take the
# user-agent from the one constant, and the constant must be the IAEA string.
gui = open(os.path.join(ROOT, "gui", "src", "ExforData.cpp")).read()
constant = re.search(r'kExforUserAgent\s*=\s*"([^"]*)"', gui)
if not constant or constant.group(1) != EXPECTED:
    failures.append(f"GUI: kExforUserAgent is {constant.group(1) if constant else None!r}")
setters = re.findall(r"UserAgentHeader\s*,\s*([^)]*)\)", gui)
if len(setters) != 2 or any(s.strip() != "kExforUserAgent" for s in setters):
    failures.append(f"GUI: request user-agents {setters}")
else:
    print(f"ok  GUI: {len(setters)} requests use kExforUserAgent")

for f in failures:
    print("FAIL", f)
sys.exit(1 if failures else 0)
