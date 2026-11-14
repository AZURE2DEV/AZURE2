#!/usr/bin/env python3
"""Small defects of pyazr, each with the check that failed before its repair (docs/upgrade/HANDOFF.md, the open
items of 2026-11-14):

  (a) Session.save_fit to another directory writes a project that opens there and gives the fit's chi-squared (its
      data paths rewritten relative to the new directory); (a') a project whose data covariance file names its data
      files relative to the project's directory (T2) is refused there before anything is written, with the three
      parts, and written in the project's own directory;
  (b) `python -m pyazr --help` runs without the engine module (exit 0, the usage), and `import pyazr` without it
      still raises the ImportError that says how to build it;
  (c) `pyazr ui FILE --json` prints its URL document and keeps serving (it exited with status 2 at once: uvicorn's
      log set-up asked the captured standard output whether it is a terminal); an interrupt ends it with status 0;
  (d) ChannelWidth.jpi, Parameter.jpi and ScatteringLength.jpi are the engine's J^pi text ("1/2+", not "0.5+");
  (e) a pair record's is_entrance: on T1 the proton pair is an entrance, the photon pairs are not.
The message of Project.set_pair(1, radius=-1) (the value and what to give) is a case of errors_test.py.

Run from anywhere:  python3 tests/pyazr/small_defects_test.py
"""
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import pyazr  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

T1 = os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")
tmp = tempfile.mkdtemp(prefix="small_defects_")
try:
    d = os.path.join(tmp, "T1")
    shutil.copytree(os.path.dirname(T1), d, ignore=shutil.ignore_patterns("output*", "checks*"))
    t1 = os.path.join(d, "T1.azr")

    print("(a) save_fit to another directory")
    other = os.path.join(tmp, "elsewhere", "deeper")
    os.makedirs(other)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(t1) as s:
            x = s.parameters.x0
            chi2 = s.evaluate(x).chi2
            path, _ = s.save_fit(os.path.join(other, "fit.azr"), x)
        try:
            with pyazr.Session(path) as s2:
                again = s2.evaluate(s2.parameters.x0).chi2
            check(f"the snapshot opens in its own directory with the fit's chi2 ({again!r})",
                  abs(again - chi2) <= 1e-9 * chi2, f"{again} vs {chi2}")
        except Exception as e:  # noqa: BLE001
            check("the snapshot opens in its own directory", False, f"{type(e).__name__}: {str(e)[:300]}")
    text = open(path).read()
    check("its data files are named from there (../../T1/data/...)", "../../T1/data/" in text)

    print("(a') save_fit to another directory of a project whose covariance file names data files (T2)")
    T2 = os.path.join(ROOT, "tests", "upgrade", "T2", "T2.azr")
    d2 = os.path.join(tmp, "T2")
    shutil.copytree(os.path.dirname(T2), d2, ignore=shutil.ignore_patterns("output*", "checks*"))
    t2 = os.path.join(d2, "T2.azr")
    snap = os.path.join(tmp, "elsewhere", "t2", "snap.azr")
    os.makedirs(os.path.dirname(snap))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(t2) as s:
            try:
                s.save_fit(snap, param_sav=True)
                check("refused before anything is written", False, "no exception")
            except pyazr.ParameterError as e:
                parts = [e.what, e.allowed, e.fix]
                check("refused with its three parts (what, allowed, fix), each in the message",
                      all(parts) and all(q in str(e) for q in parts), f"{parts}")
                check("the refusal names the covariance file's data files and the project's directory",
                      "data/" in e.what and d2 in e.fix, str(e)[:300])
            check("nothing is left in the other directory", os.listdir(os.path.dirname(snap)) == [],
                  str(os.listdir(os.path.dirname(snap))))
            chi2 = s.evaluate().chi2
            kept, _ = s.save_fit(os.path.join(d2, "snap.azr"))
        with pyazr.Session(kept) as s2:
            again = s2.evaluate().chi2
        check(f"in the project's own directory it is written and reopens with the fit's chi2 ({again!r})",
              abs(again - chi2) <= 1e-9 * abs(chi2), f"{again} vs {chi2}")

    print("(b) the command line's help without the engine")
    env = {k: v for k, v in os.environ.items() if k != "PYAZR_ENGINE_DIR"}
    env["PYTHONPATH"] = ROOT
    engine_in_tree = any(f.startswith("_azure2") for f in os.listdir(os.path.join(ROOT, "pyazr")))
    if engine_in_tree:
        print("  skip (b): the engine module is in pyazr/ itself")
    else:
        r = subprocess.run([sys.executable, "-m", "pyazr", "--help"], env=env, capture_output=True, text=True,
                           cwd=tmp, timeout=60)
        check("python -m pyazr --help: exit 0 and the usage", r.returncode == 0 and "usage:" in r.stdout,
              f"{r.returncode}: {r.stderr[-300:]}")
        r = subprocess.run([sys.executable, "-c", "import pyazr"], env=env, capture_output=True, text=True, cwd=tmp,
                           timeout=60)
        check("import pyazr without the engine: the ImportError that says how to build it",
              r.returncode != 0 and "ImportError" in r.stderr and "not built" in r.stderr, r.stderr[-300:])

    print("(c) pyazr ui FILE --json")
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401
        web = True
    except ImportError:
        print("  skip (c): the extra 'web' (fastapi, uvicorn) is not installed")
        web = False
    if web:
        p = subprocess.Popen([sys.executable, "-m", "pyazr", "ui", t1, "--json", "--threads", "1"], cwd=tmp,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                             env=dict(os.environ, PYTHONPATH=ROOT))
        try:
            doc = json.loads(p.stdout.readline())
            check("the first line is the URL document", doc.get("ok") and doc["result"]["url"].startswith("http://"),
                  str(doc))
            time.sleep(3)
            check("the server still runs 3 s after the URL (it exited with status 2 at once)", p.poll() is None,
                  f"exit {p.poll()}: {p.stderr.read()[-300:] if p.poll() is not None else ''}")
            import urllib.request
            url = doc["result"]["url"]
            token = url.split("token=")[1]
            req = urllib.request.Request(url.split("/?")[0] + "/api/describe",
                                         headers={"Authorization": f"Bearer {token}"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                check("it answers /api/describe", resp.status == 200)
        finally:
            p.send_signal(signal.SIGINT)
            try:
                code = p.wait(timeout=30)
            except subprocess.TimeoutExpired:
                p.kill()
                code = p.wait()
        check("an interrupt ends it with status 0", code == 0, f"exit {code}")

    print("(d), (e) on a session of T1")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(t1) as s:
            rows = s.dimensionless_widths()
            jpis = {r.jpi for r in rows}
            check(f"ChannelWidth.jpi as the engine writes J^pi ({sorted(jpis)})",
                  jpis and all("/" in j or "." not in j for j in jpis) and "1/2+" in jpis, str(jpis))
            pjpis = {q.jpi for q in s.parameters if q.jpi is not None}
            check(f"Parameter.jpi as the engine writes J^pi, the text of ChannelWidth.jpi ({sorted(pjpis)})",
                  pjpis and all("." not in j for j in pjpis) and "1/2+" in pjpis and jpis <= pjpis, str(pjpis))
            sl = {v.jpi for v in s.scattering_lengths()}
            check(f"ScatteringLength.jpi as the engine writes J^pi ({sorted(sl)})",
                  sl and all("." not in j for j in sl), str(sl))
            entrance = {r["key"]: r["is_entrance"] for r in s.pairs.records()} if hasattr(s.pairs, "records") else \
                {q.key: q.is_entrance for q in s.pairs}
            check("is_entrance: the proton pair yes, the photon pairs no", entrance == {1: True, 2: False, 3: False},
                  str(entrance))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
