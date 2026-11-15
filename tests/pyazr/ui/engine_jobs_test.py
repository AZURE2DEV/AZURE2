#!/usr/bin/env python3
"""The engine's fitters as jobs of the web UI (docs/upgrade/web_ui.md section 11.6), on copies of T3 and T2 at the
environment's threads (one under ctest):

  1. every fit row of the engine's table runs as a job to its end, at the objective and end point of
     ``Session.fit_with_engine`` called directly on a session of the same text from the same start (the same code:
     equal to the last bit), with its covariance and errors by name where the fitter gives them; the start objective
     is the objective at the start point;
  2. the engine's MINOS as a job that follows an engine fit, from its best point: asymmetric intervals of every free
     parameter; MINOS as a setting of the fit job gives the same kind of intervals;
  3. a call limit ends the fit with the engine's status said ("ended", the status "call_limit" in its words, the
     result's outcome "call_limit"; on T2 with a limit of 40 too, where a fit without a limit is "converged"); a
     fitter that ignores a call limit says so in the result's notices, as the engine's result does; MINOS's
     validity per side is the engine's; the fitter rows carry the table's "call_limit" and "fallback" columns;
  4. Accept writes the job's end point into the project (one history entry; the χ² of the project as written equals
     the job's χ² at its best point to rounding of the written digits) and Undo restores the start χ² exactly;
  5. on T2, a stop request ends a running engine fit (MIGRAD and LM) within one objective call of the request, with
     the best point seen as its best point, which is then accepted.

The stop times are printed (measured on the product).
Run: python3 tests/pyazr/ui/engine_jobs_test.py
"""
import os
import sys
import time

import numpy as np
try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client, wait_job, wait_until  # noqa: E402


@pytest.fixture(scope="module")
def t3(tmp_path_factory):
    path = copy_project("T3", tmp_path_factory.mktemp("engine_t3"))
    with make_client(path) as client:
        client.path = path
        yield client


@pytest.fixture(scope="module")
def t2(tmp_path_factory):
    path = copy_project("T2", tmp_path_factory.mktemp("engine_t2"))
    with make_client(path) as client:
        client.path = path
        yield client


def rev(client):
    return client.get("/api/project").json()["rev"]


def start(client, **body):
    r = client.post("/api/jobs", json=dict({"kind": "fit", "x0": "text"}, **body))
    assert r.status_code == 200, r.text
    return r.json()


def direct(client, method, minos=False, max_calls=None):
    """Session.fit_with_engine called in this process on a session of the text the job read, from the text's x0."""
    import pyazr
    state = client.app.state.server.state
    path = state.session_file("direct")
    with open(path, "wb") as f:
        f.write(state.session_text().encode())
    with pyazr.Session(path, cwd=state.directory) as s:
        r = s.fit_with_engine(method, minos=minos, max_calls=max_calls)
        start = float(s.objective()(s.parameters.x0))
    os.remove(path)
    return r, start


def engine_fit_rows(client):
    return [r for r in client.get("/api/jobs/kinds").json() if r["kind"] == "fit" and r["method"].startswith("engine_")]


def test_every_engine_fitter_ends_where_the_direct_call_ends(t3):
    rows = engine_fit_rows(t3)
    import pyazr
    assert sorted(r["method"] for r in rows) == sorted(f"engine_{n}" for n in pyazr.kinds()["fitters"])
    for row in rows:
        name = row["method"][len("engine_"):]
        done = wait_job(t3, start(t3, method=row["method"])["id"])
        r, start_objective = direct(t3, name)
        result = done["result"]
        print(f"\n  {name}: job {done['best']['objective']!r}, direct {r.objective!r}, {result['calls']} calls "
              f"(fit {result['fit_calls']}), {result['elapsed']:.2f} s, {result['message']}")
        assert done["state"] == ("converged" if r.status == "converged" else "ended"), done["state"]
        assert result["engine_status"] == r.status and r.status in result["message"]
        assert done["best"]["objective"] == r.objective
        assert [done["best"]["x"][n] for n in r.names] == r.x.tolist()
        assert result["start_objective"] == start_objective
        if row["facts"]["covariance"]:
            assert result["covariance"]["names"] == r.names
            assert np.array_equal(np.array(result["covariance"]["matrix"]), r.covariance)
            assert [result["errors"][n] for n in r.names] == r.errors.tolist()
            assert result["errors_from"]
        calls = [p[0] for p in t3.get(f"/api/jobs/{done['id']}").json()["trace"]]
        assert calls == sorted(calls) and len(calls) >= 1


def test_engine_minos_follows_an_engine_fit(t3):
    table = {(r["kind"], r["method"]): r for r in t3.get("/api/jobs/kinds").json()}
    assert table[("minos", "engine_migrad")]["follows"] == ["fit", "engine_migrad"]
    fit = wait_job(t3, start(t3, method="engine_migrad")["id"])
    job = start(t3, kind="minos", method="engine_migrad", x0={"job": fit["id"]})
    done = wait_job(t3, job["id"])
    assert done["state"] == "converged", done
    minos = done["result"]["minos"]
    assert sorted(minos) == sorted(fit["best"]["x"])
    for n, m in minos.items():
        assert m["lower"] < 0 < m["upper"], (n, m)
    assert done["result"]["start_x"] == fit["best"]["x"]
    with_minos = wait_job(t3, start(t3, method="engine_migrad", settings={"minos": "true"})["id"])
    r, _ = direct(t3, "migrad", minos=True)
    assert {n: [m["lower"], m["upper"]] for n, m in with_minos["result"]["minos"].items()} == \
        {n: list(v) for n, v in r.minos.items()}
    assert {n: m["sides"] for n, m in with_minos["result"]["minos"].items()} == r.minos_valid
    assert all(m["valid"] == all(v["valid"] for v in m["sides"].values())
               for m in with_minos["result"]["minos"].values())
    r = t3.post("/api/jobs", json={"kind": "fit", "method": "engine_lm", "settings": {"minos": True}})
    assert r.status_code == 422 and r.json()["allowed"] == ["max_calls"], r.text


def test_call_limit_is_said_in_the_engine_words(t3):
    done = wait_job(t3, start(t3, method="engine_migrad", settings={"max_calls": "20"})["id"])
    r, _ = direct(t3, "migrad", max_calls=20)
    assert r.status == "call_limit" and done["state"] == "ended" and "call_limit" in done["result"]["message"]
    assert done["best"]["objective"] == r.objective
    assert done["result"]["outcome"] == r.outcome == "call_limit" and done["result"]["notices"] == r.notices
    lm = wait_job(t3, start(t3, method="engine_lm", settings={"max_calls": "20"})["id"])
    r, _ = direct(t3, "lm", max_calls=20)
    assert lm["result"]["notices"] == r.notices and r.notices, (lm["result"], r.notices)
    import pyazr
    fitters = pyazr.kinds()["fitters"]
    for row in engine_fit_rows(t3):
        f = fitters[row["method"].removeprefix("engine_")]
        assert row["facts"]["honours a call limit"] == f["call_limit"], row
        assert row["facts"]["falls back to"] == (f["fallback"] or "none"), row


def test_t2_call_limit_of_40_is_not_a_convergence(t2):
    limited = wait_job(t2, start(t2, method="engine_migrad", settings={"max_calls": "40"})["id"])
    print(f"\n  T2, max_calls 40: state {limited['state']}, outcome {limited['result']['outcome']}, best "
          f"{limited['best']['objective']!r}, notices {limited['result']['notices']}")
    assert limited["state"] == "ended" and limited["result"]["outcome"] == "call_limit", limited["result"]
    assert limited["result"]["valid"] is False


def test_accept_then_undo(t3):
    chi2_start = t3.post("/api/evaluate").json()["chi2"]
    text = t3.get("/api/project/text").text
    done = wait_job(t3, start(t3, method="engine_lm")["id"])
    r = t3.post(f"/api/jobs/{done['id']}/accept", json={"rev": rev(t3)})
    assert r.status_code == 200, r.text
    at_best = t3.get(f"/api/jobs/{done['id']}").json()["accepted"]["at_best"]
    written = t3.post("/api/evaluate").json()
    print(f"\n  accept: chi2 at the best point {at_best['chi2']!r}, of the project as written {written['chi2']!r}")
    assert written["chi2"] == pytest.approx(at_best["chi2"], rel=1e-9)
    t3.post("/api/project/undo", json={"rev": rev(t3)})
    assert t3.get("/api/project/text").text == text
    assert t3.post("/api/evaluate").json()["chi2"] == chi2_start


@pytest.mark.parametrize("method", ["engine_migrad", "engine_lm"])
def test_stop_ends_a_running_engine_fit_within_one_call(t2, method):
    job = start(t2, method=method)
    path = f"/api/jobs/{job['id']}"
    assert wait_until(lambda: ((t2.get(path).json().get("progress") or {}).get("calls") or 0) >= 3, 300)
    before = t2.get(path).json()["progress"]
    t0 = time.perf_counter()
    assert t2.post(f"{path}/stop").status_code == 200
    stopped = wait_job(t2, job["id"])
    dt = time.perf_counter() - t0
    per_call = before["elapsed"] / before["calls"]
    print(f"\n  {method} on T2: stop request to 'stopped' {dt:.3f} s; {before['calls']} calls before it, "
          f"{stopped['progress']['calls']} at the end; {stopped['result']['reports_after_stop']} report(s) after the "
          f"request; {per_call:.3f} s per call before it")
    assert stopped["state"] == "stopped" and stopped["result"]["engine_status"] == "stopped", stopped["result"]
    # T2's THM data sets have no analytic Jacobian: LM runs MIGRAD in its place, and the result says so
    fitter, notices = stopped["result"]["fitter"], stopped["result"]["notices"]
    if method == "engine_lm":
        assert fitter == {"asked": "lm", "ran": "migrad"} and any("MIGRAD" in n for n in notices), stopped["result"]
        assert "migrad in place of lm" in stopped["result"]["message"], stopped["result"]["message"]
    else:
        assert fitter == {"asked": "migrad", "ran": "migrad"} and notices == [], stopped["result"]
    assert stopped["best"]["objective"] == stopped["progress"]["best"] or method == "engine_lm"
    assert stopped["best"]["objective"] <= stopped["result"]["start_objective"] * (1 + 1e-12)  # see 11.6
    # the evaluation (MIGRAD) or iteration (LM) running when the request came is the last one reported
    assert stopped["result"]["reports_after_stop"] <= uitest.THREADS, stopped["result"]
    text = t2.get("/api/project/text").text
    r = t2.post(f"{path}/accept", json={"rev": rev(t2)})
    assert r.status_code == 200, r.text
    t2.post("/api/project/undo", json={"rev": rev(t2)})
    assert t2.get("/api/project/text").text == text


if __name__ == "__main__":
    uitest.run(__file__)
