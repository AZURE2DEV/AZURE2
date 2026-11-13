#!/usr/bin/env python3
"""Jobs of the web UI on a copy of T1 (docs/upgrade/web_ui.md sections 3.2, 4.4, 9):

  1. a recipe fit (scipy L-BFGS-B) runs to its end in a job worker and ends at the same objective, after the same
     number of calls, as the same recipe called directly on a session of the same text; its progress reports
     increasing call numbers; its result can be accepted into the project text (one history entry), after which the
     evaluation is the job's best;
  2. a cooperative stop ends it before the next objective call with the best point (state "stopped"); a kill ends it
     at once (state "killed", the worker gone by its PID); a new job resumes from the killed job's checkpoint (its
     first objective is the checkpoint's best);
  3. a request over the compute budget is refused (409) with the running job's id and its stop route;
  4. the iminuit recipe (MIGRAD and HESSE) as a row of the same table (skipped without iminuit);
  5. a job of the table that waits for pyazr answers 501 naming the call;
  6. the job table lists every setting with its type, default and doc; a request's value is converted by its type
     or refused with what the type allows;
  7. accepting a fit writes into the project text exactly what pyazr's own write-back (``save_fit``) writes for the
     best point, and Undo gives the text back byte for byte; a stopped job's best point is accepted the same way;
  8. a job's trace (objective against call, thinned past TRACE_POINTS) and log tail; a GET with ``wait`` answers when
     the state changes; the parameter table and x stay readable while the job holds the budget;
  9. iminuit's MINOS as a job from an iminuit fit's best point: asymmetric intervals, covariance and correlation by
     name.

The stop and kill times are printed (measured on the product).
Run: python3 tests/pyazr/ui/jobs_test.py
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import THREADS, copy_project, make_client, wait_job, wait_until  # noqa: E402

LONG = {"max_iterations": 100000}


@pytest.fixture(scope="module")
def t1(tmp_path_factory):
    path = copy_project("T1", tmp_path_factory.mktemp("jobs"))
    with make_client(path) as client:
        client.path = path
        yield client


def rev(client):
    return client.get("/api/project").json()["rev"]


def start(client, **body):
    r = client.post("/api/jobs", json=dict({"kind": "fit", "method": "scipy_lbfgsb", "x0": "text"}, **body))
    assert r.status_code == 200, r.text
    return r.json()


def direct(client, recipe, settings):
    """The recipe called directly, in this process, on a session of the text the job read."""
    import pyazr
    from pyazr.ui.jobs import Tracker
    state = client.app.state.server.state
    path = state.session_file("direct")  # the session's text of the current revision, beside the project
    with open(path, "wb") as f:
        f.write(state.session_text().encode())
    s = pyazr.azure2(path, cwd=state.directory)
    names = [p.name for p in s.parameters if p.free_index is not None and not p.fixed]
    tracker = Tracker(s.objective_and_grad, names)
    status = recipe(tracker, np.asarray(s.params_rwa, float), settings, names)
    s.close()
    os.remove(path)
    return tracker, status


def test_fit_to_end_equals_the_recipe_called_directly(t1):
    seen = []
    job = start(t1)
    wait_until(lambda: seen.append((t1.get(f"/api/jobs/{job['id']}").json().get("progress") or {}).get("calls"))
               or seen[-1] is not None and seen[-1] > 3, 60)
    done = wait_job(t1, job["id"])
    assert done["state"] == "converged", done
    from pyazr.ui.jobs import fit_scipy_lbfgsb
    tracker, status = direct(t1, fit_scipy_lbfgsb, done["settings"])
    assert status[0] == "converged"
    assert done["best"]["objective"] == tracker.best
    assert done["result"]["calls"] == tracker.calls
    calls = [c for c in seen if c is not None]
    assert calls == sorted(calls) and calls[-1] > 3
    assert uitest.gone(done["pid"]) if done.get("pid") else True
    checkpoint = os.path.join(t1.app.state.server.state.ui_dir, "jobs", f"{job['id']}.json")
    with open(checkpoint) as f:
        assert json.load(f)["best"]["objective"] == tracker.best
    # accept: the best point becomes the project's text
    r0 = rev(t1)
    r = t1.post(f"/api/jobs/{job['id']}/accept", json={"rev": r0})
    assert r.status_code == 200, r.text
    assert r.json()["rev"] == r0 + 1
    assert t1.post("/api/evaluate").json()["objective"] == pytest.approx(tracker.best, rel=1e-6)
    t1.post("/api/project/undo", json={"rev": rev(t1)})


def _running(client, job_id, calls=3):
    return wait_until(lambda: ((client.get(f"/api/jobs/{job_id}").json().get("progress") or {}).get("calls") or 0)
                      >= calls, 60)


def test_stop_kill_and_resume_from_the_checkpoint(t1):
    job = start(t1, settings=LONG)
    assert _running(t1, job["id"])
    t0 = time.perf_counter()
    assert t1.post(f"/api/jobs/{job['id']}/stop").status_code == 200
    stopped = wait_job(t1, job["id"])
    print(f"\n  stop: request to 'stopped' {time.perf_counter() - t0:.3f} s")
    assert stopped["state"] == "stopped" and stopped["best"]["x"], stopped

    job = start(t1, settings=LONG)
    assert _running(t1, job["id"])
    pid = t1.get(f"/api/jobs/{job['id']}").json()["pid"]
    t0 = time.perf_counter()
    assert t1.post(f"/api/jobs/{job['id']}/kill").status_code == 200
    killed = wait_job(t1, job["id"])
    print(f"  kill: request to 'killed' {time.perf_counter() - t0:.3f} s")
    assert killed["state"] == "killed"
    assert wait_until(lambda: uitest.gone(pid), 10)
    checkpoint = os.path.join(t1.app.state.server.state.ui_dir, "jobs", f"{job['id']}.json")
    with open(checkpoint) as f:
        saved = json.load(f)
    assert saved["state"] == "killed" and saved["best"]["x"]

    resumed = start(t1, x0={"job": job["id"]}, settings=LONG)
    assert _running(t1, resumed["id"], 1)
    t1.post(f"/api/jobs/{resumed['id']}/stop")
    again = wait_job(t1, resumed["id"])
    assert again["state"] == "stopped"
    assert again["result"]["start_objective"] == saved["best"]["objective"]
    assert again["best"]["objective"] <= saved["best"]["objective"]


def test_budget_refuses_with_the_job(t1):
    job = start(t1, settings=LONG)
    assert _running(t1, job["id"], 1)
    try:
        if THREADS == 1:
            r = t1.post("/api/evaluate")
            assert r.status_code == 409, r.text
            assert r.json()["job"] == job["id"] and r.json()["stop"] == f"/api/jobs/{job['id']}/stop"
        others = [start(t1, settings=LONG)["id"] for _ in range(max(0, THREADS - 2))]
        r = t1.post("/api/jobs", json={"kind": "fit", "method": "scipy_lbfgsb"})
        assert r.status_code == 409 and r.json()["job"] in [job["id"]] + others, r.text
    finally:
        for j in t1.get("/api/jobs").json():
            t1.post(f"/api/jobs/{j['id']}/kill")
            wait_job(t1, j["id"])


@pytest.mark.skipif(importlib.util.find_spec("iminuit") is None, reason="iminuit is not installed")
def test_iminuit_recipe(t1):
    job = start(t1, method="iminuit_migrad")
    done = wait_job(t1, job["id"])
    assert done["state"] == "converged" and done["result"]["errors"], done
    assert done["result"]["message"] == "MIGRAD: valid minimum; HESSE: still valid", done["result"]
    from pyazr.ui.jobs import fit_iminuit_migrad
    tracker, status = direct(t1, fit_iminuit_migrad, done["settings"])
    assert done["best"]["objective"] == tracker.best and status[0] == "converged"


def test_missing_job_kinds_answer_501(t1):
    table = t1.get("/api/jobs/kinds").json()
    for row in table:
        if row["missing"]:
            r = t1.post("/api/jobs", json={"kind": row["kind"], "method": row["method"]})
            assert r.status_code == 501 and r.json()["missing"] == row["missing"]


def test_settings_listed_with_types_and_converted(t1):
    from pyazr.ui.jobs import SETTING_TYPES
    for row in t1.get("/api/jobs/kinds").json():
        for setting in row["settings"]:
            assert set(setting) == {"name", "type", "default", "doc", "allowed"}, setting
            assert setting["type"] in SETTING_TYPES and setting["doc"], setting
    r = t1.post("/api/jobs", json={"kind": "fit", "method": "scipy_lbfgsb", "settings": {"max_iterations": "1.5"}})
    assert r.status_code == 422 and "integer" in r.json()["allowed"], r.text
    r = t1.post("/api/jobs", json={"kind": "fit", "method": "scipy_lbfgsb", "settings": {"max_calls": 3}})
    assert r.status_code == 422 and r.json()["allowed"] == ["max_iterations"], r.text
    job = start(t1, settings={"max_iterations": "2"})
    done = wait_job(t1, job["id"])
    assert done["settings"]["max_iterations"] == 2 and done["state"] in ("converged", "ended"), done


def written_back(client, best_x):
    """What pyazr's own write-back gives for ``best_x`` ({name: value}): save_fit on a session of the user's file."""
    import pyazr
    path = client.path
    s = pyazr.azure2(path, cwd=os.path.dirname(path))
    names = [p.name for p in s.parameters if p.free_index is not None and not p.fixed]
    out = os.path.join(os.path.dirname(path), "written_back.azr")
    s.save_fit(out, np.array([best_x[n] for n in names]), param_sav=False)
    s.close()
    with open(out, "rb") as f:
        text = f.read().decode("utf-8", "surrogateescape")
    os.remove(out)
    return text


def accept_and_compare(client, job_id):
    before = client.get("/api/project/text").text
    best = client.get(f"/api/jobs/{job_id}").json()["best"]["x"]
    r = client.post(f"/api/jobs/{job_id}/accept", json={"rev": rev(client)})
    assert r.status_code == 200, r.text
    assert client.get("/api/project/text").text == written_back(client, best)
    accepted = client.get(f"/api/jobs/{job_id}").json()["accepted"]
    assert accepted["rev"] == r.json()["rev"] == rev(client)
    best_objective = client.get(f"/api/jobs/{job_id}").json()["best"]["objective"]
    assert accepted["at_best"]["objective"] == pytest.approx(best_objective, rel=1e-9)
    client.post("/api/project/undo", json={"rev": rev(client)})
    assert client.get("/api/project/text").text == before


def test_accept_writes_what_pyazr_writes_back_and_undo_restores(t1):
    job = start(t1, method="iminuit_migrad")
    done = wait_job(t1, job["id"])
    assert done["state"] == "converged", done
    accept_and_compare(t1, job["id"])


def test_stopped_job_leaves_its_best_point_which_is_accepted(t1):
    job = start(t1, settings=LONG)
    assert _running(t1, job["id"], 3)
    t1.post(f"/api/jobs/{job['id']}/stop")
    stopped = wait_job(t1, job["id"])
    assert stopped["state"] == "stopped" and stopped["best"]["objective"] <= stopped["result"]["start_objective"]
    assert set(stopped["best"]["x"]) == set(stopped["result"]["start_x"])
    accept_and_compare(t1, job["id"])


def test_trace_log_wait_and_parameters_during_a_job(t1):
    from pyazr.ui.jobs import TRACE_POINTS, Job
    assert t1.get("/api/parameters").status_code == 200  # the session of the current text is built
    job = start(t1, settings=LONG)
    try:
        assert _running(t1, job["id"], 5)
        rows = t1.get("/api/parameters")
        assert rows.status_code == 200 and any(r["free"] for r in rows.json()), rows.text
        assert t1.get("/api/x").status_code == 200
        if THREADS == 1:
            assert t1.post("/api/evaluate").status_code == 409
        record = t1.get(f"/api/jobs/{job['id']}").json()
        # an answer's "rev" is the project's current revision (the page re-reads on a new one); a job names the
        # revision of the text it runs on otherwise
        assert "rev" not in record and record["text_rev"] == rev(t1), record.keys()
        calls = [p[0] for p in record["trace"]]
        assert calls == sorted(calls) and len(calls) >= 5, record["trace"][:10]
        assert all(p[2] <= p[1] for p in record["trace"])  # best so far <= the call's objective
        assert isinstance(record["log"], list)
        t0 = time.perf_counter()
        t1.post(f"/api/jobs/{job['id']}/stop")
        waited = t1.get(f"/api/jobs/{job['id']}", params={"wait": 8}).json()
        assert waited["state"] != "running" and time.perf_counter() - t0 < 6, waited["state"]
    finally:
        wait_job(t1, job["id"])
    t0 = time.perf_counter()
    t1.get(f"/api/jobs/{job['id']}", params={"wait": 8})  # an ended job answers at once
    assert time.perf_counter() - t0 < 1
    lister = t1.get("/api/jobs").json()
    assert all("trace" not in j and "best" not in j for j in lister)
    j = Job({"id": "x"})  # the thinning of a long trace keeps the first call and the shape of the whole run
    with j.lock:
        for c in range(1, 5 * TRACE_POINTS):
            j.add_to_trace(c, 1000.0 / c, 1000.0 / c)
    kept = [p[0] for p in j.record["trace"]]
    assert len(kept) <= TRACE_POINTS and kept[-1] > 4 * TRACE_POINTS and kept == sorted(kept), kept[-5:]


@pytest.mark.skipif(importlib.util.find_spec("iminuit") is None, reason="iminuit is not installed")
def test_minos_after_an_iminuit_fit(t1):
    kinds = {(r["kind"], r["method"]): r for r in t1.get("/api/jobs/kinds").json()}
    assert kinds[("minos", "iminuit")]["follows"] == ["fit", "iminuit_migrad"]
    fit = wait_job(t1, start(t1, method="iminuit_migrad")["id"])
    names = list(fit["best"]["x"])[:2]
    job = start(t1, kind="minos", method="iminuit", x0={"job": fit["id"]}, settings={"parameters": " ".join(names)})
    done = wait_job(t1, job["id"])
    assert done["state"] == "converged", done
    result = done["result"]
    assert sorted(result["minos"]) == sorted(names)
    for n in names:
        assert result["minos"][n]["lower"] < 0 < result["minos"][n]["upper"], result["minos"][n]
    cov, cor = result["covariance"], result["correlation"]
    assert cov["names"] == cor["names"] == list(result["errors"])
    k = len(cov["names"])
    assert all(abs(cor["matrix"][i][i] - 1) < 1e-12 for i in range(k))
    assert all(abs(cov["matrix"][i][i] ** 0.5 - result["errors"][cov["names"][i]]) < 1e-9 * max(1, cov["matrix"][i][i])
               for i in range(k))


if __name__ == "__main__":
    uitest.run(__file__)
