#!/usr/bin/env python3
"""The web UI's API through FastAPI's test client on copies of T1 and T2 (docs/upgrade/web_ui.md sections 4, 5, 8, 9).

  1. evaluate equals a direct session of the same text bit for bit (objective, data chi2, every data set's model);
     the chi-squared per data set equals the command-line program's chiSquared.out on T1 and T2 (T2's s02/s03 are
     joined by a data covariance: the program's split, not evaluate().data[i].chi2's);
  2. an edit is validated by the engine before it enters the history: a refused edit returns the diagnostics and
     leaves the text and the revision as they were; undo and redo restore the text byte for byte and the objective
     exactly;
  3. a stale revision is refused (412) for a second client, a missing one (428);
  4. the user's file is written only by save (edits, evaluations and a job leave it as it was);
  5. no path outside the project's directory is read or written (403): save, open, a record's file, a whole text,
     a symbolic link leading out;
  6. the token is required on every route of the table, the Host header must be a loopback name;
  7. the error classes give their statuses with their fields; a route waiting for pyazr gives 501 naming the call;
  8. the working x moves without a rebuild; apply writes it into the text (one history entry).

Run: python3 tests/pyazr/ui/api_test.py (needs fastapi, httpx, the engine; AZURE2_BIN for the program's numbers).
"""
import hashlib
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import AUTH, copy_project, make_client  # noqa: E402


@pytest.fixture(scope="module")
def t1(tmp_path_factory):
    path = copy_project("T1", tmp_path_factory.mktemp("api"))
    with make_client(path) as client:
        client.path = path
        yield client
        workers = [p["pid"] for p in client.get("/api/server").json() if p["role"] != "server"]
    for pid in workers:
        assert uitest.wait_until(lambda: uitest.gone(pid), 10), f"worker {pid} still runs after the app ended"


def rev(client):
    return client.get("/api/project").json()["rev"]


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def direct_session(client):
    """A session of today's pyazr on the text the session worker read, in the project's directory."""
    import pyazr
    state = client.app.state.server.state
    return pyazr.Session(state.session_file(), cwd=state.directory)


# -- 1. numbers -----------------------------------------------------------------------------------------------------

def test_evaluate_equals_direct_session(t1):
    e = t1.post("/api/evaluate").json()
    s = direct_session(t1)
    x = np.array(s.params_rwa)
    assert e["objective"] == s.objective()(x)
    assert e["chi2"] == float(s._engine.calculate_chi2_rwa(x))
    models = [d.model for d in s.evaluate(x).data]
    data = t1.get("/api/evaluation").json()["datasets"]
    assert [d["name"] for d in data] == [d.name for d in s.active_datasets]
    for d, m in zip(data, models):
        assert d["model"] == np.asarray(m, float).tolist()
        assert len(d["e_cm"]) == len(d["value"]) == len(d["error"]) == len(d["model"]) == len(d["residual"])
    one = t1.get(f"/api/evaluation/data/{data[0]['name']}").json()
    assert one["model"] == data[0]["model"]
    s.close()


def _program_rows(text):
    rows = [line.split(",") for line in text.splitlines() if line and line[0].isdigit()]
    return [(float(r[1]), int(r[2]), float(r[3]), float(r[4])) for r in rows]


@pytest.mark.skipif(not os.environ.get("AZURE2_BIN"), reason="AZURE2_BIN (the command-line program) is not set")
@pytest.mark.parametrize("name", ["T1", "T2"])
def test_chi2_per_dataset_is_the_programs(name, tmp_path):
    program = _program_rows(uitest.run_program(copy_project(name, tmp_path / "cli")))
    with make_client(copy_project(name, tmp_path / "api")) as client:
        e = client.post("/api/evaluate").json()
    api = [(d["chi2"], d["n"], d["norm"], d["norm_chi2"]) for d in e["datasets"]]
    assert api == program
    total = sum(r[0] for r in program)
    assert abs(e["chi2"] - total) <= 1e-5 * total  # the file's rows have six significant digits
    if name == "T2":  # s02 and s03 share a data covariance: the program's split (web_ui.md section 2.2)
        by_name = {d["name"]: d["chi2"] for d in e["datasets"]}
        assert (by_name["s02_diff090"], by_name["s03_diff150"]) == (99.1682, 414.728)
        assert abs(e["objective"] - (e["chi2"] + e["terms"]["prior covariance"])) < 1e-9 * e["objective"]


# -- 2. edits, refusals, undo / redo -----------------------------------------------------------------------------------

def _first_channel(client):
    level = client.get("/api/project/levels").json()[0]
    ch = next(c for c in level["channels"] if not c["fixed"])
    return level["label"], ch


def test_edit_undo_redo_restore_text_and_objective(t1):
    text0 = t1.get("/api/project/text").text
    obj0 = t1.post("/api/evaluate").json()["objective"]
    label, ch = _first_channel(t1)
    r0 = rev(t1)
    r = t1.patch(f"/api/project/levels/{label}/channels/{ch['pair']}/{ch['l']}/{ch['s']}",
                 json={"width": ch["width"] * 1.5, "rev": r0})
    assert r.status_code == 200, r.text
    assert r.json()["rev"] == r0 + 1
    text1 = t1.get("/api/project/text").text
    changed = [a for a, b in zip(text0.splitlines(), text1.splitlines()) if a != b]
    assert len(changed) == 1 and len(text0.splitlines()) == len(text1.splitlines())
    obj1 = t1.post("/api/evaluate").json()["objective"]
    assert obj1 != obj0
    assert t1.post("/api/project/undo", json={"rev": r0 + 1}).status_code == 200
    assert t1.get("/api/project/text").text == text0
    assert t1.post("/api/evaluate").json()["objective"] == obj0
    assert t1.post("/api/project/redo", json={"rev": rev(t1)}).status_code == 200
    assert t1.get("/api/project/text").text == text1
    assert t1.post("/api/evaluate").json()["objective"] == obj1
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    history = t1.get("/api/project/history").json()
    assert [h["current"] for h in history].index(True) == 0


def test_refused_edit_leaves_state(t1):
    text0, r0 = t1.get("/api/project/text").text, rev(t1)
    name = next(r["name"] for r in t1.get("/api/project/records").json() if r["kind"] == "dataset")
    r = t1.patch(f"/api/project/records/{name}", json={"fields": {"norm": "abc"}, "rev": r0})
    assert r.status_code == 422, r.text
    body = r.json()
    assert body["error"] == "FormatError" and body["diagnostics"] and body["diagnostics"][0]["code"]
    r = t1.patch(f"/api/project/records/{name}", json={"fields": {"normalization": "no_such_record"}, "rev": r0})
    assert r.status_code == 422 and "E_" in r.text
    r = t1.patch("/api/project/levels/9/2+ 99.0/channels/1/0/0.5", json={"width": 1.0, "rev": r0})
    assert r.status_code == 422, r.text
    lines = text0.splitlines(keepends=True)
    header = t1.get(f"/api/project/records/{name}").json()["line"]
    k = next(i for i in range(header, len(lines)) if lines[i].split()[:1] == ["norm"])
    unknown = "".join(lines[:k] + ["  norm  no_such_normalization\n"] + lines[k + 1:])
    r = t1.put("/api/project/text", json={"text": unknown, "rev": r0})
    assert r.status_code == 422 and "E_UNKNOWN_REFERENCE" in r.text, r.text
    bad = text0.replace("<levels>", "<levels>\nnot a level row", 1)
    r = t1.put("/api/project/text", json={"text": bad, "rev": r0})
    assert r.status_code == 422, r.text
    assert t1.get("/api/project/text").text == text0 and rev(t1) == r0


def test_stale_revision_refused(t1):
    other = t1.app  # a second client on the same server: the same state
    from fastapi.testclient import TestClient
    second = TestClient(other, base_url="http://127.0.0.1", headers=AUTH)
    label, ch = _first_channel(t1)
    url = f"/api/project/levels/{label}/channels/{ch['pair']}/{ch['l']}/{ch['s']}"
    r0 = rev(t1)
    assert t1.patch(url, json={"width": ch["width"] * 1.1, "rev": r0}).status_code == 200
    r = second.patch(url, json={"width": ch["width"] * 1.2, "rev": r0})
    assert r.status_code == 412 and r.json()["rev"] == r0 + 1, r.text
    r = second.patch(url, json={"width": ch["width"] * 1.2})
    assert r.status_code == 428, r.text
    r = second.patch(url, json={"width": ch["width"]}, headers={"If-Match": str(r0 + 1)})
    assert r.status_code == 200, r.text


# -- 4. the user's file ------------------------------------------------------------------------------------------------

def test_user_file_written_only_by_save(t1):
    path = t1.path
    before, mtime = sha(path), os.stat(path).st_mtime_ns
    label, ch = _first_channel(t1)
    r = t1.patch(f"/api/project/levels/{label}/channels/{ch['pair']}/{ch['l']}/{ch['s']}",
                 json={"width": ch["width"] * 1.3, "rev": rev(t1)})
    assert r.status_code == 200
    t1.post("/api/evaluate")
    assert sha(path) == before and os.stat(path).st_mtime_ns == mtime
    state = t1.get("/api/project").json()
    assert state["unsaved"] is True
    with open(state["working_file"]) as f:
        assert f.read() == t1.get("/api/project/text").text
    r = t1.post("/api/project/save", json={"rev": state["rev"]})
    assert r.status_code == 200
    with open(path) as f:
        assert f.read() == t1.get("/api/project/text").text
    assert t1.get("/api/project").json()["unsaved"] is False
    assert not os.path.exists(state["working_file"])


def test_crash_copy_offered_back(tmp_path):
    import pyazr
    path = copy_project("T1", tmp_path)
    p = pyazr.Project.open(path)
    level = p.levels()[0]
    ch = level.channels[0]
    p.set_channel(level.label, ch.pair, ch.l, ch.s, width=2 * ch.width)
    with open(os.path.join(os.path.dirname(path), ".T1.ui.azr"), "w") as f:
        f.write(p.text)
    before = sha(path)
    with make_client(path) as client:
        state = client.get("/api/project").json()
        assert state["crash_copy"] and state["unsaved"] is False
        r = client.post("/api/project/restore", json={"rev": state["rev"]})
        assert r.status_code == 200, r.text
        assert client.get("/api/project/text").text == p.text
        assert client.get("/api/project").json()["unsaved"] is True
    assert sha(path) == before


# -- 5. paths ----------------------------------------------------------------------------------------------------------

def test_paths_outside_the_root_refused(t1, tmp_path):
    r0 = rev(t1)
    outside = os.path.join(os.path.dirname(os.path.dirname(t1.path)), "escape.azr")
    r = t1.post("/api/project/save", json={"rev": r0, "path": "../escape.azr"})
    assert r.status_code == 403 and not os.path.exists(outside), r.text
    assert t1.post("/api/project/open", json={"path": "/etc/hostname"}).status_code == 403
    name = next(r["name"] for r in t1.get("/api/project/records").json() if r["kind"] == "dataset")
    r = t1.patch(f"/api/project/records/{name}", json={"fields": {"file": "../../secret.dat"}, "rev": r0})
    assert r.status_code == 403, r.text
    link = os.path.join(os.path.dirname(t1.path), "link_out")
    os.symlink(str(tmp_path), link)
    r = t1.patch(f"/api/project/records/{name}", json={"fields": {"file": "link_out/x.dat"}, "rev": r0})
    assert r.status_code == 403, r.text
    os.remove(link)
    text = t1.get("/api/project/text").text
    data_file = t1.get(f"/api/project/records/{name}").json()["fields"]["file"].split()[0]
    r = t1.put("/api/project/text", json={"text": text.replace(data_file, "/etc/hostname", 1), "rev": r0})
    assert r.status_code == 403, r.text
    r = t1.patch("/api/project/config", json={"entries": {"priorCovariance": "/etc/hostname"}, "rev": r0})
    assert r.status_code == 403, r.text
    assert rev(t1) == r0


# -- 6. token and host -------------------------------------------------------------------------------------------------

def _sample_path(path):
    return (path.replace("{label:path}", "1/2- 0.0").replace("{name}", "x").replace("{pair}", "1")
            .replace("{l}", "1").replace("{s}", "0.5").replace("{key}", "1").replace("{job_id}", "x")
            .replace("{fn}", "x").replace("{nuclide}", "x"))


def test_token_required_on_every_route(t1):
    from fastapi.testclient import TestClient
    from pyazr.ui.routes import ROUTES
    bare = TestClient(t1.app, base_url="http://127.0.0.1")
    wrong = TestClient(t1.app, base_url="http://127.0.0.1", headers={"Authorization": "Bearer nope"})
    for r in ROUTES:
        for c in (bare, wrong):
            got = c.request(r.method, _sample_path(r.path), json={})
            assert got.status_code == 401, f"{r.method} {r.path}: {got.status_code}"
    for path in ("/openapi.json", "/docs", "/api/nonexistent"):
        assert bare.get(path).status_code == 401
    ok = bare.get(f"/api/project?token={uitest.TOKEN}")
    assert ok.status_code == 200 and "azure2_ui_token" in ok.headers.get("set-cookie", "")


def test_host_header_must_be_loopback(t1):
    r = t1.get("/api/project", headers={"Host": "evil.example:80"})
    assert r.status_code == 403
    for host in ("localhost:1234", "127.0.0.1", "[::1]", "[::1]:28400", "localhost"):
        assert t1.get("/api/project", headers={"Host": host}).status_code == 200, host
    for host in ("[::2]:80", "127.0.0.2", "localhost.evil.example"):
        assert t1.get("/api/project", headers={"Host": host}).status_code == 403, host


def test_origin_must_be_the_servers_own(t1):
    label, ch = _first_channel(t1)
    url = f"/api/project/levels/{label}/channels/{ch['pair']}/{ch['l']}/{ch['s']}"
    for origin in ("http://127.0.0.1:8888", "http://localhost", "null", "https://127.0.0.1"):
        r = t1.patch(url, json={"width": 1.0, "rev": rev(t1)}, headers={"Origin": origin})
        assert r.status_code == 403 and "Origin" in r.json()["message"], (origin, r.text)
    r = t1.get("/api/project", headers={"Origin": "http://127.0.0.1"})  # TestClient's own origin
    assert r.status_code == 200
    r = t1.get("/api/project", headers={"Host": "localhost:28400", "Origin": "http://localhost:28400"})
    assert r.status_code == 200


# -- 7. errors ---------------------------------------------------------------------------------------------------------

def test_error_statuses_and_fields(t1):
    r = t1.get("/api/evaluation/data/no_such_set")  # raised in the session worker: its classes travel
    assert r.status_code == 404 and "no_such_set" in r.json()["message"], r.text
    assert r.json()["classes"][:2] == ["ParameterKeyError", "ParameterError"]
    r = t1.get("/api/project/records/no_such_record")
    assert r.status_code == 404 and "no_such_record" in r.json()["message"], r.text
    r = t1.get("/api/jobs/nojob")
    assert r.status_code == 404 and r.json()["error"] == "NotFound"
    x = t1.get("/api/x").json()
    fixed = next(p["name"] for p in t1.get("/api/parameters").json() if not p["free"])
    r = t1.patch("/api/x", json={"values": {fixed: 1.0}})
    assert r.status_code == 422 and "fixed" in r.json()["message"], r.text
    r = t1.patch("/api/x", json={"values": {"no_such_parameter": 1.0}})
    assert r.status_code == 404 and set(r.json()["known"]) == set(x["x"])
    r = t1.post("/api/evaluate/grid", json={})
    assert r.status_code == 501 and r.json()["missing"] == "pyazr.Session.evaluate.energies"
    waiting = next(row for row in t1.get("/api/jobs/kinds").json() if row["missing"])
    r = t1.post("/api/jobs", json={"kind": waiting["kind"], "method": waiting["method"]})
    assert r.status_code == 501 and r.json()["missing"] == waiting["missing"]
    r = t1.post("/api/evaluate", content=b"[1, 2]", headers={"Content-Type": "application/json"})
    assert r.status_code == 422


def test_every_pyazr_error_class_has_a_status():
    """Every public exception class of pyazr.errors maps to a status other than 500, by the table's rule; a lookup
    that found nothing is 404, a closed session 409, a missing file 404, any other refusal 422."""
    import inspect
    import pyazr.errors as errors
    from pyazr.ui.errors import STATUS, status_of, to_http
    expected = {"SessionClosed": 409, "ParameterKeyError": 404, "ParameterIndexError": 404, "MissingFile": 404}
    classes = [c for _, c in inspect.getmembers(errors, inspect.isclass)
               if issubclass(c, Exception) and not issubclass(c, Warning)  # EngineWarning is warned, not raised
               and c.__module__ == errors.__name__ and not c.__name__.startswith("_")]
    assert classes
    for c in classes:
        names = [b.__name__ for b in c.__mro__ if b not in (BaseException, object)]
        status = status_of(names)
        assert status != 500 and any(n in dict(STATUS) for n in names), c
        assert status == expected.get(c.__name__, 422), (c, status)
        e = c("/x", "no file /x") if c.__name__ == "MissingFile" else c("a message")
        assert to_http(e)[0] == status


def test_tables_of_the_engine(t1):
    import pyazr
    kinds = t1.get("/api/kinds").json()
    assert kinds == pyazr.kinds()
    assert t1.get("/api/physics").json() == kinds["physics"]
    d = t1.get("/api/describe").json()
    assert d["package"]["pyazr"] == pyazr.__version__ and d["package"]["coulomb_table"]
    assert set(d["session_options"]) == set(kinds["options"])
    assert all(o["value"] == o["default"] and o["source"] == "default" for o in d["session_options"].values())


def test_session_options_are_config_keys(t1):
    """An option of the calculation is a <config> key (project_format3.md section 3.12): set through the config route
    (one history entry, a new session), reported with source "project" by /api/describe and /api/session, and the
    session computes what a direct session with the keyword computes; a contradiction is refused, the text unchanged."""
    r0 = rev(t1)
    r = t1.patch("/api/project/config", json={"entries": {"use_brune": "0"}, "rev": r0})
    assert r.status_code == 200, r.text
    assert rev(t1) != r0
    d = t1.get("/api/describe").json()["session_options"]
    assert d["use_brune"]["value"] is False and d["use_brune"]["source"] == "project"
    assert d["transform"]["source"] == "default"
    session = t1.get("/api/session").json()
    assert session["session_options"]["use_brune"]["source"] == "project"
    assert session["options"]["use_brune"] is False
    e = t1.post("/api/evaluate").json()
    import pyazr
    with pyazr.Session(t1.path, use_brune=False) as s:
        assert e["objective"] == s.objective()(s.params_rwa)
    r1 = rev(t1)
    r = t1.patch("/api/project/config", json={"entries": {"use_park": "1"}, "rev": r1})
    assert r.status_code == 422 and "E_CONFIG_OPTIONS_CONTRADICT" in r.text, r.text
    assert rev(t1) == r1
    r = t1.patch("/api/project/config", json={"entries": {"use_brune": None}, "rev": r1})  # back to T1 for the others
    assert r.status_code == 200 and t1.get("/api/describe").json()["session_options"]["use_brune"]["source"] == "default"


def test_evaluation_carries_the_engines_notices(tmp_path):
    """An engine notice (pyazr.EngineWarning) of the session worker reaches the client in the answer's warnings: the
    UI's output directory made read-only, the session build cannot write its external-capture cache."""
    path = copy_project("T1", tmp_path)
    engine_dir = os.path.join(os.path.dirname(path), ".T1.ui", "engine")
    os.makedirs(engine_dir)
    os.chmod(engine_dir, 0o555)
    try:
        with make_client(path) as client:
            name = next(r["name"] for r in client.get("/api/project/records").json() if r["kind"] == "dataset")
            r = client.get(f"/api/evaluation/data/{name}")  # the first request: the session is built in it
    finally:
        os.chmod(engine_dir, 0o755)
    assert r.status_code == 200, r.text
    notices = [w for w in r.json()["warnings"] if w["category"] == "EngineWarning"]
    assert notices and "intEC.dat" in notices[0]["message"], r.json()["warnings"]


def test_openapi_lists_every_route(t1):
    from pyazr.ui.routes import ROUTES
    spec = t1.get("/openapi.json").json()
    listed = {(m.upper(), p) for p, ops in spec["paths"].items() for m in ops}
    table = {(r.method, r.path.replace(":path", "")) for r in ROUTES}
    assert table == {(m, p) for m, p in listed if p.startswith("/api")}


# -- 8. the working x -----------------------------------------------------------------------------------------------

def test_x_moves_without_rebuild_and_apply_writes_it(t1):
    session = t1.get("/api/session").json()
    x = t1.get("/api/x").json()["x"]
    name = next(iter(x))
    r = t1.patch("/api/x", json={"values": {name: x[name] * 1.01}})
    assert r.status_code == 200
    e = t1.post("/api/evaluate").json()
    assert t1.get("/api/session").json()["pid"] == session["pid"]
    s = direct_session(t1)
    names = [p.name for p in s.parameters if p.free_index is not None and not p.fixed]
    xs = np.array([r.json()["x"][n] for n in names])
    assert e["objective"] == s.objective()(xs)
    s.close()
    r0 = rev(t1)
    r = t1.post("/api/project/apply", json={"rev": r0})
    assert r.status_code == 200, r.text
    assert r.json()["rev"] == r0 + 1
    after = t1.post("/api/evaluate").json()
    assert after["objective"] == pytest.approx(e["objective"], rel=1e-6)
    assert t1.get("/api/x").json()["x0"][name] == pytest.approx(x[name] * 1.01, rel=1e-6)
    t1.post("/api/project/undo", json={"rev": rev(t1)})


if __name__ == "__main__":
    uitest.run(__file__)
