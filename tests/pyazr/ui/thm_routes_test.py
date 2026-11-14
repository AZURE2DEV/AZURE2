#!/usr/bin/env python3
"""The THM tab's routes (docs/upgrade/web_ui.md sections 6.7 and 11.8) through FastAPI's test client:

  1. on T2 (one plane-wave experiment) and on a copy of tests/18O_p_a_thm with a coherent background, a Coulomb
     distortion factor and the line shape: every diagnostic of /api/session/functions and every raw route
     /api/session/thm_* equals the direct Session call on the text the worker read, bit for bit;
  2. an option edit and an experiment edit through the project routes: one history entry each, the chi-squared of
     the edited text equal to the program's chiSquared.out on the same text, Undo back to the start chi-squared
     exactly; a new experiment with the data sets that join it, its removal refused while they name it;
  3. a project without THM content (T1): GET /api/project/thm is null, the diagnostics answer with a refusal
     (what / allowed / fix), never a 500;
  4. the keys served (GET /api/thm/keys) are the engine's table (pyazr.kinds()["thm"]), a key added to it is served
     with no line of the page; each diagnostic row is pyazr's row (Session.describe()["functions"]) with the physics
     rows' fields, and a diagnostic added to pyazr's table is drawn with no line of the page.

Run: python3 tests/pyazr/ui/thm_routes_test.py (needs fastapi, httpx, the engine; AZURE2_BIN for the program).
"""
import os
import re
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import uitest  # noqa: E402
from uitest import ROOT, copy_project, make_client  # noqa: E402

E_AXIS = {"from": 0.1, "to": 1.6, "points": 31, "spacing": "linear"}
ENERGIES = np.linspace(0.1, 1.6, 31)


# The two THM copies of tests/18O_p_a_thm (the engine refuses the line shape beside a coherent background).
VARIANTS = {"dist": dict(beam="18O", target="3He", spectator="d", Ebeam="115", distortion="coulomb", lineshape=True),
            "cbkg": dict(background="linear", cbackground="1/2+:2:0.5,0,0.5,1=0.25,-0.5")}


def thm_copy(dest, variant):
    """A copy of tests/18O_p_a_thm whose two THM data sets form experiment E with the keys of ``variant``: the
    kinematics of a 3He Trojan horse with a Coulomb distortion factor and the line shape (dist), or a linear and a
    coherent background (cbkg)."""
    from thm_records import set_experiment
    from pyazr import Project
    target = os.path.join(str(dest), "o18thm")
    shutil.copytree(os.path.join(ROOT, "tests", "18O_p_a_thm"), target,
                    ignore=shutil.ignore_patterns("output", "checks", "expected", ".*.ui*"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(target, d), exist_ok=True)
    p = Project.open(os.path.join(target, "18O_p_a_thm.azr"))
    set_experiment(p, "E", [1, 2], **VARIANTS[variant])
    path = os.path.join(target, "o18thm.azr")
    p.save(path)
    return path


@pytest.fixture(scope="module")
def t2(tmp_path_factory):
    path = copy_project("T2", tmp_path_factory.mktemp("thm_t2"))
    with make_client(path) as client:
        client.path = path
        yield client


@pytest.fixture(scope="module", params=sorted(VARIANTS))
def o18(request, tmp_path_factory):
    path = thm_copy(tmp_path_factory.mktemp("thm_o18"), request.param)
    with make_client(path) as client:
        client.variant = request.param
        client.path = path
        yield client


def direct_session(client):
    """A session of pyazr on the text the session worker read, in the project's directory."""
    import pyazr
    state = client.app.state.server.state
    return pyazr.Session(state.session_file(), cwd=state.directory)


def call(client, fn, **arguments):
    r = client.post(f"/api/session/functions/{fn}", json={"arguments": arguments})
    assert r.status_code == 200, r.text
    return r.json()


def plain(v):
    """A direct call's value as the JSON the route sends (numpy arrays as lists)."""
    if isinstance(v, dict):
        return {k: plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    return v.tolist() if hasattr(v, "tolist") else v


# -- 1. the diagnostics equal the session's calls --------------------------------------------------------------------

def check_diagnostics(client, experiment, calls):
    client.post("/api/evaluate")  # the worker's session of the current text
    s = direct_session(client)
    try:
        x = np.array(s.parameters.x0)
        direct = s.thm_experiments(x)
        answer = call(client, "thm_experiments")
        rows = {r["experiment"]: r for r in answer["tables"][0]["rows"]}
        assert set(rows) == set(direct)
        for name, d in direct.items():
            r = rows[name]
            assert (r["norm"], r["chi2"], r["points"], r["status"]) == (d["norm"], d["chi2"], d["points"], d["status"])
            assert [r[f"b{k}"] for k in range(len(d["b"]))] == d["b"].tolist()
            assert [r["norm_sigma"]] + [r[f"b{k}_sigma"] for k in range(len(d["b"]))] == d["sigma"].tolist()
        raw = client.post("/api/session/thm_experiments", json={}).json()
        assert {k: v for k, v in raw.items() if k not in ("rev", "warnings")} == plain(direct)
        b = call(client, "thm_background", experiment=experiment)
        db = s.thm_background(experiment, x)
        assert [r["value"] for r in b["tables"][0]["rows"]] == [db["norm"], *db["b"].tolist()]
        assert [list(r.values())[1:] for r in b["tables"][1]["rows"]] == db["cov"].tolist()
        for fn in calls:
            answer = call(client, fn, experiment=experiment, energies=E_AXIS)
            assert answer["arguments"]["energies"] == ENERGIES.tolist()
            d = getattr(s, fn)(experiment, ENERGIES, *([x] if fn != "thm_distortion" else []))
            raw = client.post(f"/api/session/{fn}", json={"experiment": experiment,
                                                           "energies": ENERGIES.tolist()}).json()
            assert {k: v for k, v in raw.items() if k not in ("rev", "warnings")} == plain(d), fn
            if fn == "thm_vertex":
                levels = [lev for c in d["channels"] for lev in c["levels"]]
                assert answer["returns"]["M2"] == [lev["M2"].tolist() for lev in levels]
                assert answer["returns"]["M2_qf"] == [lev["M2_qf"].tolist() for lev in levels]
                assert len(answer["lines"]["M2"]) == len(levels)
                assert answer["status"] == [0 if r else -3 for r in d["reached"]]
            elif fn == "thm_lineshape":
                levels = [lev for e in d["exits"] for lev in e["levels"]]
                assert answer["returns"]["NC2"] == [lev["NC2"].tolist() for lev in levels]
                assert answer["returns"]["zeta"] == [e["zeta"].tolist() for e in d["exits"]]
                assert answer["returns"]["E_sF"] == d["E_sF"].tolist()
            else:
                for k in ("R", "R_model", "M2", "M2_PW", "q"):
                    assert answer["returns"][k] == d[k].tolist(), k
            facts = {f["name"]: f["value"] for f in answer["facts"]}
            assert facts["B"] == d["B"], fn
            assert f"session.{fn}(experiment={experiment!r}, energies=energies)" in answer["python"]
    finally:
        s.close()


def test_t2_diagnostics_equal_the_session(t2):
    assert t2.get("/api/project/thm").json()["experiments"] == {"LC09": {"beam": "15N", "target": "d",
                                                                         "spectator": "n", "Ebeam": "60"}}
    check_diagnostics(t2, "LC09", ["thm_vertex"])
    r = t2.post("/api/session/functions/thm_lineshape", json={"arguments": {"experiment": "LC09",
                                                                             "energies": E_AXIS}})
    # the engine's refusal, its words: a message of the engine's own text has no parts (what, allowed, fix are None)
    assert r.status_code == 422 and "line shape" in r.json()["message"], r.text


def test_distortion_lineshape_and_coherent_background_equal_the_session(o18):
    if o18.variant == "dist":
        check_diagnostics(o18, "E", ["thm_vertex", "thm_lineshape", "thm_distortion"])
        # the distortion factor's own bound (about 17 ms per energy with an optical potential)
        from pyazr.ui.thm import DISTORTION_ELEMENTS
        r = o18.post("/api/session/functions/thm_distortion", json={"arguments": {"experiment": "E", "energies": {
            "from": 0.1, "to": 1.0, "points": DISTORTION_ELEMENTS + 1, "spacing": "linear"}}})
        assert r.status_code == 422 and str(DISTORTION_ELEMENTS) in r.json()["allowed"], r.text
        return
    check_diagnostics(o18, "E", ["thm_vertex"])
    table = call(o18, "thm_experiments")["tables"]
    assert table[1]["title"].startswith("coherent background")
    assert [r["value"] for r in table[1]["rows"]] == [0.25, -0.5]


# -- 2. edits: one history entry, the program's chi-squared, undo exact -------------------------------------------

def _program_rows(text):
    rows = [line.split(",") for line in text.splitlines() if line and line[0].isdigit()]
    return [(float(r[1]), int(r[2]), float(r[3]), float(r[4])) for r in rows]


def _program_of_text(client, tmp_path):
    """The program's chiSquared.out rows on the project's current text (a copy of its directory)."""
    target = tmp_path / f"cli{client.get('/api/project').json()['rev']}"
    source = os.path.dirname(client.path)
    shutil.copytree(source, target, ignore=shutil.ignore_patterns("output", "checks", ".*.ui*"))
    for d in ("output", "checks"):
        os.makedirs(target / d, exist_ok=True)
    azr = target / os.path.basename(client.path)
    azr.write_text(client.get("/api/project/text").text)
    return _program_rows(uitest.run_program(str(azr)))


def _rows(e):
    return [(d["chi2"], d["n"], d["norm"], d["norm_chi2"]) for d in e["datasets"]]


@pytest.mark.skipif(not os.environ.get("AZURE2_BIN"), reason="AZURE2_BIN (the command-line program) is not set")
def test_option_and_experiment_edits_move_chi2_to_the_programs_and_undo(t2, tmp_path):
    rev0 = t2.get("/api/project").json()["rev"]
    e0 = t2.post("/api/evaluate").json()
    for what, request in (("option", lambda r: t2.patch("/api/project/thm", json={"options": {"kinematics": "triple"},
                                                                                   "rev": r})),
                          ("experiment", lambda r: t2.put("/api/project/thm/experiments/LC09",
                                                          json={"keys": {"background": "const"}, "rev": r}))):
        r = request(rev0)
        assert r.status_code == 200 and r.json()["rev"] == rev0 + 1, (what, r.text)
        e1 = t2.post("/api/evaluate").json()
        assert e1["chi2"] != e0["chi2"], what
        assert _rows(e1) == _program_of_text(t2, tmp_path), what
        assert t2.post("/api/project/undo", json={"rev": rev0 + 1}).status_code == 200
        back = t2.post("/api/evaluate").json()
        assert (back["chi2"], back["objective"]) == (e0["chi2"], e0["objective"]), what
        rev0 = t2.get("/api/project").json()["rev"]


def test_new_experiment_joined_and_removal_refused(t2):
    rev0 = t2.get("/api/project").json()["rev"]
    sets = [r["name"] for r in t2.get("/api/project/records").json() if "thm" in r["fields"]]
    r = t2.put("/api/project/thm/experiments/NEW", json={"keys": {"background": "const"}, "datasets": sets[:1],
                                                          "rev": rev0})
    assert r.status_code == 200, r.text
    assert t2.get(f"/api/project/records/{sets[0]}").json()["fields"]["thm"] == "NEW"
    assert t2.get("/api/project/thm").json()["experiments"]["NEW"] == {"background": "const"}
    refused = t2.delete("/api/project/thm/experiments/NEW", params={"rev": rev0 + 1})
    assert refused.status_code == 422 and refused.json()["what"], refused.text
    assert t2.get("/api/project").json()["rev"] == rev0 + 1
    assert t2.post("/api/project/undo", json={"rev": rev0 + 1}).status_code == 200
    assert "NEW" not in t2.get("/api/project/thm").json()["experiments"]


# -- 3. a project without THM content ---------------------------------------------------------------------------------

def test_no_thm_content_refuses(tmp_path):
    with make_client(copy_project("T1", tmp_path)) as client:
        assert client.get("/api/project/thm").json() is None
        assert call(client, "thm_experiments")["tables"][0]["rows"] == []
        # a refusal of pyazr in the session worker carries its parts (e.what, e.allowed, e.fix) as the body's fields
        import pyazr
        client.post("/api/evaluate")
        s = direct_session(client)
        try:
            with pytest.raises(pyazr.Error) as direct:
                s.thm_background("LC09")
        finally:
            s.close()
        body = client.post("/api/session/functions/thm_background", json={"arguments": {"experiment": "LC09"}}).json()
        e = direct.value
        assert e.what and (body["what"], body["allowed"], body["fix"], body["message"]) == (
            e.what, e.allowed, e.fix, str(e).strip("'\"")), body
        for fn in ("thm_background", "thm_vertex", "thm_lineshape", "thm_distortion"):
            r = client.post(f"/api/session/functions/{fn}", json={"arguments": {"experiment": "LC09",
                                                                                 "energies": E_AXIS}})
            assert 400 <= r.status_code < 500 and r.json()["message"], (fn, r.status_code, r.text)
            r = client.post(f"/api/session/functions/{fn}", json={"arguments": {"energies": E_AXIS}})
            body = r.json()
            assert r.status_code == 422 and "experiment" in body["what"], (fn, r.text)
            assert body["allowed"] and body["fix"], (fn, r.text)


# -- 4. the tables -----------------------------------------------------------------------------------------------

def test_keys_are_the_engines_table(t2, monkeypatch):
    import pyazr
    served = t2.get("/api/thm/keys").json()
    assert served == pyazr.kinds()["thm"]
    for kind in served.values():
        for key, field in kind["fields"].items():
            assert field["type"] in ("enum", "number", "integer", "text"), key
            assert field["type"] != "enum" or field.get("values"), key
    # a key the engine's table gains (here: pyazr's answer with one row more) is served with no line of the page
    table = pyazr.kinds()
    table["thm"]["experiment"]["fields"]["newKey"] = {"type": "number", "doc": "a key added", "required": False}
    monkeypatch.setattr(pyazr, "kinds", lambda: table)
    assert t2.get("/api/thm/keys").json()["experiment"]["fields"]["newKey"]["doc"] == "a key added"


def test_diagnostic_rows_are_pyazrs_rows(t2):
    from pyazr.session_functions import session_functions
    rows = t2.get("/api/session/functions").json()
    physics = next(iter(t2.get("/api/physics").json().values()))
    pyazr_rows = session_functions()
    assert sorted(rows) == sorted(pyazr_rows)
    thm = {n: r for n, r in rows.items() if r.get("tab") == "thm"}
    assert sorted(thm) == ["thm_background", "thm_distortion", "thm_experiments", "thm_lineshape", "thm_vertex"]
    for name, row in rows.items():
        assert set(row) >= set(physics), name
        assert all(r["name"] and "unit" in r for r in row["returns"]), name
        if row["axis"]:
            assert next(a for a in row["arguments"] if a["name"] == row["axis"])["sweep"], name
        mine = pyazr_rows[name]
        assert (row["summary"], row["doc"], row["axis"]) == (mine["summary"], mine["doc"], mine["axis"]), name
        # every argument of pyazr's row but the fit vector (the server gives the working x), with its own fields
        assert [a["name"] for a in row["arguments"]] == [a["name"] for a in mine["arguments"] if a["name"] != "params"]
        for a in row["arguments"]:
            assert {k: a[k] for k in next(m for m in mine["arguments"] if m["name"] == a["name"])} == next(
                m for m in mine["arguments"] if m["name"] == a["name"]), (name, a["name"])
        units = {r["name"]: r["unit"] for r in mine["returns"]}
        assert all(units[r["name"]] == r["unit"] for r in row["returns"]), name
    assert next(a for a in rows["channel_functions"]["arguments"] if a["name"] == "J")["required"] is False
    assert next(a for a in rows["thm_vertex"]["arguments"] if a["name"] == "experiment")["required"] is True


def test_a_new_diagnostic_needs_no_line_of_the_page(t2, monkeypatch):
    """A row added to pyazr's table (here a copy of thm_vertex under a new name, its method on the session class) is
    listed by the Physics tab and answered with the numbers of the call it makes."""
    import pyazr
    import pyazr.session_functions as sf
    rows = sf.session_functions()
    monkeypatch.setattr(pyazr.Session, "thm_vertex_copy", pyazr.Session.thm_vertex, raising=False)
    new = dict(rows["thm_vertex"], name="thm_vertex_copy", library="pyazr: Session.thm_vertex_copy")
    monkeypatch.setitem(sf.SESSION_FUNCTIONS, "thm_vertex_copy", new)
    served = t2.get("/api/session/functions").json()["thm_vertex_copy"]
    assert served["tab"] is None and served["group"] == "this project"
    a = call(t2, "thm_vertex", experiment="LC09", energies=E_AXIS)
    # the worker is another process: its session class gains the method by the row's name only there, so the
    # answer is made here from the session's own answer, walked by the new row
    from pyazr.ui.physics import answer_parts
    s = direct_session(t2)
    try:
        out = plain(s.thm_vertex("LC09", ENERGIES, np.array(s.parameters.x0)))
    finally:
        s.close()
    parts = answer_parts("thm_vertex_copy", new, out, len(ENERGIES))
    assert parts["returns"]["M2"] == a["returns"]["M2"] and parts["lines"]["M2"] == a["lines"]["M2"]
    assert len(parts["status"]) == len(ENERGIES)  # the default rule: -1 where a field is not finite
    assert {f["name"]: f["value"] for f in parts["facts"]} == {f["name"]: f["value"] for f in a["facts"]}


if __name__ == "__main__":
    uitest.run(__file__)
