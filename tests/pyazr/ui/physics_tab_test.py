#!/usr/bin/env python3
"""What the web UI's Physics tab asks of the server (docs/upgrade/web_ui.md section 6.6):

  1. every row of /api/physics plots without typing: its defaults and examples, the axis from its sweep made by the
     server, every element computed;
  2. the answer's Python text ("Copy as Python"), executed, makes the route's arrays bit for bit, NaN where the
     route's status is not 0 (every row; a masked axis; several values as lines; a list of choices by name);
  3. the axis {from, to, points, spacing} is made with numpy and refused with what / allowed / fix when it cannot be;
  4. the session's functions: rows of the same shape, the answer equal to the session's own call bit for bit, its
     Python text (a session opened on the project file) too; one value per argument but the axis.

Run as a script (python3 tests/pyazr/ui/physics_tab_test.py) or by ctest (pyazr_ui_physics_tab).
"""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402


@pytest.fixture(scope="module")
def t1(tmp_path_factory):
    path = copy_project("T1", tmp_path_factory.mktemp("physics_tab"))
    with make_client(path) as client:
        client.path = path
        yield client


def first_request(row, points=50):
    """The request of the tab's first plot: defaults left out, examples given, the axis from its sweep."""
    args = {}
    for a in row["arguments"]:
        if a["name"] == row["axis"]:
            w = a["sweep"]
            args[a["name"]] = ({"from": w["lower"], "to": w["upper"], "spacing": "integer"} if a["kind"] == "integer"
                               else {"from": w["lower"], "to": w["upper"], "points": points,
                                     "spacing": w["spacing"]})
        elif a["default"] is None:
            args[a["name"]] = a["example"]
    return args


def run_python(text, cwd=None):
    """Execute an answer's Python text; its namespace."""
    namespace = {}
    here = os.getcwd()
    try:
        if cwd:
            os.chdir(cwd)
        exec(compile(text, "<copy as python>", "exec"), namespace)
    finally:
        os.chdir(here)
    return namespace


def same(route_value, direct, status):
    """The route's list equals the direct array bit for bit where status is 0, and is "nan" (NaN there) elsewhere."""
    import numpy as np
    direct = np.asarray(direct, float)
    flat_route = np.ravel(np.asarray(route_value, dtype=object))
    flat_status = np.ravel(np.asarray(status))
    for v, d, s in zip(flat_route, np.ravel(direct), flat_status):
        if s == 0:
            if not (isinstance(v, float) and v == d and math.copysign(1, v) == math.copysign(1, d)):
                return False
        elif not ((v == "nan" or (isinstance(v, float) and math.isnan(v))) and math.isnan(d)):
            return False
    return flat_route.size == direct.size


def check_answer(body):
    import numpy as np
    ns = run_python(body["python"])
    result, status = ns["result"], ns["status"]
    assert np.array_equal(np.asarray(status), np.asarray(body["status"])), body["function"]
    fields = result if isinstance(result, dict) else {next(iter(body["returns"])): result}
    for name, direct in fields.items():
        assert same(body["returns"][name], direct, body["status"]), (body["function"], name)


def test_every_row_plots_without_typing_and_its_python_is_the_route(t1):
    rows = t1.get("/api/physics").json()
    for name, row in rows.items():
        r = t1.post(f"/api/physics/{name}", json={"arguments": first_request(row)})
        assert r.status_code == 200, (name, r.text)
        body = r.json()
        codes = sorted({int(c) for c in __import__("numpy").ravel(body["status"])})
        assert codes == [0], (name, body["status_words"])
        assert len(body["arguments"][row["axis"]]) == body["shape"][-1], name
        assert set(body["arguments"]) == {a["name"] for a in row["arguments"]}, name
        check_answer(body)


def test_masked_elements_lines_and_choices(t1):
    import numpy as np
    import pyazr
    r = t1.post("/api/physics/coulomb_waves", json={"arguments": {
        "l": [[0], [2]], "eta": 0.5, "method": [["coulsm"], ["coulcc"]],
        "rho": {"from": -2, "to": 20, "points": 23, "spacing": "linear"}}})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["shape"] == [2, 23]
    assert body["arguments"]["rho"] == np.linspace(-2, 20, 23).tolist()
    status = np.asarray(body["status"])
    assert (status[:, :3] == 3).all() and (status[:, 3:] == 0).all(), status  # rho = -2, -1, 0 outside rho > 0
    assert body["status_words"]["3"] == pyazr.physics.status_text("coulomb_waves", 3)
    assert body["returns"]["F"][0][0] == "nan"
    assert "np.linspace(-2, 20, 23)" in body["python"] and "coulsm, coulcc" in body["python"]
    check_answer(body)
    direct, st = pyazr.physics.coulomb_waves(np.array([[0], [2]]), 0.5, np.linspace(-2, 20, 23),
                                             np.array([[0], [3]]), errors="mask")
    assert same(body["returns"]["G"], direct["G"], st)
    # an integer axis: every whole number from to to
    r = t1.post("/api/physics/statistical_factor", json={"arguments": {
        "two_J": {"from": 0, "to": 6, "spacing": "integer"}, "two_ja": 1, "two_jA": 0}})
    assert r.status_code == 200 and r.json()["arguments"]["two_J"] == [0, 1, 2, 3, 4, 5, 6], r.text
    assert "np.arange(0, 7)" in r.json()["python"]
    check_answer(r.json())


def test_axis_refusals(t1):
    from pyazr.ui.physics import MAX_ELEMENTS
    base = {"l": 0, "eta": 0.5}
    for axis, word in (({"from": -1, "to": 10, "points": 20, "spacing": "log"}, "log"),
                       ({"from": 0.1, "to": 10, "points": MAX_ELEMENTS + 1, "spacing": "linear"}, "points"),
                       ({"from": 0.1, "to": 10, "points": 20, "spacing": "cubic"}, "spacing"),
                       ({"from": 0.1, "points": 20}, "to"), ({"from": 0.5, "to": 3, "spacing": "integer"}, "integer")):
        r = t1.post("/api/physics/coulomb_waves", json={"arguments": dict(base, rho=axis)})
        assert r.status_code == 422, (axis, r.text)
        body = r.json()
        assert body["what"] and body["allowed"] and body["fix"] and word in body["allowed"], (axis, body)
    r = t1.post("/api/physics/coulomb_waves", json={"arguments": dict(
        base, l=[[0], [1], [2]], rho={"from": 0.1, "to": 10, "points": MAX_ELEMENTS // 2, "spacing": "log"})})
    assert r.status_code == 422 and str(MAX_ELEMENTS) in r.json()["allowed"], r.text


def test_session_functions(t1):
    import numpy as np
    import pyazr
    rows = t1.get("/api/session/functions").json()
    physics = t1.get("/api/physics").json()
    row = rows["channel_functions"]
    assert set(row) >= set(next(iter(physics.values()))), "the session's rows have the physics rows' fields"
    assert row["axis"] == "e_cm" and row["group"] == "this project"
    r = t1.post("/api/session/functions/channel_functions", json={"arguments": {
        "pair": 1, "l": 1, "e_cm": {"from": -0.5, "to": 3, "points": 15, "spacing": "linear"}}})
    assert r.status_code == 200, r.text
    body = r.json()
    e = np.linspace(-0.5, 3, 15)
    direct = t1.post("/api/session/channel_functions", json={"pair": 1, "e_cm": e.tolist(), "l": 1}).json()
    for k in ("F", "Fp", "G", "Gp", "P", "S", "delta_hs"):
        assert body["returns"][k] == direct[k], k
    assert body["status"] == [0 if o else -2 for o in direct["open"]] and -2 in body["status"]
    assert repr(t1.path) in body["python"]
    ns = run_python(body["python"], cwd=os.path.dirname(t1.path))
    try:
        for k in ("F", "P", "S", "delta_hs"):
            assert np.array_equal(np.asarray(ns["result"][k]), np.asarray(body["returns"][k], float)), k
    finally:
        ns["session"].close()
    assert pyazr  # the import is the test's requirement
    r = t1.post("/api/session/functions/channel_functions", json={"arguments": {
        "pair": 1, "l": [[0], [1]], "e_cm": {"from": 0.1, "to": 3, "points": 5, "spacing": "linear"}}})
    assert r.status_code == 422 and "one value" in r.json()["allowed"], r.text
    assert t1.post("/api/session/functions/nothing", json={"arguments": {}}).status_code == 404


if __name__ == "__main__":
    uitest.run(__file__)
