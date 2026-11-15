#!/usr/bin/env python3
"""The model as a curve over a data set on the page (GET /api/evaluation/curve/{name}) and an ENDF-6 file's values
on its grid (POST /api/endf/curves), against the pyazr calls made directly:

  1. T2 and T1: every data set's curve over HTTP equals pyazr's curves made as the server makes them at the text's
     x0 (``refine_plan`` with the page's session budget, then one ``CurveSession`` of every data set) bit for bit: grid, c.m. and lab energies,
     angles, model, the data set's own points, the branches; the S factor where offered is the model times the
     engine's conversion (T1's capture), none for an analyzing power;
  2. at the data set's own points the curve equals GET /api/evaluation/data/{name}'s model (T1: to rounding);
  3. the working x moved (PATCH /api/x): the curve follows (the grid kept, or refined again where a seed moved,
     evaluated at that x as the direct calls do);
  4. refused with the three parts: T2's beam-profile data set (422), an unknown data set (404 with the names);
  5. the ENDF curves route: a file's values on the grid of a data set's curve (``on_curve``) equal
     ``Material.curves`` of that curve; the values at the points and the ratios are unchanged.

Run: python3 tests/pyazr/ui/curves_test.py
"""
import os
import shutil
import sys
import warnings

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402

np = pytest.importorskip("numpy")
DATA = os.path.join(uitest.ROOT, "tests", "pyazr", "data", "endf")
HAND = "p-C012_hand.endf"


def direct(path, x=None):
    import pyazr
    warnings.simplefilter("ignore", pyazr.EngineWarning)
    s = pyazr.Session(path)
    return s, page_curves(s, x)


def page_curves(s, x, refined_at=None):
    """The curves as the page's server makes them (``pyazr/ui/engine.py _model_curves``), by pyazr's calls: the grid
    refined at ``refined_at`` (default x) with the page's session budget, refined again at x where a seed moved."""
    from pyazr.model_curve import CurveSession, grid_seeds, plan_curves, refine_plan, seeds_moved
    from pyazr.ui.engine import CURVE_SESSIONS
    first = x if refined_at is None else refined_at
    plan = refine_plan(s, plan_curves(s, None, first), first, max_sessions=CURVE_SESSIONS)
    if seeds_moved(plan, grid_seeds(s, plan, first), grid_seeds(s, plan, x)):
        plan = refine_plan(s, plan_curves(s, None, x), x, max_sessions=CURVE_SESSIONS)
    with CurveSession(s, plan) as c:
        return c.evaluate(x)


def same_curve(answer, c):
    assert answer["variable"] == c.variable and answer["points"] == c.n
    for key, values in (("e_cm", c.energies), ("e_lab", c.energy_lab), ("e_x", c.excitation_energies),
                        ("angle_cm", c.angles), ("model", c.model)):
        assert answer[key] == [float(v) for v in values], key
    assert answer["at_data"] == [bool(v) for v in c.at_data]
    assert [(b["start"], b["stop"], b["fixed"], b["value"]) for b in answer["branches"]] == \
        [(b.start, b.stop, b.fixed, b.value) for b in c.branches]


@pytest.mark.parametrize("name", ["T1", "T2"])
def test_curves_equal_the_direct_call(tmp_path, name):
    path = copy_project(name, tmp_path)
    s, curves = direct(path)
    with make_client(path) as client:
        for n, c in curves.items():
            r = client.get(f"/api/evaluation/curve/{n}")
            if isinstance(c, Exception):
                assert r.status_code == 422, r.text
                body = r.json()
                assert body["what"] == c.what and body["allowed"] == c.allowed and body["fix"] == c.fix
                continue
            assert r.status_code == 200, r.text
            answer = r.json()
            same_curve(answer, c)
            data = client.get(f"/api/evaluation/data/{n}").json()
            k = np.flatnonzero(c.at_data)
            if name == "T1":  # no fold: at its own points the model of the data
                for i, e in enumerate(data["e_cm"]):
                    hit = [j for j in k if abs(answer["e_cm"][j] - e) <= 1e-12 * e
                           and abs(answer["angle_cm"][j] - data["angle_cm"][i]) <= 1e-9 * max(1, data["angle_cm"][i])]
                    assert hit, (n, i)
                    assert abs(answer["model"][hit[0]] - data["model"][i]) <= 1e-12 * abs(data["model"][i]), (n, i)
            if "s_factor" in data:
                assert answer["s_factor"]["model"] == [float(v) for v in c.model * c.sfactor_conversion]
            else:
                assert "s_factor" not in answer
        if name == "T1":
            assert "s_factor" in client.get("/api/evaluation/curve/cap_g0").json()
            assert "s_factor" not in client.get("/api/evaluation/curve/ay_cm090").json()
        else:
            r = client.get("/api/evaluation/curve/s07_profile")
            assert r.status_code == 422 and "window_lo" in r.json()["what"] and r.json()["fix"]
        r = client.get("/api/evaluation/curve/no_such_set")
        assert r.status_code == 404, r.text
    s.close()


def test_curve_follows_the_working_x(tmp_path):
    path = copy_project("T1", tmp_path)
    s, _ = direct(path)
    with make_client(path) as client:
        before = client.get("/api/evaluation/curve/cap_g0").json()
        names = client.get("/api/x").json()
        free = list(names["x"]) if isinstance(names.get("x"), dict) else list(names)
        first = free[0]
        value = client.get("/api/x").json()["x"][first]
        client.patch("/api/x", json={"values": {first: value * 1.01}})
        after = client.get("/api/evaluation/curve/cap_g0").json()
        x = np.array(s.parameters.x0, float)
        x[[p.name for p in s.parameters.free].index(first)] = value * 1.01
        same_curve(after, page_curves(s, x, refined_at=np.array(s.parameters.x0, float))["cap_g0"])
        assert after["model"] != before["model"]
    s.close()


def test_endf_file_on_the_curve(tmp_path):
    path = copy_project("T1", tmp_path)
    shutil.copyfile(os.path.join(DATA, HAND), os.path.join(os.path.dirname(path), HAND))
    import pyazr.endf as E
    s, curves = direct(path)
    m = E.read(os.path.join(os.path.dirname(path), HAND))
    with make_client(path) as client:
        answer = client.post("/api/endf/curves", json={"files": [HAND]}).json()
    seen = 0
    for entry in answer["datasets"]:
        c = curves[entry["name"]]
        if isinstance(c, Exception):
            assert entry["curve"]["refused"]["what"] == c.what
            continue
        assert entry["curve"]["model"] == [float(v) for v in c.model]
        file_entry = entry["curves"][HAND]
        if "refused" in file_entry:
            assert "on_curve" not in file_entry
            continue
        expected = m.curves(c, session=s)
        assert file_entry["on_curve"]["values"] == [float(v) if np.isfinite(v) else None for v in expected.values]
        seen += 1
    assert seen >= 1, "no data set of T1 has the hand file's values on its curve"
    s.close()


if __name__ == "__main__":
    uitest.run(__file__)
