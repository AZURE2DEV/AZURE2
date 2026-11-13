#!/usr/bin/env python3
"""The web UI's physics routes (docs/upgrade/web_ui.md section 4.2):

  1. GET /api/physics is pyazr.kinds()["physics"];
  2. POST /api/physics/{fn} gives the numbers of a direct pyazr.physics call bit for bit (arrays broadcast, scalars,
     several returned fields, a choice by name), NaN with the status and its words outside the domain;
  3. its refusals carry the three parts (what / allowed / fix): an unknown function (404), an unknown or missing
     argument, a wrong type, shapes that do not broadcast, more elements than the bound (422);
  4. POST /api/session/channel_functions gives the session's channel_functions bit for bit, and the same bound.

Run as a script (python3 tests/pyazr/ui/physics_routes_test.py) or by ctest (pyazr_ui_physics_routes).
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402


@pytest.fixture(scope="module")
def t1(tmp_path_factory):
    path = copy_project("T1", tmp_path_factory.mktemp("physics"))
    with make_client(path) as client:
        client.path = path
        yield client


def test_rows_are_the_kinds(t1):
    import pyazr
    assert t1.get("/api/physics").json() == pyazr.kinds()["physics"]


def test_values_are_the_direct_call(t1):
    import numpy as np
    import pyazr
    e = [0.1, 0.5, 1.0, 2.5, 7.0]
    r = t1.post("/api/physics/penetrability", json={"arguments": {"l": [[0], [1], [2]], "z1": 1, "z2": 6,
                                                                   "m1_u": 1.00728, "m2_u": 12.0, "radius": 5.0,
                                                                   "e_cm": e}})
    assert r.status_code == 200, r.text
    body = r.json()
    direct = pyazr.physics.penetrability(np.array([[0], [1], [2]]), 1, 6, 1.00728, 12.0, 5.0, e)
    assert body["shape"] == [3, 5] and body["axis"] == "e_cm" and body["units"] == {"P": ""}
    assert np.array_equal(np.asarray(body["returns"]["P"], float), direct)
    assert body["status"] == [[0] * 5] * 3 and body["status_words"] == {"0": "computed"}
    # several fields, a choice by name, scalars stay scalars
    r = t1.post("/api/physics/coulomb_waves", json={"arguments": {"l": 1, "eta": 0.5, "rho": 3.0,
                                                                   "method": "coulcc"}})
    w = pyazr.physics.coulomb_waves(1, 0.5, 3.0, method="coulcc")
    assert r.status_code == 200, r.text
    assert r.json()["returns"] == w and r.json()["status"] == 0 and r.json()["shape"] == []


def test_outside_the_domain_is_nan_with_status(t1):
    import pyazr
    r = t1.post("/api/physics/penetrability", json={"arguments": {"l": 0, "z1": 1, "z2": 6, "m1_u": 1.00728,
                                                                   "m2_u": 12.0, "radius": 5.0,
                                                                   "e_cm": [1.0, -1.0]}})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["returns"]["P"][0] == pyazr.physics.penetrability(0, 1, 6, 1.00728, 12.0, 5.0, 1.0)
    assert body["returns"]["P"][1] == "nan" and body["status"] == [0, 7]
    assert body["status_words"]["7"] == "e_cm outside its domain (finite > 0)"


def _refused(r, status, *words):
    assert r.status_code == status, r.text
    body = r.json()
    assert body["what"] and body["allowed"] and body["fix"], body
    for w in words:
        assert w in body["message"] + body["allowed"] + body["fix"], (w, body)
    return body


def test_refusals(t1):
    args = {"l": 0, "z1": 1, "z2": 6, "m1_u": 1.00728, "m2_u": 12.0, "radius": 5.0}
    body = _refused(t1.post("/api/physics/no_such_function", json={"arguments": {}}), 404, "penetrability")
    assert "penetrability" in body["known"]
    _refused(t1.post("/api/physics/penetrability", json={"arguments": args}), 422, "missing", "e_cm")
    _refused(t1.post("/api/physics/penetrability", json={"arguments": dict(args, e_cm=1.0, energy=1.0)}), 422,
             "unknown", "energy")
    _refused(t1.post("/api/physics/penetrability", json={"arguments": dict(args, e_cm="one")}), 422, "e_cm",
             "a number or a (nested) list")
    _refused(t1.post("/api/physics/penetrability", json={"arguments": dict(args, l=[0, 1], e_cm=[1.0, 2.0, 3.0])}),
             422, "broadcast")
    _refused(t1.post("/api/physics/coulomb_waves", json={"arguments": {"l": 0, "eta": 1.0, "rho": 1.0,
                                                                      "method": "gsl"}}), 422, "coulsm")


def test_the_bound(t1):
    from pyazr.ui.physics import MAX_ELEMENTS
    args = {"l": 0, "z1": 1, "z2": 6, "m1_u": 1.00728, "m2_u": 12.0, "radius": 5.0}
    r = t1.post("/api/physics/sommerfeld_eta", json={"arguments": {"z1": 1, "z2": 6, "m1_u": 1.0, "m2_u": 12.0,
                                                                    "e_cm": [1.0] * MAX_ELEMENTS}})
    assert r.status_code == 200 and len(r.json()["returns"]["eta"]) == MAX_ELEMENTS
    _refused(t1.post("/api/physics/penetrability", json={"arguments": dict(args, e_cm=[1.0] * (MAX_ELEMENTS + 1))}),
             422, str(MAX_ELEMENTS))
    # refused in the session worker: the class and the message travel (the three parts are in the message)
    r = t1.post("/api/session/channel_functions", json={"pair": 1, "e_cm": [1.0] * (MAX_ELEMENTS + 1)})
    assert r.status_code == 422 and r.json()["error"] == "Refused", r.text
    assert f"above the bound of {MAX_ELEMENTS}" in r.json()["message"] and "split" in r.json()["message"], r.text


def test_session_channel_functions_are_the_sessions(t1):
    import numpy as np
    import pyazr
    e = [0.5, -0.2, 1.0, 2.0]
    r = t1.post("/api/session/channel_functions", json={"pair": 1, "e_cm": e, "l": 1})
    assert r.status_code == 200, r.text
    s = pyazr.Session(os.path.join(os.path.dirname(t1.path), os.path.basename(t1.path)))
    try:
        direct = s.channel_functions(1, np.asarray(e), 1)
    finally:
        s.close()
    body = r.json()
    for k in ("F", "Fp", "G", "Gp", "P", "S", "delta_hs"):
        assert np.array_equal(np.asarray(body[k], float), direct[k]), k
    assert body["open"] == direct["open"].tolist() and body["e_cm"] == e


if __name__ == "__main__":
    uitest.run(__file__)
