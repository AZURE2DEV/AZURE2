#!/usr/bin/env python3
"""The Nucleus tab's pair and channel edits through the project routes (docs/upgrade/web_ui.md section 6.2), each
against the same pyazr call made directly on a copy of the text:

  1. a pair added (POST /api/project/pairs) and removed (DELETE) on T1: the text equals Project.add_pair's and
     Project.remove_pair's; an unused pair leaves the objective exactly as it was; Undo restores the text byte for
     byte;
  2. removing a pair a level names is refused (422) with the message of Project.remove_pair, the project unchanged;
     with remove_channels the channels go too (the text equals the direct call's) and Undo restores it;
  3. a pair's THM binding energy on T2 (PATCH /api/project/pairs/{key}): the text and the objective equal
     Project.set_pair's and a session's on it;
  4. a channel's "width is an amplitude" flag (PATCH .../channels/{pair}/{l}/{s}): the text equals
     Project.set_channel's; a photon channel's refusal comes with the engine's message;
  5. the channels a J^pi allows (GET /api/project/allowed_channels) equal Project.allowed_channels; a new level with
     them ("allowed") gives the text of Project.add_level.

Run: python3 tests/pyazr/ui/nucleus_routes_test.py
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402


@pytest.fixture()
def t1(tmp_path):
    path = copy_project("T1", tmp_path)
    with make_client(path) as client:
        client.path = path
        yield client


@pytest.fixture()
def t2(tmp_path):
    path = copy_project("T2", tmp_path)
    with make_client(path) as client:
        client.path = path
        yield client


def rev(client):
    return client.get("/api/project").json()["rev"]


def project_of(client):
    """A Project of the text the server holds (a copy beside the project, for the direct call)."""
    import pyazr
    path = client.path + ".direct.azr"
    with open(path, "w") as f:
        f.write(client.get("/api/project/text").text)
    return pyazr.Project.open(path)


def objective(client):
    return client.post("/api/evaluate").json()["objective"]


# d + 11C, a pair of T1's compound nucleus 13N that no data set names (a second p + 12C would make the data sets'
# reactions ambiguous: the engine refuses that)
NEW_PAIR = {"kind": "particle", "light": {"J": 1, "parity": 1, "Z": 1, "mass": 2.014102},
            "heavy": {"J": 1.5, "parity": -1, "Z": 6, "mass": 11.011433}, "separation_energy": 18.68, "radius": 4.5}


def test_pair_added_and_removed(t1):
    text, before = t1.get("/api/project/text").text, objective(t1)
    p = project_of(t1)
    key = p.add_pair(**NEW_PAIR)
    r = t1.post("/api/project/pairs", json=dict(NEW_PAIR, rev=rev(t1)))
    assert r.status_code == 200 and r.json()["result"] == key, r.text
    assert t1.get("/api/project/text").text == p.text
    assert objective(t1) == before  # a pair no level names changes nothing the engine computes
    p.remove_pair(key)
    r = t1.request("DELETE", f"/api/project/pairs/{key}", json={"rev": rev(t1)})
    assert r.status_code == 200, r.text
    assert t1.get("/api/project/text").text == p.text
    for _ in range(2):
        t1.post("/api/project/undo", json={"rev": rev(t1)})
    assert t1.get("/api/project/text").text == text and objective(t1) == before


def test_pair_in_use_refused_then_removed_with_its_channels(t1):
    import pyazr
    text = t1.get("/api/project/text").text
    pairs = t1.get("/api/project/pairs").json()
    used = {c["pair"] for lv in t1.get("/api/project/levels").json() for c in lv["channels"]}
    key = next(p["key"] for p in pairs if p["key"] in used)
    p = project_of(t1)
    with pytest.raises(pyazr.Error) as direct:
        p.remove_pair(key)
    r = t1.request("DELETE", f"/api/project/pairs/{key}", json={"rev": rev(t1)})
    assert r.status_code == 422 and r.json()["message"] == str(direct.value), r.text
    assert t1.get("/api/project/text").text == text
    try:
        p.remove_pair(key, remove_channels=True)
        expected = p.text
    except pyazr.Error as e:  # a level would be left without a channel: the route refuses with the same words
        r = t1.request("DELETE", f"/api/project/pairs/{key}", json={"rev": rev(t1), "remove_channels": True})
        assert r.status_code == 422 and r.json()["message"] == str(e), r.text
        return
    r = t1.request("DELETE", f"/api/project/pairs/{key}", json={"rev": rev(t1), "remove_channels": True})
    assert r.status_code == 200, r.text
    assert t1.get("/api/project/text").text == expected
    t1.post("/api/project/undo", json={"rev": rev(t1)})
    assert t1.get("/api/project/text").text == text


def test_binding_energy_of_a_pair(t2):
    import pyazr
    pair = next(p for p in t2.get("/api/project/pairs").json() if p["binding_energy"] > 0)
    B = pair["binding_energy"] + 0.15
    p = project_of(t2)
    p.set_pair(pair["key"], binding_energy=B)
    r = t2.patch(f"/api/project/pairs/{pair['key']}", json={"binding_energy": B, "rev": rev(t2)})
    assert r.status_code == 200, r.text
    assert t2.get("/api/project/text").text == p.text
    assert next(q for q in t2.get("/api/project/pairs").json() if q["key"] == pair["key"])["binding_energy"] == B
    path = t2.path + ".b.azr"
    p.save(path)
    with pyazr.Session(path) as s:
        direct = float(s.objective()(s.parameters.x0))
    assert objective(t2) == direct
    print(f"\n  T2: B {pair['binding_energy']} -> {B} MeV: objective {direct!r}")


def test_width_is_amplitude_flag(t2):
    import pyazr
    levels = t2.get("/api/project/levels").json()
    pairs = {p["key"]: p for p in t2.get("/api/project/pairs").json()}
    level, c = next((lv, c) for lv in levels for c in lv["channels"] if pairs[c["pair"]]["kind"] == 0)
    p = project_of(t2)
    p.set_channel(level["label"], c["pair"], c["l"], c["s"], width_is_amplitude=not c["amplitude"])
    path = f"/api/project/levels/{level['label']}/channels/{c['pair']}/{c['l']}/{c['s']}"
    r = t2.patch(path, json={"width_is_amplitude": not c["amplitude"], "rev": rev(t2)})
    assert r.status_code == 200, r.text
    assert t2.get("/api/project/text").text == p.text
    photon = next(((lv, ch) for lv in levels for ch in lv["channels"] if pairs[ch["pair"]]["kind"] == 10), None)
    if photon:
        lv, ch = photon
        with pytest.raises((pyazr.Error, ValueError, TypeError)) as direct:
            p.set_channel(lv["label"], ch["pair"], ch["l"], ch["s"], width_is_amplitude=True)
        r = t2.patch(f"/api/project/levels/{lv['label']}/channels/{ch['pair']}/{ch['l']}/{ch['s']}",
                     json={"width_is_amplitude": True, "rev": rev(t2)})
        assert r.status_code == 422 and r.json()["message"] == str(direct.value), r.text


def test_allowed_channels_and_a_new_level(t1):
    p = project_of(t1)
    expected = p.allowed_channels(1.5, -1)
    r = t1.get("/api/project/allowed_channels", params={"J": 1.5, "parity": -1})
    assert r.status_code == 200 and r.json() == expected, r.text
    # a limit given in the query is the call's: 5/2+ of p + 12C needs l = 2, so l_max 1 leaves fewer channels
    narrow = t1.get("/api/project/allowed_channels", params={"J": 2.5, "parity": 1, "l_max": 1}).json()
    assert narrow == p.allowed_channels(2.5, 1, l_max=1) != p.allowed_channels(2.5, 1)
    p.add_level(1.5, -1, 9.5, "allowed")
    r = t1.post("/api/project/levels", json={"J": 1.5, "parity": -1, "energy": 9.5, "channels": "allowed",
                                             "rev": rev(t1)})
    assert r.status_code == 200, r.text
    assert t1.get("/api/project/text").text == p.text
    r = t1.get("/api/project/allowed_channels", params={"J": -1, "parity": 1})
    assert r.status_code == 422 and r.json()["message"], r.text


if __name__ == "__main__":
    uitest.run(__file__)
