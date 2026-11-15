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
     them ("allowed") gives the text of Project.add_level;
  6. the lists in words on T2: the pairs carry their partners by name, J^pi and kind in words (the engine's label
     split at its "+"), every field of a pair has its words in pyazr.kinds()["pairs"]; the references
     (GET /api/project/references) equal Project.references; a record the engine does not describe (a target) is
     refused by GET /api/project/describe/{name} with the records that use it, not answered with another entry;
  7. a prior and limits on a channel's width (PUT /api/project/parameter_settings): the text equals
     Project.set_parameter_setting's, the rows Project.parameter_settings', the prior in the objective only, a fixed
     parameter refused with the three parts, Undo.

Run: python3 tests/pyazr/ui/nucleus_routes_test.py
"""
import os
import sys

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

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


def test_parameter_settings_prior_and_limits(t1):
    """A prior and limits on a level parameter (PUT /api/project/parameter_settings): the text equals
    Project.set_parameter_setting's, GET lists the rows as Project.parameter_settings, the prior's term is in the
    objective and not in the data's chi2, a fixed parameter is refused with the three parts, Undo restores the text."""
    before_text = t1.get("/api/project/text").text
    before = t1.get("/api/evaluation").json()
    level, channel = "1/2+ 2.365", [1, 0, 0.5]
    width = 33000.0
    r = t1.put("/api/project/parameter_settings", json={"rev": rev(t1), "level": level, "channel": channel,
                                                         "prior": [width * 1.1, width * 0.05],
                                                         "limits": [1000.0, 1e5]})
    assert r.status_code == 200, r.text
    direct = project_of_text(t1, before_text)
    direct.set_parameter_setting(level, tuple(channel), prior=(width * 1.1, width * 0.05), limits=(1000.0, 1e5))
    assert t1.get("/api/project/text").text == direct.text
    rows = t1.get("/api/project/parameter_settings").json()
    assert rows == [dataclass_plain(row) for row in direct.parameter_settings()], rows
    assert rows[0]["level"] == level and rows[0]["channel"] == channel and rows[0]["nuisance"], rows
    after = t1.get("/api/evaluation").json()
    assert after["chi2"] == before["chi2"] and after["objective"] > before["objective"], (before, after)
    fixed = t1.put("/api/project/parameter_settings", json={"rev": rev(t1), "level": "1/2- 0.0",
                                                             "prior": [0.0, 0.01]})
    assert fixed.status_code == 422 and all(fixed.json()[k] for k in ("what", "allowed", "fix")), fixed.text
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    assert t1.get("/api/project/text").text == before_text
    assert t1.get("/api/project/parameter_settings").json() == []


def test_masses_by_name_and_an_empty_project(tmp_path):
    """A new .azr in an existing folder: the server starts an empty project there (Project.new(path), unsaved, the
    output and checks directories made); GET /api/masses gives pyazr.sources.masses.pair_masses by the nuclides'
    names (p + 12C: S_p(13N)), the heavy nucleus alone without light, an unknown nuclide refused (404)."""
    from pyazr.sources import masses
    path = str(tmp_path / "first.azr")
    with make_client(path) as client:
        state = client.get("/api/project").json()
        assert state["unsaved"] and client.get("/api/project/levels").json() == [], state
        assert client.get("/api/project/pairs").json() == []
        # no record yet: a list (the page's lists filter it), not the empty mapping Project.records() is
        assert client.get("/api/project/records").json() == []
        assert os.path.isdir(tmp_path / "output") and os.path.isdir(tmp_path / "checks")
        assert not os.path.exists(path)
        m = client.get("/api/masses", params={"light": "p", "heavy": "12C"}).json()
        assert {k: m[k] for k in ("Z1", "M1", "Z2", "M2", "S")} == {
            k: v for k, v in masses.pair_masses("p", "12C")._asdict().items() if k != "compound"}, m
        assert m["compound"] == [7, 13] and m["table"] == "AME2020"
        alone = client.get("/api/masses", params={"heavy": "13N"}).json()
        assert alone["Z2"] == 7 and alone["M2"] == masses.atomic_mass("13N") and "S" not in alone
        assert client.get("/api/masses", params={"heavy": "300Xx"}).status_code in (404, 422)
        assert client.post("/api/project/save", json={"rev": client.get("/api/project").json()["rev"]}) \
            .status_code == 200
        assert os.path.isfile(path)


def project_of_text(client, text):
    import pyazr
    path = client.path + ".direct.azr"
    with open(path, "w") as f:
        f.write(text)
    return pyazr.Project.open(path)


def dataclass_plain(value):
    import dataclasses
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in dataclasses.asdict(value).items()}


def test_lists_in_words_pairs_references_description(t2):
    import pyazr
    pairs = t2.get("/api/project/pairs").json()
    assert [(p["light_name"], p["heavy_name"], p["light_jpi"], p["heavy_jpi"], p["kind_name"]) for p in pairs][:2] == [
        ("p", "15N", "1/2+", "1/2-", "particle"), ("a", "12C", "0+", "0+", "particle")]
    words = pyazr.kinds()["pairs"]
    assert set(pairs[0]) == set(words["fields"]), "every field of a pair has its words"
    assert all(f["title"] and "_" not in f["title"] for f in words["fields"].values())
    assert set(words["partner_fields"]) == {"J", "parity", "Z", "mass", "g"}
    direct = pyazr.Project.open(t2.path)
    references = t2.get("/api/project/references").json()
    assert references == direct.references()
    assert {"record": "s01_int_thin", "field": "target", "name": "tgt_s01_int_thin", "epoch": ""} in references
    refused = t2.get("/api/project/describe/tgt_s01_int_thin")
    assert refused.status_code >= 400 and "s01_int_thin.target" in refused.text, refused.text
    described = t2.get("/api/project/describe/s01_int_thin").json()
    assert described["name"] == "s01_int_thin"
    assert set(pyazr.kinds()["description"]["fields"]) <= set(described) | {"angle_window_deg"}
