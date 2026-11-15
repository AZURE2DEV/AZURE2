#!/usr/bin/env python3
"""A first project made from nothing through the page's routes (docs/upgrade/web_ui.md section 11.18): T1 again, as
the page's controls send it, from a folder that holds only T1's data files.

  1. ``python -m pyazr.ui first.azr`` on a file that does not exist: an empty project (``Project.new(path)``), unsaved;
  2. the pairs by naming the nuclei: p + 12C and the two photon pairs of 13N, the charges, masses and the separation
     energy from the mass table (``GET /api/masses``), J, parity and the radius chosen as the form has them;
  3. the levels of T1 with their channels (``POST /api/project/levels`` with the channels the J^pi allows, then each
     channel's width and fixed flag as the level's screen sends them);
  4. the data sets from the folder's files (``GET /api/files``, the first rows ``GET /api/files/{path}``), each added
     with its reaction, quantity, unit, columns, windows and normalization (``POST /api/project/records``);
  5. the calculation (``GET /api/evaluation``): with T1's own numbers typed the data's chi-squared equals T1's bit
     for bit; with the mass table's masses and separation energy (T1's are older ones) to 1e-4 of itself;
  6. the model curve of a data set (``GET /api/evaluation/curve/{name}``), a short fit as a job, accepted, and Save:
     the file written, a session on it gives the objective the page showed.

Prints the steps with the requests sent and the fields typed (the walk's table of section 11.18 is built from it).

Run: python3 tests/pyazr/ui/first_project_test.py
"""
import os
import shutil
import sys
import time

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import make_client  # noqa: E402

T1 = os.path.join(uitest.ROOT, "tests", "upgrade", "T1", "T1.azr")
# T1's objective and data chi-squared (the program's, tests/upgrade/T1)
T1_CHI2 = 482.27523961726735
# What T1's <levels> holds that Project.pairs() does not list: the proton's g-factor, the photon pairs' external
# capture (multipolarities E1, M1, E2: 7)
PROTON_G = 5.5856946893
EXTERNAL_CAPTURE = 7


class Walk:
    """The page's client with a count of what the walk sent and typed."""

    def __init__(self, client):
        self.client, self.steps = client, []

    def step(self, name):
        self.steps.append({"step": name, "requests": 0, "typed": 0})

    def send(self, method, path, typed=0, **kwargs):
        self.steps[-1]["requests"] += 1
        self.steps[-1]["typed"] += typed
        if method in ("POST", "PUT", "PATCH", "DELETE") and "json" in kwargs:
            kwargs["json"] = dict(kwargs["json"], rev=self.client.get("/api/project").json()["rev"])
        r = self.client.request(method, path, **kwargs)
        assert r.status_code == 200, (method, path, r.text)
        return r.json() if r.headers.get("content-type", "").startswith("application/json") else r.text


@pytest.fixture()
def folder(tmp_path):
    """A folder with T1's data files only (what a user brings), no project."""
    shutil.copytree(os.path.join(os.path.dirname(T1), "data"), tmp_path / "data")
    return tmp_path


def walk(client, use_table_masses):
    import pyazr
    t1 = pyazr.Project.open(T1)
    w = Walk(client)

    w.step("the empty project")
    state = w.send("GET", "/api/project")
    assert state["unsaved"] and w.send("GET", "/api/project/levels") == []

    w.step("the pairs by name")
    names = {1: ("p", "12C"), 2: (None, "13N"), 3: (None, "13N")}
    for pair in t1.pairs():
        light, heavy = names[pair.key]
        m = w.send("GET", "/api/masses", typed=2 if light else 1,
                   params={"heavy": heavy, **({"light": light} if light else {})})
        photon = pair.kind_name == "photon"
        # the proton's g-factor is typed (T1's 5.5856946893: no table of the page holds it); Project.pairs() does
        # not list the g-factors nor the external capture, which a user reads from the old file or knows
        light_part = None if photon else {"J": pair.two_j_light / 2, "parity": pair.light_parity, "Z": m["Z1"],
                                          "mass": m["M1"] if use_table_masses else pair.light_mass,
                                          "g": PROTON_G}
        heavy_part = {"J": pair.two_j_heavy / 2, "parity": pair.heavy_parity, "Z": m["Z2"],
                      "mass": m["M2"] if use_table_masses else pair.heavy_mass}
        body = {"kind": pair.kind_name, "light": light_part, "heavy": heavy_part,
                "separation_energy": (m["S"] if use_table_masses else pair.separation_energy) if not photon else 0.0,
                "radius": pair.radius if not photon else None, "heavy_excitation": pair.heavy_excitation,
                "external_capture": EXTERNAL_CAPTURE if photon else None}
        body = {k: v for k, v in body.items() if v is not None}
        # typed: the J and parity of each partner (choices), the radius, the excitation of an excited final state
        w.send("POST", "/api/project/pairs", typed=2 * (1 + (not photon)) + 2 * (not photon) + photon
               + (pair.heavy_excitation > 0), json=body)
    assert [(p["kind"], p["light"], p["heavy"]) for p in w.send("GET", "/api/project/pairs")] == \
        [(p.kind, list(p.light), list(p.heavy)) for p in t1.pairs()]

    w.step("the levels and their channels")
    for level in t1.levels():
        allowed = w.send("GET", "/api/project/allowed_channels",
                         params={"J": level.J, "parity": level.parity})
        wanted = {(c.pair, c.l, c.s) for c in level.channels}
        chosen = [c for c in allowed if (c["pair"], c["l"], c["s"]) in wanted]
        assert len(chosen) == len(wanted), (level.label, allowed)
        w.send("POST", "/api/project/levels", typed=3 + len(chosen),
               json={"J": level.J, "parity": level.parity, "energy": level.energy, "fixed": level.fixed,
                     "channels": [{"pair": c["pair"], "l": c["l"], "s": c["s"], "width": 0, "fixed": True}
                                  for c in chosen]})
        label = f"{level.jpi} {level.energy_text}"
        for c in level.channels:
            path = f"/api/project/levels/{label}/channels/{c.pair}/{c.l}/{c.s}"
            if c.width != 0:
                w.send("PATCH", path, typed=1, json={"width": c.width})
            if not c.fixed:
                w.send("PATCH", path, typed=1, json={"fixed": False})

    w.step("the data sets from the folder's files")
    files = w.send("GET", "/api/files")
    listed = {f["path"] for f in files["files"]} if isinstance(files, dict) else {f["path"] for f in files}
    for record in t1.records().values():
        if record.record_kind != "dataset":
            continue
        assert record.fields["file"] in listed, (record.fields["file"], sorted(listed)[:5])
        w.send("GET", f"/api/files/{record.fields['file']}", params={"rows": 5})
        # typed or chosen: the name and each field (the reaction, quantity, unit, columns, windows, norm)
        w.send("POST", "/api/project/records", typed=1 + len(record.fields),
               json={"kind": "dataset", "name": record.record_name, "fields": dict(record.fields)})

    w.step("the calculation")
    ev = w.send("GET", "/api/evaluation")
    w.step("the curve, a fit, Save")
    curve = w.send("GET", f"/api/evaluation/curve/{ev['datasets'][0]['name']}")
    assert curve["points"] > len(ev["datasets"][0]["e_cm"])
    job = w.send("POST", "/api/jobs", json={"kind": "fit", "method": "engine_migrad", "x0": "text",
                                            "settings": {"max_calls": "30"}})
    for _ in range(600):
        rec = w.send("GET", f"/api/jobs/{job['id']}", params={"wait": 1})
        if rec["state"] not in ("building", "running", "stopping"):
            break
    w.send("POST", f"/api/jobs/{job['id']}/accept", json={})
    after = w.send("GET", "/api/evaluation")
    saved = w.send("POST", "/api/project/save", json={})
    return w, ev, after, saved


def table(w):
    return "\n".join(f"  {s['step']:<40} requests {s['requests']:>3}  typed or chosen {s['typed']:>3}" for s in w.steps)


def test_t1_from_nothing_with_its_own_numbers(folder):
    import pyazr
    path = str(folder / "first.azr")
    with make_client(path) as client:
        w, ev, after, saved = walk(client, use_table_masses=False)
        print("\nthe walk with T1's numbers typed:\n" + table(w))
        # to 1e-8: T1_CHI2 is computed with the Coulomb table; without one (the exact functions: the CI's machines)
        # the same calculation gives 482.2752395734501, 9e-11 away
        assert abs(ev["chi2"] / T1_CHI2 - 1) < 1e-8, ev["chi2"]
        assert after["objective"] <= ev["objective"]
        s = pyazr.Session(os.path.join(str(folder), saved["path"]) if not os.path.isabs(saved["path"])
                          else saved["path"])
        assert s.evaluate().objective == after["objective"]
        s.close()


def test_t1_from_nothing_with_the_tables_masses(folder):
    path = str(folder / "first.azr")
    with make_client(path) as client:
        t0 = time.perf_counter()
        w, ev, after, _ = walk(client, use_table_masses=True)
        print(f"\nthe walk with the mass table's numbers ({time.perf_counter() - t0:.1f} s):\n" + table(w))
        # T1's numbers are not the table's: the proton 1.00782503223 u (AME2020: 1.007825031898), S_p(13N)
        # 1.9434900024 MeV (AME2020: 1.943490064, 62 meV more), 13N 13.00573861 u (13.005738609); near the narrow
        # 1/2+ level the chi-squared moves by 2.2e-5 of itself (482.26440355 against 482.27523962)
        assert abs(ev["chi2"] / T1_CHI2 - 1) < 1e-4, ev["chi2"]


if __name__ == "__main__":
    uitest.run(__file__)
