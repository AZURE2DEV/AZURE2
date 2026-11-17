#!/usr/bin/env python3
"""The keys of <thm> as the engine's table (src/ThmKeys.cpp): pyazr.kinds()["thm"].

  (a) the two record kinds "<thm>" and "experiment" with the shape of a schema's fields; their keys are the ones the
      engine's refusals name: the unknown experiment key lists the table's keys in order, the refused option lists
      the table's options;
  (b) the web page serves this table as it is (pyazr/ui/thm.py keys(), GET /api/thm/keys): it holds no table of
      its own;
  (c) every key reads a value the table offers: each enum word, and a sample of each other key, is read by the
      engine's reader in an edit of examples/o18_lacognata2010 (Project.set_thm / set_thm_experiment; a refusal by a
      rule between keys, as psNodes without ps, is not the reader's), and a word outside an enum is refused;
  (d) docs/source/reference/thm_keys.rst is what tools/gen_thm_keys_docs.py writes from the tables;
  (e) kinds()["nuclides"] and kinds()["optical_potentials"] are the engine's lists, and each name is read by the
      key that takes it (spectator=, opticalAA=).

Run from anywhere:  python3 tests/pyazr/thm_keys_test.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import pyazr  # noqa: E402
    from pyazr import Project, FormatError  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

THM = os.path.join(ROOT, "examples", "o18_lacognata2010", "o18_lacognata2010.azr")
thm = pyazr.kinds()["thm"]

print("(a) the kinds and the engine's refusals")
check("two kinds, '<thm>' and 'experiment'", list(thm) == ["<thm>", "experiment"], str(list(thm)))
for kind in thm.values():
    for key, f in kind["fields"].items():
        if f["type"] not in ("enum", "number", "integer", "text") or (f["type"] == "enum") != bool(f.get("values")):
            check(f"{key}: a schema type, values exactly for an enum", False, str(f))
p = Project.open(THM)
experiment = next(iter(p.thm()["experiments"]))
try:
    p.set_thm_experiment(experiment, foo="1")
    check("an unknown experiment key refused", False)
except FormatError as e:
    listed = str(e).split("(keys: ")[1].split(")")[0].split(", ")
    check("the unknown key's refusal lists the table's experiment keys, in order",
          listed == list(thm["experiment"]["fields"]), f"{listed}")
try:
    p.set_thm(vertex="nonsense")
    check("a refused option", False)
except FormatError as e:
    text = str(e)
    check("the refused option lists every option of the table",
          all(f"{k}=" in text for k in thm["<thm>"]["fields"]) and "spectatorEnergy[<pair key>]=" in text, text)

print("(b) the web page serves the table")
try:
    from pyazr.ui import thm as ui_thm
except Exception as err:  # the page's extra (fastapi) not installed
    print(f"  skip (b): pyazr.ui not importable ({type(err).__name__}: {err})")
    ui_thm = None
if ui_thm is not None:
    check("pyazr.ui.thm.keys() == kinds()['thm'] (kinds, keys in order, every field of every key)",
          ui_thm.keys() == thm)
    check("the page has no key table of its own", not hasattr(ui_thm, "KEYS"))

print("(c) every key reads what the table offers")
SAMPLES = {"spectatorEnergy": "0.1", "beam": "18O", "target": "d", "spectator": "n", "Ebeam": "54",
           "ps": "hulthen:0-40", "psNodes": "8", "distortion": "coulomb", "opticalAA": "plane", "opticalSF": "plane",
           "spectatorAngle": "qf", "distortionRef": "0.5", "boundState": "whittaker", "theta": "all",
           "spectatorAngles": "0-10", "spectatorAngleNodes": "4", "cbackground": "1-:2",
           "ffWell": "1.2,0.6", "ffNodes": "1/2+:1"}
for kind, setter in (("<thm>", lambda q, k, v: q.set_thm(**{k: v})),
                     ("experiment", lambda q, k, v: q.set_thm_experiment(experiment, **{k: v}))):
    for key, f in thm[kind]["fields"].items():
        values = f.get("values") or [SAMPLES[key]]
        bad = []
        for v in values:
            q = Project.open(THM)
            try:
                setter(q, key, v)
            except FormatError as e:  # a refusal of the reader, not of a rule between keys (CheckThmExperiments)
                if any(w in str(e) for w in (f"{key}='{v}'", "unknown key", "line not understood", f"{key}: ")):
                    bad.append(f"{v}: {str(e)[:160]}")
        check(f"{kind} {key}: {', '.join(values)} read by the engine's reader", not bad, "; ".join(bad))
        if f.get("values"):
            q = Project.open(THM)
            try:
                setter(q, key, "not_a_word")
                check(f"{kind} {key}: a word outside the enum refused", False)
            except FormatError:
                pass

print("(e) the built-in nuclei and the global optical potentials")
from pyazr import _azure2  # noqa: E402  (the binding itself, to compare the public listing with)
kinds = pyazr.kinds()
check("kinds()['nuclides'] is the engine's table (name, Z, A, nuclear mass in u)",
      kinds["nuclides"] == {n: {"Z": z, "A": a, "mass": m} for n, z, a, m in _azure2.thm_nuclides()}
      and kinds["nuclides"]["d"]["Z"] == 1 and kinds["nuclides"]["18O"]["A"] == 18)
check("kinds()['optical_potentials'] names the engine's sets", sorted(kinds["optical_potentials"])
      == sorted(g["name"] for g in _azure2.thm_global_opticals()) and "kd03" in kinds["optical_potentials"])
for name in kinds["optical_potentials"]:
    q = Project.open(THM)
    try:
        q.set_thm_experiment(experiment, opticalAA=name)
        ok = True
    except FormatError as e:
        ok = "opticalAA=" not in str(e)
    check(f"optical potential {name}: read by opticalAA=", ok)
for name in ("d", "3He", "24Mg"):
    q = Project.open(THM)
    try:
        q.set_thm_experiment(experiment, spectator=name)
        ok = True
    except FormatError as e:
        ok = "spectator: " not in str(e)
    check(f"nuclide {name}: read by spectator=", ok)

print("(d) the reference page")
spec = importlib.util.spec_from_file_location("gen_thm_keys_docs", os.path.join(ROOT, "tools", "gen_thm_keys_docs.py"))
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)
with open(gen.PAGE) as f:
    page = f.read()
check("docs/source/reference/thm_keys.rst is the table's (python3 tools/gen_thm_keys_docs.py)",
      page == gen.generate(thm, kinds["nuclides"], kinds["optical_potentials"]))

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
