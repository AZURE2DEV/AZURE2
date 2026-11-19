#!/usr/bin/env python3
"""The session's diagnostic calls described as rows (pyazr.session_functions, Session.describe()["functions"]).

  (a) Session.describe()["functions"] has a row for channel_functions, jpi_block and the five THM calls, each with the
      shape of a physics row (arguments with kind and unit, returns with unit, path and shape, the axis an argument);
  (b) the rows describe the answers as they are: on T1 (channel_functions, jpi_block of two J groups), on two copies of tests/18O_p_a_thm whose
      THM data sets form an experiment with a Coulomb distortion factor and the line shape, or with a coherent
      background and a window of spectator directions, and on
      examples/o18_lacognata2010 (the DW vertex, with a window of spectator directions), every key of every answer is a
      path of its row, and every path of a row is in some answer; an "axis" field has the length of the energies;
  (c) the answers are not changed by the description (the call's dict is what it was: keys of thm_vertex on o18).

Run from anywhere:  python3 tests/pyazr/session_functions_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import pyazr  # noqa: E402
    from pyazr import Project  # noqa: E402
    from pyazr.session_functions import SESSION_FUNCTIONS  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

ENERGIES = np.linspace(0.2, 1.2, 6)


def paths(answer, prefix="", names=False):
    """Every key path of an answer: a list of dicts is 'key[]', the experiments' names '{experiment}'."""
    out = set()
    for k, v in answer.items():
        here = "{experiment}" if names else f"{prefix}{k}"
        if names:
            out |= paths(v, "{experiment}.")
            continue
        if isinstance(v, list) and v and isinstance(v[0], dict):
            for item in v:
                out |= paths(item, f"{here}[].")
        elif isinstance(v, dict) and k != "cbkg":
            out |= paths(v, f"{here}.")
        else:
            out.add(here)
    return out


def leaf(answer, path):
    """The values of a path in an answer (one per element of its lists)."""
    parts = path.replace("[]", "").split(".")
    values = [answer]
    for p in parts:
        if p == "{experiment}":
            values = [v for d in values for v in d.values()]
            continue
        nxt = []
        for d in values:
            v = d.get(p) if isinstance(d, dict) else None
            if v is None:
                continue
            nxt += v if (isinstance(v, list) and v and isinstance(v[0], dict)) else [v]
        values = nxt
    return values


print("(a) the rows in Session.describe()")
T1 = os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    with pyazr.Session(T1) as s:
        described = s.describe()["functions"]
        answers = {"channel_functions": [s.channel_functions(1, ENERGIES, l=1)],
                   "jpi_block": [s.jpi_block(1, ENERGIES + 1.944), s.jpi_block(3, ENERGIES + 1.944)]}
check("seven rows", sorted(described) == sorted(["channel_functions", "jpi_block", "thm_experiments", "thm_background",
                                                   "thm_lineshape", "thm_vertex", "thm_distortion"]),
      str(sorted(described)))
for name, row in described.items():
    args = [a["name"] for a in row["arguments"]]
    ok = (row["axis"] in args or row["axis"] is None) and all("unit" in a and "kind" in a for a in row["arguments"]) \
        and all({"name", "unit", "doc", "path", "shape"} <= set(r) for r in row["returns"])
    check(f"{name}: arguments with kind and unit, returns with unit, path and shape, axis an argument", ok)

print("(b) the rows describe the answers")
tmp = tempfile.mkdtemp(prefix="session_functions_")
try:
    from thm_records import set_experiment
    target = os.path.join(tmp, "o18thm")
    shutil.copytree(os.path.join(ROOT, "tests", "18O_p_a_thm"), target,
                    ignore=shutil.ignore_patterns("output", "checks", "expected"))
    p = Project.open(os.path.join(target, "18O_p_a_thm.azr"))
    set_experiment(p, "E", [1, 2], beam="18O", target="3He", spectator="d", Ebeam="115", distortion="coulomb",
                   lineshape=True, ps="hulthen:0-40")
    dist = p.save(os.path.join(target, "o18thm.azr"))
    # A coherent background (it excludes the line shape) and a window of spectator directions.
    p = Project.open(os.path.join(target, "18O_p_a_thm.azr"))
    set_experiment(p, "C", [1, 2], beam="18O", target="3He", spectator="d", Ebeam="115", distortion="coulomb",
                   background="linear", cbackground="1/2+:2:0.5,0,0.5,1=0.25,-0.5", spectatorAngles="0-10")
    cbkg = p.save(os.path.join(target, "o18cbkg.azr"))
    o18 = os.path.join(tmp, "o18_2010")
    shutil.copytree(os.path.dirname(os.path.join(ROOT, "examples", "o18_lacognata2010", "x")), o18,
                    ignore=shutil.ignore_patterns("output*", "checks*"))
    q = Project.open(os.path.join(o18, "o18_lacognata2010.azr"))
    q.set_thm_experiment("LC10", spectatorAngles="0-10")  # a window of spectator directions for the DW vertex
    q.save()
    for path, experiment, calls in ((dist, "E", ("thm_experiments", "thm_background", "thm_lineshape", "thm_vertex",
                                                 "thm_distortion")),
                                    (cbkg, "C", ("thm_experiments", "thm_background", "thm_vertex",
                                                 "thm_distortion")),
                                    (os.path.join(o18, "o18_lacognata2010.azr"), "LC10",
                                     ("thm_experiments", "thm_background", "thm_vertex"))):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pyazr.Session(path) as s:
                for fn in calls:
                    call = getattr(s, fn)
                    out = call() if fn == "thm_experiments" else call(experiment) if fn == "thm_background" \
                        else call(experiment, ENERGIES)
                    answers.setdefault(fn, []).append(out)
    for fn, row in SESSION_FUNCTIONS.items():
        described_paths = {r["path"] for r in row["returns"]}
        seen = set()
        for out in answers.get(fn, []):
            got = paths(out, names=(fn == "thm_experiments"))
            seen |= got
            check(f"{fn}: every key of the answer is described", got <= described_paths,
                  str(sorted(got - described_paths)))
            for r in row["returns"]:
                if r["shape"] == "axis":
                    bad = [v for v in leaf(out, r["path"]) if np.size(v) != ENERGIES.size]
                    if bad:
                        check(f"{fn}: {r['path']} has the length of the energies", False, str(bad[:1]))
        check(f"{fn}: every path of the row is in an answer", described_paths <= seen,
              str(sorted(described_paths - seen)))
    print("(c) the answers as they were")
    check("thm_vertex on o18_lacognata2010 keeps its keys",
          set(answers["thm_vertex"][-1]) >= {"E", "channels", "reached", "window", "model", "B", "pair"})
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
