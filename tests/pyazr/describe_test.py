#!/usr/bin/env python3
"""The package describes itself from the engine's tables (docs/upgrade/pyazr3.md section 2.5).

1. ``pyazr.describe()`` (no path): version, engine build, Coulomb table (path and content identity, or "exact"),
   stopping and mass data, threads; plain data (``json.dumps`` takes it).  ``pyazr.describe(path)`` is still the
   engine's description of a project.
2. ``pyazr.kinds()``: every table non-empty and JSON-clean, and equal to the engine's: the record kinds and the data
   kinds are those of ``AZURE2 --schema``; the parameter kinds those of the engine's table and they cover every kind
   ``AZURE2 --describe --json`` names on tests/upgrade/T1, T2, T3; the options are the option keys of ``<config>`` and
   every other field of the binding's ``RuntimeOptions`` (a field without a row fails here, a row without a reason
   too), each a keyword of ``pyazr.Session`` and each flag one of
   ``AZURE2 --help``; the physics functions are the rows of the binding's table.  A kind added to the engine thus
   appears in ``kinds()`` with no edit of Python.
3. ``Session.describe()`` on tests/upgrade/T1: plain data, its parameters, pairs, levels and data sets those of the
   session and of ``AZURE2 --describe --json``.
4. The session of the package docstring (``help(pyazr)``) runs, in a copy of tests/upgrade/T1.

Needs the compiled engine; without it exits 77 (ctest: Skipped).  The program is AZURE2_BIN (ctest sets it).
Threads: as the caller's environment sets them.

Run from anywhere:  AZURE2_BIN=<build>/src/AZURE2 python3 tests/pyazr/describe_test.py
"""
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


try:
    import pyazr
    from pyazr import _azure2
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)

BIN = os.environ.get("AZURE2_BIN", "")
if not (BIN and os.access(BIN, os.X_OK)):
    print("skip: AZURE2_BIN does not name the program (ctest sets it)")
    sys.exit(77)


def json_clean(label, value):
    try:
        json.dumps(value, allow_nan=True)
        check(f"{label}: json.dumps takes it", True)
    except (TypeError, ValueError) as e:
        check(f"{label}: json.dumps takes it", False, e)


def cli(*args, cwd=None):
    run = subprocess.run([BIN, *args], capture_output=True, text=True, timeout=600, cwd=cwd)
    return run.returncode, run.stdout


# -- 1. the package ---------------------------------------------------------------------------------------------------
d = pyazr.describe()
json_clean("describe()", d)
for key in ("pyazr", "engine", "build", "threads", "coulomb_table", "stopping", "masses", "python"):
    check(f"describe() has {key}", key in d and d[key] not in (None, "", {}), d.get(key))
check("describe(): the pyazr version is the package's", d["pyazr"] == pyazr.__version__)
check("describe(): engine version and module", d["engine"]["version"] not in ("", "unknown")
      and os.path.isfile(d["engine"]["module"]), d["engine"])
check("describe(): build options are booleans", all(isinstance(v, bool) for v in d["build"].values()), d["build"])
check("describe(): the thread count is a positive integer", isinstance(d["threads"]["max_threads"], int)
      and d["threads"]["max_threads"] >= 1, d["threads"])
table = os.environ.get("AZURE2_COULOMB_TABLE", "")
if table and os.path.isfile(table):
    ct = d["coulomb_table"]
    check("describe(): the Coulomb table in use is $AZURE2_COULOMB_TABLE, with its identity",
          ct["in_use"] == "table" and ct["path"] == table and len(ct["identity"]) > 16, ct)
check("describe(): the Coulomb table is 'exact' with the option off",
      pyazr.described.coulomb_table(False)["in_use"] == "exact")
check("describe(): stopping cites SRIM", "SRIM" in d["stopping"]["citation"], d["stopping"])
check("describe(): the mass table", d["masses"]["version"].startswith("AME") and len(d["masses"]["sha256"]) == 64,
      d["masses"])

# -- 2. the kinds -----------------------------------------------------------------------------------------------------
k = pyazr.kinds()
json_clean("kinds()", k)
for key in ("records", "data", "parameters", "options", "physics", "stopping", "errors", "units", "dimensions"):
    check(f"kinds()[{key!r}] is not empty", bool(k.get(key)), type(k.get(key)))
check("kinds() returns a fresh copy", pyazr.kinds() is not k and pyazr.kinds() == k)

rc, out = cli("--schema")
schema = json.loads(out)
check("data kinds == the quantities of AZURE2 --schema", list(k["data"]) == list(schema["quantities"]),
      set(k["data"]) ^ set(schema["quantities"]))
check("record kinds == the records of AZURE2 --schema", list(k["records"]) == list(schema["records"]),
      set(k["records"]) ^ set(schema["records"]))
check("every record kind's fields are the schema's",
      all(set(k["records"][r]["fields"]) == set(schema["records"][r]["fields"]) for r in schema["records"]))
check("every data kind has units and frames", all(q["units"] and q["engine_unit"] for q in k["data"].values()))
check("parameter kinds == the engine's table", list(k["parameters"]) == [r["key"] for r in _azure2.parameter_kinds()])
check("every parameter kind has a unit and owner fields",
      all(p["unit"] and p["owner"] for p in k["parameters"].values()), k["parameters"])
# The program prints no table of parameter kinds; --describe names each data set's parameters (with words of its
# own for the kinds, "energy shift" for "shift": Describe.cpp).  Every parameter it names is a parameter of a session
# on the same project, whose kind (the engine's record) is one of kinds().
for name in ("T1", "T2", "T3"):
    folder = os.path.join(ROOT, "tests", "upgrade", name)
    rc, out = cli("--describe", f"{name}.azr", "--json", cwd=folder)
    described = json.loads(out) if rc == 0 else {}
    check(f"AZURE2 --describe --json {name}", rc == 0 and "datasets" in described, rc)
    named = {p["name"] for ds in described.get("datasets", []) for p in ds.get("parameters", []) if p.get("name")}
    with tempfile.TemporaryDirectory(prefix="pyazr_describe_kinds_") as copy_dir:
        copy = os.path.join(copy_dir, name)
        shutil.copytree(folder, copy, ignore=shutil.ignore_patterns("output", "checks"))
        os.makedirs(os.path.join(copy, "output"), exist_ok=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", pyazr.EngineWarning)
            with pyazr.Session(os.path.join(copy, f"{name}.azr")) as s:
                kind_of = {p.name: p.kind for p in s.parameters}
    check(f"{name}: the parameters --describe names are the session's ({len(named)})",
          named and named <= set(kind_of), sorted(named - set(kind_of))[:5])
    check(f"{name}: every kind of the session's parameters is a kind of kinds()",
          set(kind_of.values()) <= set(k["parameters"]), set(kind_of.values()) - set(k["parameters"]))

# The one table (docs/upgrade/project_format3.md section 3.12): the options of the calculation are the option keys of
# <config> (given to RuntimeOptions.options by name); every other field of RuntimeOptions is a row with its reason.
fields = sorted(f for f in dir(_azure2.RuntimeOptions()) if not f.startswith("_") and f != "options")
check("options: one per field of the binding's RuntimeOptions that is not project state, each with its reason",
      sorted(o["engine_name"] for o in k["options"].values() if not o["key"]) == fields
      and all(o["reason"] for o in k["options"].values() if not o["key"]),
      (fields, sorted(o["engine_name"] for o in k["options"].values() if not o["key"])))
_, config_keys = _azure2.config_keys()
check("options: the others are bool keys of <config>, under their keyword",
      all(o["key"] == kw and (kw, "bool") in [(n, t) for n, t, _ in config_keys]
          for kw, o in k["options"].items() if o["key"]) and sum(1 for o in k["options"].values() if o["key"]) == 10)
init = inspect.signature(pyazr.Session.__init__).parameters
check("options: each is a keyword of pyazr.Session", all(kw in init for kw in k["options"]),
      [kw for kw in k["options"] if kw not in init])
check("options: the defaults are pyazr.Session's",
      all(init[kw].default == o["default"] or (kw in ("coulomb_table", "cwd") and init[kw].default is None
                                               and o["default"] == "")
          for kw, o in k["options"].items() if kw in init),
      {kw: (init[kw].default, o["default"]) for kw, o in k["options"].items() if kw in init})
rc, helptext = cli("--help")
flags = [f.split("=")[0] for o in k["options"].values() for f in (o["flag"], o["negation"]) if f]
check("options: every flag is one of AZURE2 --help", flags and all(re.search(r"\t" + re.escape(f) + r"[\[:=]",
                                                                             helptext) for f in flags),
      [f for f in flags if not re.search(r"\t" + re.escape(f) + r"[\[:=]", helptext)])
check("options: every one has a type and a doc", all(o["type"] and o["doc"] for o in k["options"].values()))

check("physics: the rows of the binding's physics table, each a function of pyazr.physics",
      list(k["physics"]) == [r["name"] for r in _azure2.physics_functions()]
      and set(k["physics"]) == set(pyazr.physics.functions) and "penetrability" in k["physics"],
      list(k["physics"])[:10])
check("physics: every row has units, a domain per argument and an axis among its arguments",
      all(r["axis"] in [a["name"] for a in r["arguments"]] and all(a["domain"] for a in r["arguments"])
          and r["returns"] for r in k["physics"].values()))
check("stopping: the stopping functions are listed", "stopping_power" in k["stopping"] or len(k["stopping"]) > 5,
      list(k["stopping"])[:10])
check("errors: the classes of pyazr.errors", {"Error", "FormatError", "ParameterError", "AZURE2Error",
                                              "EngineError", "DomainError", "SessionClosed"} <= set(k["errors"]),
      list(k["errors"]))

# -- 3. a session -----------------------------------------------------------------------------------------------------
tmp = tempfile.mkdtemp(prefix="pyazr_describe_")
try:
    t1 = os.path.join(tmp, "tests", "upgrade", "T1")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), t1,
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(t1, "output"), exist_ok=True)
    with pyazr.Session(os.path.join(t1, "T1.azr")) as s:
        sd = s.describe()
        json_clean("Session.describe()", sd)
        check("Session.describe(): one entry per parameter", len(sd["parameters"]) == len(s.parameters))
        check("Session.describe(): the parameters are the records of Session.parameters",
              sd["parameters"] == json.loads(json.dumps([dict(p.to_dict(), bounds=[b if abs(b) != float("inf") else None
                                                                                   for b in p.bounds])
                                                         for p in s.parameters], default=float)))
        check("Session.describe(): the free ones start at parameters.x0",
              [p["start"] for p in sd["parameters"] if not p["fixed"]] == [float(v) for v in s.parameters.x0])
        check("Session.describe(): every parameter's unit is its kind's",
              all(p["unit"] == k["parameters"][p["kind"]]["unit"] for p in sd["parameters"]))
        check("Session.describe(): pairs and levels", len(sd["pairs"]) == len(s.pairs)
              and len(sd["levels"]) == len(s.parameters.of_kind("energy")))
        rc, out = cli("--describe", "T1.azr", "--json", cwd=t1)
        names = [x["name"] for x in json.loads(out)["datasets"]]
        check("Session.describe(): the data sets of AZURE2 --describe", [x["name"] for x in sd["datasets"]] == names,
              names)
        check("Session.describe(): data sets have kind, unit, frame and points",
              all(x["quantity"] and x["unit"] and x["frame"] for x in sd["datasets"])
              and all(x["points"] for x in sd["datasets"] if x["active"]), sd["datasets"][:2])
        check("Session.describe(): the options in effect", sd["options"]["use_brune"] is True
              and sd["options"]["cwd"] == s.cwd and sd["mode"] == "data")

    # -- 4. the package docstring's session ------------------------------------------------------------------------
    doc = pyazr.__doc__ or ""
    check("pyazr.__doc__ is about thirty lines", 20 <= len(doc.strip().splitlines()) <= 45,
          len(doc.strip().splitlines()))
    block = doc.split("::", 1)[1] if "::" in doc else ""
    code = textwrap.dedent("\n".join(line for line in block.splitlines() if line.startswith("    ") or not line))
    check("the docstring's session is about ten lines", 8 <= len([c for c in code.splitlines() if c.strip()]) <= 14,
          code)
    run = subprocess.run([sys.executable, "-c", code], cwd=tmp, capture_output=True, text=True, timeout=900,
                         env=dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")))
    check("the docstring's session runs", run.returncode == 0, (run.stdout[-500:], run.stderr[-800:]))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
