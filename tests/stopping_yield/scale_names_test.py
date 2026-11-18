#!/usr/bin/env python3
"""The names of the stopping scales (docs/source/reference/data_formats.rst, the stopping field): one parameter per name.

A target that writes ``phase solid`` and one of the same composition that leaves the phase out have two scale keys
(experiment/src/BuildEffects.cpp ScaleOf), and both scales were named ``stopping_scale.<formula>.<projectile>``.
Minuit2 refuses a second parameter of a name with a warning only, so the vector held one parameter fewer than its
records and every session failed with "N parameters but N+1 records" (2026-11-18: a 10B target with ``phase solid``
beside another 10B target; the 13C(p,n) data set of Ramstroem 1979 on a C target beside the C targets of Kashy 1961,
which write ``phase solid``).  This test, on models/stopping_yield.azr (three TiN layers, protons):

  1. TS with ``phase solid``: the evaluations' scales are stopping_scale.TiN.1H and stopping_scale.TiN_solid.1H
     (the table's is stopping_scale.TT.1H); ``AZURE2 calculate`` runs and its model values are those of the project
     as shipped (solid is TiN's default phase: the same stopping); with pyazr, Session(...).objective() is made and
     evaluates to the session's chi^2;
  2. the table target TT renamed TiN (its scale's label is its name): two scales of one name, refused by --validate
     and by the run with E_STOPPING_SCALE_CONFLICT naming stopping_scale.TiN.1H.

    AZURE2_BIN=<build>/src/AZURE2 [PYAZR_ENGINE_DIR=<build>/pyazr-out] python3 tests/stopping_yield/scale_names_test.py
Exit 77 without AZURE2_BIN, and without the compiled module once the CLI's checks have passed.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
AZURE2 = os.environ.get("AZURE2_BIN")
if not AZURE2 or not os.access(AZURE2, os.X_OK):
    print("skip: AZURE2_BIN is not an executable")
    sys.exit(77)

failures = []


def ok(name, good, detail=""):
    print(f"  {'ok  ' if good else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not good:
        failures.append(name)


def project(tmp, name, *edits):
    """A copy of models/stopping_yield.azr with its data in tmp/name, each (old, new) replaced once."""
    work = os.path.join(tmp, name)
    shutil.copytree(os.path.join(HERE, "models", "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    text = open(os.path.join(HERE, "models", "stopping_yield.azr")).read()
    for old, new in edits:
        assert old in text, old
        text = text.replace(old, new, 1)
    open(os.path.join(work, "m.azr"), "w").write(text)
    return work


def cli(work, *args):
    r = subprocess.run([AZURE2, *args, "m.azr"], cwd=work, capture_output=True, text=True, timeout=600)
    return r.returncode, r.stdout + r.stderr


def model_values(work):
    """Every AZUREOut file's columns but the last (the data's error) of the output directory, by file name."""
    out = {}
    for f in sorted(glob.glob(os.path.join(work, "output", "AZUREOut_*.out"))):
        out[os.path.basename(f)] = [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
    return out


TS = ("target TS\n    kind           layer\n    active         15N\n    composition    Ti 1 N 1\n",
      "target TS\n    kind           layer\n    active         15N\n    composition    Ti 1 N 1\n"
      "    phase          solid\n")

with tempfile.TemporaryDirectory() as tmp:
    print("1. a target with phase solid beside one without a phase")
    base = project(tmp, "base")
    solid = project(tmp, "solid", TS)
    code, log = cli(solid, "--validate")
    ok("--validate: the project is valid", code == 0 and "valid." in log, log.strip().splitlines()[-1:])
    code_b, log_b = cli(base, "calculate", "--no-gui")
    code_s, log_s = cli(solid, "calculate", "--no-gui")
    ok("the calculation runs", code_s == 0 and "records (a parameter added" not in log_s,
       [l for l in log_s.splitlines() if "rror" in l][:2])
    vb, vs = model_values(base), model_values(solid)
    ok("its model values are those of the phase left out (TiN's default phase is solid)",
       code_b == 0 and vb and vb == vs, f"{len(vs)} output files")
    try:
        import numpy as np
        from pyazr import Session
    except Exception as e:  # pragma: no cover
        print(f"  skip  pyazr not importable ({e})")
        Session = None
    else:
        with Session(os.path.join(solid, "m.azr"), cwd=solid) as s:
            params = list(s.parameters)
            names = [q.name for q in params]
            # what a parameter is, from its record (check 23: never from its name)
            scales = sorted(q.name for q in params if q.kind == "apparatus" and q.apparatus == "stopping_scale")
            ok("three scales: stopping_scale.TiN.1H (TP), stopping_scale.TiN_solid.1H (TS), stopping_scale.TT.1H (the table)",
               scales == ["stopping_scale.TT.1H", "stopping_scale.TiN.1H", "stopping_scale.TiN_solid.1H"], scales)
            ok("each parameter name once", len(names) == len(set(names)), f"{len(names)} names")
            f = s.objective()
            x = np.asarray(s.parameters.x0, float)
            value, chi2 = float(f(x)), s.evaluate(x).chi2
            ok("Session(...).objective() evaluates: the chi^2 (no prior on a fixed scale)",
               np.isfinite(value) and np.isclose(value, chi2, rtol=1e-12), f"{value} vs {chi2}")

    print("2. two scales of one name")
    clash = project(tmp, "clash", ("target TT\n", "target TiN\n"), ("    target         TT\n", "    target         TiN\n"))
    for what in ("--validate", "calculate"):
        code, log = cli(clash, *([what] if what == "--validate" else [what, "--no-gui"]))
        ok(f"{what}: refused, E_STOPPING_SCALE_CONFLICT naming stopping_scale.TiN.1H",
           code != 0 and "E_STOPPING_SCALE_CONFLICT" in log and "stopping_scale.TiN.1H" in log,
           [l for l in log.splitlines() if "E_" in l or "rror" in l][:2])

print()
if not failures and Session is None:
    print("skip: the CLI's checks passed; pyazr's (Session(...).objective()) not run")
    sys.exit(77)
if failures:
    print(f"scale_names: {len(failures)} failure(s)")
    sys.exit(1)
print("scale_names: all checks passed")
