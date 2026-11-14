#!/usr/bin/env python3
"""A thick target with a positive energy shift: the translated stopping stays inside its domain.

A shifted folded point reads the stopping at E_n - Delta at its nodes (include/ShiftedFold.h, question 3).  For a
thick target the nominal support ends at the 1 keV floor; with a positive shift of more than 1 keV inside the reach the
enclosing node below the moved support lies at the floor, and E_n - Delta < 0 there, where a power-law stopping formula
has no real value (src/Equation.cpp printed "Exponent results in unsupported imaginary number" for every such call,
thousands per fit of T2).  The translated energy is held at the nominal support's lower edge.

Model: tests/upgrade/T2 with only s04_thick active (15N(p,a) thick-target yield, stopping a0 x^-0.45), its sqrt(E)
shift set to 0.002 MeV^1/2 (about 2 keV at 1 MeV) in two forms: (A) varied with prior 0.003, so the mesh is planned
around 0 with a reach and the moved support's lower edge falls inside an interval above the floor node; (B) fixed,
so the mesh is planned for the moved support itself (its lower edge a hard node, no node below it).  Checked:
  (a) run A prints no warning of the stopping formula (17 on the engine before the fix, one per point whose shift
      exceeds 1 keV);
  (b) A's yields equal B's to 2e-4 relative (the two meshes differ: the quadrature difference; measured 2.2e-5);
  (c) the shift is active: A differs from the unshifted run.

Run from anywhere:  AZURE2_BIN=path/to/AZURE2 python3 tests/pyazr/shift_fold_floor_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    from pyazr import Project                  # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)
BIN = os.environ.get("AZURE2_BIN", "")
if not BIN or not os.access(BIN, os.X_OK):
    print("skip: AZURE2_BIN not set")
    sys.exit(77)

T2 = os.path.join(ROOT, "tests", "upgrade", "T2")


def model(work, shift):
    p = Project.open(os.path.join(T2, "T2.azr"))
    for name, r in p.records().items():
        if r.record_kind == "dataset" and name != "s04_thick":
            p.set(name, active="no")
    p.set("s04_thick", shift=shift)
    path = p.save(os.path.join(work, "model.azr"))
    # the prior covariance names parameters of the data sets switched off
    with open(path) as f:
        text = f.read()
    edited = text.replace("priorCovariance=prior_cov.f3.txt\n", "")
    if edited == text:
        raise SystemExit("the model's priorCovariance line was not found")
    with open(path, "w") as f:
        f.write(edited)
    return path


def run(label, shift):
    work = os.path.join(TMP, label)
    shutil.copytree(T2, work)
    for d in ("output", "checks"):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
        os.makedirs(os.path.join(work, d))
    azr = model(work, shift)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"))
    out = subprocess.run([BIN, "calculate", "--no-gui", azr], stdin=subprocess.DEVNULL, cwd=work, env=env,
                         capture_output=True, text=True)
    yields = []
    with open(os.path.join(work, "output", "AZUREOut_aa=1_R=2.out")) as f:
        for line in f:
            t = line.split()
            if len(t) >= 4 and not line.lstrip().startswith("#"):
                try:
                    yields.append(float(t[3]))
                except ValueError:
                    pass
    return out.returncode, out.stdout + out.stderr, yields


TMP = tempfile.mkdtemp(prefix="pyazr_shift_floor_")
try:
    rcA, logA, yA = run("varied", "sqrt 0.002 MeV^1/2 fit prior 0 0.003")
    rcB, logB, yB = run("fixed", "sqrt 0.002 MeV^1/2 fixed")
    rc0, log0, y0 = run("unshifted", "sqrt 0 MeV^1/2 fixed")
    check("the three runs end normally", rcA == 0 and rcB == 0 and rc0 == 0 and len(yA) > 0, (rcA, rcB, rc0, len(yA)))
    warnings = logA.count("unsupported imaginary")
    check("(a) no warning of the stopping formula in the varied run", warnings == 0, f"{warnings} warnings")
    same = len(yA) == len(yB) and len(yA) > 0
    rel = max(abs(a - b) / abs(b) for a, b in zip(yA, yB) if b != 0) if same else float("inf")
    print(f"    largest relative difference varied / fixed: {rel:.3e} over {len(yA)} yields")
    check("(b) varied = fixed shift to the quadrature difference (2e-4)", same and rel < 2.0e-4, rel)
    moved = max(abs(a - b) / abs(b) for a, b in zip(yA, y0) if b != 0) if len(yA) == len(y0) else 0.0
    check("(c) the shift moves the yields", moved > 1.0e-3, moved)
finally:
    shutil.rmtree(TMP, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
