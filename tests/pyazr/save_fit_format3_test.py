#!/usr/bin/env python3
"""save_fit and the hybrid potential on format-3 projects (the pyazr debts of HANDOFF 2026-11-10):

  1. save_fit(norms="fitted") of a format-3 project (T2's conversion, with one data set's shift in keV): the free
     norms (a data set's own, a shared normalization), shift terms (an offset in keV, a sqrt term, a shared energy
     scale's offset) and the shift function's parameter are written into the records' nuisances in their units, the
     priors untouched; verify=True reopens the snapshot and finds the fit's chi^2 and priors (it raised before);
     with a free abundance and a fitted areal density moved off their prior centres too: the snapshot keeps the
     original centres and pays the same priors (the engine took the start value as the centre until 2026-11-13, and
     verify raised), and Parameter.prior_centre is the centre written;
  2. a record's matching_radius and step (fm), stated in the file or set by Project.set: the view reports them, the
     channel functions of the edited project equal those of the record stating them, bit for bit, and save_fit keeps
     them in the format-3 record.

Needs the compiled module.  Run from anywhere:  python3 tests/pyazr/save_fit_format3_test.py
"""
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T2 = os.path.join(ROOT, "tests", "upgrade", "T2")
HYBRID = os.path.join(ROOT, "tests", "hybrid_potential")
FIXTURES = os.path.join(ROOT, "convert", "tests", "fixtures")
ORIGINALS = os.path.join(ROOT, "tests", "convert", "fixtures", "tests")  # the originals of converted projects

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np  # noqa: E402
from pyazr import Session, convert_azr, Project  # noqa: E402


def functions(A, pair=1):
    out = []
    for L in (0, 1, 2):
        cf = A.channel_functions(pair, [0.5, 1.0, 2.0], l=L)
        out += [float(v) for v in cf["P"]] + [float(v) for v in cf["S"]]
    return out


def main():
    work = tempfile.mkdtemp(prefix="save_fit_f3.")
    try:
        fitted(work)
        potentials(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(failures)} failure(s)" if failures else "PASS")
    return 1 if failures else 0


def fitted(work):
    print("1. save_fit(norms='fitted') of a format-3 project")
    d = os.path.join(work, "T2")
    os.makedirs(os.path.join(d, "output"))
    os.symlink(os.path.join(T2, "data"), os.path.join(d, "data"))
    shutil.copy(os.path.join(T2, "cov_diff.txt"), d)
    shutil.copy(os.path.join(FIXTURES, "T2.prior_cov.f3.txt"), os.path.join(d, "prior_cov.f3.txt"))
    shutil.copy(os.path.join(FIXTURES, "T2.f3.azr"), os.path.join(d, "T2.azr"))
    p = Project.open(os.path.join(d, "T2.azr"))
    p["s01_int_thin"].set(shift="offset 0 keV fit prior 0 3")  # the offset in keV: written back in keV
    p.save()
    with Session(p.path, cwd=d) as m:
        x = np.array(m.parameters.x0, float)
        by_name = {q.name: q for q in m.parameters}
        moved = {"s01_int_thin.norm": 1.0234, "D.value": 0.9871, "s01_int_thin.shift.offset": 0.0012,
                 "G.offset": -0.0007, "s04_thick.shift.sqrt": 0.0013, "G.a0": 0.0042,
                 "natN.isotope.15N": 0.0036, "tgt_s05_implanted.areal_density": 2.1e17}
        iso = by_name["natN.isotope.15N"]
        check("the abundance's prior_centre is the centre written (0.0037)", iso.prior_centre == 0.0037, iso.prior_centre)
        for name, v in moved.items():
            x[by_name[name].free_index] = v
        out = os.path.join(d, "T2_fit.azr")
        try:
            m.save_fit(out, x, norms="fitted", param_sav=False)
            check("save_fit(norms='fitted', verify=True) of a format-3 project: written, read back as the fit", True)
        except Exception as e:  # noqa: BLE001 -- what failed before the fix
            check("save_fit(norms='fitted', verify=True) of a format-3 project: written, read back as the fit", False,
                  f"{type(e).__name__}: {e}")
            return
    q = Project.open(out)
    got = {"s01_int_thin.norm": q["s01_int_thin"].norm, "D.value": q["D"].value,
           "s01_int_thin.shift.offset": q["s01_int_thin"].shift, "G.offset": q["G"].offset,
           "s04_thick.shift.sqrt": [s for s in q["s04_thick"].get("shift") if s.startswith("sqrt")]
           if isinstance(q["s04_thick"].get("shift"), list) else q["s04_thick"].shift, "G.a0": q["G"].param,
           "natN.isotope.15N": q["natN"].isotope, "tgt_s05_implanted.areal_density": q["tgt_s05_implanted"].areal_density}
    want = {"s01_int_thin.norm": "1.0234 fit prior 1 5%", "D.value": "0.9871 fit prior 1 8%",
            "s01_int_thin.shift.offset": "offset 1.2 keV fit prior 0 3", "G.offset": "-0.0007 MeV fit prior 0 0.003",
            "G.a0": "a0 0.0042 fit prior 0 0.01", "natN.isotope.15N": "15N 0.0036 fit prior 0.0037 0.0002",
            "tgt_s05_implanted.areal_density": "2.1e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17"}
    def same(a, b):  # equal words, numbers equal as numbers (the writer's shortest form: -7e-04 for -0.0007)
        a, b = a.split(), b.split()
        return len(a) == len(b) and all(x == y or (re.match(r"^[-+.\de]+$", x) and float(x) == float(y))
                                        for x, y in zip(a, b))
    for k, w in want.items():
        check(f"{k}: written as '{w}' (its unit, prior untouched)", same(got[k], w), got[k])
    check("s04_thick.shift.sqrt: written as 'sqrt 0.0013 MeV^1/2 fit ...'", "sqrt 0.0013 MeV^1/2 fit" in str(got["s04_thick.shift.sqrt"]),
          got["s04_thick.shift.sqrt"])


def potentials(work):
    print("2. matching radius and step")
    d = os.path.join(work, "hybrid")
    os.makedirs(os.path.join(d, "output"))
    shutil.copytree(os.path.join(HYBRID, "data"), os.path.join(d, "data"))
    old = os.path.join(d, "hybrid_potential.azr")
    shutil.copy(os.path.join(ORIGINALS, "hybrid_potential", "hybrid_potential.azr"), old)  # the original
    f3 = os.path.join(d, "f3.azr")
    convert_azr(old, f3, to=3)
    text = open(f3).read()
    stated = os.path.join(d, "stated.azr")
    with open(stated, "w") as f:
        f.write(re.sub(r"(step +)0\.0275 fm", r"\g<1>0.02 fm", re.sub(r"(matching_radius +)15 fm", r"\g<1>12 fm", text)))
    with Session(stated, cwd=d) as s:
        want = functions(s)
        view = s.nuclear_potential(1)
    check("the record's matching radius and step in the view (12 fm, 0.02 fm)",
          view.matching_radius == 12.0 and view.step == 0.02, repr((view.matching_radius, view.step)))
    with Session(f3, cwd=d) as m:
        base = functions(m)
    edited = Project.open(f3)
    edited.set("p12C", matching_radius="12 fm", step="0.02 fm")
    with Session(edited.save(os.path.join(d, "edited.azr")), cwd=d) as m:
        p = m.nuclear_potential(1)
        got = functions(m)
        check("Project.set(p12C, matching_radius=12 fm, step=0.02 fm): the view", p.matching_radius == 12.0
              and p.step == 0.02, repr(p))
        check("  the channel functions equal those of the record stating them, bit for bit (and differ from 15 fm, "
              "0.0275 fm)", got == want and got != base, max(abs(a - b) for a, b in zip(got, want)))
        m.save_fit(os.path.join(d, "written.azr"), verify=False, norms="nominal")
    written = open(os.path.join(d, "written.azr")).read()
    check("save_fit keeps the matching radius and step in the record",
          re.search(r"matching_radius +12 fm", written) and re.search(r"step +0\.02 fm", written), written[-600:])
    with Session(os.path.join(d, "written.azr"), cwd=d) as w:
        check("  and the written project gives the same channel functions", functions(w) == want)


if __name__ == "__main__":
    sys.exit(main())
