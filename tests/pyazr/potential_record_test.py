#!/usr/bin/env python3
"""The hybrid potential as a format-3 record (`potential` in <potentials>, docs/upgrade/project_format3.md section
3.10a): its Coulomb radius reaches the channel functions, and pyazr's set_nuclear_potential, its view and its writer
agree with the record.

Checks, on the old original of tests/hybrid_potential (tests/convert/fixtures) converted to format 3 (pyazr.convert_azr,
to=3), pair 1 (p + 12C, Woods-Saxon
20 / 3.6 / 0.6 MeV, fm, a_c = 3.6 fm):
  (1) the conversion writes one record, coulomb_radius 3.5 fm, and no <potential> block;
  (2) with `coulomb_radius 5 fm` in the record, P and S from channel_functions (l = 0-2, E = 0.5-3 MeV, and S at the
      13N ground state, -1.94351 MeV) against an independent DOP853 integration with R_C = 5 fm
      (data/coulomb_radius.csv, data/make_coulomb_radius_reference.py) to 1e-6 (P relative, S absolute); the R_C = 3.5
      fm rows differ from those by more than 1e-4, so the check sees the radius;
  (3) the record as it was converted (3.5 fm) gives the 3.5 fm rows;
  (4) the session view reads coulomb_radius 5.0 for pair 1 and no default; set_nuclear_potential(coulomb_radius=5)
      on the conversion (3.5 fm) gives the 5 fm record's channel functions bit for bit;
  (5) save_fit of the format-3 project after set_nuclear_potential(1, coulomb_radius=5) writes the record with
      coulomb_radius 5 fm, and the written project gives the same channel functions; without a change the
      <potentials> section is written as it was read.
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/potential_record_test.py
"""
import csv
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "hybrid_potential")
sys.path.insert(0, HERE)
REF = os.path.join(HERE, "data", "coulomb_radius.csv")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session, convert_azr
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

ref = {}
with open(REF) as f:
    for row in csv.DictReader(f):
        ref[(float(row["R_C_fm"]), int(row["l"]), float(row["E_MeV"]), row["quantity"])] = float(row["value"])
ENERGIES = sorted({k[2] for k in ref if k[2] > 0})
BOUND = -1.94351


def section(text, name):
    a, b = text.index(f"<{name}>"), text.index(f"</{name}>")
    return text[a:b + len(name) + 3]


def functions(A):
    """P and S at the scattering energies and S at the bound state, l = 0-2, pair 1."""
    out = {}
    for L in (0, 1, 2):
        cf = A.channel_functions(1, ENERGIES, l=L)
        for i, E in enumerate(ENERGIES):
            out[(L, E, "P")] = float(cf["P"][i])
            out[(L, E, "S")] = float(cf["S"][i])
        out[(L, BOUND, "S")] = float(A.channel_functions(1, [BOUND], l=L)["S"][0])
    return out


def against(values, rc, label):
    worst_p = max(abs(values[(L, E, "P")] / ref[(rc, L, E, "P")] - 1) for L in (0, 1, 2) for E in ENERGIES)
    worst_s = max(abs(values[k] - ref[(rc,) + k]) for k in values if k[2] == "S")
    check(f"{label}: P to 1e-6 relative ({worst_p:.1e}), S to 1e-6 ({worst_s:.1e}) against DOP853 with R_C = {rc} fm",
          worst_p < 1e-6 and worst_s < 1e-6)


work = tempfile.mkdtemp(prefix="potential_record.")
try:
    os.makedirs(os.path.join(work, "output"))
    shutil.copytree(os.path.join(SRC, "data"), os.path.join(work, "data"))
    old = os.path.join(work, "hybrid_potential.azr")
    shutil.copy(os.path.join(ROOT, "tests", "convert", "fixtures", "tests", "hybrid_potential", "hybrid_potential.azr"),
                old)  # the original, kept for the converter
    f35 = os.path.join(work, "f3_35.azr")
    convert_azr(old, f35, to=3)
    text = open(f35).read()
    pot = section(text, "potentials")
    check("(1) one record, coulomb_radius 3.5 fm, no <potential> block",
          pot.count("\npotential ") == 1 and re.search(r"coulomb_radius +3\.5 fm", pot) and "<potential>" not in text,
          pot)
    f50 = os.path.join(work, "f3_50.azr")
    with open(f50, "w") as f:
        f.write(re.sub(r"(coulomb_radius +)3\.5 fm", r"\g<1>5 fm", text))

    contrast = min(abs(ref[(5.0, L, E, "P")] / ref[(3.5, L, E, "P")] - 1) for L in (0, 1, 2) for E in ENERGIES)
    check(f"the table sees the radius: P(5 fm) and P(3.5 fm) differ by >= {contrast:.1e} > 1e-4", contrast > 1e-4)

    A50 = Session(f50, cwd=work)
    v50 = functions(A50)
    against(v50, 5.0, "(2) record with coulomb_radius 5 fm")
    A35 = Session(f35, cwd=work)
    against(functions(A35), 3.5, "(3) record as converted (3.5 fm)")

    p1 = A50.nuclear_potential(1)
    check("(4) the view: pair 1 enabled, own, coulomb_radius 5.0; no default",
          p1.enabled and p1.own and p1.coulomb_radius == 5.0 and not A50.nuclear_potential(0).enabled, repr(p1))
    B = Session(f35, cwd=work)
    check("(4) the converted record's Coulomb radius is 3.5 fm", B.nuclear_potential(1).coulomb_radius == 3.5)
    B.set_nuclear_potential(1, coulomb_radius=5.0)
    vB = functions(B)
    check("(4) set_nuclear_potential(1, coulomb_radius=5) = the record, bit for bit", vB == v50,
          str(max(abs(vB[k] - v50[k]) for k in v50)))

    # (5) the writer
    C = Session(f35, cwd=work)
    same = os.path.join(work, "same.azr")
    C.save_fit(same, verify=False, norms="nominal")
    check("(5) save_fit without a change keeps <potentials> as read", section(open(same).read(), "potentials") == pot)
    C.set_nuclear_potential(1, coulomb_radius=5.0)
    moved = os.path.join(work, "moved.azr")
    C.save_fit(moved, verify=False, norms="nominal")
    written = section(open(moved).read(), "potentials")
    check("(5) save_fit writes the record p12C with coulomb_radius 5 fm",
          "potential p12C" in written and re.search(r"coulomb_radius +5 fm", written) and written.count("\npotential ") == 1,
          written)
    D = Session(moved, cwd=work)
    check("(5) the written project gives the record's channel functions", functions(D) == v50)
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"FAILED: {len(failures)} check(s)")
    sys.exit(1)
print("potential_record: all checks passed")
