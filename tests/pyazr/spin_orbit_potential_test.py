#!/usr/bin/env python3
"""The spin-orbit term of a hybrid pair through pyazr (docs/upgrade/spin_orbit_potential.md section 4.2, tests 3, 5
and 9; the project tests/spin_orbit_potential).

  (1) channel_functions(pair, E, L, J) of the format-3 project tests/spin_orbit_potential (p + 12C, Woods-Saxon
      20 / 3.6 / 0.6 with V_so = 6 MeV on the proton) against the independent DOP853 integration
      physics/tests/data/spin_orbit.csv (make_spin_orbit_reference.py): S, P and delta_hs = -phi for l = 0-2, both j,
      E = 0.5-3 MeV, and S below threshold, to 1e-6 (P relative);
  (2) the same term set through set_nuclear_potential on tests/hybrid_potential (a project without it) gives the
      record's values bit for bit, and nuclear_potential() reads it back; the record writer (potential_records)
      writes the spin_orbit word;
  (3) the refusals: J missing for the spin-orbit pair (E_SO_NEEDS_J), J given for a pair without the term
      (E_SO_NO_TERM), a j that l = 1 and spin 1/2 cannot make, the term on the spin-0 particle (E_SO_ZERO_SPIN), the
      term on the default pair 0;
  (4) external capture, entrance wave and bound state per J (models/capture.azr: 12C(p,g)13N into the 1/2- ground
      state): pyazr.ec_integrals with the term divided by the same without it, for the s- and d-wave E1 pathways
      (direct and channel capture), against ec_ratio.csv (ec_reference.py: DOP853 + mpmath, the engine's formula of
      the integrals' combination only) to 1e-5.
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/spin_orbit_potential_test.py
"""
import csv
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SO = os.path.join(ROOT, "tests", "spin_orbit_potential")
HP = os.path.join(ROOT, "tests", "hybrid_potential")
REF = os.path.join(ROOT, "physics", "tests", "data", "spin_orbit.csv")
EC_REF = os.path.join(SO, "expected", "ec_ratio.csv")
TERM = {"depth": 6.0, "radius": 3.6, "diffuseness": 0.6, "on": "light"}
SO_LINE = "    spin_orbit     woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s\n"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)


def refused(code, call):
    try:
        call()
    except ValueError as e:
        return code in str(e)
    return False


ref = {}
with open(REF) as f:
    for row in csv.DictReader(f):
        ref[(row["kind"], int(row["l"]), int(row["two_j"]), float(row["E_MeV"]))] = float(row["value"])

work = tempfile.mkdtemp(prefix="spin_orbit_potential.")
try:
    for d in ("so", "hp", "cap"):
        os.makedirs(os.path.join(work, d, "output"))
    shutil.copytree(os.path.join(SO, "data"), os.path.join(work, "so", "data"))
    shutil.copy(os.path.join(SO, "spin_orbit_potential.azr"), os.path.join(work, "so"))
    shutil.copytree(os.path.join(HP, "data"), os.path.join(work, "hp", "data"))
    shutil.copy(os.path.join(HP, "hybrid_potential.azr"), os.path.join(work, "hp"))
    A = Session(os.path.join(work, "so", "spin_orbit_potential.azr"), cwd=os.path.join(work, "so"))

    # (1) against the independent integration
    above = np.array([0.5, 1.0, 2.0, 3.0])
    below = np.array([-1.94351, -0.5])
    worst = 0.0
    for l in (0, 1, 2):
        for two_j in sorted({abs(2 * l - 1), 2 * l + 1}):
            c = A.channel_functions(1, above, l=l, J=two_j / 2)
            b = A.channel_functions(1, below, l=l, J=two_j / 2)
            for i, E in enumerate(above):
                worst = max(worst, abs(c["S"][i] - ref[("S", l, two_j, E)]), abs(c["P"][i] / ref[("P", l, two_j, E)] - 1),
                            abs(c["delta_hs"][i] + ref[("phi", l, two_j, E)]))
            for i, E in enumerate(below):
                worst = max(worst, abs(b["S"][i] - ref[("S", l, two_j, E)]))
    check(f"channel_functions with J = the independent integration (S, P, delta_hs; 60 values; worst {worst:.1e})",
          worst <= 1e-6, f"{worst:.3e} > 1e-6")

    # (2) the API's term = the record's term
    B = Session(os.path.join(work, "hp", "hybrid_potential.azr"), cwd=os.path.join(work, "hp"))
    B.set_nuclear_potential(1, spin_orbit=TERM)
    same = all(np.array_equal(A.channel_functions(1, above, l=1, J=j)[k], B.channel_functions(1, above, l=1, J=j)[k])
               for j in (0.5, 1.5) for k in ("S", "P", "delta_hs"))
    check("set_nuclear_potential(spin_orbit=...) = the record's term (bits)", same)
    check("nuclear_potential() reads the term back", B.nuclear_potential(1).spin_orbit == TERM,
          repr(B.nuclear_potential(1).spin_orbit))
    lines = A.potential_records()[1]
    check("potential_records writes the spin_orbit word",
          any(l.split()[:1] == ["spin_orbit"] and "operator 2l.s" in l and "on light" in l for l in lines), repr(lines))

    # (3) refusals
    check("E_SO_NEEDS_J: no J for the spin-orbit pair", refused("E_SO_NEEDS_J", lambda: A.channel_functions(1, above, l=1)))
    C = Session(os.path.join(work, "hp", "hybrid_potential.azr"), cwd=os.path.join(work, "hp"))
    check("E_SO_NO_TERM: a J for a pair without the term",
          refused("E_SO_NO_TERM", lambda: C.channel_functions(1, above, l=1, J=0.5)))
    check("a j that l = 1 and spin 1/2 cannot make", refused("TwoLDotI", lambda: A.channel_functions(1, above, l=1, J=2.5)))
    check("E_SO_ZERO_SPIN: the term on 12C",
          refused("E_SO_ZERO_SPIN", lambda: C.set_nuclear_potential(1, spin_orbit=dict(TERM, on="heavy"))))
    check("a term on the default pair 0 is refused",
          refused("name the pair", lambda: C.set_nuclear_potential(0, spin_orbit=TERM)))

    # (4) external capture per J: with / without the term against ec_reference.py
    cap = os.path.join(work, "cap")
    shutil.copytree(os.path.join(SO, "models", "data"), os.path.join(cap, "data"))
    text = open(os.path.join(SO, "models", "capture.azr")).read()
    with open(os.path.join(cap, "with.azr"), "w") as f:
        f.write(text)
    with open(os.path.join(cap, "without.azr"), "w") as f:
        f.write(text.replace(SO_LINE, ""))
    check("models/capture.azr has the term", SO_LINE in text)
    E = np.array([0.5, 1.0])
    paths = {}
    for name in ("with", "without"):
        M = Session(os.path.join(cap, name + ".azr"), cwd=cap)
        for p in M.ec_integrals(1, E):
            if p["multipolarity"] == 1 and p["lf"] == 1 and p["li"] in (0, 2):
                paths.setdefault((name, p["li"]), []).append(np.asarray(p["value"]))
    worst, n = 0.0, 0
    with open(EC_REF) as f:
        for row in csv.DictReader(f):
            li, e = int(row["li"]), float(row["E_MeV"])
            i = int(np.argmin(np.abs(E - e)))
            # the two pathways of (l_i, E1, l_f): direct (purely imaginary here) and channel capture
            for w, wo in zip(paths[("with", li)], paths[("without", li)]):
                r = w[i] / wo[i]
                want = (float(row["direct_ratio"]) if abs(w[i].real) == 0.0
                        else complex(float(row["channel_ratio_re"]), float(row["channel_ratio_im"])))
                worst = max(worst, abs(r - want))
                n += 1
    check(f"EC integrals with/without the term = ec_reference.py ({n} values, worst {worst:.1e})", n == 8 and worst <= 1e-5,
          f"n = {n}, {worst:.3e}")
except Exception as e:  # pragma: no cover
    check("session", False, repr(e))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("all checks passed")
