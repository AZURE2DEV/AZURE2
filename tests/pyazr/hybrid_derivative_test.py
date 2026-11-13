#!/usr/bin/env python3
"""The hybrid (nuclear-potential) Coulomb functions with a vanishing nuclear potential are the plain
Coulomb functions: P and S agree with the exact routine's.

`CoulFunc::computeHybrid` integrates H+ = G + iF inward from the matching radius by Numerov and
reads F, G and their derivatives at the channel radius.  Up to dev 2d48d61 the value was a linear
interpolation between the bracketing grid nodes and the derivative the chord slope of that cell, an
O(dr) quantity: with a nuclear potential of 1e-9 MeV (the parser refuses 0) the shift function of
p + 12C at the channel radius was off by 3.8e-3 (l = 0) to 4.2e-2 (l = 3), while F and G were good
to 1e-5 (audit 2026-10-25, item 16).  The grid is now laid so that the requested radius is a node,
and the derivative is a five-point central difference.  The ~1e-5 left in P and S until brief O
(2026-11-10) was not the integration but the Coulomb term of the potential, e^2 = 1.4399764 MeV fm
against the alpha hbar c of the boundary values; with alpha hbar c the integration agrees with the
Coulomb functions to ~1e-8 (physics_hybrid_potential), and this test asks for 1e-6.

Model: tests/hybrid_potential/hybrid_potential.azr (p + 12C, pair 1) with its <potentials> section
replaced: once empty (no potential), once one record for p + 12C, a Woods-Saxon of depth 1e-9 MeV.  Compared through
`Session.coulomb_functions` at l = 0..3, 0.3-4 MeV, at the channel radius and at an off-node radius.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/hybrid_derivative_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "hybrid_potential", "hybrid_potential.azr")

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

work = tempfile.mkdtemp(prefix="hybrid_derivative.")
try:
    os.makedirs(os.path.join(work, "output"))
    shutil.copytree(os.path.join(os.path.dirname(SRC), "data"), os.path.join(work, "data"))
    text = open(SRC).read()
    s = text.index("<potentials>")
    e = text.index("</potentials>") + len("</potentials>")

    def variant(name, records):
        path = os.path.join(work, name + ".azr")
        open(path, "w").write(text[:s] + "<potentials>\n" + records + "</potentials>" + text[e:])
        return path

    plain = Session(variant("plain", ""), cwd=work)
    tiny = Session(variant("tiny", "potential p12C\n    pair            p+12C\n    coulomb_radius  3.5 fm\n"
                                  "    central         woods_saxon depth 1e-9 MeV radius 3.6 fm diffuseness 0.6 fm\n"
                                  "    matching_radius 15 fm\n    step            0.0275 fm\n"), cwd=work)
    energies = np.array([0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])
    worst_s = worst_p = 0.0
    for L in (0, 1, 2, 3):
        for radius, label in ((0.0, "channel radius"), (5.4321, "5.4321 fm")):
            a = plain.coulomb_functions(1, energies, L=L, radius=radius)
            b = tiny.coulomb_functions(1, energies, L=L, radius=radius)
            ds = float(np.max(np.abs(b["S"] - a["S"])))
            dp = float(np.max(np.abs(b["P"] / a["P"] - 1.0)))
            worst_s = max(worst_s, ds)
            worst_p = max(worst_p, dp)
            check(f"l = {L}, {label}: S of the 1e-9 MeV potential = plain S (max |dS| {ds:.1e}, S {a['S'][2]:.4f} at 1 MeV)",
                  ds <= 1e-6, f"max |dS| {ds:.2e} > 1e-6")
            check(f"l = {L}, {label}: P agrees to 1e-6 (max rel {dp:.1e})", dp <= 1e-6, f"max rel {dp:.2e}")
    # not vacuous: the hybrid path really ran (its values are not bit-identical to the plain ones)
    a = plain.coulomb_functions(1, energies, L=0)
    b = tiny.coulomb_functions(1, energies, L=0)
    check("the hybrid path ran (values differ from the plain ones at the 1e-6 level)",
          0.0 < float(np.max(np.abs(b["G"] / a["G"] - 1.0))) < 1e-3)
    print(f"worst: |dS| {worst_s:.2e}, P rel {worst_p:.2e}")
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
