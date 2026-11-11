#!/usr/bin/env python3
"""External-capture integrals of pyazr in the formal basis are those the engine computes.

The engine evaluates the external-capture radial integrals at the final state's capture energy (ALevel::GetECEnergy):
its input (Brune) energy in every basis, since the formal parameters change the level energies of a block of several
levels (docs/upgrade/cnuc_split.md section 8).  ``azure2.ec_integrals`` must use the same energy.  Model: tests/13N
with a background 1/2- level at E_x = -3 MeV (proton ANC 2 fm^-1/2, no photon width) below the final state (the 1/2-
ground state), so that the boundary condition is set at that level and the final state's formal energy is not its
input energy.  Builds before 2026-11-10 evaluated the integrals of the formal basis at the
formal energy (cnuc_split.md section 11.9, item 4): the two bases differed by up to a factor 3 on this grid.

Check: on an energy grid, every pathway's integral in the formal basis equals the Brune basis's (relative 1e-12), and
the pathways are the same, in the same order.

Run from anywhere:  python3 tests/pyazr/ec_integrals_formal_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

SRC = os.path.join(ROOT, "tests", "13N")


def integrals(work, use_brune):
    from pyazr import Project, azure2
    d = os.path.join(work, "brune" if use_brune else "formal")
    shutil.copytree(SRC, d, ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(d, "output"), exist_ok=True)
    os.makedirs(os.path.join(d, "checks"), exist_ok=True)
    m = Project.open(os.path.join(d, "13N.azr"))
    m.add_level(J=0.5, parity=-1, energy=-3.0, fixed=True,
                channels=[{"pair": 1, "l": 1, "s": 0.5, "width": 2.0, "fixed": True},
                          {"pair": 2, "l": 1, "s": 0.5, "width": 0.0, "fixed": True}])
    path = m.save(os.path.join(d, "model.azr"))
    s = azure2(path, use_brune=use_brune)
    return s.ec_integrals(1, [0.2, 0.5, 1.0, 2.0])


def main():
    try:
        import numpy as np
        import pyazr  # noqa: F401
    except ImportError as e:
        print(f"skipping: {e}")
        return 77
    work = tempfile.mkdtemp(prefix="ec_integrals_formal_")
    failures = []
    try:
        brune, formal = integrals(work, True), integrals(work, False)
        keys = ("li", "lf", "si", "sf", "multipolarity", "radiation")
        if not brune:
            failures.append("no external-capture pathways in the Brune basis")
        if [tuple(p[k] for k in keys) for p in brune] != [tuple(p[k] for k in keys) for p in formal]:
            failures.append("the pathways differ between the bases")
        worst = 0.0
        for b, f in zip(brune, formal):
            vb, vf = np.asarray(b["value"]), np.asarray(f["value"])
            scale = np.maximum(np.abs(vb), 1e-300)
            worst = max(worst, float(np.max(np.abs(vf - vb) / scale)))
        print(f"  {len(brune)} pathways, largest relative difference formal / Brune {worst:.3e}")
        if worst > 1e-12:
            failures.append(f"formal-basis integrals differ from the Brune basis's by {worst:.3e} (relative)")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    for f in failures:
        print("  FAIL ", f)
    print("ec_integrals_formal_test:", "PASS" if not failures else f"FAIL ({len(failures)})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
