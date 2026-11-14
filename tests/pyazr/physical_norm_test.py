#!/usr/bin/env python3
"""The engine's chi-squared of a physical parameter vector must be that of the fit vector (a test of the binding).

pyazr's API takes the fit vector only; the engine's physical path (``CalculateChi2Physical``, the binding's
``calculate_chi2_physical``) is what its fitters' physical basis uses, so it is tested here through the binding
(``Session._engine``), against the binding's ``calculate_chi2_rwa`` (the data chi-squared plus the model's terms).


The physical vector shares its tail (normalizations, energy shifts, screening) with the RWA
vector, but CalculateChi2Physical never copied the norms out of it: the segments kept whatever
norm the previous call had set, so the physical chi-squared ignored the norm entries of x and the value
depended on the call history (audit item 6, fixed 2026-10-26).  It also left out the Park,
abundance and scattering-length terms, and returned 0 for a vector the transform rejects.

Checks, on tests/13N/13N.azr with one normalization freed (5 % systematic error):
  1. a moved norm changes calculate_chi2_physical exactly as it changes calculate_chi2_rwa;
  2. no hysteresis: calculate_chi2_physical at the nominal vector gives the nominal value after a call
     with the moved norm, and calculate_chi2_rwa afterwards is unchanged;
  3. with a measured scattering length (a datum of <modelData>) the two bases still agree (the term is in both).
(A vector the transform rejects now gives +inf instead of 0; CNuc::TransformIn only warns on a
negative denominator, so no vector exercises that branch here.)

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/physical_norm_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import Session, Project
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 77

    work = tempfile.mkdtemp(prefix="pyazr_physical_norm_")
    try:
        src = os.path.join(ROOT, "tests", "13N")
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        model = Project.open(os.path.join(src, "13N.azr"))
        model.set("meyer_84", norm="1 fit prior 1 5%")
        azr = model.save(os.path.join(work, "free_norm.azr"))

        print("1. a moved normalization moves calculate_chi2_physical as it moves calculate_chi2_rwa")
        with Session(azr, cwd=work) as m:
            x = np.asarray(m.parameters.x0, float)
            kinds = [p.kind for p in m.parameters if not p.fixed]
            inorm = kinds.index("norm")
            ref = float(m._engine.calculate_chi2_rwa(x))
            x8 = x.copy()
            x8[inorm] = 0.8
            ref8 = float(m._engine.calculate_chi2_rwa(x8))
            check("the norm matters in the RWA basis", rel(ref, ref8) > 1e-3, f"{ref} vs {ref8}")
            # transform maps the R-matrix block; the tail (norms, shifts, screening) is the same in both bases,
            # so the physical free vector is transform's free entries
            ph = m.transform(x)[m.parameters.free_mask]
            check("transform keeps the tail after the level block", np.array_equal(ph[inorm:], x[inorm:]),
                  f"{ph[inorm:]} vs {x[inorm:]}")
            ph8 = ph.copy()
            ph8[inorm] = 0.8
            p1 = float(m._engine.calculate_chi2_physical(ph))
            check("physical == RWA at the nominal norm", rel(p1, ref) <= 1e-6, f"{p1} vs {ref}")
            p8 = float(m._engine.calculate_chi2_physical(ph8))
            check("physical == RWA at norm 0.8", rel(p8, ref8) <= 1e-6, f"{p8} vs {ref8}")

            print("\n2. no hysteresis between calls")
            p1b = float(m._engine.calculate_chi2_physical(ph))
            check("nominal again after the moved call", rel(p1b, ref) <= 1e-6, f"{p1b} vs {ref}")
            again = float(m._engine.calculate_chi2_rwa(x))
            check("calculate_chi2_rwa unchanged afterwards", rel(again, ref) <= 1e-12, f"{again} vs {ref}")

        print("\n3. with a scattering-length constraint the two bases still agree")
        model.add_datum("a_p12C", quantity="scattering_length", pair="p+12C", s="1/2", value="1.0 fm", error="0.2 fm")
        azr2 = model.save(os.path.join(work, "free_norm_sl.azr"))
        with Session(azr2, cwd=work) as m:
            x = np.asarray(m.parameters.x0, float)
            kinds = [p.kind for p in m.parameters if not p.fixed]
            inorm = kinds.index("norm")
            x8 = x.copy()
            x8[inorm] = 0.9
            r = float(m._engine.calculate_chi2_rwa(x8))
            sl = m.evaluate(x8).terms["model data"]
            check("the constraint term is nonzero", sl > 1e-3, f"{sl}")
            ph8 = m.transform(x8)[m.parameters.free_mask]
            p = float(m._engine.calculate_chi2_physical(ph8))
            check("physical == RWA with the constraint", rel(p, r) <= 1e-6, f"{p} vs {r}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
