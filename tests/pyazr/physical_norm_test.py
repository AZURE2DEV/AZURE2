#!/usr/bin/env python3
"""calculate_chi2 (physical parameters) must be the same objective as calculate_chi2_rwa.

The physical vector shares its tail (normalizations, energy shifts, screening) with the RWA
vector, but CalculateChi2Physical never copied the norms out of it: the segments kept whatever
norm the previous call had set, so calculate_chi2(x) ignored the norm entries of x and the value
depended on the call history (audit item 6, fixed 2026-10-26).  It also left out the Park,
abundance and scattering-length terms, and returned 0 for a vector the transform rejects.

Checks, on tests/13N/13N.azr with one normalization freed (5 % systematic error):
  1. a moved norm changes calculate_chi2 exactly as it changes calculate_chi2_rwa;
  2. no hysteresis: calculate_chi2 at the nominal vector gives the nominal value after a call
     with the moved norm, and calculate_chi2_rwa afterwards is unchanged;
  3. with a <scatteringLengths> constraint the two bases still agree (the term is in both).
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
        from pyazr import azure2, AzrModel
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 0

    work = tempfile.mkdtemp(prefix="pyazr_physical_norm_")
    try:
        src = os.path.join(ROOT, "tests", "13N")
        shutil.copy(os.path.join(src, "13N.azr"), work)
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        model = AzrModel.from_file(os.path.join(work, "13N.azr"))
        model.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
        azr = os.path.join(work, "free_norm.azr")
        model.write(azr)

        print("1. a moved normalization moves calculate_chi2 as it moves calculate_chi2_rwa")
        with azure2(azr, cwd=work) as m:
            x = np.asarray(m.params_rwa, float)
            kinds = [p.kind for p in m.parameters if not p.fixed]
            inorm = kinds.index("norm")
            ref = float(m.calculate_chi2_rwa(x)[0])
            x8 = x.copy()
            x8[inorm] = 0.8
            ref8 = float(m.calculate_chi2_rwa(x8)[0])
            check("the norm matters in the RWA basis", rel(ref, ref8) > 1e-3, f"{ref} vs {ref8}")
            # transform_rwa maps the free R-matrix block; the tail (norms, shifts, screening)
            # is the same in both bases, so the physical free vector is the two joined
            lev = np.asarray(m.transform_rwa(x), float)
            check("transform_rwa covers the level block only", lev.size == inorm, f"{lev.size} vs {inorm}")
            ph = np.concatenate([lev, x[lev.size:]])
            ph8 = ph.copy()
            ph8[inorm] = 0.8
            p1 = float(m.calculate_chi2(ph)[0])
            check("physical == RWA at the nominal norm", rel(p1, ref) <= 1e-6, f"{p1} vs {ref}")
            p8 = float(m.calculate_chi2(ph8)[0])
            check("physical == RWA at norm 0.8", rel(p8, ref8) <= 1e-6, f"{p8} vs {ref8}")

            print("\n2. no hysteresis between calls")
            p1b = float(m.calculate_chi2(ph)[0])
            check("nominal again after the moved call", rel(p1b, ref) <= 1e-6, f"{p1b} vs {ref}")
            again = float(m.calculate_chi2_rwa(x)[0])
            check("calculate_chi2_rwa unchanged afterwards", rel(again, ref) <= 1e-12, f"{again} vs {ref}")

        print("\n3. with a scattering-length constraint the two bases still agree")
        model.add_scattering_length(1, 0.5, 1.0, 0.2)
        azr2 = os.path.join(work, "free_norm_sl.azr")
        model.write(azr2)
        with azure2(azr2, cwd=work) as m:
            x = np.asarray(m.params_rwa, float)
            kinds = [p.kind for p in m.parameters if not p.fixed]
            inorm = kinds.index("norm")
            x8 = x.copy()
            x8[inorm] = 0.9
            r = float(m.calculate_chi2_rwa(x8)[0])
            sl = float(m.scattering_length_chi2(x8))
            check("the constraint term is nonzero", sl > 1e-3, f"{sl}")
            lev8 = np.asarray(m.transform_rwa(x8), float)
            ph8 = np.concatenate([lev8, x8[lev8.size:]])
            p = float(m.calculate_chi2(ph8)[0])
            check("physical == RWA with the constraint", rel(p, r) <= 1e-6, f"{p} vs {r}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
