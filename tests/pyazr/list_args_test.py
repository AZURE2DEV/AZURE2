#!/usr/bin/env python3
"""A parameter vector passed as a Python list must work like an ndarray.

The binding's array arguments were taken by value.  pybind11 moves the cast's
result into such a parameter, so when the caller passes a list (or any
non-float64 array) the temporary ndarray ``forcecast`` makes is owned by the
parameter and freed at function exit -- inside the ``call_guard<gil_scoped_release>``
of ``calculate_chi2_rwa``, ``update_segments_rwa`` & co., i.e. with the GIL
released: a segfault in ``PyObject_Free`` (seen 2026-10-13 on the 22Ne+n member
of the natural-neon joint fit; rmfit never hit it because it passes arrays).
With the arguments taken by const reference the temporary lives in pybind's
argument loader, which outlives the guard.

Checks, on ``tests/identical_pp_res``: every array-taking call gives the same
result for a list, an int32 array and a float64 array.

Run from anywhere:  python3 tests/pyazr/list_args_test.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PROJECT = os.path.join(ROOT, "tests", "identical_pp_res")
AZR = os.path.join(PROJECT, "identical_pp_res.azr")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    try:
        from pyazr import azure2
        import pyazr._azure2  # noqa: F401
    except ImportError as e:
        print(f"pyazr engine module not built ({e}); nothing to test")
        return 0
    with azure2(AZR, cwd=PROJECT, use_brune=True) as s:
        x = np.asarray(s.params_rwa, float)
        xa = np.asarray(s.sess.params_all_rwa(), float)
        print(f"{os.path.basename(AZR)}: {x.size} free of {xa.size} parameters")
        forms = {"list": list, "int-ish list": lambda v: [float(f"{t:.6g}") for t in v],
                 "float32 array": lambda v: np.asarray(v, np.float32),
                 "float64 array": lambda v: np.asarray(v, float)}
        ref = float(s.calculate_chi2_rwa(x)[0])
        for name, f in forms.items():
            v = f(x)
            c = float(s.calculate_chi2_rwa(v)[0])
            tol = 1e-3 * abs(ref) if "float32" in name or "int-ish" in name else 1e-12 * abs(ref)
            check(f"calculate_chi2_rwa from a {name}", abs(c - ref) <= tol, f"{c} vs {ref}")
        # the other GIL-releasing array calls, from a list
        xl = list(x)
        check("update_segments_rwa from a list", s.sess.update_segments_rwa(xl) >= 0)
        s.sess.set_params_all_rwa(list(xa))
        check("set_params_all_rwa from a list", True)
        check("update_segments_all_rwa from a list", s.sess.update_segments_all_rwa(list(xa)) >= 0)
        phys = np.asarray(s.transform_rwa(x), float)
        check("transform_rwa from a list",
              np.allclose(np.asarray(s.transform_rwa(xl), float), phys))
        cp = float(s.sess.calculate_chi2_physical(list(phys)))
        check("calculate_chi2_physical from a list", abs(cp - ref) <= 1e-6 * abs(ref), f"{cp} vs {ref}")
        g = np.asarray(s.sess.calculate_chi2_grad_rwa(xl), float)
        check("calculate_chi2_grad_rwa from a list", g.size == x.size + 1 and np.all(np.isfinite(g)),
              f"size {g.size} vs {x.size + 1}")       # [chi2, grad...]
        r, J = s.residual_jacobian(xl)
        check("residual_jacobian from a list", len(r) > 0 and np.all(np.isfinite(J)))
        check("chi2 unchanged afterwards", abs(float(s.calculate_chi2_rwa(x)[0]) - ref) <= 1e-12 * abs(ref))
    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
