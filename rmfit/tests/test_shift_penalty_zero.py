#!/usr/bin/env python3
"""rmfit's shift penalty rows follow the engine's convention (rmfit.engine._shift_penalty_from_zero).

AZURE2 dev 2026-10-25 penalises every energy-shift term from 0, (dE/sigma)^2 (pyazr azure2.SHIFT_PENALTY_FROM_ZERO);
older builds from the .azr value.  Fixture tests/fixtures/shift_fn_groups/sfg.azr with the free terms' widths cut to
2 keV / 2e-3 so the convention matters (F: constant shifts 0.001 and -0.002 MeV, a0 = 0.001, a1 = 0.0002, all free):
SingleSession's and Sharded's objective equal pyazr's objective() (chi2 + penalties) at the start and at a perturbed
vector, standard and Brune.  Run with the pyazr of either convention:
    AZURE2_ROOT=<build> python3 rmfit/tests/test_shift_penalty_zero.py"""
import os
import re
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
FIX = os.path.join(HERE, "fixtures", "shift_fn_groups")
failures = []


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        failures.append(what)


def main():
    import rmfit  # noqa: F401
    import pyazr.datasets as ds
    if "shiftfngroup" not in getattr(ds, "_KEYWORDS", {}):
        print("SKIP: pyazr without shiftfngroup")
        return 0
    from pyazr import azure2
    from rmfit.engine import SingleSession, Sharded
    try:
        from rmfit.engine import _shift_penalty_from_zero
    except ImportError:   # an rmfit from before 2026-10-25: report the engine's convention from pyazr itself
        def _shift_penalty_from_zero():
            return bool(getattr(azure2, "SHIFT_PENALTY_FROM_ZERO", False))
    from rmfit.keys import take
    print(f"pyazr {os.path.dirname(ds.__file__)}: penalties from zero = {_shift_penalty_from_zero()}")
    work = tempfile.mkdtemp(prefix="rmfit_spz_")
    try:
        fix = os.path.join(work, "fix")
        shutil.copytree(FIX, fix)
        for d in ("output", "checks"):
            os.makedirs(os.path.join(fix, d), exist_ok=True)
        azr = os.path.join(fix, "sfg.azr")
        s = open(azr).read()
        s = s.replace("0 0.001 1 1 data/xs_40_lin_mislabelled.dat", "0 0.001 0.002 1 data/xs_40_lin_mislabelled.dat")
        s = s.replace('2 0.001 1 1 0.0002 1 1', '2 0.001 0.002 1 0.0002 0.002 1')
        s = s.replace("0 -0.002 1 1 data/xs_90_lin_mislabelled_p3.dat", "0 -0.002 0.002 1 data/xs_90_lin_mislabelled_p3.dat")
        open(azr, "w").write(s)
        for brune in (False, True):
            flags = dict(use_brune=brune)
            S = SingleSession(azr, cwd=fix, flags=flags, name=f"s_{int(brune)}")
            E = azure2(azr, cwd=fix, **flags)
            H = Sharded(azr, os.path.join(work, f"sh_{int(brune)}"), nshards=3, flags=flags)
            try:
                pen0 = E.penalties(np.asarray(S.x0, float))
                shiftpen = float(np.sum(pen0["shift"]) + np.sum(pen0["shift_fn"]))
                check(shiftpen > 0.1 if _shift_penalty_from_zero() else shiftpen < 1e-12,
                      f"use_brune={brune}: engine shift penalty at the start {shiftpen:.6g} "
                      f"({'from zero' if _shift_penalty_from_zero() else 'from the .azr value'})")
                for tag, x in (("start", np.asarray(S.x0, float)),
                               ("perturbed", np.asarray(S.x0, float) + np.where(
                                   [k.kind in ("shift", "shift_fn") for k in S.keys], 0.0015, 0.0))):
                    oE = float(E.objective(x))
                    oS = float(S.objective(x))
                    oH = float(H.objective(take(x, S.keys, H.keys, H.x0)))
                    check(abs(oS - oE) < 1e-9 * oE and abs(oH - oE) < 1e-9 * oE,
                          f"use_brune={brune} {tag}: Single {oS:.9f} = Sharded {oH:.9f} = engine objective {oE:.9f}")
            finally:
                H.close()
                S.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
