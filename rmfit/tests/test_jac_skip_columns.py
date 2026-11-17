#!/usr/bin/env python3
"""rmfit asks the engine only for the Jacobian columns it fits (2026-10-30, the evaluator "go ahead with the column mask").

Energy-shift, sqrt(E)-shift, shift-function and density columns are finite-differenced by the engine (two residual
evaluations each).  Before this, every free one was differenced on every Jacobian even when polish's bounds held it:
a polish with several compound nuclei (8 held shift/shift_fn terms) spent ~16 of ~18 evaluations per Jacobian (41 min each) on
columns it threw away and ran 7.5 h without finishing one.  Now ``Evaluator.residuals(x, jac=True, cols=...)``
passes the complement to the workers (pyazr ``residual_jacobian(skip=)``), and polish / gn_step ask for theirs.

Fixture tests/fixtures/shift_fn_groups/sfg.azr (two free constant shifts and two free shift-function terms) with the
0.8 and 1.0 MeV levels' energies and widths freed here, so there is something to fit with the shifts held.
  (a) SingleSession: with cols = everything but the four shift columns, those columns are exactly zero on the data
      rows (the analytic shift-penalty rows keep theirs) and the residuals and every other column are bit-identical
      to the full call (whose shift columns are not zero);
  (b) Sharded (4 shards): the same;
  (c) polish with the four shifts held by its bounds moves, and gives the same result as the old behaviour (an
      evaluator without ``cols``); the times are printed.
Needs a pyazr with residual_jacobian(skip=) (AZURE2 dev, 2026-10-30):
    AZURE2_ROOT=<build> python3 rmfit/tests/test_jac_skip_columns.py      (about two minutes)"""
import os
import shutil
import sys
import tempfile
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")   # bit-for-bit comparisons: no thread-order rounding

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_jac_skip_columns.py)")
    sys.exit(77)
from _fixtures import fixture  # noqa: E402  (fixtures/ or fixtures_p3/ by backend)
FIX = fixture("shift_fn_groups")
failures = []
SHIFT_KINDS = ("shift", "shift_sqrt", "shift_fn", "density")


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        failures.append(what)


class NoCols:
    """The evaluator as before the change: residuals() without ``cols``."""
    def __init__(self, ev):
        self._ev = ev

    def __getattr__(self, name):
        return getattr(self._ev, name)

    def residuals(self, x, jac=False):
        return self._ev.residuals(x, jac=jac)


def compare(ev, x, label):
    sh = [j for j, k in enumerate(ev.keys) if k.kind in SHIFT_KINDS]
    other = [j for j in range(x.size) if j not in sh]
    r, J = ev.residuals(x, jac=True)
    n = ev.npts                                   # data rows; the penalty rows follow
    check(len(sh) == 4 and len(other) > 0 and all(np.any(J[:n, j] != 0) for j in sh),
          f"({label}) full call: {len(sh)} shift columns, all nonzero on the data rows; {len(other)} other columns")
    try:
        r1, J1 = ev.residuals(x, jac=True, cols=other)
    except TypeError as err:
        check(False, f"({label}) residuals accepts cols ({err})")
        return
    check(not np.any(J1[:n, sh]), f"({label}) the held shift columns are zero on the data rows (not differenced)")
    check(np.array_equal(J1[n:], J[n:]), f"({label}) the penalty rows are unchanged")
    check(np.array_equal(r1, r), f"({label}) residuals bit-identical")
    check(np.array_equal(J1[:, other], J[:, other]), f"({label}) every other column bit-identical")
    mask = np.zeros(x.size, bool); mask[other] = True
    r2, J2 = ev.residuals(x, jac=True, cols=mask)
    check(np.array_equal(J2, J1) and np.array_equal(r2, r1), f"({label}) a bool mask gives the same")


def main():
    import rmfit  # noqa: F401  (puts AZURE2_ROOT on sys.path)
    import inspect
    from pyazr import azure2
    if "skip" not in inspect.signature(azure2.residual_jacobian).parameters:
        print(f"FAIL: the pyazr at {os.path.dirname(inspect.getfile(azure2))} has no residual_jacobian(skip=)")
        return 1
    from rmfit.engine import SingleSession, Sharded
    from rmfit.local import polish, default_bounds, BoundsPolicy
    work = tempfile.mkdtemp(prefix="rmfit_jacskip_")
    try:
        fix = os.path.join(work, "fix")
        shutil.copytree(FIX, fix)
        os.makedirs(os.path.join(fix, "output"), exist_ok=True); os.makedirs(os.path.join(fix, "checks"), exist_ok=True)
        azr = os.path.join(fix, "sfg.azr")
        lines = open(azr).read().split("\n")
        i0, i1 = lines.index("<levels>"), lines.index("</levels>")
        for i in range(i0 + 1, i1):
            f = lines[i].split()
            if len(f) > 10 and float(f[2]) in (0.8, 1.0):  # free these levels' energies and widths
                f[3], f[10] = "0", "0"
                lines[i] = " ".join(f)
        open(azr, "w").write("\n".join(lines))
        S = SingleSession(azr, cwd=fix, name="js")
        try:
            x1 = np.array(S.x0, float)
            for j, k in enumerate(S.keys):
                if k.kind == "shift":
                    x1[j] += 0.004 if k.seg == 1 else -0.003
                elif k.kind == "shift_fn":
                    x1[j] += 0.002 * (k.term + 1)
            compare(S, x1, "a")
            H = Sharded(azr, os.path.join(work, "sh"), nshards=4)
            try:
                from rmfit.keys import take
                compare(H, take(x1, S.keys, H.keys, H.x0), "b")
            finally:
                H.close()
            # (c)
            b = default_bounds(S, x1, BoundsPolicy())
            sh = [j for j, k in enumerate(S.keys) if k.kind in SHIFT_KINDS]
            b.fixed[sh] = True
            t = time.time(); new = polish(S, x1, bounds=b, max_nfev=6); t_new = time.time() - t
            t = time.time(); old = polish(NoCols(S), x1, bounds=b, max_nfev=6); t_old = time.time() - t
            check(new.nfev > 1 and new.objective < new.start_objective,
                  f"(c) the polish moved: {new.start_objective:.4f} -> {new.objective:.4f} in {new.nfev} evaluations")
            check(new.status == old.status and abs(new.objective - old.objective) <= 1e-12 * abs(old.objective)
                  and np.allclose(new.x, old.x, rtol=1e-12, atol=0),
                  f"(c) polish with held shifts: {new.objective:.10f} = old behaviour {old.objective:.10f}")
            print(f"  info  (c) time {t_new:.2f} s with the mask, {t_old:.2f} s without", flush=True)
        finally:
            S.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
