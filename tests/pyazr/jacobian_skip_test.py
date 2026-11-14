#!/usr/bin/env python3
"""objective().residuals_and_jacobian(x, skip=cols) leaves the finite differences of the skipped columns out, and nothing else.

Energy-shift, sqrt(E)-shift, shift-function and target-density columns are central differences, two residual
evaluations each.  Up to 2026-10-30 every free one was differenced on every call, including the columns a fit
holds through its own bounds: a 17O+a rmfit polish with 8 such terms spent ~16 of its ~18 evaluations per Jacobian
(41 min each) on columns it then threw away.  ``skip`` names free columns not to difference; they come back zero.

Model: tests/13N with only the data sets meyer_84 and meyer_114 (12C(p,p) at 84.3 and 114.5 deg) active, each with a free
constant energy shift (10 keV).  Checks: (a) without skip both shift columns are nonzero; (b) skipping the first
shift column leaves it exactly zero in the data rows while the residuals, the other shift column and every adjoint
column are bit-identical to the full call, and the rows after the data rows (the shifts' priors, analytic) are the
full call's; (c) skipping both gives two zero columns, the rest identical; (d) a skip index
outside the free vector is refused.  A build without the option fails (a TypeError on ``skip``).

Run from anywhere:  python3 tests/pyazr/jacobian_skip_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np                         # noqa: E402
    from pyazr import Session, Project          # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")


def make_model(work):
    p = Project.open(SRC)
    for name, r in p.records().items():
        if r.record_kind != "dataset":
            continue
        if name in ("meyer_84", "meyer_114"):  # segments 3 and 4 of the old file: a free shift each (the first
            # with the prior of its old line's error, 5 keV)
            p.set(name, shift="offset 0.01 MeV fit" + (" prior 0 0.005" if name == "meyer_84" else ""))
        else:
            p.set(name, active="no")
    return p.save(os.path.join(work, "skip.azr"))


work = tempfile.mkdtemp(prefix="pyazr_jac_skip_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    a = Session(make_model(work), cwd=work)
    x = np.array(a.parameters.x0, float)
    free = [p for p in a.parameters if not p.fixed]
    sh = [i for i, p in enumerate(free) if p.kind == "shift"]
    check("model has two free shifts", len(sh) == 2, f"{len(sh)}")
    # skip concerns the finite-difference columns of the data rows (the engine's residual Jacobian); the objective's
    # rows after them (here the shifts' priors) carry analytic derivatives, the same with and without skip
    n = a.evaluate(x).data[-1].rows.stop
    r_all, J_all = a.objective().residuals_and_jacobian(x)
    r, J = np.asarray(r_all, float)[:n], np.asarray(J_all, float)[:n]
    check("(a) without skip both shift columns are nonzero", all(np.any(J[:, k] != 0) for k in sh))
    others = [c for c in range(x.size) if c not in sh]
    try:
        r1_all, J1_all = a.objective().residuals_and_jacobian(x, skip=[sh[0]])
        check("(b) the rows after the data rows (the priors) are those of the full call",
              np.array_equal(r1_all[n:], r_all[n:]) and np.array_equal(J1_all[n:], J_all[n:]) and len(r_all) > n)
        r1, J1 = np.asarray(r1_all, float)[:n], np.asarray(J1_all, float)[:n]
        check("(b) the skipped column is zero", not np.any(J1[:, sh[0]]))
        check("(b) residuals bit-identical", np.array_equal(r1, r))
        check("(b) the other shift column bit-identical", np.array_equal(J1[:, sh[1]], J[:, sh[1]]))
        check("(b) every adjoint column bit-identical", np.array_equal(J1[:, others], J[:, others]))
        r2, J2 = a.objective().residuals_and_jacobian(x, skip=sh)
        r2, J2 = np.asarray(r2, float)[:n], np.asarray(J2, float)[:n]
        check("(c) both skipped: two zero columns, the rest identical",
              not np.any(J2[:, sh]) and np.array_equal(J2[:, others], J[:, others]) and np.array_equal(r2, r))
        try:
            a.objective().residuals_and_jacobian(x, skip=[x.size])
            check("(d) a skip index outside the free vector is refused", False, "no error")
        except Exception as err:
            check("(d) a skip index outside the free vector is refused", "skip column" in str(err), str(err))
    except TypeError as err:
        check("objective().residuals_and_jacobian accepts skip", False, str(err))
    a.close() if hasattr(a, "close") else None
finally:
    shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("all checks passed")
