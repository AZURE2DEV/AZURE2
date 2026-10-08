#!/usr/bin/env python3
"""Two engine sessions in one process must not see each other's state.

Reported from the 8Be+alpha rate work (2026-10-07): a data-mode session that is
already open when an extrapolation-mode session is created in the same process
computed wrong residuals from then on -- capture segments only, by up to 28 %,
silently, and still after the second session was closed.  Cause: the external-
capture amplitude cache (`g_ecAmplitudeCache`) was one process-global object,
keyed by (k-group, EC m-group, entrance, exit, segment key) with nothing that
identifies the session; `EPoint::GetECAmplitudeWithShift`, which every capture
point goes through, interpolates from it at the point's energy.  Creating the
second session re-created the global (deleting the first session's entries)
and then filled it with the second model's grid under colliding keys, so the
first session interpolated its capture amplitudes from the wrong grid.  The
binding now gives every session its own cache and swaps it in per call, like
the Config.

Pinned on tests/13N (12C+p capture with external capture):

  * a data session's cross sections and per-segment chi2 are unchanged after
    an extrapolation session on the same model is opened, used, and closed;
  * and after a second DATA session on the same model is opened and used;
  * the extrapolation session's own cross sections are unaffected too.

Needs the compiled engine; skips cleanly without it.
Run from anywhere:  python3 tests/pyazr/shared_state_test.py
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
    sys.stdout.flush()
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2, AzrModel
except Exception as err:
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)


def stage(tmp, tag):
    work = os.path.join(tmp, "13N" + tag)
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    return os.path.join(work, "13N.azr")


def flat(seqs):
    return np.concatenate([np.asarray(s, float).ravel() for s in seqs])


with tempfile.TemporaryDirectory() as tmp:
    azr_a = stage(tmp, ".data")
    # a capture extrapolation grid for the second session (lab energies)
    azr_b = stage(tmp, ".extrap")
    mdl = AzrModel.from_file(azr_b)
    mdl.clear_extrapolations()
    mdl.add_extrapolation(entrance=1, exit=2, e_min=0.2, e_max=2.5, e_step=0.02,
                          observable="angle-integrated")
    mdl.write(azr_b)

    A = azure2(azr_a)
    A.__enter__()
    x = np.asarray(A.params_rwa, float)

    def probe(m, z):
        # cross sections and per-segment chi2: cheap, and exactly what the bug moved
        return flat(m.calculate_rwa(z)), np.asarray(m.segment_chi2(z), float)

    xs0, seg0 = probe(A, x)
    r0 = xs0
    check("data session evaluates (reference cross sections)", xs0.size > 0 and np.all(np.isfinite(xs0)))

    # --- an extrapolation session on the same model, used, then closed
    with azure2(azr_b, data_mode=False) as B:
        xsB = flat(B.calculate_rwa(np.asarray(B.params_rwa, float)))
        check("extrapolation session evaluates", xsB.size > 0 and np.all(np.isfinite(xsB)))
        # the data session, while B is open
        r1, seg1 = probe(A, x)
        d = float(np.max(np.abs(r1 - r0))) if r1.shape == r0.shape else float("inf")
        check(f"data cross sections unchanged while an extrapolation session is open (max |d sigma| {d:.1e})",
              r1.shape == r0.shape and np.allclose(r1, r0, rtol=1e-10, atol=0.0),
              f"per-segment chi2 {seg0} -> {seg1}")
        # and B's own numbers are not disturbed by A's evaluation
        xsB2 = flat(B.calculate_rwa(np.asarray(B.params_rwa, float)))
        check("extrapolation cross sections stable while the data session evaluates",
              np.allclose(xsB2, xsB, rtol=1e-12, atol=0.0))

    # --- after B is closed
    r2, seg2 = probe(A, x)
    check("data cross sections unchanged after the extrapolation session closed",
          r2.shape == r0.shape and np.allclose(r2, r0, rtol=1e-10, atol=0.0),
          f"per-segment chi2 {seg0} -> {seg2}")
    check("per-segment chi2 unchanged", np.allclose(seg2, seg0, rtol=1e-10, atol=0.0))

    # --- a second DATA session on the same model
    azr_c = stage(tmp, ".data2")
    with azure2(azr_c) as C:
        xc = np.asarray(C.params_rwa, float)
        rc, _ = probe(C, xc)
        check("second data session reproduces the first's cross sections",
              rc.shape == r0.shape and np.allclose(rc, r0, rtol=1e-10, atol=0.0))
        r3, _ = probe(A, x)
        check("first data session unchanged while a second data session is open",
              r3.shape == r0.shape and np.allclose(r3, r0, rtol=1e-10, atol=0.0))
    # one Jacobian pass: the adjoint reads the same per-session cache
    rj = np.asarray(A.residual_jacobian(x)[0], float)
    check("residuals consistent with the per-segment chi2 after all of the above",
          np.isclose(float(np.sum(rj ** 2)), float(np.sum(seg0)), rtol=1e-8))
    A.__exit__(None, None, None)

print()
if failures:
    print(f"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
