#!/usr/bin/env python3
"""Park formalism: analytic derivatives, the J > 0 penalty, and Brune equivalence.

`--use-park` (pyazr `use_park=True`) fits the observed reduced width amplitudes.
Its level matrix is Brune's with every basis state rescaled by sqrt(J),
J = 1 - sum_c gamma^2 dS_c/dE, so three things have to hold and are pinned here:

  * the cross sections and chi-squared equal the default (Brune) run's, from the
    same .azr, to rounding;
  * `chi2_and_grad` and `residual_jacobian` agree with central finite
    differences of `objective` / the residuals -- the adjoint differentiates
    Park's level matrix, including the width- and energy-dependent overlap J and
    the J-dependent bound-state normalization of external capture.  The
    tolerance is 1e-4 (relative to the largest component), not the 1e-7 Brune
    mode reaches: Park's level matrix contains dS/dE at the level energy, which
    AZURE2 differentiates numerically (gsl_deriv_central), and a finite
    difference in E_lambda picks up that noise divided by its step.  On p+p the
    finite differences themselves scatter by 1e-5 between steps of 1e-4 and
    1e-3 MeV, bracketing the analytic value;
  * a level pushed past the bound (J < 0) raises the objective by
    (J/kParkNormScale)^2, with the matching gradient, and is reported by
    `park_norms`.

Models: tests/13N (capture with external capture and a bound final state) and
tests/identical_pp_res (two-channel 3P2-3F2 group, several levels per J^pi).

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/park_gradient_test.py
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


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "4")
    sys.stdout.reconfigure(line_buffering=True)
    import numpy as np
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)


def stage(tmp, name, tag=""):
    # one copy per session: the API reuses an output/intEC.dat it finds, and
    # a second session reloading another's cache is a known (separate) bug
    work = os.path.join(tmp, name + tag)
    shutil.copytree(os.path.join(ROOT, "tests", name), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    return os.path.join(work, name + ".azr")


def fd_gradient(f, x, h):
    """Five-point central differences (error O(h^4))."""
    g = np.zeros_like(x)
    for i in range(x.size):
        def at(d):
            z = x.copy(); z[i] += d
            return f(z)
        g[i] = (8.0 * (at(h[i]) - at(-h[i])) - (at(2 * h[i]) - at(-2 * h[i]))) / (12.0 * h[i])
    return g


def fd_jacobian(f, x, h):
    """Three-point central differences: each column costs two residual passes."""
    r0 = f(x)
    J = np.zeros((r0.size, x.size))
    for i in range(x.size):
        xp = x.copy(); xp[i] += h[i]
        xm = x.copy(); xm[i] -= h[i]
        J[:, i] = (f(xp) - f(xm)) / (2.0 * h[i])
    return J


def rel(a, b):
    scale = max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300)
    return float(np.max(np.abs(a - b)) / scale)


with tempfile.TemporaryDirectory() as tmp:
    # optional model names on the command line, for a quick single-model run
    for name in (sys.argv[1:] or ("identical_pp_res", "13N")):
        print(f"\n{name}")
        azr_b = stage(tmp, name, ".brune")
        azr = stage(tmp, name, ".park")

        # --- Brune and Park give the same cross sections from the same file.
        #     Brune's analytic gradient is established; its finite-difference
        #     agreement calibrates what the same check can resolve for Park.
        with azure2(azr_b, use_brune=True) as mb:
            xb = np.asarray(mb.params_rwa, float)
            xs_b = np.concatenate([np.asarray(s, float).ravel() for s in mb.calculate_rwa(xb)])
            chi_b = mb.objective(xb)
            hb = 1e-4 * (np.abs(xb) + 1.0)
            _, gb = mb.chi2_and_grad(xb)
            gfd_b = fd_gradient(lambda z: float(np.sum(mb.calculate_chi2_rwa(z))), xb, hb)
            ref = rel(gb, gfd_b)
            print(f"  (Brune mode: analytic gradient vs finite differences, rel {ref:.1e})")
        with azure2(azr, use_park=True) as m:
            x = np.asarray(m.params_rwa, float)
            xs_p = np.concatenate([np.asarray(s, float).ravel() for s in m.calculate_rwa(x)])
            chi_p = m.objective(x)
            check(f"{name}: Park cross sections = Brune's", rel(xs_p, xs_b) < 1e-9,
                  f"max rel diff {rel(xs_p, xs_b):.2e}")
            check(f"{name}: Park objective = Brune's", abs(chi_p - chi_b) <= 1e-7 * abs(chi_b),
                  f"{chi_p} vs {chi_b}")

            # Widths differ by sqrt(J) per level, and every J is in (0, 1].
            J = np.asarray(m.park_norms(x), float)
            check(f"{name}: every level has 0 < J <= 1", bool(np.all((J > 0) & (J <= 1 + 1e-12))),
                  f"J = {J}")
            pen = m.penalties(x)["park"]
            check(f"{name}: penalty is zero at the file's parameters", float(pen) == 0.0, f"{pen}")

            # --- analytic gradient of the objective vs central differences.
            # Steps scaled to the parameter: energies in MeV, amplitudes MeV^1/2.
            h = 1e-4 * (np.abs(x) + 1.0)
            val, grad = m.chi2_and_grad(x)
            check(f"{name}: chi2_and_grad value = calculate_chi2", abs(val - np.sum(m.calculate_chi2_rwa(x))) < 1e-6 * val)
            gfd = fd_gradient(lambda z: float(np.sum(m.calculate_chi2_rwa(z))), x, h)
            d = rel(grad, gfd)
            check(f"{name}: analytic gradient = finite differences (rel {d:.1e})", d < max(1e-4, 3.0 * ref),
                  "\n      analytic " + np.array2string(grad, precision=3) +
                  "\n      numeric  " + np.array2string(gfd, precision=3))

            # --- residual Jacobian vs central differences.
            r, Jac = m.residual_jacobian(x)
            Jac = np.asarray(Jac, float)
            Jfd = fd_jacobian(lambda z: np.asarray(m.residual_jacobian(z)[0], float), x, h)
            d = rel(Jac, Jfd)
            check(f"{name}: residual Jacobian = finite differences (rel {d:.1e})", d < max(1e-4, 3.0 * ref))
            check(f"{name}: sum(r^2) = chi2", abs(np.sum(np.asarray(r) ** 2) - val) < 1e-6 * val)

            # --- the J > 0 wall: scale one level's free widths up until J < 0 and
            #     check the penalty, its gradient, and that the fit objective sees it.
            widths = [p for p in m.parameters.widths if not p.fixed and p.free_index is not None]
            target = widths[0]
            same_level = [p.free_index for p in widths if (p.jgroup, p.level) == (target.jgroup, target.level)]
            for scale in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
                xw = x.copy()
                for fi in same_level:
                    xw[fi] *= scale
                Jw = np.asarray(m.park_norms(xw), float)
                if np.any(Jw < 0):
                    break
            check(f"{name}: a level can be driven to J < 0 (scale {scale:g}, min J {Jw.min():.3f})", bool(np.any(Jw < 0)))
            penw = float(m.penalties(xw)["park"])
            expect = float(np.sum((Jw[Jw < 0] / 1.0e-3) ** 2))
            check(f"{name}: penalty = sum (J/1e-3)^2 over J<0", abs(penw - expect) < 1e-6 * expect, f"{penw} vs {expect}")
            # The engine's chi2 carries the wall; the per-segment values do not.
            chi_wall = float(np.sum(m.calculate_chi2_rwa(xw)))
            chi_data = float(np.sum(m.segment_chi2(xw)))
            check(f"{name}: engine chi2 at J<0 = data chi2 + penalty",
                  abs((chi_wall - chi_data) - penw) < 1e-6 * penw, f"{chi_wall - chi_data} vs {penw}")
            valw, gradw = m.chi2_and_grad(xw)
            hw = 1e-5 * (np.abs(xw) + 1.0)
            gfdw = fd_gradient(lambda z: float(np.sum(m.calculate_chi2_rwa(z))), xw, hw)
            d = rel(gradw, gfdw)
            # At J = -1.3 the penalty gradient is 2J/1e-6 times dJ/dE, so the noise
            # of S' and S'' in dJ/dE is amplified a thousandfold: 1 % is what a
            # finite difference can confirm here, and enough for a minimizer to
            # be walked out of the wall.
            check(f"{name}: gradient with the wall active = finite differences (rel {d:.1e})", d < 1e-2,
                  "\n      analytic " + np.array2string(gradw, precision=3) +
                  "\n      numeric  " + np.array2string(gfdw, precision=3))

print()
if failures:
    print(f"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
