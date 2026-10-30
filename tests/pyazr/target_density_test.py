#!/usr/bin/env python3
"""Fitted target density (``densityfit`` / ``densitygroup`` on a <targetInt> line).

The areal density N of a thin target is a fit parameter ``density_<effect>``.  It moves the energy-loss window
dE (int_{E-dE}^{E} dE'/epsilon = N) and nothing else: the yield is the integral over the window divided by the same
N, so the computed value is the depth average of the cross section, with no factor of N.

Model (tests/13N, 12C(p,p) at 84.3 deg lab, Meyer, 0.3-0.9 MeV across the 0.457 MeV resonance), a constant stopping
cross section eps = 1e-20 MeV cm^2/atom, N0 = 1e18 /cm^2 (10 keV lab, 9.2 keV c.m.):
  key 3   the data, target integration, ``densityfit 1 5%``            (effect 1: density_1, prior 5e16, max 1.15e18)
  key 17  a 0.25 keV grid at the same angle, no experimental effect    (the bare cross section)
  key 18  the data, a plain line at N = 1.1 N0                         (effect 2: no parameter)
  key 19  the data, ``densityfit 1 1e16 densitygroup g`` at N0          (effect 3: density_3, max 1.1e18)
  key 20  the data, ``densitygroup g`` at 2 N0 (takes the group's N)    (effect 4: follows density_3)
Elastic scattering keeps the c.m. angle of a fixed lab angle at every energy, so key 17 is the bare c.m. cross
section the sub-points of keys 3, 18-20 sample.

1. the parameters: two densities, their priors (5 % converted), their maxima, the follower has none
2. the window, and no factor of N: key 3 at N = 0.9, 1.0, 1.1 N0 equals the average of key 17 over
   [E_cm - F eps N, E_cm] (F = m2/(m1+m2)), 3e-4 relative (adaptive sub-point grid against a 0.25 keV trapezoid)
3. a fitted density at 1.1 N0 = the plain line at 1.1 N0 (key 18), 3e-4
4. the group: keys 19 and 20 identical, both follow density_3, and equal key 3 at the same N (3e-4)
5. penalties()["density"] = sum ((N - N0)/sigma)^2 and objective() includes it
6. the residual-Jacobian column of density_1 = central differences of the residuals, 1e-3
7. an energy shift on key 3 moves the window with the point: equal to the bare average below the shifted energy
8. beyond its maximum a density is computed at the maximum (warned once), the penalty with it
9. refused at load: densityfit on a line without target integration, a group whose first line has no densityfit

Run from anywhere:  python3 tests/pyazr/target_density_test.py
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
    from pyazr import azure2                   # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(0)

# numpy >= 2.0 renamed trapz to trapezoid; the linux CI job runs Ubuntu 22.04's numpy 1.21, which has only trapz.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")
N0, EPS = 1.0e18, 1.0e-20
STOP = f'"a0" 1 {EPS:g}'
TARGET = {
    3: f'1 "3" 50 0 0 1 {N0:g} {STOP} 0 0 0 "" 0 0 densityfit 1 5%',
    18: f'1 "18" 50 0 0 1 {1.1 * N0:g} {STOP} 0 0 0 "" 0 0',
    19: f'1 "19" 50 0 0 1 {N0:g} {STOP} 0 0 0 "" 0 0 densityfit 1 1e16 densitygroup g',
    20: f'1 "20" 50 0 0 1 {2 * N0:g} {STOP} 0 0 0 "" 0 0 densitygroup g',
}


def make_model(work, shift=None, target=None):
    lines = open(SRC).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    seg = [l.split() for l in lines[i0 + 1:i1] if l.split()]
    meyer = list(seg[2])
    meyer[3], meyer[4] = "0.3", "0.9"
    for f in seg:
        f[0] = "0"
    a = list(meyer)
    if shift is not None:
        a[11], a[13] = f"{shift:g}", "1"
    seg[2] = a
    fine = list(meyer)
    fine[3], fine[4], fine[14] = "0.28", "0.92", "data/fine_84.dat"
    seg += [fine, list(meyer), list(meyer), list(meyer)]          # keys 17-20
    lines[i0 + 1:i1] = [" ".join(f) for f in seg]
    j0, j1 = lines.index("<targetInt>"), lines.index("</targetInt>")
    lines[j0 + 1:j1] = target if target is not None else [TARGET[k] for k in (3, 18, 19, 20)]
    path = os.path.join(work, "density.azr")
    open(path, "w").write("\n".join(lines))
    e = np.arange(0.28, 0.92 + 1e-9, 0.00025)
    np.savetxt(os.path.join(work, "data", "fine_84.dat"), np.column_stack([e, np.full_like(e, 84.3), np.ones_like(e),
                                                                           np.full_like(e, 0.01)]), fmt="%.6f")
    return path


def bare_average(e_fine, s_fine, lo, hi):
    """Trapezoid average of the bare cross section over [lo, hi] (linear interpolation at the ends)."""
    grid = np.concatenate([[lo], e_fine[(e_fine > lo) & (e_fine < hi)], [hi]])
    vals = np.interp(grid, e_fine, s_fine)
    return _trapezoid(vals, grid) / (hi - lo)


def session(work, **kw):
    return azure2(make_model(work, **kw), cwd=work)


work = tempfile.mkdtemp(prefix="pyazr_density_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    print("1. parameters")
    a = session(work)
    dens = list(a.parameters.densities)
    names = [p.name for p in dens]
    check("two density parameters, density_1 and density_3", names == ["density_1", "density_3"], names)
    if len(dens) == 2:
        d1, d3 = dens
        check("density_1: free, N0, prior 5 % = 5e16, max 1.15 N0",
              not d1.fixed and d1.value == N0 and d1.prior_nominal == N0 and abs(d1.prior_sigma - 5e16) < 1e4
              and abs(d1.density_max - 1.15e18) < 1e6, d1)
        check("density_3 (group g): prior 1e16, max 1.1 N0", abs(d3.prior_sigma - 1e16) < 1e4
              and abs(d3.density_max - 1.1e18) < 1e6, d3)
    free = [p for p in a.parameters if not p.fixed]
    k1 = [i for i, p in enumerate(free) if p.name == "density_1"][0]
    k3 = [i for i, p in enumerate(free) if p.name == "density_3"][0]
    x0 = np.array(a.params_rwa, float)

    # result segments are the active ones in order: 3, 17, 18, 19, 20; F = m2/(m1+m2) of p + 12C from the
    # session's own conversion, E_cm / E_lab of the fine grid
    seg = {key: i for i, key in enumerate((3, 17, 18, 19, 20))}
    en = a.calculate_energies(a.params)
    e_fine = np.asarray(en[seg[17]], float)
    lab = np.arange(0.28, 0.92 + 1e-9, 0.00025)
    F = float(np.median(e_fine / lab))

    def values(x):
        v = a.calculate_rwa(x)
        return {key: np.asarray(v[i], float) for key, i in seg.items()}

    v0 = values(x0)
    s_fine = v0[17]
    e_pts = np.asarray(en[seg[3]], float)

    print("\n2. the window, and no factor of N")
    for f in (0.9, 1.0, 1.1):
        x = x0.copy()
        x[k1] = f * N0
        v = values(x)[3]
        ref = np.array([bare_average(e_fine, s_fine, e - F * EPS * f * N0, e) for e in e_pts])
        rel = np.max(np.abs(v / ref - 1))
        check(f"N = {f:.1f} N0: depth average of the bare cross section (worst {rel:.1e})", rel < 3e-4, rel)

    print("\n3. a fitted density = the plain line at the same density")
    x = x0.copy()
    x[k1] = 1.1 * N0
    v = values(x)
    rel = np.max(np.abs(v[3] / v[18] - 1))
    check(f"key 3 at 1.1 N0 = key 18 (plain, 1.1 N0) (worst {rel:.1e})", rel < 3e-4, rel)

    print("\n4. the density group")
    x = x0.copy()
    x[k1] = 0.95 * N0
    x[k3] = 0.95 * N0
    v = values(x)
    check("keys 19 and 20 identical (20's own 2 N0 replaced by the group's)", np.array_equal(v[19], v[20]),
          np.max(np.abs(v[19] - v[20])))
    rel = np.max(np.abs(v[19] / v[3] - 1))
    check(f"the group at 0.95 N0 = key 3 at 0.95 N0 (worst {rel:.1e})", rel < 3e-4, rel)

    print("\n5. penalties")
    x = x0.copy()
    x[k1] = 1.08 * N0
    x[k3] = 0.97 * N0
    want = ((0.08 * N0) / 5e16) ** 2 + ((0.03 * N0) / 1e16) ** 2
    pen = a.penalties(x)["density"]
    check(f"penalties()['density'] = {pen:.6f} (want {want:.6f})", abs(pen - want) < 1e-9 * want, pen)
    obj = a.objective(x)
    chi = float(np.sum(a.calculate_chi2_rwa(x)))
    check("objective() = chi2 + the density penalty", abs(obj - chi - want) < 1e-6 * obj, (obj, chi, want))

    print("\n6. the residual-Jacobian column")
    x = x0.copy()
    x[k1] = 1.03 * N0
    r, J = a.residual_jacobian(x)
    h = 1e-4 * N0
    xp, xm = x.copy(), x.copy()
    xp[k1] += h
    xm[k1] -= h
    rp, _ = a.residual_jacobian(xp)
    rm, _ = a.residual_jacobian(xm)
    fd = (rp - rm) / (2 * h)
    col = J[:, k1]
    scale = np.max(np.abs(fd))
    err = np.max(np.abs(col - fd)) / scale if scale > 0 else np.inf
    check(f"Jacobian column of density_1 = central differences (worst {err:.1e} of the largest entry)",
          scale > 0 and err < 1e-3, err)

    print("\n7. an energy shift moves the window with the point")
    b = session(work, shift=0.005)
    xb = np.array(b.params_rwa, float)
    enb = b.calculate_energies(b.params)
    vb = b.calculate_rwa(xb)
    e_shift = np.asarray(enb[0], float)
    rel = np.max(np.abs(np.asarray(vb[0]) / np.array([bare_average(e_fine, s_fine, e - F * EPS * N0, e)
                                                       for e in e_shift]) - 1))
    check(f"key 3 shifted 5 keV: bare average below the shifted energy (worst {rel:.1e})", rel < 3e-4, rel)
    check("its energies moved by F x 5 keV", abs(np.median(e_shift - e_pts) - F * 0.005) < 1e-9,
          np.median(e_shift - e_pts))

    print("\n8. beyond the maximum")
    x = x0.copy()
    x[k1] = 1.3 * N0
    v_hi = values(x)[3]
    x[k1] = 1.15 * N0
    v_max = values(x)[3]
    rel = np.max(np.abs(v_hi / v_max - 1))
    check(f"1.3 N0 is computed at the maximum 1.15 N0 ({rel:.1e})", rel < 1e-12, rel)
    x[k1] = 1.3 * N0
    check("and penalised there", abs(a.penalties(x)["density"] - (0.15 * N0 / 5e16) ** 2) < 1e-9)

    print("\n9. refused at load")
    for label, bad in (("densityfit on a Gaussian-only line", ['1 "3" 50 1 0.002 0 0 "" 0 0 0 0 "" 0 densityfit 1 5%']),
                       ("densityfit on a thick target (density 0)", [f'1 "3" 50 0 0 1 0 {STOP} 0 0 0 "" 0 0 densityfit 1 5%']),
                       ("a group whose first line has no densityfit",
                        [f'1 "3" 50 0 0 1 {N0:g} {STOP} 0 0 0 "" 0 0 densitygroup g',
                         f'1 "18" 50 0 0 1 {N0:g} {STOP} 0 0 0 "" 0 0 densityfit 1 5% densitygroup g'])):
        try:
            session(work, target=bad)
            check(f"{label}: refused", False, "loaded")
        except Exception as e:
            check(f"{label}: refused", True)
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
