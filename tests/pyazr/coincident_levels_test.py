#!/usr/bin/env python3
"""Two levels of one J^pi at one energy: the limits of the divided differences (LevelPairShift.h).

The Brune level matrix, Park's overlap and the Brune -> formal transformation hold divided differences
(S_a - S_b)/(E_a - E_b) of the shift function between two levels of one J^pi: 0/0 when the levels coincide.  Below
kCoincidentLevelGap_MeV (1e-6 MeV) they are replaced by their Taylor forms about the pair.  Model:
tests/identical_pp_res with a second 2- level (two channels), its energy put on the first one's through the
parameter vector of `Session.objective()`.  Pinned here, in the Brune and Park bases (the formal basis fits
formal energies; its transformation at a coincidence is section 8 of physics_level_parameters_unit):

  * chi-squared is continuous: at the gaps 0, 5e-10, 2e-9, 1e-7, 0.99e-6 (inside the limit), 1.01e-6 and 1e-5 MeV
    (outside) each value equals chi2(0) + g (dchi2/dE_b(0) + dchi2/dE_b(g)) / 2 from the analytic gradients.
    Inside the limit the Taylor forms are exact to gap^2 S'''/12: bound 1e-10 relative (measured 3e-12 to 9e-12,
    which also bounds the trapezoid's own error).  Across the threshold chi2 steps by what the numerical S' (Richardson
    differences, about 1e-8 relative) differs from the divided difference of S: bound 1e-7 relative (measured 7e-9
    Brune, 5e-9 Park with the default Coulomb table, 2.3e-8 and 2.0e-8 with the exact Coulomb functions, which is
    what a machine without a table file computes; the step is inferred to be that of S', not shown term by term).
    Without the limits chi2(0)
    is 134,108 instead of 53,300.75 (no error), and the gaps 5e-10 and 2e-9 are off by 1e-5 relative;
  * the analytic gradient (energies and widths of the two levels) equals five-point central differences at gap 0
    (steps 2e-7 MeV, inside the limit on both sides) to 1e-5 relative to its largest component (measured 4e-8
    with the table, 2e-7 without) and at gap 1e-5 (steps 1e-6 MeV) to 1e-4 (measured 2e-6 with the table, 2.7e-5
    without: the differences of chi2 see the rounding of the exact functions in the divided difference of S);
  * the same with a third 2- level, so that the group takes the transformed (channel-space) route: before, its K and
    N were not finite at gap 0 and the group fell back to the level route (chi2 198,494, gradient NaN);
  * a path that drives E_b through E_a (41 points over +-2e-6 MeV) has no jump beyond the step of S' at the
    threshold: every second difference of chi2 along it is below 1e-7 chi2 (measured 1.2e-8 with the table, 4.1e-8
    without, at +-1e-6).

Needs numpy and the compiled engine; exits 77 (skipped) without them.  Run: python3 tests/pyazr/coincident_levels_test.py
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
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    sys.stdout.reconfigure(line_buffering=True)
    import numpy as np
    from pyazr import Session
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

SECOND_LEVEL = ("2 -1 1.25 0 1 1 2 2 7 1 0 50000.0 0.5 1 0.5 1 0 1.00727646688 1.00727646688 1 1 -1.0 -1.0 0 0 0.0 0 "
                "3.0 0 0 0\n"
                "2 -1 1.25 0 1 1 2 6 7 1 0 -3000.0 0.5 1 0.5 1 0 1.00727646688 1.00727646688 1 1 -1.0 -1.0 0 0 0.0 0 "
                "3.0 0 0 0\n\n")


# A third 2- level (1.4 MeV): three levels on two channels take the transformed (channel-space) route, whose K and N
# hold the same divided differences; with two levels the group takes the level route.
THIRD_LEVEL = SECOND_LEVEL.replace("1.25", "1.4").replace(" 7 1 0 50000.0", " 8 1 0 30000.0").replace(
    " 7 1 0 -3000.0", " 8 1 0 2000.0")


def stage(tmp, tag, levels):
    work = os.path.join(tmp, tag)
    shutil.copytree(os.path.join(ROOT, "tests", "identical_pp_res"), work)
    for junk in ("output", "checks", "expected"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    for need in ("output", "checks"):
        os.makedirs(os.path.join(work, need))
    path = os.path.join(work, "identical_pp_res.azr")
    with open(path) as f:
        text = f.read()
    anchor = "2 1 1.5 0 1 1 0 4 5"
    assert anchor in text
    with open(path, "w") as f:
        f.write(text.replace(anchor, levels + anchor, 1))
    return path


def fd_gradient(f, x, h, which):
    g = np.zeros(len(which))
    for k, i in enumerate(which):
        def at(d):
            z = x.copy()
            z[i] += d
            return f(z)
        g[k] = (8.0 * (at(h) - at(-h)) - (at(2 * h) - at(-2 * h))) / (12.0 * h)
    return g


# The formal basis fits formal energies (no divided differences); its transformation at a coincidence is
# physics_level_parameters_unit's section 8.
BASES = (("Brune", dict(use_brune=True), SECOND_LEVEL), ("Park", dict(use_park=True), SECOND_LEVEL),
         ("Brune, transformed route", dict(use_brune=True), SECOND_LEVEL + THIRD_LEVEL),
         ("Park, transformed route", dict(use_park=True), SECOND_LEVEL + THIRD_LEVEL))

with tempfile.TemporaryDirectory() as tmp:
    for basis, options, levels in BASES:
        print(f"\n{basis}")
        with Session(stage(tmp, basis.replace(",", "").replace(" ", "_"), levels), **options) as m:
            names = list(m.parameters.names)
            x = np.asarray(m.parameters.x0, float)
            # Levels 4 (2-, 1.2 MeV) and 5 (the added 2-, 1.25 MeV) in the order of the file.
            a, b = names.index("energy_4"), names.index("energy_5")
            f = m.objective()

            def at_gap(gap):
                z = x.copy()
                z[b] = z[a] + gap
                return z

            x0 = at_gap(0.0)
            chi0, grad0 = f.value_and_grad(x0)
            check(f"{basis}: chi2 at gap 0 finite", np.isfinite(chi0) and np.all(np.isfinite(grad0)),
                  f"{chi0}, {grad0}")
            for gap in (5e-10, 2e-9, 1e-7, 0.99e-6, 1.01e-6, 1e-5):
                chi, grad = f.value_and_grad(at_gap(gap))
                trapezoid = chi0 + 0.5 * gap * (grad0[b] + grad[b])
                d = abs(chi - trapezoid) / chi0
                bound = 1e-10 if gap < 1e-6 else 1e-7
                check(f"{basis}: chi2 at gap {gap:.2e} = chi2(0) + its slopes' trapezoid (rel {d:.1e} < {bound:.0e})",
                      d < bound,
                      f"{chi:.10f} vs {trapezoid:.10f}")

            # The analytic gradient against central differences, at the coincidence and beside it.
            for gap, h in ((0.0, 2e-7), (1e-5, 1e-6)):
                z = at_gap(gap)
                _, g = f.value_and_grad(z)
                which = [a, b] + [i for i, n in enumerate(names) if n.startswith(("width_4_", "width_5_"))]
                num = fd_gradient(f, z, h, which)
                scale = max(np.max(np.abs(g[which])), 1e-300)
                d = float(np.max(np.abs(g[which] - num)) / scale)
                check(f"{basis}: gradient at gap {gap:.0e} = central differences (rel {d:.1e})",
                      d < (1e-5 if gap == 0.0 else 1e-4),
                      "\n      analytic " + np.array2string(g[which], precision=6) +
                      "\n      numeric  " + np.array2string(num, precision=6))

            # A path that drives E_b through E_a: no jump.
            path = np.array([f(at_gap(g)) for g in np.linspace(-2e-6, 2e-6, 41)])
            second = np.abs(path[2:] - 2 * path[1:-1] + path[:-2])
            check(f"{basis}: no jump along E_b through E_a (max second difference {second.max() / chi0:.1e} chi2)",
                  bool(np.all(np.isfinite(path))) and second.max() < 1e-7 * chi0)

print("\nall passed" if not failures else f"\n{len(failures)} failed: {failures}")
sys.exit(1 if failures else 0)
