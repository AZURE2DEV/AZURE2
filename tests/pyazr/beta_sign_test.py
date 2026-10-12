#!/usr/bin/env python3
"""A beta-decay feeding amplitude keeps its sign through the physical transform.

``CNuc::GetTransformParams`` used to return ``sign(B) * B = |B|`` for a Fermi
or Gamow-Teller channel, so every consumer of the physical vector -- the GUI's
level write-back after a fit, the MCMC's physical parameters and pyazr's
``save_fit`` -- turned a negative amplitude positive.  On the 12C+a evaluation
(2026-10-12) a snapshot written that way re-read with the beta-delayed alpha
chi-squared 750 -> 29724, and ``save_fit(verify=True)`` could not see it: the
verification compares physical values, and both sides lose the sign alike.

The model is ``tests/beta_delayed_alpha/16N.azr``: 16N -> 16O* + beta, three
1- levels feeding the 12C+alpha channel, the 20 MeV background pole with a
negative beta amplitude, and the model's own spectrum as data (chi2 ~ 0).

Checks, in the Brune, Park and standard parametrizations:
  * transform_rwa carries the sign of every beta amplitude (Brune: identity);
  * a save_fit snapshot re-reads to the same chi-squared;
  * flipping that sign changes the chi-squared, so the two checks above
    are not vacuous.

Run from anywhere:  python3 tests/pyazr/beta_sign_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PROJECT = os.path.join(ROOT, "tests", "beta_delayed_alpha")
AZR = os.path.join(PROJECT, "16N.azr")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except Exception as err:                                  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

COMMON = dict(ignore_externals=True, use_long_wavelength=False)
MODES = (("brune", dict(use_brune=True)),
         ("park", dict(use_brune=False, use_park=True)),
         ("standard", dict(use_brune=False)))

for label, opts in MODES:
    print(f"\n{label} parametrization")
    snap = os.path.join(PROJECT, f"_beta_sign_{label}.azr")
    sav = os.path.splitext(snap)[0] + ".sav"
    try:
        with azure2(AZR, cwd=PROJECT, **COMMON, **opts) as m:
            x = np.asarray(m.params_rwa, float)
            ph = np.asarray(m.transform_rwa(x), float)
            beta = [p for p in m.parameters
                    if not p.fixed and p.kind == "width"
                    and p.radiation_type in ("F", "G")]
            idx = [p.free_index for p in beta]
            # Brune and Park amplitudes are the file's; the standard
            # parametrization's formal amplitudes mix the levels, so only
            # the count and "at least one negative" hold there.
            check("three free beta amplitudes, one negative",
                  len(idx) == 3 and (sum(x[i] < 0 for i in idx) == 1
                                     if label != "standard"
                                     else any(x[i] < 0 for i in idx)),
                  f"{[round(float(x[i]), 4) for i in idx]}")
            check("transform keeps the sign of every beta amplitude",
                  all(np.sign(ph[i]) == np.sign(x[i]) for i in idx),
                  f"x {[round(float(x[i]), 4) for i in idx]} -> "
                  f"{[round(float(ph[i]), 4) for i in idx]}")
            if label == "brune":
                check("Brune: the transform is the identity on them",
                      all(abs(ph[i] - x[i]) < 1e-12 for i in idx))
            chi2 = float(m.calculate_chi2_rwa(x)[0])
            flipped = x.copy()
            neg = [i for i in idx if x[i] < 0][0]
            flipped[neg] = -flipped[neg]
            chi2_flipped = float(m.calculate_chi2_rwa(flipped)[0])
            check("the sign matters to the data",
                  chi2_flipped > chi2 + 1.0, f"{chi2:.4g} vs {chi2_flipped:.4g}")
            m.save_fit(snap, x=x, param_sav=True, verify=True)
        with azure2(snap, cwd=PROJECT, **COMMON, **opts) as back:
            xb = np.asarray(back.params_rwa, float)
            chi2_back = float(back.calculate_chi2_rwa(xb)[0])
        # Brune and Park read the file's amplitudes back exactly; the standard
        # parametrization inverts its transform iteratively (energies to
        # 1e-6 MeV), which leaves ~1e-6 in the amplitudes and ~1e-5 in the
        # chi-squared of this steep spectrum.
        rtol = 1e-6 if label != "standard" else 1e-4
        check("snapshot re-reads to the same chi-squared",
              abs(chi2_back - chi2) <= rtol * max(1.0, chi2),
              f"{chi2:.6g} -> {chi2_back:.6g}")
        # In the standard parametrization the transform fixes a level's
        # amplitudes only up to a common sign (an eigenvector convention, no
        # change of model), so there the magnitudes are compared.
        same = (np.allclose(xb, x, rtol=1e-6, atol=1e-9) if label != "standard"
                else np.allclose(np.abs(xb), np.abs(x), rtol=1e-5, atol=1e-9))
        check("snapshot re-reads to the same level parameters"
              + (" (up to a level's sign)" if label == "standard" else ""),
              xb.shape == x.shape and same,
              f"max rel diff {np.max(np.abs(np.abs(xb) - np.abs(x)) / np.maximum(np.abs(x), 1e-12)):.3g}"
              if xb.shape == x.shape else f"{xb.shape} vs {x.shape}")
    finally:
        for f in (snap, sav):
            if os.path.exists(f):
                os.remove(f)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all beta-sign checks passed")
