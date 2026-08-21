#!/usr/bin/env python3
"""The hybrid Coulomb model, switched on for one particle pair at a time.

AZURE2 can replace the pure Coulomb radial wave functions of a channel with the
solution through a nuclear potential -- Woods-Saxon or Gaussian -- integrated
outward by Numerov and matched to the Coulomb functions outside.  A potential
belongs to a *pair*: it bends the wave functions of that channel and no other.
So each pair carries its own setting, and ``pair=0`` is the default that a pair
without one inherits.

This example runs on ``tests/hybrid_potential`` -- ``13N`` cut down to its
three elastic p + 12C segments, so there is no external capture to recompute
and a potential change costs a second rather than a minute:

  1. reads back what every pair currently resolves to,
  2. scans the Woods-Saxon depth on the entrance pair and prints chi-squared,
  3. shows a per-pair setting overriding the default,
  4. shows the one thing that will bite -- the parameter vector is re-derived
     by the change, so it has to be re-read.

Run it from anywhere:

    python3 pyazr/examples/nuclear_potential.py
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "4")   # before numpy

import numpy as np

from pyazr import azure2

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, "..", "..", "tests",
                                        "hybrid_potential"))
AZR = os.path.join(PROJECT, "hybrid_potential.azr")


def chi2(model):
    """Chi-squared at the model's *current* parameters.

    Always re-read ``params_rwa`` rather than caching it across a potential
    change: the potential moves the penetrabilities and shift functions, which
    are what map physical widths to reduced-width amplitudes, so the vector is
    re-derived by the rebuild.
    """
    return float(np.sum(model.calculate_chi2_rwa(np.asarray(model.params_rwa, float))))


with azure2(AZR, cwd=PROJECT) as m:

    print("what each pair resolves to as the file is written")
    print("-" * 62)
    for key, setting in m.nuclear_potentials().items():
        who = "default" if key == 0 else f"pair {key}"
        print(f"  {who:<10} {setting}")

    # The committed project already switches the model on, so turn it off to
    # get something to measure against.
    m.clear_nuclear_potential(0)
    baseline = chi2(m)
    print(f"\n  chi-squared with no hybrid potential: {baseline:.3f}")

    # ------------------------------------------------------------------ 2 ---
    # Pair 1 is p + 12C, the entrance channel.  Pair 2 is the photon pair: it
    # has no Coulomb functions, so a potential there would be meaningless.
    print("\nWoods-Saxon depth scan on pair 1 (R = 3.6 fm, a = 0.6 fm)")
    print("-" * 62)
    print(f"  {'V0 (MeV)':>10}  {'chi-squared':>16}  {'vs baseline':>12}")
    for V0 in (0.0, 5.0, 10.0, 20.0, 30.0):
        if V0 == 0.0:
            m.set_nuclear_potential(1, enabled=False)
        else:
            m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True,
                                    V0=V0, R=3.6, a=0.6)
        c = chi2(m)
        print(f"  {V0:>10.1f}  {c:>16.3f}  {c / baseline:>11.4f}x")

    # ------------------------------------------------------------------ 3 ---
    # A pair's own setting wins over the default.  Here the default asks for a
    # Gaussian and pair 1 keeps its Woods-Saxon, so pair 1 is unmoved.
    print("\na per-pair setting overrides the default")
    print("-" * 62)
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True,
                            V0=20.0, R=3.6, a=0.6)
    own = chi2(m)
    m.set_nuclear_potential(0, type="Gaussian", enabled=True, V0=60.0, r0=4.0)
    after = chi2(m)
    print(f"  pair 1 Woods-Saxon V0 = 20        : {own:.3f}")
    print(f"  default switched to a Gaussian    : {after:.3f}")
    print(f"  pair 1 unaffected                 : {abs(after - own) < 1e-6}")
    print(f"  pair 1 still reports              : {m.nuclear_potential(1)}")

    # Dropping a pair's own setting makes it follow the default again.
    m.clear_nuclear_potential(1)
    print(f"  after clear_nuclear_potential(1)  : {m.nuclear_potential(1)}")
    print(f"  chi-squared now                   : {chi2(m):.3f}")

    # ------------------------------------------------------------------ 4 ---
    print("\nthe parameter vector is re-derived -- re-read it")
    print("-" * 62)
    m.clear_nuclear_potential(0)                       # back to no potential
    stale = np.asarray(m.params_rwa, float)            # captured *before*
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True,
                            V0=20.0, R=3.6, a=0.6)
    with_stale = float(np.sum(m.calculate_chi2_rwa(stale)))
    with_fresh = chi2(m)
    print(f"  reusing the pre-change vector : {with_stale:.3f}   <- wrong")
    print(f"  re-reading params_rwa         : {with_fresh:.3f}   <- right")
    print("\n  The two differ because a reduced-width amplitude means something")
    print("  different once the penetrabilities have moved.  Same caveat as a")
    print("  channel-radius change, and the same fix: re-read, then refit.")
