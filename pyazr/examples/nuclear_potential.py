#!/usr/bin/env python3
"""The hybrid Coulomb model, switched on for one particle pair at a time.

AZURE2 can replace the pure Coulomb radial wave functions of a channel with the
solution through a nuclear potential -- Woods-Saxon or Gaussian -- integrated
outward by Numerov and matched to the Coulomb functions outside.  A potential
belongs to a *pair*: it bends the wave functions of that channel and no other.
In a project it is a ``potential`` record of the ``<potentials>`` section, one
per pair; a change of it is an edit of the :class:`pyazr.Project` and a new
session, as a change of a channel radius is.

This example runs on ``tests/hybrid_potential`` -- ``13N`` cut down to its
three elastic p + 12C segments, so there is no external capture to recompute
and a new session costs about a second:

  1. reads back what every pair resolves to,
  2. scans the Woods-Saxon depth on the entrance pair and prints chi-squared,
     one edited copy of the project and one session per depth,
  3. shows the one thing that will bite -- the fit vector of a session without
     the potential does not describe the model with it.

Run it from anywhere:

    python3 pyazr/examples/nuclear_potential.py
"""
import os
import shutil
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "4")   # before numpy

import numpy as np

from pyazr import Project, Session

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, "..", "..", "tests",
                                        "hybrid_potential"))
AZR = os.path.join(PROJECT, "hybrid_potential.azr")
RECORD = "p12C"                                   # the project's record of pair 1 (p + 12C)


def chi2(model, x=None):
    """Chi-squared at ``x``, by default the session's start ``parameters.x0``."""
    return model.evaluate(np.asarray(model.parameters.x0 if x is None else x, float)).chi2


def variant(work, name, depth):
    """A copy of the project with pair 1's potential at ``depth`` MeV (no record at all for ``None``)."""
    p = Project.open(AZR)
    if depth is None:
        p.remove(RECORD)
    else:
        p.set(RECORD, central=f"woods_saxon depth {depth} MeV radius 3.6 fm diffuseness 0.6 fm")
    return p.save(os.path.join(work, name))


work = tempfile.mkdtemp(prefix="nuclear_potential.")
try:
    shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))

    with Session(AZR, cwd=work) as m:
        print("what each pair resolves to as the file is written")
        print("-" * 62)
        for key, setting in m.nuclear_potentials().items():
            who = "default" if key == 0 else f"pair {key}"
            print(f"  {who:<10} {setting}")

    with Session(variant(work, "_none.azr", None), cwd=work) as m:
        baseline = chi2(m)
        stale = np.asarray(m.parameters.x0, float)  # the fit vector of the model without a potential
    print(f"\n  chi-squared with no hybrid potential: {baseline:.3f}")

    # ------------------------------------------------------------------ 2 ---
    # Pair 1 is p + 12C, the entrance channel.  Pair 2 is the photon pair: it
    # has no Coulomb functions, so a potential there would be meaningless.
    print("\nWoods-Saxon depth scan on pair 1 (R = 3.6 fm, a = 0.6 fm)")
    print("-" * 62)
    print(f"  {'V0 (MeV)':>10}  {'chi-squared':>16}  {'vs baseline':>12}")
    for V0 in (5.0, 10.0, 20.0, 30.0):
        with Session(variant(work, f"_V{V0:g}.azr", V0), cwd=work) as m:
            c = chi2(m)
        print(f"  {V0:>10.1f}  {c:>16.3f}  {c / baseline:>11.4f}x")

    # ------------------------------------------------------------------ 3 ---
    print("\nthe fit vector belongs to its session -- read it from the new one")
    print("-" * 62)
    with Session(variant(work, "_V20.azr", 20.0), cwd=work) as m:
        with_stale = chi2(m, stale)
        with_fresh = chi2(m)
    print(f"  the vector of the session without it : {with_stale:.3f}   <- wrong")
    print(f"  this session's parameters.x0         : {with_fresh:.3f}   <- right")
    print("\n  The two differ because the project holds physical widths and a")
    print("  reduced-width amplitude means something different once the")
    print("  penetrabilities have moved.  Same caveat as a channel-radius change,")
    print("  and the same fix: take the new session's vector, then refit.")
finally:
    shutil.rmtree(work, ignore_errors=True)
