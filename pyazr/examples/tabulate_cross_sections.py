"""Tabulate cross sections on an adaptive energy grid.

pyazr.tabulate() produces, for any (entrance, exit) particle-pair
combination, a non-uniform grid on which sigma(E) is tabulated such that
log-linear interpolation between the knots reproduces the engine to the
requested relative tolerance.  Knots pile up across resonances and thin out
on smooth stretches, so a table that would need ~1e5 uniform points for the
same fidelity typically needs a few hundred to a few thousand.

Run from the repository root against any project, e.g. the 13N test:

    python3 pyazr/examples/tabulate_cross_sections.py tests/13N/13N.azr
"""

import sys

from pyazr import tabulate

azr = sys.argv[1] if len(sys.argv) > 1 else "tests/13N/13N.azr"

# Pair keys are the file pair keys (the `ir` column of the level lines).
# Here: entrance pair 1 into exit pair 2, tabulated over the c.m. range
# 0.1 - 2.0 MeV at 0.5% interpolation fidelity.
tabs = tabulate(azr, pairs=[(1, 2)], e_min=0.1, e_max=2.0, rel_tol=5e-3)

tab = tabs[(1, 2)]
print(tab)
print(f"{len(tab.e_cm)} knots; finest spacing "
      f"{1e3 * min(b - a for a, b in zip(tab.e_cm, tab.e_cm[1:])):.2f} keV")

# The object interpolates directly ...
print("sigma(0.424 MeV) =", tab(0.424), "b")

# ... and writes a plain three-column table (E_cm, E_lab, sigma).
tab.save("sigma_1_to_2.dat")
print("wrote sigma_1_to_2.dat")

# Central values only, by design: to add uncertainty bands, evaluate the
# analytic parameter Jacobian on these knots and fold in your own fit
# covariance -- pyazr does not own fitting.
