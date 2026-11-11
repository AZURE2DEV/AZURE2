#!/usr/bin/env python3
"""pyazr.stopping (the binding of stopping/) against LunaSRIM's values in stopping/tests/data/srim_reference.csv.

The reference rows were written by stopping/tools/make_srim_reference.py with LunaSRIM 1.0.0 (energies in keV there):
  se        electronic stopping of an ion in one element, solid (0) or gas (1)    tolerance 1e-10 relative
  compound  electronic and nuclear stopping of a compound (C applied once)        tolerance 1e-10 relative
  layer     energy after a layer and the TRIM and Bohr straggling widths          tolerance 1e-5 relative (LunaSRIM's RK4)
the same tolerances as the C++ tests (stopping/tests).  Also: N(E_in, E_out(E_in, N)) = N, the provenance cites SRIM,
the per-atom-of-an-element rule of the module's docstring, and refusals as Python exceptions.

    PYAZR_ENGINE_DIR=<build>/pyazr-out python3 tests/pyazr/stopping_test.py
"""
import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
REFERENCE = os.path.join(ROOT, "stopping", "tests", "data", "srim_reference.csv")

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    from pyazr import stopping
    from pyazr.azure2 import _azure2  # noqa: F401
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)

SYMBOL = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg",
          13: "Al", 14: "Si", 20: "Ca", 22: "Ti", 26: "Fe", 28: "Ni", 29: "Cu", 40: "Zr", 47: "Ag", 73: "Ta", 74: "W",
          79: "Au", 82: "Pb"}
PHASE = {"0": "solid", "1": "gas"}

rows = defaultdict(list)
with open(REFERENCE) as f:
    for r in csv.reader(line for line in f if not line.startswith("#")):
        rows[r[0]].append(r)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.max(np.abs(a / b - 1.0)) if a.size else 0.0


# electronic stopping of ions in elements, grouped by (Z1, M1, Z2, phase)
groups = defaultdict(list)
for _, z1, m1, z2, phase, e, se in rows["se"]:
    groups[(int(z1), float(m1), int(z2), phase)].append((float(e) / 1000.0, float(se)))
worst = 0.0
for (z1, m1, z2, phase), pts in groups.items():
    e, se = np.array(pts).T
    out = stopping.stopping(stopping.Material(SYMBOL[z2], phase=PHASE[phase]), (z1, m1), e)
    worst = max(worst, rel(out["electronic"], se))
check(f"electronic stopping, {len(rows['se'])} rows of {len(groups)} ion-element-phase pairs, to 1e-10", worst <= 1e-10,
      f"largest {worst:.1e}")

# compounds: electronic and nuclear with the compound correction applied once
groups = defaultdict(list)
for _, formula, phase, c, z1, m1, e, se, sn in rows["compound"]:
    groups[(formula, phase, float(c), int(z1), float(m1))].append((float(e) / 1000.0, float(se), float(sn)))
we = wn = 0.0
for (formula, phase, c, z1, m1), pts in groups.items():
    e, se, sn = np.array(pts).T
    out = stopping.stopping(stopping.Material(formula, phase=PHASE[phase], compound_correction=c), (z1, m1), e)
    we, wn = max(we, rel(out["electronic"], se)), max(wn, rel(out["nuclear"], sn))
check(f"compounds, {len(rows['compound'])} rows: electronic and nuclear to 1e-10", we <= 1e-10 and wn <= 1e-10,
      f"largest {we:.1e}, {wn:.1e}")

# layers: energy loss and straggling against lunasrim.straggling.energy_loss
wd = wt = wb = 0.0
for _, formula, phase, z1, m1, e, n, eout, s_trim, s_bohr in rows["layer"]:
    mat, p, e, n = stopping.Material(formula, phase=PHASE[phase]), (int(z1), float(m1)), float(e) / 1000.0, float(n)
    loss = stopping.energy_loss(mat, p, [e], [n])[0]
    wd = max(wd, abs(loss / (e - float(eout) / 1000.0) - 1.0))
    wt = max(wt, rel(stopping.straggling(mat, p, [e], [n], model="trim")["straggling"], float(s_trim) / 1000.0))
    wb = max(wb, rel(stopping.straggling(mat, p, [e], [n], model="bohr")["straggling"], float(s_bohr) / 1000.0))
check(f"{len(rows['layer'])} layers: energy loss, TRIM and Bohr straggling to 1e-5", max(wd, wt, wb) <= 1e-5,
      f"largest {wd:.1e}, {wt:.1e}, {wb:.1e}")

# the integrals agree with each other; the range bounds the energy after
ta2o5 = stopping.Material("Ta2O5")
e_in = np.array([0.3, 1.0, 3.0])
n = np.array([1e17, 1e18, 5e18])
back = stopping.areal_density_between(ta2o5, "p", e_in, stopping.energy_after(ta2o5, "p", e_in, n))
check("N(E_in, E_out(E_in, N)) = N for p in Ta2O5, to 1e-8", rel(back, n) <= 1e-8, f"{rel(back, n):.1e}")
r = stopping.range_of(ta2o5, "a", [1.0, 5.0])
check("the range of alpha in Ta2O5 grows with E", r[0] > 0 and r[1] > r[0], str(r))

# a stopping per atom of one element is eps / its atom fraction (the engine's stopping evaluation, TiN with N active)
tin = stopping.Material("TiN")
eps = stopping.stopping(tin, (1, 1.007825), [0.5])["total"][0]
check("TiN: two atoms, N's fraction 1/2; p (1.007825 u) at 0.5 MeV 12.00133921 eV/(1e15 atoms/cm2) per atom",
      tin.atom_fraction(7) == 0.5 and abs(eps / 12.00133921 - 1) < 1e-9, str(eps))
check("the phase of TiN is solid, of N gas, of CH4 solid unless said",
      tin.phase_used == "solid" and stopping.Material("N").phase_used == "gas" and stopping.Material("CH4").phase_used == "solid"
      and stopping.Material("CH4", phase="gas").phase_used == "gas")
text = stopping.provenance(ta2o5, "p")
check("provenance cites SRIM and names LunaSRIM", stopping.CITATION in text and "LunaSRIM" in text and
      "Ziegler" in stopping.CITATION)
check("projectiles by name: p and a carry AME2020 masses", stopping.projectile("p")[0] == 1 and
      abs(stopping.projectile("a")[1] - 4.002603) < 1e-6)

# refusals
for what, call in (("an unknown element", lambda: stopping.Material("Xx2O")),
                   ("an energy <= 0", lambda: stopping.stopping(ta2o5, "p", [0.0])),
                   ("a layer beyond the range", lambda: stopping.energy_after(ta2o5, "p", [0.1], [1e21])),
                   ("an unknown straggling model", lambda: stopping.straggling(ta2o5, "p", [1.0], [1e18], model="chu"))):
    try:
        call()
        check(f"{what} is refused", False, "no exception")
    except Exception as err:  # the C++ exceptions arrive as ValueError / RuntimeError
        check(f"{what} is refused ({type(err).__name__})", True)

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("pyazr.stopping: ok")
