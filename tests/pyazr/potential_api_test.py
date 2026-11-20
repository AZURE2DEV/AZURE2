#!/usr/bin/env python3
"""The session's view of a pair's potential (docs/upgrade/potential_background.md; review of 2026-11-20, finding 3).

On tests/potential_background/models/identity_a35.azr (3He + 4He, the Gaussian potential VD, a `wave` with its own J and
an interior background on every partial wave, a_c = 3.5 fm):
  (a) nuclear_potential(1).waves lists the seven `wave` fields with their l, J, depth and background;
  (b) channel_functions needs J (the waves name J) and refuses it for a pair whose functions do not depend on J;
  (c) channel_functions gives the dressed functions together: P = rho/(F^2 + G^2) and delta_hs = -atan2(F, G) of the
      same F, G (before, F, G and delta_hs were undressed while P and S were dressed), and they differ from the
      undressed ones (the same model without its backgrounds);
  (d) potential_bound_state(1, l=1, J=3/2) finds the p3/2 state at -1.4548867 MeV with ANC 4.09316 fm^-1/2
      (tests/potential_capture/make_reference.py, scipy and mpmath: -1.4548867432 MeV, 4.0931604450);
  (e) potential_box_states(1, l=0, J=1/2) numbers its states 0, 1, 2, ... ascending, with positive residues;
  (f) potential_scattering_length(1, l=0, J=1/2) gives the s-wave a0 and r0 of VD (two bound s states) as
      physics/tests/data/make_potential_tuning_reference.py (scipy, mpmath, numpy: a0 30.953078 fm, r0 0.929216 fm) to
      3e-5 and 1e-3, and refuses the charged p wave (the Coulomb-modified expansion of l > 0 is not built).
"""
import math
import os
import sys
import tempfile
import warnings

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok or not detail else f"  -- {detail}"))
    if not ok:
        failures.append(name)


# pyazr from this checkout (ctest sets no PYTHONPATH; before 2026-11-21 the test was skipped under ctest for that)
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
try:
    import numpy as np
    import pyazr
except Exception as e:  # pragma: no cover
    print(f"skip: pyazr not importable ({e})")
    sys.exit(77)

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "potential_background", "models", "identity_a35.azr")
# a copy in a temporary directory with its output directory (a session of the extrapolations writes there)
WORK = tempfile.mkdtemp(prefix="potential_api_")
os.makedirs(os.path.join(WORK, "output"))
os.makedirs(os.path.join(WORK, "checks"))
MODEL = os.path.join(WORK, "identity_a35.azr")
with open(SOURCE) as f, open(MODEL, "w") as g:
    g.write(f.read())

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    with pyazr.Session(MODEL, extrapolations=True) as s:
        pot = s.nuclear_potential(1)
        w = pot.waves
        check("(a) seven waves", len(w) == 7, str(len(w)))
        check("(a) the s wave: l 0, J 1/2, depth 67.5 MeV, background B 0, no capture",
              w and w[0]["l"] == 0 and w[0]["J"] == 0.5 and w[0]["depth"] == 67.5 and w[0]["background"] is not None
              and w[0]["background"]["boundary"] == 0.0 and not w[0]["background"]["capture_interior"]
              and not w[0]["only_bound_state"], str(w[:1]))
        try:
            s.channel_functions(1, 1.0, l=0)
            check("(b) J needed for a wave with its own J", False, "no error")
        except ValueError as e:
            check("(b) J needed for a wave with its own J", "give J" in str(e), str(e))
        E = np.array([0.5, 1.0, 2.0, 4.0])
        cf = s.channel_functions(1, E, l=0, J=0.5)
        rho = math.sqrt(2.0 * 3.0160293 * 4.0026032 / (3.0160293 + 4.0026032) * 931.49410242) / 197.3269804 * 3.5 * np.sqrt(E)
        p_ok = np.allclose(cf["P"], rho / (cf["F"] ** 2 + cf["G"] ** 2), rtol=1e-6)
        d_ok = np.allclose(cf["delta_hs"], -np.arctan2(cf["F"], cf["G"]), rtol=1e-12, atol=1e-14)
        check("(c) P = rho/(F^2 + G^2) and delta_hs = -atan2(F, G): one set of dressed functions", p_ok and d_ok,
              f"P {cf['P']} vs {rho / (cf['F'] ** 2 + cf['G'] ** 2)}")
        b = s.potential_bound_state(1, l=1, J=1.5, energy=-1.5)
        check("(d) the p3/2 bound state of the potential: E and ANC as scipy's",
              b["found"] and abs(b["E"] + 1.4548867432) < 2e-6 and abs(b["anc"] - 4.0931604450) < 2e-4, str(b))
        bs = s.potential_box_states(1, l=0, J=0.5, e_min=-50.0, e_max=30.0)
        n = bs["nodes"]
        check("(e) box states numbered 0, 1, ... ascending, positive residues",
              len(n) >= 3 and list(n) == list(range(len(n))) and np.all(np.diff(bs["E"]) > 0) and np.all(bs["gamma2"] > 0),
              str(bs))
        sl = s.potential_scattering_length(1, l=0, J=0.5)
        ref = {}
        with open(os.path.join(ROOT, "physics", "tests", "data", "potential_tuning.csv")) as f:
            for line in f:
                kind, system, l, x, value = line.strip().split(",")
                if system == "he3" and l == "0" and x == "67.5" and kind in ("a", "r"):
                    ref[kind] = float(value)
        check("(f) a0 and r0 of the potential's s wave as scipy's, two bound s states",
              abs(sl["a"] / ref["a"] - 1) < 3e-5 and abs(sl["r"] / ref["r"] - 1) < 1e-3 and sl["bound_states"] == 2
              and sl["e_window"][0] < sl["e_window"][1], f"{sl} vs {ref}")
        try:
            s.potential_scattering_length(1, l=1, J=1.5)
            check("(f) the charged p wave refused", False, "no error")
        except ValueError as e:
            check("(f) the charged p wave refused", "l = 0 only" in str(e), str(e))

# (c) continued: without the backgrounds the functions are the undressed ones, and differ
try:
    plain = os.path.join(WORK, "plain.azr")
    with open(MODEL) as f, open(plain, "w") as g:
        for line in f:
            g.write(line.replace(" background interior boundary 0.0", ""))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(plain, extrapolations=True) as s:
            cp = s.channel_functions(1, E, l=0, J=0.5)
    check("(c) the dressed phase differs from the undressed one", np.max(np.abs(cp["delta_hs"] - cf["delta_hs"])) > 1e-3,
          f"{cp['delta_hs']} {cf['delta_hs']}")
finally:
    import shutil
    shutil.rmtree(WORK, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
