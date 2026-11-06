#!/usr/bin/env python3
"""save_fit reproduces a fit in which a width is fixed while its level's other
widths move.

A width marked fixed holds its reduced-width amplitude during a fit (in pyazr
as in the CLI's MIGRAD and LM), so its physical value -- the partial width a
<levels> line holds -- follows the level's other widths through the Brune
denominator 1 + sum_c gamma_c^2 dS_c/dE.  save_fit wrote only the free
parameters, so the fixed width kept the value it was read with and the
snapshot reloaded as a different model: 4 % off in chi2 here.

Model: tests/15N_p_a with the alpha width of the 1- level at 13.1 MeV fixed
and its proton width free.  The proton width amplitude is scaled by 1.5; then

  1. save_fit writes a snapshot that reloads with the same chi2 (rel 1e-9)
     and the same amplitude for the fixed width (rel 1e-9);
  2. the fixed width in the snapshot is the physical value the fit had,
     transform_all_rwa(..., include_fixed=True), not the one read in (which
     it differs from by several per cent here);
  3. transform_all_rwa without include_fixed still returns the free
     R-matrix parameters only.

Run from anywhere:  python3 tests/pyazr/fixed_width_roundtrip_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "15N_p_a")
LEVEL_E = "13.104982826821646"

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
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "15N")
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    # The 13.1 MeV level: alpha channel (pair 2) fixed, proton channels free;
    # no <parameterSettings>, so the flags in <levels> are what count.
    out, skip = [], False
    for line in open(os.path.join(SOURCE, "15N_p_a.azr")).read().split("\n"):
        tag = line.strip()
        if tag == "<parameterSettings>":
            skip = True
            out.append(line)
            continue
        if tag == "</parameterSettings>":
            skip = False
        if skip:
            continue
        f = line.split()
        if len(f) > 30 and f[2] == LEVEL_E:
            f[10] = "1" if f[5] == "2" else "0"
            line = " ".join(f)
        out.append(line)
    azr = os.path.join(work, "fx.azr")
    with open(azr, "w") as fh:
        fh.write("\n".join(out))

    with azure2(azr) as m:
        ps = m.parameters
        free = [p.name for p in ps if not p.fixed]
        names = [p.name for p in ps]
        fixed_name = "width_2_3"          # the 13.1 MeV level's alpha width
        k = names.index(fixed_name)
        check("the alpha width is fixed and the proton width free",
              ps[k].fixed and "width_2_1" in free, [(p.name, p.fixed) for p in ps[4:8]])
        x = np.asarray(m.params_rwa, float).copy()
        x[free.index("width_2_1")] *= 1.5
        chi2 = float(np.sum(m.calculate_chi2_rwa(x)))
        allx = m._all_rwa(x)
        phys_all = np.asarray(m.transform_all_rwa(allx, include_fixed=True), float)
        phys_free = np.asarray(m.transform_all_rwa(allx), float)
        read_in = float(np.asarray(m.transform_all_rwa(
            np.asarray(m.sess.params_all_rwa(), float), include_fixed=True))[k])
        check("3. transform_all_rwa returns the free R-matrix parameters by default",
              phys_free.size == m.n_rmatrix, (phys_free.size, m.n_rmatrix))
        print(f"  {fixed_name}: {read_in:.6g} eV as read, {phys_all[k]:.6g} eV "
              f"with width_2_1 x 1.5 (fixed amplitude {allx[k]:.6g})")
        saved, _ = m.save_fit(os.path.join(work, "saved.azr"), x)

    with azure2(saved) as s:
        chi2_back = float(np.sum(s.calculate_chi2_rwa(np.asarray(s.params_rwa, float))))
        amp_back = float(np.asarray(s.sess.params_all_rwa(), float)[k])
        written = float(np.asarray(s.transform_all_rwa(
            np.asarray(s.sess.params_all_rwa(), float), include_fixed=True))[k])
    print(f"  chi2 {chi2:.10g} before save_fit, {chi2_back:.10g} reloaded")
    check("1. the snapshot reloads with the fit's chi2 (rel 1e-9)",
          rel(chi2_back, chi2) < 1e-9, f"{chi2_back} vs {chi2}")
    check("1. the fixed width keeps its amplitude (rel 1e-9)",
          rel(amp_back, allx[k]) < 1e-9, f"{amp_back} vs {allx[k]}")
    check("2. the snapshot holds the fit's physical value of the fixed width",
          rel(written, phys_all[k]) < 1e-9 and rel(phys_all[k], read_in) > 1e-3,
          f"written {written}, fit {phys_all[k]}, read in {read_in}")

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: save_fit reproduces a fit with fixed widths")
