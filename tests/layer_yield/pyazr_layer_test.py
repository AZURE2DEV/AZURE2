#!/usr/bin/env python3
"""pyazr: a yield on a layer built as records and run (docs/upgrade/energy_grid.md section 5; README.md of this directory).

  1. Project.new on the nuclear part of models/layer_yield.azr, then add_target (kind layer, the areal density fitted
     with a 10 % prior, a constant stopping) and add_dataset (quantity yield, 1/particle) across the narrow resonance,
     and a bare sigma data set at the back of the layer; the project validates and its session runs;
  2. the yield at full precision equals the arctan form of the Breit-Wigner integral (README.md eq. R3, computed here)
     to 3e-4 of the plateau (the R-matrix line shape; check.py has the same limit from the CLI's output);
  3. dY/dN, the density column of objective().residuals_and_jacobian (the engine's central differences of the fitted areal density,
     data 0 and errors 1 so that the column is dY/dN), equals 1e-24 sigma(E_out) with sigma the engine's bare cross
     section at E_out (README.md eq. D1): to 1e-4 within 10 widths of the resonance, 3e-3 beyond (the kernel's linear
     interpolation of sigma between its sub-points);
  4. the same records written by the CLI (menu 1) give the same model values as the session to 1e-10.

Needs the compiled module (PYAZR_ENGINE_DIR) and AZURE2_BIN for item 4.  Exit 77 (skip) without the module.
    PYAZR_ENGINE_DIR=<build>/pyazr-out AZURE2_BIN=<build>/src/AZURE2 python3 tests/layer_yield/pyazr_layer_test.py
"""
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    import numpy as np
    from pyazr import Project
    import pyazr.described  # noqa: F401
    Project  # the compiled module is loaded on first use
except Exception as e:  # pragma: no cover
    print(f"skip: pyazr not importable ({e})")
    sys.exit(77)

sys.path.insert(0, HERE)
import check as closed  # noqa: E402  (the closed forms of check.py, which read nothing from the engine)

failures = []


def ok(name, good, detail=""):
    print(f"  {'ok  ' if good else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not good:
        failures.append(name)


def main():
    truth = json.load(open(os.path.join(HERE, "models", "truth.json")))
    lv = truth["levels"]
    pN, L = lv["p15N"], lv["level"]
    eps, N = truth["eps_const"], 1e18
    Delta = N * eps
    Er_lab = truth["Er_lab"]
    # across the front and the back of the layer, every 0.1 keV within 0.5 keV of E_out = E_r
    E0 = sorted(set([Er_lab - 0.001 + 0.0005 * i for i in range(17)] +
                    [Er_lab + Delta + 0.0001 * i for i in range(-5, 6)]))
    Eout = [e - Delta for e in E0]
    tmp = tempfile.mkdtemp(prefix="layer_pyazr.")
    try:
        os.makedirs(os.path.join(tmp, "data"))
        with open(os.path.join(tmp, "data", "scan.dat"), "w") as f:
            f.writelines(f"{e:.12g} 0 0 1\n" for e in E0)  # data 0: r_i = m_i, so the differences keep their digits
        with open(os.path.join(tmp, "data", "back.dat"), "w") as f:
            f.writelines(f"{e:.12g} 0 1 1\n" for e in Eout)
        p = Project.new(os.path.join(tmp, "p.azr"), nuclear_part=os.path.join(HERE, "models", "layer_yield.azr"))
        p.add_target("TP", kind="layer", active="15N", areal_density="1e18 atoms/cm2 fit prior 1e18 10%",
                     stopping=f'formula "a0" a0={eps:g} MeV*cm2/atom', straggling="none")
        p.add_dataset("scan", reaction="15N(p,a)12C", quantity="yield", unit="1/particle", file="data/scan.dat",
                      columns="E_lab[MeV] _ value error", norm="1 fixed", target="TP", quadrature="points 40")
        p.add_dataset("back", reaction="15N(p,a)12C", quantity="sigma", unit="b", file="data/back.dat",
                      columns="E_lab[MeV] _ value error", norm="1 fixed")
        ok("1. the records validate", p.validate() == [], str(p.validate()))
        p.save()
        for d in ("output", "checks"):
            os.makedirs(os.path.join(tmp, d), exist_ok=True)
        print("1. session")
        s = p.session()
        x = np.asarray(s.parameters.x0, float)
        vals = [d.model for d in s.evaluate(x).data]
        y, back = np.asarray(vals[0], float), np.asarray(vals[1], float)
        ok("1. a session of the project runs (two data sets)", len(vals) == 2 and len(y) == len(E0))

        print("2. the arctan form")
        worst, plateau = 0.0, None
        for e, v in zip(E0, y):
            c, plateau = closed.resonance_yield(pN, L, eps, e, Delta)
            worst = max(worst, abs(v - c) / plateau)
        ok("2. Y = arctan form to 3e-4 of the plateau", worst <= 3e-4, f"{worst:.2e} of {plateau:.6e}")

        print("3. dY/dN")
        dens = [q for q in s.parameters if q.kind == "density" and not q.fixed]
        ok("3. one free areal density (TP.areal_density)", len(dens) == 1 and dens[0].name == "TP.areal_density",
           str([q.name for q in dens]))
        r, J = s.objective().residuals_and_jacobian(x)
        col = dens[0].free_index
        # r_i = (m_i - d_i)/e_i with n = 1 and e_i = 1: dY_i/dN = J_i,col
        dydn = J[:len(E0), col]
        # The kernel interpolates sigma linearly between its sub-points, so the derivative of the yield it computes
        # is 1e-24 sigma_lin(E_out): within 20 widths of the resonance the nodes are Gamma/50 apart and
        # h^2 sigma''/8 <= 1e-4 of sigma; in the geometric tail of the grid beyond, up to 3e-3.
        G = (L["Gp_eV"] + L["Ga_eV"]) * 1e-6
        worst_core, worst_tail = 0.0, 0.0
        for d, sg, e in zip(dydn, back, Eout):
            r = abs(d / (1e-24 * sg) - 1)
            if abs(closed.f_cm(pN) * e - L["Er_cm"]) <= 10 * G:
                worst_core = max(worst_core, r)
            else:
                worst_tail = max(worst_tail, r)
        ok("3. dY/dN (engine, central differences) = 1e-24 sigma(E_out), |E_out - E_r| <= 10 Gamma, to 1e-4",
           worst_core <= 1e-4, f"worst {worst_core:.2e}")
        ok("3. the same beyond 10 Gamma (the grid's geometric tail), to 3e-3", worst_tail <= 3e-3,
           f"worst {worst_tail:.2e}")
        k = int(np.argmax(back))
        bw = closed.sigma_bw(pN, L, closed.f_cm(pN) * Eout[k])
        print(f"      at the peak: dY/dN = {dydn[k]:.10e}, 1e-24 sigma = {1e-24 * back[k]:.10e}, "
              f"1e-24 sigma_BW = {1e-24 * bw:.10e}")
        s.close()

        print("4. the CLI on the same records")
        bin_ = os.environ.get("AZURE2_BIN")
        if not bin_:
            print("  skip  4. (AZURE2_BIN not set)")
        else:
            subprocess.run([bin_, "calculate", "--no-gui", "p.azr"], stdin=subprocess.DEVNULL, text=True, cwd=tmp,
                           capture_output=True, timeout=600)
            rows = closed.blocks(os.path.join(tmp, "output", "AZUREOut_aa=1_R=2.out"))
            cli = [row[3] for row in rows[0]] if rows else []
            worst = max((abs(a / b - 1) for a, b in zip(cli, y)), default=1.0) if len(cli) == len(y) else 1.0
            ok("4. the CLI's yield = the session's to 1e-10 (the output's 11 digits)", worst <= 1e-10, f"{worst:.2e}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if failures:
        print(f"FAIL pyazr_layer_yield ({len(failures)})")
        return 1
    print("PASS pyazr_layer_yield")
    return 0


if __name__ == "__main__":
    sys.exit(main())
