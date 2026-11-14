#!/usr/bin/env python3
"""The independent integral of the thick layer (README.md, "(K') thick layer against an independent integral"):
models/thick_reference.json.

    python3 -I tests/layer_yield/make_thick_reference.py <AZURE2 binary> [--write]

Y(E_0) = 1e-24 int_{E_out}^{E_0} sigma(f E) / eps(E) dE for the data set `thick` of models/layer_yield.azr (target TK:
N = 1.1e19 atoms/cm^2, eps = a0 x^a1 lab MeV cm^2 per 15N, x = E/MeV), with
  - E_out from the closed form (C2) of README.md (no quadrature),
  - sigma the engine's own bare 15N(p,a)12C cross section (b, c.m.), computed point by point -- no fold, no sub-point
    grid -- by a menu-1 calculation of one `quantity sigma` data set on the integration nodes (the levels of
    models/layer_yield.azr),
  - the trapezoid rule on nodes every 50 eV (lab) from below the lowest E_out to the highest E_0, and every Gamma_lab/50
    (2.1 eV) within 200 Gamma of the resonance, with E_out and E_0 of each point added as nodes.
The trapezoid rule on these nodes is tested by the same integral on every second node: the two differ by the size
printed as `trapezoid_check` (Richardson: the error of the fine rule is about a third of it).

With --write, models/thick_reference.json is written: the three E_0, the reference yields, the node count, the
trapezoid check, and the engine's values and relative deviation today (from the calculation of models/layer_yield.azr
by the binary given), whose largest value is the bound check.py holds the engine to.  Without it the numbers are
printed only.  Only the standard library.
"""
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
sys.path.insert(0, HERE)
import check  # noqa: E402  (blocks, f_cm, exit_pow: the closed forms of check.py)


def ceil3(x):
    """x rounded up to three significant digits (the bound is never below today's deviation)."""
    scale = 10.0 ** (math.floor(math.log10(x)) - 2)
    return math.ceil(x / scale) * scale


def project_with(azr_text, data_block):
    """The model file with its <data> section replaced by `data_block` and no extrapolation."""
    head = azr_text[:azr_text.index("<data>")]
    tail = azr_text[azr_text.index("</data>") + len("</data>"):]
    return head + "<data>\n" + data_block + "</data>" + tail


def run_menu1(binary, workdir, azr_name):
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    os.makedirs(os.path.join(workdir, "checks"), exist_ok=True)
    p = subprocess.run([binary, "calculate", "--no-gui", azr_name], stdin=subprocess.DEVNULL, cwd=workdir,
                       capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        sys.exit(f"the calculation in {workdir} failed:\n{p.stdout[-2000:]}{p.stderr[-2000:]}")


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    binary = os.path.abspath(sys.argv[1])
    write = "--write" in sys.argv[2:]
    truth = json.load(open(os.path.join(MODELS, "truth.json")))
    pN = truth["levels"]["p15N"]
    L = truth["levels"]["level"]
    a0, a1 = truth["eps_pow"]["a0"], truth["eps_pow"]["a1"]
    S = truth["sets"]["thick"]
    N, E0s = S["N"], S["E0"]
    f = check.f_cm(pN)
    E_outs = [check.exit_pow(E0, N, a0, a1) for E0 in E0s]

    # The nodes (lab MeV): 50 eV everywhere, Gamma_lab / 50 within 200 Gamma of the resonance, the limits.
    gamma_lab = (L["Gp_eV"] + L["Ga_eV"]) * 1e-6 / f
    er_lab = L["Er_cm"] / f
    lo, hi = min(E_outs) - 1e-4, max(E0s) + 1e-4
    nodes = set()
    n = int(round((hi - lo) / 50e-6))
    nodes.update(round(lo + k * (hi - lo) / n, 12) for k in range(n + 1))
    h = gamma_lab / 50.0
    m = int(200 * gamma_lab / h)
    nodes.update(round(er_lab + k * h, 12) for k in range(-m, m + 1))
    nodes.update(E_outs)
    nodes.update(E0s)
    nodes = sorted(nodes)

    azr = open(os.path.join(MODELS, "layer_yield.azr")).read()
    work = tempfile.mkdtemp(prefix="thick_reference.")
    try:
        os.makedirs(os.path.join(work, "data"))
        with open(os.path.join(work, "data", "nodes.dat"), "w") as fh:
            fh.write("# E_lab(MeV) unused value error: the integration nodes (make_thick_reference.py)\n")
            for e in nodes:
                fh.write(f"{e!r} 0 1 1\n")
        block = ("dataset nodes\n    reaction       15N(p,a)12C\n    quantity       sigma\n    unit           b\n"
                 "    file           data/nodes.dat\n    columns        E_lab[MeV] _ value error\n    norm           1 fixed\n")
        with open(os.path.join(work, "nodes.azr"), "w") as fh:
            fh.write(project_with(azr, block))
        run_menu1(binary, work, "nodes.azr")
        rows = check.blocks(os.path.join(work, "output", "AZUREOut_aa=1_R=2.out"))[0]
        if len(rows) != len(nodes):
            sys.exit(f"{len(rows)} values for {len(nodes)} nodes")
        sigma = [r[3] for r in rows]

        def integral(E_out, E0, step=1):
            idx = [i for i, e in enumerate(nodes) if E_out <= e <= E0]
            keep = idx[::step]
            if keep[-1] != idx[-1]:
                keep.append(idx[-1])
            g = [sigma[i] / (a0 * nodes[i] ** a1) for i in keep]
            return 1e-24 * sum(0.5 * (g[k] + g[k + 1]) * (nodes[keep[k + 1]] - nodes[keep[k]]) for k in range(len(keep) - 1))

        ref = [integral(eo, e0) for eo, e0 in zip(E_outs, E0s)]
        coarse = [integral(eo, e0, 2) for eo, e0 in zip(E_outs, E0s)]
        trap = max(abs(c / r - 1) for c, r in zip(coarse, ref))

        # The engine today: the data set `thick` of the model as it is.
        shutil.copytree(os.path.join(MODELS, "data"), os.path.join(work, "model", "data"))
        shutil.copy(os.path.join(MODELS, "layer_yield.azr"), os.path.join(work, "model"))
        run_menu1(binary, os.path.join(work, "model"), "layer_yield.azr")
        pb = check.blocks(os.path.join(work, "model", "output", "AZUREOut_aa=1_R=2.out"))
        engine = [r[3] for r in pb[3]]  # thin, thin_sigma, thin_bare, thick, ...
    finally:
        shutil.rmtree(work, ignore_errors=True)
    dev = [e / r - 1 for e, r in zip(engine, ref)]
    for E0, Eo, r, e, d in zip(E0s, E_outs, ref, engine, dev):
        print(f"E0 {E0:.4f} MeV  E_out {Eo:.10f}  reference {r:.10e}  engine {e:.10e}  engine/reference - 1 {d:+.3e}")
    print(f"{len(nodes)} nodes; trapezoid on every second node differs by {trap:.2e}")
    if write:
        out = {
            "doc": "make_thick_reference.py: Y = 1e-24 int sigma/eps dE of the data set `thick` (target TK), the trapezoid "
                   "rule on the engine's bare sigma at the nodes, E_out in closed form; engine_today and bound: the "
                   "engine's deviation when this file was written (check.py fails above the bound)",
            "E0": E0s, "E_out": E_outs, "Y_reference": ref, "nodes": len(nodes), "trapezoid_check": trap,
            "engine_today": engine, "deviation_today": dev,
            "bound": ceil3(max(abs(d) for d in dev)),
        }
        with open(os.path.join(MODELS, "thick_reference.json"), "w") as fh:
            json.dump(out, fh, indent=1)
            fh.write("\n")
        print("wrote models/thick_reference.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
