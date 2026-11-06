#!/usr/bin/env python3
"""THM (HOES) derivatives: the adjoint against central differences.

The model Jacobian of a THM point comes from the reverse-mode adjoint of the
HOES cross section (THMMatrixFunc::THMPointAdjoint, used by ComputeTHMRows);
AZURE2_THM_FD=1 restores the central differences it replaced.  Checked here on
tests/18O_p_a_thm (two 1/2+ levels, Gaussian folding, two segments of one THM
experiment), every level energy and amplitude free, for <thm> blocks that
between them take every branch of the adjoint:

  * lineshape=on (N_C and the width derivatives), vertex=constant;
  * vertex=perlevel (B = S_c(E_lambda), Brune) with entranceL=coherent, a
    spectator-momentum window and the line shape;
  * vertex=onshell;
  * the DW vertex (vertexModel=dw) with a spectator window;
  * an angular window with a spectator window and a coherent background
    (cbackground=, its two columns);
  * a linear coherent background with a fitted background;
  * the Park formalism with the line shape.

For each, the residuals bit for bit, and every residual_jacobian column
analytic vs differences, max |diff| / max |column| <= 1e-5; a column beyond
that (the difference's own truncation error: a coherent-background column, a
narrow level's energy) against a Richardson-extrapolated difference of the
residuals (h, h/2, h/4) to 1e-6.  model_gradients to 1e-4 of the differences.

Needs numpy and the compiled engine; exits 77 (skipped) without them.

Run from anywhere:  python3 tests/pyazr/thm_adjoint_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
KIN = "beam=18O target=3He spectator=d Ebeam=115"

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

CASES = {
    "lineshape": (f"experiment[A] segments=1,2 {KIN} lineshape=on", {}),
    "perlevel+coherentL+window+lineshape": (
        f"entranceL=coherent\nvertex=perlevel\n"
        f"experiment[A] segments=1,2 {KIN} ps=hulthen:0-40 psNodes=6 lineshape=on", {}),
    "onshell": (f"vertex=onshell\nexperiment[A] segments=1,2 {KIN}", {}),
    "dw+window": (f"experiment[A] segments=1,2 {KIN} ps=hulthen:0-30 psNodes=8 "
                  f"distortion=coulomb vertexModel=dw", {}),
    "angle+window+cbackground": (
        f"experiment[A] segments=1,2 {KIN} ps=hulthen:0-40 psNodes=6 theta=50-70 "
        f"cbackground=1/2+:2=0.3,-0.2", {}),
    "cbackground linear+background": (
        "experiment[A] segments=1,2 background=linear cbackground=1/2+:2:linear=0.25,-0.5,0.1,1e-3", {}),
    "park+lineshape": (f"experiment[A] segments=1,2 {KIN} lineshape=on", {"use_park": True}),
}


def stage(tmp, name, block):
    """The project with every level energy and amplitude free and BLOCK as
    its <thm> block (the folding kept)."""
    work = os.path.join(tmp, name.replace(" ", "_").replace("+", "_"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    out, levels = [], False
    with open(os.path.join(SOURCE, "18O_p_a_thm.azr")) as f:
        for line in f.read().replace("\r", "").split("\n"):
            if line.startswith("<levels>"):
                levels = True
            elif line.startswith("</levels>"):
                levels = False
            t = line.split()
            if levels and len(t) > 30:
                t[3], t[10] = "0", "0"
                line = " ".join(t)
            out.append(line)
    out.append(f"<thm>\n{block}\n</thm>\n")
    path = os.path.join(work, "run.azr")
    with open(path, "w") as f:
        f.write("\n".join(out))
    return path


def jacobians(m, x, fd):
    if fd:
        os.environ["AZURE2_THM_FD"] = "1"
    else:
        os.environ.pop("AZURE2_THM_FD", None)
    try:
        r, J = m.residual_jacobian(x)
        G = np.vstack(m.model_gradients(x))
    finally:
        os.environ.pop("AZURE2_THM_FD", None)
    return r, J, G


def colrel(a, b):
    scale = np.maximum(np.max(np.abs(a), axis=0), np.max(np.abs(b), axis=0))
    return np.max(np.abs(a - b), axis=0) / (scale + 1e-300)


tmp = tempfile.mkdtemp(prefix="thm_adjoint_")
try:
    for name, (block, kw) in CASES.items():
        print(name)
        path = stage(tmp, name, block)
        m = azure2(path, **kw)
        x = np.asarray(m.params_rwa, float)
        rA, JA, GA = jacobians(m, x, False)
        rF, JF, GF = jacobians(m, x, True)
        check("residuals identical in both modes", np.array_equal(rA, rF))
        check("free columns", JA.shape[1] >= 6, f"{JA.shape}")
        rel = colrel(JA, JF)
        relg = colrel(GA, GF)
        check("model_gradients = central differences (every column, 1e-4)",
              relg.max() <= 1e-4, f"max rel {relg.max():.2e} at column {relg.argmax()}")
        # A column off the differences by more than 1e-5 is re-checked against
        # a Richardson-extrapolated difference of the residuals (h, h/2, h/4):
        # the step 1e-6 (|x| + 1) leaves a truncation error of that size on a
        # quadratic coherent-background column or a narrow level's energy.
        worst = 0.0
        for c in np.where(rel > 1e-5)[0]:
            def d(h):
                xp, xm = x.copy(), x.copy()
                xp[c] += h
                xm[c] -= h
                return (m.residuals(xp) - m.residuals(xm)) / (2.0 * h)
            h = 1e-6 * (abs(x[c]) + 1.0)
            d1, d2, d4 = d(h), d(h / 2), d(h / 4)
            R = (16.0 * (4.0 * d4 - d2) / 3.0 - (4.0 * d2 - d1) / 3.0) / 15.0
            worst = max(worst, np.max(np.abs(JA[:, c] - R)) / (np.max(np.abs(R)) + 1e-300))
            rel[c] = 0.0
        check("residual_jacobian = central differences (every column, 1e-5)",
              rel.max() <= 1e-5, f"max rel {rel.max():.2e} at column {rel.argmax()}")
        check("columns beyond: = Richardson differences (1e-6)", worst <= 1e-6, f"max rel {worst:.2e}")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\nFAIL: " + ", ".join(failures) if failures else "\nPASS")
sys.exit(1 if failures else 0)
