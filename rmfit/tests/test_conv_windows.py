#!/usr/bin/env python3
"""<targetInt> restricted to lab-energy windows (AZURE2 TargetEffect ranges / blend width / auto tolerance) in rmfit's sharded
convolution (model.parse_target_effect, TargetEffectSpec.blend_weight, build_conv_spec -> _windowed_conv_spec), 2026-10-15.
No engine.  Run:  python3 rmfit/tests/test_conv_windows.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3: pyazr 2 only")
    sys.exit(77)
from rmfit.model import build_conv_spec, parse_target_effect  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


HEAD = '1  "11,12"  10  1  0.00000  0  0  ""  0  0  0  1  "sqrt(a0+a1*x^3)"  2  2.0e-06  1.0e-07  0  0.04  20  50'


def smoothstep_weight(x, lo, hi, bw):
    up = min(1.0, max(0.0, (x - (lo - 0.5 * bw)) / bw)); dn = min(1.0, max(0.0, ((hi + 0.5 * bw) - x) / bw))
    return up * up * (3 - 2 * up) * dn * dn * (3 - 2 * dn)


def main():
    # -- parsing ------------------------------------------------------------------------------------------------
    old = parse_target_effect(HEAD + "  sigmalab")
    check("old-format line: no ranges, blend 0, tolerance 0", old.ranges == [] and old.blend_width == 0 and old.auto_tolerance == 0
          and old.sigma_lab and old.conv_coeffs == [2.0e-06, 1.0e-07])
    w = parse_target_effect(HEAD + '  "2.0-2.1,3.0000-3.2500"  0.002  sigmalab')
    check("ranges + blend parsed, sigmalab still seen", w.ranges == [(2.0, 2.1), (3.0, 3.25)] and w.blend_width == 0.002
          and w.auto_tolerance == 0 and w.sigma_lab, (w.ranges, w.blend_width, w.auto_tolerance, w.sigma_lab))
    t = parse_target_effect(HEAD + '  "1.0-2.0"  0  0.001')
    check("blend 0 then a tolerance (engine reads both)", t.ranges == [(1.0, 2.0)] and t.blend_width == 0 and t.auto_tolerance == 0.001,
          (t.ranges, t.blend_width, t.auto_tolerance))
    e = parse_target_effect(HEAD + '  ""  sigmalab')
    check('empty "" = unrestricted', e.ranges == [] and e.sigma_lab)
    real = open(os.path.join(HERE, "..", "..", "..", "w_new.azr")).read() if os.path.exists(os.path.join(HERE, "..", "..", "..", "w_new.azr")) else None
    if real:
        lines = [l for l in real.split("\n") if "sigmalab" in l and '"11,12,' in l]
        r = parse_target_effect(lines[0])
        check("archive targetInt line: 10 windows, blend 2 keV, 6 coefficients", len(r.ranges) == 10 and r.blend_width == 0.002
              and len(r.conv_coeffs) == 6, (len(r.ranges), r.blend_width))
    # -- blend weight = TargetEffect::BlendWeight ---------------------------------------------------------------
    xs = np.array([1.95, 1.999, 2.0, 2.001, 2.05, 2.1, 2.2, 3.1, 3.25, 3.2505])
    ref = np.array([max(smoothstep_weight(x, 2.0, 2.1, 0.002), smoothstep_weight(x, 3.0, 3.25, 0.002)) for x in xs])
    check("blend weight = engine smoothstep (max over windows)", np.allclose(w.blend_weight(xs), ref, atol=1e-15), (w.blend_weight(xs), ref))
    check("edge weight 0.5, inside 1, outside 0", np.isclose(w.blend_weight(2.0)[0], 0.5) and w.blend_weight(2.05)[0] == 1
          and w.blend_weight(1.95)[0] == 0)
    hard = parse_target_effect(HEAD + '  "2.0-2.1"  0')
    check("blend 0: hard window, edges inclusive", list(hard.blend_weight([1.9999, 2.0, 2.1, 2.1001])) == [0, 1, 1, 0])
    # -- windowed kernel -------------------------------------------------------------------------------------------
    cmc = 13.0 / 17.0
    e_lab = np.round(np.arange(1.90, 3.40, 0.005), 6)
    rows = np.column_stack([e_lab, np.zeros_like(e_lab), np.ones_like(e_lab), 0.1 * np.ones_like(e_lab)])
    sp = build_conv_spec(None, w, rows, cmc)
    b = w.blend_weight(e_lab); e_cm = e_lab * cmc
    out, full, edge = b == 0, b == 1, (b > 0) & (b < 1)
    f = lambda g: np.sin(40 * g) + 2.0 + 0.3 * g
    folded = sp.W @ f(sp.grid_cm)
    check("every row normalised (a flat cross section folds to itself)", np.allclose(sp.W.sum(axis=1), 1.0, atol=1e-12))
    check("points outside every window = bare value at the point, exactly", np.allclose(folded[out], f(e_cm[out]), rtol=0, atol=1e-13),
          np.max(np.abs(folded[out] - f(e_cm[out]))))
    check("every data energy is a grid node", np.all(np.isin(e_cm, sp.grid_cm)))
    s = np.sqrt(2.0e-06 + 1.0e-07 * e_lab ** 3) * cmc
    trunc = np.sqrt(1 - 2 * 3 * np.exp(-4.5) / np.sqrt(2 * np.pi) / 0.9973002)
    rms = np.sqrt(np.array([np.sum(sp.W[i] * (sp.grid_cm - e_cm[i]) ** 2) for i in np.nonzero(full)[0]]))
    check("inside: Gaussian of sigma_lab(E_lab) x cm_conv (truncated at 3 sigma)", np.allclose(rms, s[full] * trunc, rtol=3e-3),
          np.max(np.abs(rms / (s[full] * trunc) - 1)))
    lin = sp.W @ (3 * sp.grid_cm)
    check("inside: a linear cross section folds to its value", np.allclose(lin[full], 3 * e_cm[full], rtol=0, atol=3e-6))
    gauss_only = []
    for i in np.nonzero(edge)[0]:
        gi = (sp.W[i].copy()); k = np.searchsorted(sp.grid_cm, e_cm[i]); gi[k] -= (1 - b[i]); gauss_only.append(gi / b[i])
    ok_edge = all(np.isclose(np.sum(gq), 1.0) for gq in gauss_only) and len(gauss_only) > 0
    check(f"edge points ({int(edge.sum())}): b x normalised Gaussian + (1-b) x delta at the point", ok_edge)
    # grid size: local, not one global step over the whole segment
    n_global = int((e_cm.max() - e_cm.min()) / (s.min() / 10))
    check("grid is local to the windows (much smaller than a global min(sigma)/10 grid)", sp.n_grid < 0.6 * n_global, (sp.n_grid, n_global))
    # no ranges: unchanged path (uniform grid)
    sp0 = build_conv_spec(None, old, rows[:5], cmc)
    check("no ranges: uniform-grid path unchanged", np.allclose(np.diff(sp0.grid_cm), np.diff(sp0.grid_cm)[0]))
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all windowed-convolution checks passed")


if __name__ == "__main__":
    main()
