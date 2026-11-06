#!/usr/bin/env python3
"""Checks for the angular-coefficient data segments (called by check.sh).

  verify.py prepare WORK/extrap DATADIR
      (a) the Legendre projection of the c.m. differential cross section (test segment 2, 0-180 deg
          step 5, least squares on P_0..P_8) reproduces the coefficients AZURE2 writes for test segment 1
          (isDiff 3): a_k = c_k / c_0, k = 1..6, to 1e-6 -- the convention dsigma/dOmega = A_0 [1 + sum a_k P_k];
      (b) c_0 = sigma_int / 4 pi with sigma_int from test segment 3;
      writes DATADIR/coef_norm.dat (E_lab, k, a_k, 0.01), coef_abs.dat (A_k = a_k sigma/4pi, 1 %), and the
      "_off" versions (values a few per cent off the model) for the gradient test, on the four energies
      6.0-6.3 and k = 1..4, and prints the model coefficients at 6.05-6.35 for the shift check.
  verify.py data WORK/data DATADIR
      the fitted model of both data segments equals the data (chi2 ~ 0; AZUREOut fit column = data to 1e-8).
  verify.py shift WORK/shift WORK/extrap
      the normalised segment with a 0.05 MeV shift gives the coefficients of test segment 1 at E + 0.05.
  verify.py conv WORK/conv WORK/conv_diff
      the folded coefficient segments equal the Legendre projection of the folded differential cross
      section (normalised: c_k/c_0; absolute: c_k), to 1e-6 (absolute to 1e-6 relative).
"""
import glob
import os
import sys

import numpy as np
from numpy.polynomial import legendre as L

failures = 0


def ok(msg):
    print("  ok    " + msg)


def bad(msg):
    global failures
    print("  FAIL  " + msg)
    failures += 1


def blocks(path):
    """Blank-line-separated blocks of an AZUREOut file, each an array."""
    out, cur = [], []
    for line in open(path):
        if line.strip():
            cur.append([float(x) for x in line.split()])
        elif cur:
            out.append(np.array(cur)); cur = []
    if cur:
        out.append(np.array(cur))
    return out


def project(diff_rows, kmax=8):
    """{E_cm: c_k} from rows (E_cm, E_x, theta_cm, dsigma/dOmega, ...) by least squares on P_0..P_kmax."""
    out = {}
    for e in np.unique(diff_rows[:, 0]):
        r = diff_rows[np.abs(diff_rows[:, 0] - e) < 1e-9]
        x = np.cos(np.radians(r[:, 2]))
        V = np.stack([L.Legendre.basis(k)(x) for k in range(kmax + 1)], axis=1)
        c, *_ = np.linalg.lstsq(V, r[:, 3], rcond=None)
        out[round(float(e), 9)] = c
    return out


def extrap_file(work):
    f = glob.glob(os.path.join(work, "output", "AZUREOut_aa=2_R=5.extrap"))
    return f[0] if f else None


def acoeff_file(work):
    """The angular-distribution test segment's own file: E_cm, a_0, a_1, ..."""
    f = glob.glob(os.path.join(work, "output", "AZUREOut_aa=2_R=5.acoeff"))
    return f[0] if f else None


def prepare(work, datadir):
    b = blocks(extrap_file(work))
    a = blocks(acoeff_file(work)) if acoeff_file(work) else []
    if len(b) < 2 or len(a) < 1:
        bad(f"expected the .acoeff file and 2 .extrap blocks, found {len(a)} and {len(b)}"); return
    coef, diff, integ = a[0], b[0], b[1]        # T1: E_cm, a_0..a_6 ; T2: differential ; T3: E, E_x, 0, sigma, S
    proj = project(diff)
    worst_a = worst_0 = 0.0
    for row, srow in zip(coef, integ):
        e = round(float(row[0]), 9)
        c = proj[e]
        a_engine = row[1:]                        # a_0 .. a_6 as written (a_0 should be 1)
        a_proj = c / c[0]
        worst_a = max(worst_a, np.max(np.abs(a_engine[1:7] - a_proj[1:7])))
        worst_0 = max(worst_0, abs(c[0] - srow[3] / (4 * np.pi)) / (srow[3] / (4 * np.pi)))
        if abs(a_engine[0] - 1.0) > 1e-9:
            bad(f"a_0 = {a_engine[0]} at E_cm {e}, expected 1"); return
    if worst_a <= 1e-6:
        ok(f"(a) isDiff 3 coefficients = Legendre projection of the c.m. differential cross section (max |diff| {worst_a:.1e}, {len(coef)} energies, k = 1-6)")
    else:
        bad(f"(a) isDiff 3 coefficients vs the projection of dsigma/dOmega: max |diff| {worst_a:.1e}")
    if worst_0 <= 1e-6:
        ok(f"(b) c_0 = sigma_int / 4 pi (max rel {worst_0:.1e})")
    else:
        bad(f"(b) c_0 vs sigma_int / 4 pi: max rel {worst_0:.1e}")
    big = max(abs(float(v)) for row in coef for v in row[2:5])
    if big < 0.05:
        bad(f"the test is vacuous: |a_k| < 0.05 everywhere (max {big:.3f})")
    # the data files: lab energies of the test grid are E_cm (m1+m2)/m2; T1's first column is E_cm -> use the
    # lab energies the segment was given (5.9 + 0.05 i); pick 6.0, 6.1, 6.2, 6.3 (indices 2..5)
    os.makedirs(datadir, exist_ok=True)
    e_lab = 5.9 + 0.05 * np.arange(len(coef))
    rng = np.random.default_rng(3)
    with open(os.path.join(datadir, "coef_norm.dat"), "w") as fn, open(os.path.join(datadir, "coef_abs.dat"), "w") as fa, \
         open(os.path.join(datadir, "coef_norm_off.dat"), "w") as fno, open(os.path.join(datadir, "coef_abs_off.dat"), "w") as fao:
        for i in range(2, 6):
            a0 = integ[i, 3] / (4 * np.pi)
            for k in range(1, 5):
                ak = coef[i, 1 + k]
                fn.write(f"{e_lab[i]:.6f}  {k}  {ak:.12e}  0.01\n")
                fa.write(f"{e_lab[i]:.6f}  {k}  {ak * a0:.12e}  {0.01 * abs(ak * a0) + 1e-12:.6e}\n")
                off = 1.0 + 0.03 * rng.standard_normal()
                fno.write(f"{e_lab[i]:.6f}  {k}  {ak * off:.12e}  0.01\n")
                fao.write(f"{e_lab[i]:.6f}  {k}  {ak * a0 * off:.12e}  {0.01 * abs(ak * a0) + 1e-12:.6e}\n")
    ok(f"data files written: 16 normalised and 16 absolute coefficients (|a_k| up to {big:.3f})")


def out_rows(work):
    f = glob.glob(os.path.join(work, "output", "AZUREOut_aa=2_R=5.out"))
    return blocks(f[0]) if f else []


def data(work, datadir):
    b = out_rows(work)
    if len(b) != 2:
        bad(f"data: expected 2 output blocks, found {len(b)}"); return
    for name, rows in (("normalised", b[0]), ("absolute", b[1])):
        # columns: E_cm, E_x, k, fit, fit*S, data, err, ...
        rel = np.max(np.abs(rows[:, 3] - rows[:, 5]) / np.maximum(np.abs(rows[:, 5]), 1e-300))
        if rel <= 1e-8:
            ok(f"data: {name} segment, fit = data (max rel {rel:.1e}, {len(rows)} rows, orders {sorted(set(rows[:, 2].astype(int)))})")
        else:
            bad(f"data: {name} segment, fit vs data max rel {rel:.1e}")


def shift(work, extrap_work):
    b = out_rows(work)
    coef = blocks(acoeff_file(extrap_work))[0]
    if len(b) != 1:
        bad("shift: expected 1 output block"); return
    rows = b[0]
    # the model at lab E + 0.05 = the extrapolation's next grid point: data rows i = 2..5 -> T1 rows 3..6
    worst = 0.0
    # the output E_cm (column 0) is the shifted c.m. energy; match it against T1's E_cm column
    for r in rows:
        k = int(round(r[2]))
        j = np.argmin(np.abs(coef[:, 0] - r[0]))
        if abs(coef[j, 0] - r[0]) > 1e-6:
            bad(f"shift: shifted E_cm {r[0]} not on the extrapolation grid"); return
        worst = max(worst, abs(r[3] - coef[j, 1 + k]))
    if worst <= 1e-8:
        ok(f"shift: a 0.05 MeV shift moves the coefficients to E + 0.05 (max |diff| {worst:.1e})")
    else:
        bad(f"shift: coefficients under a 0.05 MeV shift vs the extrapolation at E + 0.05: max |diff| {worst:.1e}")


def conv(work, diff_work):
    b = out_rows(work)
    d = blocks(extrap_file(diff_work))
    if len(b) != 2 or len(d) < 1:
        bad("conv: outputs missing"); return
    proj = project(d[0])
    for name, rows, absolute in (("normalised", b[0], False), ("absolute", b[1], True)):
        worst = 0.0
        for r in rows:
            e = round(float(r[0]), 9)
            if e not in proj:
                j = min(proj, key=lambda x: abs(x - e))
                if abs(j - e) > 1e-6:
                    bad(f"conv: energy {e} not in the folded extrapolation"); return
                e = j
            c = proj[e]
            k = int(round(r[2]))
            want = c[k] if absolute else c[k] / c[0]
            # absolute: relative to the largest |c_k| (the odd orders vanish by parity here)
            worst = max(worst, abs(r[3] - want) / (np.max(np.abs(c[1:5])) if absolute else 1.0))
        if worst <= 1e-6:
            ok(f"conv: folded {name} coefficients = projection of the folded differential cross section (max {'rel' if absolute else '|diff|'} {worst:.1e})")
        else:
            bad(f"conv: folded {name} coefficients vs the projection of the folded distribution: {worst:.1e}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    {"prepare": prepare, "data": data, "shift": shift, "conv": conv}[cmd](*sys.argv[2:])
    sys.exit(1 if failures else 0)
