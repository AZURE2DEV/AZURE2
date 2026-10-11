#!/usr/bin/env python3
"""Figure for DeBoer's blur idea (2026-10-11): data, full fit, broad-only model (levels < 50 keV switched off) and the same three
blurred (Gaussian kernel regression, sigma D) for three excitation functions.  Data drawn x fitted normalization:
n*d = fit*d/(d + r*e) with r = (fit - n d)/(n e) from the dump.  Broad levels (Gamma >= 150 keV) of the incumbent marked on top.
usage: blur_figure.py full.npz broad.npz parameters.out out.png [D_MeV=0.1] --panel "dataset:angle:ex_lo:ex_hi:title" ...
  angle = the detector's c.m. angle (+-1 deg) or "-" for a single-angle segment; e.g.
  --panel "Cierjacks_ntot_err4144fix_merged:-:8.4:13.25:Cierjacks n-total" --panel "ND_2021_HE_prune:0:6.95:10.65:ND 2021 0 deg"
"""
import sys, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from blur_compare import blur, detector_labels

import argparse


def levels(path, gmin=150.0):
    out, cur = [], None
    for line in open(path):
        m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
        if m:
            cur = [float(m.group(2)), 0.0]; out.append(cur); continue
        m = re.search(r"G\s+=\s+(\S+) keV", line)
        if m and cur is not None:
            cur[1] += abs(float(m.group(1)))
    return [(e, g) for e, g in out if g >= gmin]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("full"); ap.add_argument("broad"); ap.add_argument("par"); ap.add_argument("out")
    ap.add_argument("D", nargs="?", type=float, default=0.1); ap.add_argument("--panel", action="append", required=True)
    a = ap.parse_args()
    F, B = np.load(a.full), np.load(a.broad); L = levels(a.par); out, D = a.out, a.D
    PANELS = []
    for p in a.panel:
        nm, an, lo, hi, title = p.split(":", 4)
        PANELS.append((nm, None if an == "-" else float(an), title, (float(lo), float(hi))))
    fig, axes = plt.subplots(len(PANELS), 1, figsize=(10, 3.3 * len(PANELS)), squeeze=False); axes = axes[:, 0]
    for ax, (nm, ang, title, xr) in zip(axes, PANELS):
        m = np.where(F["name"] == nm)[0]
        if ang is not None:
            lab = np.full(len(F["ex"]), np.nan); lab[m] = detector_labels(F["ang"][m])
            m = m[np.abs(lab[m] - ang) < 1.0]
        m = m[np.argsort(F["ex"][m])]
        keep = np.concatenate([[True], np.diff(F["ex"][m]) > 1e-7]); m = m[keep]
        e, d, er, r, yf, yb = F["ex"][m], F["dat"][m], F["err"][m], F["res"][m], F["fit"][m], B["fit"][m]
        nd = yf * d / (d + r * er)
        bd, bf, bb = blur(e, nd, D), blur(e, yf, D), blur(e, yb, D)
        edge = (e >= e[0] + 1.5 * D) & (e <= e[-1] - 1.5 * D)
        ax.plot(e, nd, ".", ms=2, color="0.65", label="data × norm", zorder=1)
        ax.plot(e, yf, "-", lw=0.7, color="tab:blue", alpha=0.55, label="full fit", zorder=2)
        ax.plot(e, yb, "-", lw=0.9, color="tab:orange", alpha=0.75, label="broad-only model (Γ < 50 keV off)", zorder=2)
        ax.plot(e[edge], bd[edge], "-", lw=2.4, color="k", label=f"blurred data (σ = {1000*D:.0f} keV)", zorder=4)
        ax.plot(e[edge], bf[edge], "--", lw=2.0, color="tab:blue", label="blurred full fit", zorder=5)
        ax.plot(e[edge], bb[edge], ":", lw=2.2, color="tab:red", label="blurred broad-only", zorder=5)
        lo, hi = np.percentile(np.concatenate([nd, yf, yb]), [0.5, 99.7])
        ax.set_ylim(min(0, lo), hi * 1.12); ax.set_xlim(*xr)
        for el, g in L:
            if xr[0] <= el <= xr[1]:
                ax.plot([el, el], [hi * 1.04, hi * 1.11], "-", color="tab:purple", lw=1.0)
        ax.set_title(title, fontsize=10, loc="left"); ax.set_ylabel("σ (b or b/sr, c.m.)")
        ax.tick_params(labelsize=8)
    axes[-1].set_xlabel("Ex (MeV)   —   purple ticks: incumbent levels with Γ ≥ 150 keV")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, ncol=6, loc="upper center", frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.965)); fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
