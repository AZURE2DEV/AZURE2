#!/usr/bin/env python3
"""Figures of the natural-neon worked example (docs/source/examples/natural_neon.rst).

    python3 plot_natural_neon.py FIT_OUTPUT BARE_OUTPUT [FIGURE_DIR]

FIT_OUTPUT   output/ of "Calculate With Data" (mode 1) with natural_neon_fit.sav
BARE_OUTPUT  output/ of "Calculate Without Data" (mode 3) on the copy with every <targetInt> line off
FIGURE_DIR   where natural_neon_fit.png and natural_neon_decomposition.png go (default: .)
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

fit_dir, bare_dir = sys.argv[1], sys.argv[2]
fig_dir = sys.argv[3] if len(sys.argv) > 3 else "."


def blocks(path):
    """The blank-line separated blocks of an AZUREOut file, one per segment in key order."""
    out, cur = [], []
    for line in open(path):
        if line.strip() and not line.lstrip().startswith("#"):
            cur.append([float(v) for v in line.split()])
        elif cur:
            out.append(np.array(cur))
            cur = []
    if cur:
        out.append(np.array(cur))
    return out


b2 = blocks(os.path.join(fit_dir, "AZUREOut_aa=2_R=2.out"))   # keys 1-5, 11
b4 = blocks(os.path.join(fit_dir, "AZUREOut_aa=4_R=4.out"))   # keys 8-10
b6 = blocks(os.path.join(fit_dir, "AZUREOut_aa=6_R=6.out"))   # keys 6, 7
assert (len(b2), len(b4), len(b6)) == (6, 3, 2), (len(b2), len(b4), len(b6))
SETS = [(1, "Junghans, natural", b2[0]), (2, "Vaughn, natural", b2[1]),
        (3, "Hinterberger, $^{20}$Ne", b2[2]), (11, "Sikkema, $^{20}$Ne", b2[5]),
        (4, "Almeida, natural, low", b2[3]), (5, "Almeida, natural, high", b2[4]),
        (6, "Almeida, $^{21}$Ne, low", b6[0]), (7, "Almeida, $^{21}$Ne, high", b6[1]),
        (8, "Almeida, $^{22}$Ne, low", b4[0]), (9, "Almeida, $^{22}$Ne, high", b4[1]),
        (10, "Sikkema, $^{22}$Ne", b4[2])]

fig, axs = plt.subplots(3, 4, figsize=(12, 8.2))
for ax, (key, name, a) in zip(axs.flat, SETS):
    e, model, y, err = a[:, 0], a[:, 3], a[:, 5], a[:, 6]
    chi2 = (((y - model) / err) ** 2).sum()
    ax.errorbar(e, y, err, fmt=".", ms=1.5, lw=0.3, color="0.5")
    ax.plot(e, model, color="C3", lw=0.8)
    ax.set_yscale("log")
    pos = y[y > 0]
    ax.set_ylim(max(np.percentile(pos, 0.5) / 2, 1e-2), max(y.max(), model.max()) * 1.5)
    ax.set_title(f"{key}: {name}\n$\\chi^2$ = {chi2:,.0f} for {len(e)} points", fontsize=8.5)
    ax.tick_params(labelsize=7)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for ax in axs[:, 0]:
    ax.set_ylabel(r"$\sigma_\mathrm{tot}$ (b)", fontsize=8)
for ax in axs[-1, :]:
    ax.set_xlabel(r"$E_\mathrm{c.m.}$ (MeV)", fontsize=8)
axs.flat[-1].axis("off")
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "natural_neon_fit.png"), dpi=90)

d = blocks(os.path.join(bare_dir, "AZUREOut_sample_seg=1.extrap"))[0]
jung = b2[0]
fig, ax = plt.subplots(figsize=(10, 5))
ax.errorbar(jung[:, 0], jung[:, 5], jung[:, 6], fmt=".", ms=1.5, lw=0.3, color="0.65",
            label="Junghans (natural neon, with resolution)", zorder=1)
ax.plot(d[:, 0], d[:, 1], color="k", lw=0.8, label="natural sum")
for col, label, color in [(2, r"0.9048 $\sigma$(n+$^{20}$Ne) [$^{21}$Ne]", "C0"),
                          (3, r"0.0925 $\sigma$(n+$^{22}$Ne) [$^{23}$Ne]", "C3"),
                          (4, r"0.0027 $\sigma$(n+$^{21}$Ne) [$^{22}$Ne]", "C2")]:
    ax.plot(d[:, 0], d[:, col], color=color, lw=0.7, label=label)
ax.set_yscale("log")
ax.set_ylim(1e-3, 60)
ax.set_xlim(0, d[:, 0].max())
ax.set_xlabel(r"$E_\mathrm{c.m.}$ (MeV)")
ax.set_ylabel(r"$\sigma_\mathrm{tot}$ (b)")
ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "natural_neon_decomposition.png"), dpi=90)
print("wrote", os.path.join(fig_dir, "natural_neon_fit.png"), "and natural_neon_decomposition.png")
