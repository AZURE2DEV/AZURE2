#!/usr/bin/env python3
"""Blind evaluation of peakfind.py candidates against a level set, with a chance baseline (candidates shifted by random offsets
of 50-300 keV, same count and spacing).  usage: eval_peaks.py peaks.txt parameters.out"""
import sys, re
import numpy as np
sys.path.insert(0, "/users/rdeboer1/AZURE2/.claude/skills/resonance-search/scripts")
from peakfind import read_levels
rows = []
for line in open(sys.argv[1]):
    m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)k\S*\s*(?:\(res\.\))?\s+(\d+)\s+(\d+)", line)
    if m: rows.append(dict(ex=float(m.group(2)), w=float(m.group(3)), sets=int(m.group(4))))
lev = [L for L in read_levels(sys.argv[2]) if L["g"] >= 1.0]
rng = np.random.default_rng(1)
def hits(E, Ls): return sum(any(abs(e - L["e"]) <= max(0.010, 0.5e-3 * L["g"]) for e in E) for L in Ls)
def prec(E, Ls): return sum(any(abs(e - L["e"]) <= max(0.010, 0.5e-3 * L["g"]) for L in Ls) for e in E)
for smin in (1, 2, 3, 5):
    E = np.array([r["ex"] for r in rows if r["sets"] >= smin]); lo, hi = E.min(), E.max()
    Ls = [L for L in lev if lo <= L["e"] <= hi]
    out = []
    for lab, (g0, g1) in {"1-10": (1, 10), "10-100": (10, 100), "100-500": (100, 500), ">500": (500, 1e9)}.items():
        Lb = [L for L in Ls if g0 <= L["g"] < g1]
        ch = np.mean([hits(E + rng.choice([-1, 1]) * rng.uniform(0.05, 0.3), Lb) for _ in range(50)])
        out.append(f"{lab} keV {hits(E, Lb)}/{len(Lb)} (chance {ch:.1f})")
    pc = np.mean([prec(E + rng.choice([-1, 1]) * rng.uniform(0.05, 0.3), Ls) for _ in range(50)])
    print(f"sets >= {smin}: {len(E)} candidates in {lo:.2f}-{hi:.2f} MeV; at a level {prec(E, Ls)} (chance {pc:.0f}); found: " + "; ".join(out))
