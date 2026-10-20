#!/usr/bin/env python3
"""Broad-only copy of a fit: every level with total width < GMAX keV (from parameters.out) gets all its channel amplitudes zeroed
and fixed (energy fixed too), so the free parameters read from the .sav by name never touch it.
usage: make_broad_only.py fit.azr parameters.out out.azr [GMAX_keV=50]"""
import sys, os, re
sys.path.insert(0, "/users/rdeboer1/AZURE2")
from pyazr import AzrModel
azr, par, out = sys.argv[1:4]; gmax = float(sys.argv[4]) if len(sys.argv) > 4 else 50.0
lev, cur = [], None
for line in open(par):
    m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
    if m:
        J = m.group(1); cur = dict(J=float(J[:-1]), pi=J[-1], e=float(m.group(2)), g=0.0); lev.append(cur); continue
    m = re.search(r"G\s+=\s+(\S+) keV", line)
    if m and cur is not None: cur["g"] += abs(float(m.group(1)))
m = AzrModel.from_file(azr)
n = 0
for l in m.levels:
    if not l.active: continue
    c = [L for L in lev if abs(L["e"] - l.energy) < 2e-3 and l.jpi.startswith({0.5:"1/2",1.5:"3/2",2.5:"5/2",3.5:"7/2",4.5:"9/2",5.5:"11/2",6.5:"13/2"}.get(L["J"], "?")) and l.jpi.endswith(L["pi"])]
    if len(c) != 1: print("no unique match", l.jpi, l.energy, len(c)); continue
    if c[0]["g"] < gmax and l.energy > 4.2:   # keep bound / sub-threshold states (ANC only, no width)
        l.set_fixed(True)
        for ch in l.channels: ch.gamma = 0.0; ch.channel_fixed = True
        n += 1
m.set_output_dir(os.path.abspath("output") + "/"); os.makedirs("output", exist_ok=True)
m.write(out); print(f"{n} levels with Gamma < {gmax} keV zeroed and fixed -> {out}")
