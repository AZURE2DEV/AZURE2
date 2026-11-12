#!/usr/bin/env python3
"""Writes the models of the formal-parametrization external-capture final-state check (run once; the
.azr files are committed).

Source: tests/two_nuclei_ec/models/A.azr, 12C(p,g0)13N with external capture to the 1/2- ground state and
a 1/2- background pole at 20 MeV.  Only the ORDER of the <levels> lines differs between the two files of a
pair; level IDs, values and flags are the same, so the physics is the same.

  A_orig.azr  ground state first in its 1/2- group (the order of the source)
  A_swap.azr  the 20 MeV pole first
  B_orig.azr  A_orig with a second 1/2- background pole at 30 MeV: ground state, 20 MeV, 30 MeV
  B_mid.azr   20 MeV pole, ground state, 30 MeV pole

Test segments: 12C(p,g0) angle integrated, and elastic scattering at 120 deg (no capture part: it tells an
error of the internal R-matrix from one of the capture terms).
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "..", "two_nuclei_ec", "models", "A.azr")).read().split("\n")
i0, i1 = src.index("<levels>"), src.index("</levels>")
blocks, cur = [], []
for ln in src[i0 + 1:i1]:
    if ln.strip():
        cur.append(ln)
    elif cur:
        blocks.append(cur); cur = []
if cur:
    blocks.append(cur)


def key(b):
    t = b[0].split()
    return (t[0], t[1], t[2])


assert key(blocks[0]) == ("0.5", "-1", "0") and key(blocks[1]) == ("0.5", "-1", "20")
gs, bg = blocks[0], blocks[1]
last_id = max(int(b[0].split()[8]) for b in blocks)
bg30 = []
for ln in bg:                       # a copy of the 20 MeV pole at 30 MeV, with a new level ID
    t = ln.split()
    t[2], t[8] = "30", str(last_id + 1)
    if float(t[11]) != 0.0:
        t[11] = "3e+06"
    bg30.append(" " + " ".join(t))


def write(name, bl):
    out = src[:i0 + 1]
    for b in bl:
        out += b + [""]
    out += src[i1:]
    out.insert(out.index("</segmentsTest>"), "1 1 1 0.2 2.0 0.01 120 120 0 1")
    open(os.path.join(HERE, "models", name), "w").write("\n".join(out))


write("A_orig.azr", [gs, bg] + blocks[2:])
write("A_swap.azr", [bg, gs] + blocks[2:])
write("B_orig.azr", [gs, bg, bg30] + blocks[2:])
write("B_mid.azr", [bg, gs, bg30] + blocks[2:])
