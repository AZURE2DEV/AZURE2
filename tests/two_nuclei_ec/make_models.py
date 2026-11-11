#!/usr/bin/env python3
"""Writes the models of the two-nuclei external-capture check (run once; the .azr files are committed).

  A.azr   tests/13N/13N.azr (p + 12C -> 13N, capture to the ground state with external capture),
          no data, one test segment: 12C(p,g0), angle integrated
  B.azr   the same level scheme as a neutron project: n + 12C -> 13C (Z1 = 0, neutron mass,
          S_n(13C), the gamma pair's heavy particle 13C); only its role as a SECOND nucleus with
          external capture matters here, not its physics
  M.azr   A's and B's levels in one file (pyazr.merge_projects; A's pair keys 1-2, B's 3-4),
          no <abundances> block, test segments: A's, then B's

In M the n + 12C entrance pair (number 3) has two decays (pairs 3 and 4), so its capture decay
is decay number 2 while the capture pair is pair number 4.
"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("PYAZR_ROOT", os.path.join(HERE, "..", "..")))   # a checkout with a built pyazr
from pyazr import merge_projects

OUT = os.path.join(HERE, "models")
src = open(os.path.join(HERE, "..", "13N", "13N.azr")).read()


def block(text, tag, body):
    a, b = text.index(f"<{tag}>"), text.index(f"</{tag}>")
    return text[:a] + f"<{tag}>\n" + body + text[b:]


def strip(text):
    text = block(text, "segmentsData", "")
    text = block(text, "segmentsTest", "1 1 2 0.2 2.0 0.01 0 0 0 0\n")
    text = block(text, "targetInt", "")
    for tag in ("parameterSettings", "lastRun", "mcmc"):
        text = re.sub(rf"<{tag}>.*?</{tag}>\n", "", text, flags=re.S)
    lines = text.split("\n")
    i = lines.index("<config>")
    lines[i + 2] = "output/ #Full Path to Output Directory"
    lines[i + 3] = "checks/ #Full Path to Checks Directory"
    for k in range(i + 4, lines.index("</config>")):
        if "Check" in lines[k]:
            lines[k] = "none #" + lines[k].split("#", 1)[1]
    return "\n".join(lines)


A = strip(src)
open(os.path.join(OUT, "A.azr"), "w").write(A)

# B: n + 12C -> 13C.  Level-line columns (0-based): 17 light mass, 18 heavy mass, 19 Z1, 20 Z2,
# 21 entrance separation energy, 22 the pair's own separation energy (0 on a gamma pair).
MN, M13C, SN = 1.00866, 13.00335, 4.94631
out, inside = [], False
for ln in A.split("\n"):
    if ln.strip() == "<levels>": inside = True
    elif ln.strip() == "</levels>": inside = False
    elif inside and ln.strip():
        t = ln.split()
        if float(t[17]) == 0.0:                 # gamma + 13C
            t[18], t[20] = str(M13C), "6"
        else:                                   # n + 12C
            t[17], t[19], t[22] = str(MN), "0", str(SN)
        t[21] = str(SN)
        ln = " ".join(t)
    out.append(ln)
B = "\n".join(out)
open(os.path.join(OUT, "B.azr"), "w").write(B)

rep = merge_projects([os.path.join(OUT, "A.azr"), os.path.join(OUT, "B.azr")], os.path.join(OUT, "M.azr"),
                     abundances={}, natural=False)
M = open(os.path.join(OUT, "M.azr")).read()
M = re.sub(r"<abundances>\n.*?</abundances>\n", "", M, flags=re.S)
lines = M.split("\n")
i = lines.index("<config>")
lines[i + 2] = "output/ #Full Path to Output Directory"
lines[i + 3] = "checks/ #Full Path to Checks Directory"
open(os.path.join(OUT, "M.azr"), "w").write("\n".join(lines))
print("pair map", rep.pair_map)
