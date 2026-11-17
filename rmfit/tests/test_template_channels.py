#!/usr/bin/env python3
"""The template screen must only build levels the bounded refit could reach (2026-09-20).

A residual screen once reported a level worth a clear gain: its winning template was an "equal
split" of 10 keV over 15 open channels -- 667 eV each -- five of which have Wigner limits of 0.4 eV
.. 1e-15 eV (theta^2 1e3 .. 1e17).  Templates are evaluated frozen, so nothing stopped it; the
honest add gained 28.  Pure tests of the three helpers that now guard this."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.levelscan import carrying_channels, template_amplitude, align_priors

ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)

def main():
    # a level with nine channels: Gamma_W (eV) of n0 L4, a0 L3, n1 L4, a1 L4, n2 L1, n2 L3, n2 L5, a3 L6, closed n3
    gw = [5.62e5, 1.054e6, 0.43, 2.87, 2.82e5, 23.2, 9.9e-5, 6.3e-15, 0.0]
    ga = [1.85, 0.84, 1.85, 0.84, 1.85, 1.85, 1.85, 0.84, 0.0]
    o, d, c = carrying_channels(gw, ga, dead_channel_ev=1.0)
    check("channels with Gamma_W >= 1 eV carry the template", o == [0, 1, 3, 4, 5], o)
    check("sub-eV Wigner limits are dead, not open", d == [2, 6, 7], d)
    check("no Wigner width = closed", c == [8], c)
    o10, d10, _ = carrying_channels(gw, ga, dead_channel_ev=10.0)
    check("the rule follows dead_channel_ev", o10 == [0, 1, 4, 5] and 3 in d10, (o10, d10))
    a, w = template_amplitude(667.0, 5.62e5, 1.85, 1.0)
    check("a reachable width is set exactly: gamma = gamma_W sqrt(Gamma/Gamma_W)", abs(a - 1.85 * np.sqrt(667.0 / 5.62e5)) < 1e-12 and abs(w - 667.0) < 1e-9, (a, w))
    a, w = template_amplitude(667.0, 23.2, 1.85, 1.0)
    check("an unreachable width is clipped at the theta^2 cap", abs(a - 1.85) < 1e-12 and abs(w - 23.2) < 1e-9, (a, w))
    a, w = template_amplitude(667.0, 23.2, 1.85, 3.0)
    check("the cap is the policy's (3 for background levels)", abs(a - 1.85 * np.sqrt(3)) < 1e-12 and abs(w - 69.6) < 1e-9, (a, w))
    check("a closed channel gets nothing", template_amplitude(667.0, 0.0, 0.0, 1.0) == (0.0, 0.0))
    # priors: an existing level with its strength on n2 L1 (5th channel of the file order)
    order = [(1, 4, 0.5), (2, 3, 0.5), (3, 4, 0.5), (4, 4, 0.5), (5, 1, 2.5), (5, 3, 2.5)]
    pri = [dict(label="7/2+@10.42", fractions=[0.01, 0.12, 0.0, 0.0, 0.85, 0.02], signs=[-1, -1, 1, 1, -1, 1])]
    carry = [order[i] for i in (0, 1, 3, 4, 5)]                     # n1 L4 (index 2) is dead
    al = align_priors(pri, order, carry)
    check("a prior keeps its strength on the SAME channels after a dead one is removed",
          len(al) == 1 and al[0]["fractions"] == [0.01, 0.12, 0.0, 0.85, 0.02] and al[0]["signs"] == [-1, -1, 1, -1, 1], al)
    check("the old truncation [:n] would have moved 0.85 off n2 L1", pri[0]["fractions"][:5][3] == 0.0 and al[0]["fractions"][3] == 0.85)
    check("a prior with nothing on the carrying channels is dropped",
          align_priors([dict(label="x", fractions=[0, 0, 1.0, 0, 0, 0], signs=[1] * 6)], order, carry) == [])
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
