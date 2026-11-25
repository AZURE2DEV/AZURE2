#!/usr/bin/env python3
"""The level that replaces a subtracted BOUND state of the background, for check.sh section 17
(potential_background.md section 14.9).

The model: models/capture_a40.azr.in at a_c = 4.87 fm with the s1/2 depth 104 MeV, at which the s-wave potential binds
its state with 2 nodes by 0.126 MeV.  That state is subtracted from the background (`subtract 2`) and put back as a
1/2+ level given by its OBSERVED values, none of which comes from the engine:

  - the energy E_x = S - B and the ANC of the normalised bound state (make_reference.py's bound_state, scipy DOP853);
  - the photon widths to the two bound p states, Gamma = 2 k_gamma^3 g^2 with
        g = R_f I / (2 sqrt 2),   I = int_0^inf u_s(r) u_f(r) r dr,
    u_s and u_f normalised over all r, R_f the E1 factor of make_subtract_level.py, k_gamma = (B_f - B) / hbar c, and
    the sign of g on the width (the interference sign).

The engine reads the photon width as the observed one, 2 P_gamma g~^2 / N with N = 1 + gamma^2 S~'(E_b) of the dressed
channel, and splits g~ into the internal and the channel part itself.  check.sh compares the capture of this model
with the unsubtracted one.

Run: python3 make_bound_level.py     (writes models/bound_level.azr.in; no engine needed)
"""
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
V_S, A_C, NODES = 104.0, 4.87, 2
sys.argv = [sys.argv[0], f"--s-depth={V_S}"]
sys.path.insert(0, HERE)
import make_reference as m  # noqa: E402

S_SEP = 1.4548867432  # the threshold (MeV of E_x): the p3/2 binding energy


def main():
    old = "wave            l 0 J 1/2 depth 67.5 MeV background interior boundary 0.0 capture interior"
    new = f"wave            l 0 J 1/2 depth {V_S} MeV background interior boundary 0.0 subtract {NODES} capture interior"
    Bs, Cs, us_in, us_out = m.bound_state(0, 1, 0.005, 1.5)
    zeff = m.Z1 * m.M2 / (m.M1 + m.M2) - m.Z2 * m.M1 / (m.M1 + m.M2)
    width = {}
    for two_jf, pair in ((3, "2"), (1, "3")):
        Bf, _, uf_in, uf_out = m.bound_state(1, two_jf, 0.3, 4.0)
        overlap = m.overlap(uf_in, uf_out, us_in, us_out)
        W = (-1) ** 3 * m.six_j(Fraction(1), Fraction(1), Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(two_jf, 2))
        R = math.sqrt(m.ALPHA * m.HBARC) * zeff * math.sqrt(48.0) / 3.0 * math.sqrt(two_jf + 1) * W
        g = R * overlap / (2.0 * math.sqrt(2.0))
        width[pair] = math.copysign(2.0 * ((Bf - Bs) / m.HBARC) ** 3 * g * g * 1e6, g)  # eV
        print(f"# to 2J_f = {two_jf}: g = {g!r} MeV^1/2, Gamma_gamma = {width[pair]!r} eV", file=sys.stderr)
    print(f"# the s state: B = {Bs!r} MeV, ANC = {Cs!r} fm^-1/2", file=sys.stderr)
    value = {"1": Cs, **width}
    out, inl = [], False
    text = open(os.path.join(HERE, "models", "capture_a40.azr.in")).read()
    assert old in text
    for line in text.replace(old, new).split("\n"):
        if line.startswith("<levels>"):
            inl = True
        elif line.startswith("</levels>"):
            inl = False
        t = line.split()
        if inl and len(t) >= 31:
            if t[26] == "0":
                t[27] = str(A_C)
            if t[0] == "0.5" and t[1] == "1" and t[2] == "40":
                t[2] = repr(S_SEP - Bs)
                t[11] = repr(value[t[5]])
            line = " ".join(t)
        out.append(line)
    with open(os.path.join(HERE, "models", "bound_level.azr.in"), "w") as f:
        f.write("\n".join(out))


if __name__ == "__main__":
    main()
