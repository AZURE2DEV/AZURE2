#!/usr/bin/env python3
"""Models for tests/light_excitation/check.sh, written at run time.

    make_models.py <out_dir>

aC_*: alpha + 12C (16O), pair 1 the ground states, pair 2 the inelastic exit
alpha' + 12C(2+, 4.43891 MeV), four levels with widths in both pairs.
  aC_heavy       the excitation on the heavy partner, as always: e2 = E*, J2 = 2+ (31 columns)
  aC_heavy0      the same with an explicit 32nd column "0" on every line
  aC_light       the excitation moved to the light partner: e1 = E* (column 32), J1 = 2+, J2 = 0+, e2 = 0
  aC_light_noex  aC_light without column 32 -- pair 2 opens at S: the check must see the difference
  bad_gamma      aC_heavy plus a gamma pair (type 10) with e1 = 1 MeV                -> ERROR
  bad_negative   aC_light with e1 = -E*                                            -> ERROR
  bad_mixed      aC_light with a different e1 on one line of pair 2                 -> ERROR
  upos_heavy     aC_heavy plus a UPOS (secondary-gamma) segment of (a,a'g): accepted
  bad_upos       aC_light with the same segment: the gamma would come from the light particle -> ERROR

cc_*: 12C(2+) + 12C(2+) (24Mg), one 4+ level, includeAllJPi=1 with l <= 4.
  cc_new   e1 = e2 = E*, j1 = j2 = 2+: identical bosons
  cc_old   the old way, e1 = 0, e2 = 2 E*, j1 = j2 = 2+: not identical (e1 != e2)
  cc_bad   cc_new whose level also lists the forbidden channel s = 3, l = 2        -> ERROR
"""
import os
import sys

E2P = 4.43891          # 12C 2+
S_AC = 7.16191         # alpha separation energy of 16O
S_CC = 13.933          # 12C separation energy of 24Mg

CONFIG = """<config>
true #Perform A-Matrix Calculation
output/ #Full Path to Output Directory
checks/ #Full Path to Checks Directory
{cn} #Compond Nucleus Check
none #Boundary Condition Check
none #Data Check
none #Lo-Matrix and Penetrability Check
none #Legendre Polynomial Check
none #Coulomb Amplitudes Check
none #Reaction Pathway Check
none #Angular Distributions Check
{extra}</config>
"""

TAIL = """<targetInt>
</targetInt>
<lastRun>
9997
"output/param.sav"
""
0 0
0 ""
-1. -1. -1.
</lastRun>
"""


def fmt(x):
    return repr(float(x)) if not float(x).is_integer() else str(int(x)) if isinstance(x, int) else f"{x:g}"


def line(J, pi, E, lid, pair, s, l, width, p, e1=None):
    """One <levels> line; p = (j1, pi1, j2, pi2, e2, m1, m2, z1, z2, S, ptype, chRad)."""
    j1, pi1, j2, pi2, e2, m1, m2, z1, z2, S, ptype, rad = p
    cols = [J, pi, E, 1, 1, pair, int(round(2 * s)), 2 * l, lid, 1, 0, width,
            j1, pi1, j2, pi2, e2, m1, m2, z1, z2, S, S, 0, 0, 0.0, ptype, rad, 0, 0, 0]
    out = " ".join(f"{c:g}" if isinstance(c, float) else str(c) for c in cols)
    if e1 is not None:
        out += f" {e1:g}" if isinstance(e1, float) else f" {e1}"
    return out


def write(path, cfg, levels, data="", test=""):
    with open(path, "w") as f:
        f.write(cfg)
        f.write("<levels>\n" + "\n".join(levels) + "\n\n</levels>\n")
        f.write("<segmentsData>\n" + data + "</segmentsData>\n")
        f.write("<segmentsTest>\n" + test + "</segmentsTest>\n")
        f.write(TAIL)


def ac_levels(kind):
    P1 = (0.0, 1, 0.0, 1, 0.0, 4.0026, 12.0, 2, 6, S_AC, 0, 5.43)
    if kind in ("heavy", "heavy0", "bad_gamma"):
        P2, e1 = (0.0, 1, 2.0, 1, E2P, 4.0026, 12.0, 2, 6, S_AC, 0, 5.43), None
    else:
        P2, e1 = (2.0, 1, 0.0, 1, 0.0, 4.0026, 12.0, 2, 6, S_AC, 0, 5.43), E2P
        if kind == "light_noex":
            e1 = None
        elif kind == "bad_negative":
            e1 = -E2P
    G3 = (1.0, 1, 0.0, 1, 0.0, 0.0, 15.9949, 0, 8, 0.0, 10, 0.0)    # gamma to the 16O ground state
    z = "0" if kind == "heavy0" else None
    rows = []
    # (J, pi, E, id, [(pair, s, l, width)])
    levels = [
        (2.0, 1, 13.0, 1, [(1, 0, 2, 40000.0), (2, 2, 0, 20000.0), (2, 2, 2, -5000.0)]),
        (0.0, 1, 13.3, 2, [(1, 0, 0, 50000.0), (2, 2, 2, 8000.0)]),
        (1.0, -1, 13.6, 3, [(1, 0, 1, 60000.0), (2, 2, 1, 15000.0), (2, 2, 3, 3000.0)]),
        (3.0, -1, 14.2, 4, [(1, 0, 3, 30000.0), (2, 2, 1, 10000.0), (2, 2, 3, -4000.0), (2, 2, 5, 1000.0)]),
    ]
    first_p2 = True
    for J, pi, E, lid, chans in levels:
        for pair, s, l, w in chans:
            if pair == 1:
                rows.append(line(J, pi, E, lid, 1, s, l, w, P1, z))
            else:
                this_e1 = e1 if e1 is not None else z
                if kind == "bad_mixed" and not first_p2:
                    this_e1 = E2P + 0.01
                first_p2 = False
                rows.append(line(J, pi, E, lid, 2, s, l, w, P2, this_e1))
        if kind == "bad_gamma" and (J, pi) == (1.0, -1):
            rows.append(line(J, pi, E, lid, 3, 1, 1, 0.01, G3, 1.0))
        rows.append("")
    return rows[:-1]


def ac_tests():
    t = [
        "1 1 2 6.5 10.4 0.05 0 0 0 0",          # (a,a') angle integrated
        "1 1 2 6.5 10.4 0.1 20 160 20 5",       # (a,a') c.m. differential
        "1 1 2 8.5 10.4 0.1 20 160 20 1",       # (a,a') lab differential
        "1 1 1 6.5 10.4 0.1 20 160 20 5",       # elastic
        "1 2 1 1.4 4.5 0.05 0 0 0 0",           # (a',a): the excited pair as entrance
        "1 2 1 1.4 4.5 0.1 30 150 30 1",        # (a',a) lab differential
        "1 2 2 1.4 4.5 0.1 30 150 30 5",        # a' + 12C(2+) elastic
    ]
    return "".join(x + "\n" for x in t)


def ac_data(out):
    d = os.path.join(out, "data")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "inel_lab.dat"), "w") as f:
        f.write("# E_lab(MeV) theta_lab(deg) dsigma/dOmega_lab(b/sr) error -- synthetic\n")
        for e in (8.8, 9.3, 9.9, 10.3):
            for a in (30, 60, 90, 120, 150):
                f.write(f"{e} {a} 0.01 0.001\n")
    with open(os.path.join(d, "inel_int.dat"), "w") as f:
        f.write("# E_lab(MeV) dummy sigma(b) error -- synthetic\n")
        for i in range(30):
            f.write(f"{6.6 + 0.12 * i:.2f} 0 0.05 0.005\n")
    return ("1 1 2 0 100 0 180 1 1 0 0 0 0 0 data/inel_lab.dat 0 0\n"
            "1 1 2 0 100 0 0 0 1 0 0 0 0 0 data/inel_int.dat 0 0\n")


def cc_levels(kind):
    if kind == "old":
        P, e1 = (2.0, 1, 2.0, 1, 2 * E2P, 12.0, 12.0, 6, 6, S_CC, 0, 6.0), None
    else:
        P, e1 = (2.0, 1, 2.0, 1, E2P, 12.0, 12.0, 6, 6, S_CC, 0, 6.0), E2P
    chans = [(0, 4, 100.0), (2, 2, 50.0), (4, 0, 20.0)]
    if kind == "bad":
        chans.append((3, 2, 10.0))
    return [line(4.0, 1, 26.0, 1, 1, s, l, w, P, e1) for s, l, w in chans]


def main(out):
    os.makedirs(out, exist_ok=True)
    data = ac_data(out)
    for kind in ("heavy", "heavy0", "light", "light_noex", "bad_gamma", "bad_negative", "bad_mixed"):
        write(os.path.join(out, f"aC_{kind}.azr" if not kind.startswith("bad") else f"{kind}.azr"),
              CONFIG.format(cn="none", extra=""), ac_levels(kind), data, ac_tests())
    upos = "1 1 2 0 100 0 180 4 1 0 0 0 0 0 data/inel_lab.dat 0 1 2 0 0\n"   # L = 2, final spin 0, delta 0
    write(os.path.join(out, "upos_heavy.azr"), CONFIG.format(cn="none", extra=""), ac_levels("heavy"), data + upos, "")
    write(os.path.join(out, "bad_upos.azr"), CONFIG.format(cn="none", extra=""), ac_levels("light"), data + upos, "")
    cc_cfg = CONFIG.format(cn="screen", extra="includeAllJPi=1\nmaxOrbitalMomentum=4\n")
    cc_test = "1 1 1 6.0 12.0 0.5 10 170 10 5\n"
    for kind in ("new", "old", "bad"):
        write(os.path.join(out, f"cc_{kind}.azr"), cc_cfg, cc_levels(kind), "", cc_test)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
