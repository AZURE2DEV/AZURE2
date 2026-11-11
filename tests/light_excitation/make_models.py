#!/usr/bin/env python3
"""Models for tests/light_excitation/check.sh, written at run time (project format 3: the <levels> columns, data
sets and extrapolations as records, docs/upgrade/project_format3.md).

    make_models.py <out_dir>

aC_*: alpha + 12C (16O), pair 1 the ground states, pair 2 the inelastic exit
alpha' + 12C(2+, 4.43891 MeV), four levels with widths in both pairs.
  aC_heavy       the excitation on the heavy partner, as always: e2 = E*, J2 = 2+ (31 columns)
  aC_heavy0      the same with an explicit 32nd column "0" on every line
  aC_light       the excitation moved to the light partner: e1 = E* (column 32), J1 = 2+, J2 = 0+, e2 = 0
  aC_heavy_dsdo  aC_heavy with the differential (a,a') extrapolations only, the reference of aC_light_noex
  aC_light_noex  aC_light without column 32 -- pair 2 opens at S: the check must see the difference (its
                 differential (a,a') extrapolations only: pair 2 is then a second a + 12C, the keys written,
                 and an angle-integrated cross section of 12C(a,a)12C is refused as elastic)
  bad_gamma      aC_heavy plus a gamma pair (type 10) with e1 = 1 MeV                -> ERROR
  bad_negative   aC_light with e1 = -E* (the elastic extrapolation only: a reaction label
                 cannot name a negative excitation)                                -> ERROR
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

FORMAT = "format=3\n"
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
        f.write("<data>\n" + data + "</data>\n")
        f.write("<extrapolations>\n" + test + "</extrapolations>\n")
        f.write("<targets>\n</targets>\n")


def reaction(kind, entrance, exit_):
    """The reaction label of entrance -> exit pair of an aC model: pair 1 a + 12C, pair 2 with the 4.43891 MeV
    excitation on the heavy partner (12C[E*]) or on the light one (a[E*]); without column 32 (light_noex) pair 2 is a
    second a + 12C and the data set names both keys."""
    light = kind in ("light", "light_noex", "bad_negative", "bad_mixed", "bad_upos")
    ex = f"[{E2P:g}]" if kind != "light_noex" else ""
    proj = {1: "a", 2: "a" + ex if light else "a"}
    targ = {1: "12C", 2: "12C" if light else "12C" + ex}
    text = f"{targ[entrance]}({proj[entrance]},{proj[exit_]}){targ[exit_]}"
    if kind == "light_noex":   # both pairs are a + 12C: the label names each by its key
        text += f"\n    entrance_pair  {entrance}\n    exit_pair      {exit_}"
    return text


def extrapolation(name, kind, entrance, exit_, lo, hi, step, alo=0, ahi=0, astep=0, frame=None):
    """One extrapolation: E_lab lo..hi in steps (MeV); with a frame ("cm" or "lab"), dsigma/dOmega at alo..ahi deg."""
    q = "    quantity       dsigma_domega\n    unit           b/sr\n" if frame else "    quantity       sigma\n    unit           b\n"
    out = (f"extrapolation {name}\n    reaction       {reaction(kind, entrance, exit_)}\n{q}"
           f"    energies       range {lo:g} {hi:g} MeV lab step {step:g}\n")
    if frame:
        out += f"    angles         range {alo:g} {ahi:g} deg {frame} step {astep:g}\n"
    return out + "\n"


def ac_levels(kind):
    P1 = (0.0, 1, 0.0, 1, 0.0, 4.0026, 12.0, 2, 6, S_AC, 0, 5.43)
    if kind in ("heavy", "heavy0", "heavy_dsdo", "bad_gamma"):
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


def ac_tests(kind):
    if kind == "bad_negative":
        return extrapolation("el_cm", kind, 1, 1, 6.5, 10.4, 0.1, 20, 160, 20, "cm")
    t = [
        ("x_inel_int", 1, 2, 6.5, 10.4, 0.05, 0, 0, 0, None),      # (a,a') angle integrated
        ("x_inel_cm", 1, 2, 6.5, 10.4, 0.1, 20, 160, 20, "cm"),    # (a,a') c.m. differential
        ("x_inel_lab", 1, 2, 8.5, 10.4, 0.1, 20, 160, 20, "lab"),  # (a,a') lab differential
        ("x_el_cm", 1, 1, 6.5, 10.4, 0.1, 20, 160, 20, "cm"),      # elastic
        ("x_back_int", 2, 1, 1.4, 4.5, 0.05, 0, 0, 0, None),       # (a',a): the excited pair as entrance
        ("x_back_lab", 2, 1, 1.4, 4.5, 0.1, 30, 150, 30, "lab"),   # (a',a) lab differential
        ("x_exc_el_cm", 2, 2, 1.4, 4.5, 0.1, 30, 150, 30, "cm"),   # a' + 12C(2+) elastic
    ]
    if kind in ("light_noex", "heavy_dsdo"):
        t = t[1:3]   # the differential (a,a') only: sigma with the label 12C(a,a)12C is refused as elastic
    return "".join(extrapolation(*x[:1], kind, *x[1:]) for x in t)


def dataset(name, kind, quantity, upos=False):
    """A (a,a') data set over E_lab 0-100 MeV: lab dsigma/dOmega (or, upos, the secondary-gamma correlation L = 2,
    final spin 0, delta 0, at c.m. angles), or the angle-integrated cross section."""
    f, unit = ("inel_lab", "b/sr") if quantity == "dsigma_domega" else ("inel_int", "b")
    out = f"dataset {name}\n    reaction       {reaction(kind, 1, 2)}\n    quantity       {quantity}\n    unit           {unit}\n"
    if upos:
        out += "    secondary_gamma L 2 Ic 0 delta 0\n"
    if quantity == "dsigma_domega":
        out += (f"    file           data/{f}.dat\n    columns        E_lab[MeV] theta_{'cm' if upos else 'lab'}[deg] value error\n"
                "    energy_window  0 100 MeV\n    angle_window   0 180 deg\n")
    else:
        out += f"    file           data/{f}.dat\n    columns        E_lab[MeV] _ value error\n    energy_window  0 100 MeV\n"
    return out + "    norm           1 fixed\n\n"


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
    return lambda kind: ("" if kind in ("light_noex", "heavy_dsdo", "bad_negative") else
                         dataset("inel_lab", kind, "dsigma_domega") + dataset("inel_int", kind, "sigma"))


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
    for kind in ("heavy", "heavy0", "heavy_dsdo", "light", "light_noex", "bad_gamma", "bad_negative", "bad_mixed"):
        write(os.path.join(out, f"aC_{kind}.azr" if not kind.startswith("bad") else f"{kind}.azr"),
              CONFIG.format(cn="none", extra=FORMAT), ac_levels(kind), data(kind), ac_tests(kind))
    for name, kind in (("upos_heavy", "heavy"), ("bad_upos", "bad_upos")):
        write(os.path.join(out, f"{name}.azr"), CONFIG.format(cn="none", extra=FORMAT),
              ac_levels("light" if kind == "bad_upos" else kind),
              data(kind) + dataset("inel_lab_upos", kind, "dsigma_domega", upos=True), "")
    cc_cfg = CONFIG.format(cn="screen", extra="includeAllJPi=1\nmaxOrbitalMomentum=4\n" + FORMAT)
    for kind in ("new", "old", "bad"):
        pair = f"12C[{E2P:g}](12C[{E2P:g}],12C[{E2P:g}])12C[{E2P:g}]" if kind != "old" else \
            f"12C[{2 * E2P:g}](12C,12C)12C[{2 * E2P:g}]"
        cc_test = (f"extrapolation el_cm\n    reaction       {pair}\n    quantity       dsigma_domega\n"
                   "    unit           b/sr\n    energies       range 6 12 MeV lab step 0.5\n"
                   "    angles         range 10 170 deg cm step 10\n")
        write(os.path.join(out, f"cc_{kind}.azr"), cc_cfg, cc_levels(kind), "", cc_test)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
