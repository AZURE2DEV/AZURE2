#!/usr/bin/env python3
"""Writes the models of the isospin-levels check (run once; the .azr files are committed).

A toy A=4 system in Hale's charge-independent scheme (Tilley, Weller and Hale,
NPA 541 (1992); G.M. Hale, LA-UR-94-1264):

  pair 1  p + 3He  (4Li, T3 = +1)     the reference nucleus, threshold at Ex = 0
  pair 2  n + t    (4H,  T3 = -1)     threshold at Ex = 0, shift -0.86 MeV
  pair 3  p + t    (4He, T3 =  0)     Ex 19.8136, shift -0.64 MeV
  pair 4  n + 3He  (4He)              Ex 20.5776
  pair 5  d + d    (4He)              Ex 23.8465

  level 1  2-  T=1   lines in 4Li, 4H and 4He (all channels l=1, s=1)
  level 2  2-  T=0   lines in 4He only

Radii 4.9 fm (nucleon + trinucleon) and 7.0 fm (d+d), as in Hale's analysis.

  iso.azr          the model; test segments S1 p+3He elastic at 120 deg, S2 n+t
                   integrated, S3 t(p,n)3He, S4 d(d,p)t, S5 d(d,n)3He
  iso_mix.azr      the same with a fixed isospin-forbidden d+d amplitude on level 1
                   (mix = 0.05 MeV^1/2): the truth of the recovery fit
  iso_fit.azr      every level parameter fixed, shifts started at -0.5 / -0.4 and the
                   mix at 0, all free; data segments read the data written by check.sh
  standalone_4H.azr  4H alone, one T=1 level with the formal parameters check.sh copies
                   from iso.azr's isospin.out (run with --no-transform)
  bad_undeclared.azr  level 1 spans three nuclei with no <isospin> block
  bad_radius.azr      n + t at 5.0 fm
  bad_mix.azr         mix on an isospin-allowed channel (p + t of the T=1 level)
  bad_T.azr           level 2 declared T=0 with a line in 4Li
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "models")

MP, MN, MD, MT, MH = 1.007825, 1.008665, 2.014102, 3.016049, 3.016029
P1 = dict(key=1, j1=0.5, pi1=1, j2=0.5, pi2=1, m1=MP, m2=MH, z1=1, z2=2, sep=0.0, rad=4.9)
P2 = dict(key=2, j1=0.5, pi1=1, j2=0.5, pi2=1, m1=MN, m2=MT, z1=0, z2=1, sep=0.0, rad=4.9)
P3 = dict(key=3, j1=0.5, pi1=1, j2=0.5, pi2=1, m1=MP, m2=MT, z1=1, z2=1, sep=19.8136, rad=4.9)
P4 = dict(key=4, j1=0.5, pi1=1, j2=0.5, pi2=1, m1=MN, m2=MH, z1=0, z2=2, sep=20.5776, rad=4.9)
P5 = dict(key=5, j1=1.0, pi1=1, j2=1.0, pi2=1, m1=MD, m2=MD, z1=1, z2=1, sep=23.8465, rad=7.0)
ENTRANCE_SEP = {1: 0.0, 2: 0.0, 3: 19.8136, 4: 19.8136, 5: 19.8136}  # first pair of each nucleus

CONFIG = """<config>
true #Perform A-Matrix Calculation
output/ #Full Path to Output Directory
checks/ #Full Path to Checks Directory
none #Compond Nucleus Check
none #Boundary Condition Check
none #Data Check
none #Lo-Matrix and Penetrability Check
none #Legendre Polynomial Check
none #Coulomb Amplitudes Check
none #Reaction Pathway Check
none #Angular Distributions Check
</config>"""

TEST = """1 1 1 2.0 12.0 0.05 120 120 0 1
1 2 2 2.0 12.0 0.05 0 0 0 0
1 3 4 1.2 10.0 0.05 0 0 0 0
1 5 3 0.1 4.0 0.02 0 0 0 0
1 5 4 0.1 4.0 0.02 0 0 0 0"""

DATA = """1 2 2 2.0 12.0 0 0 0 1 0 0 0 0 0 nt.dat
1 3 4 1.2 10.0 0 0 0 1 0 0 0 0 0 ptn.dat
1 5 3 0.1 4.0 0 0 0 1 0 0 0 0 0 ddp.dat
1 5 4 0.1 4.0 0 0 0 1 0 0 0 0 0 ddn.dat"""


def line(P, j, pi, e, s, l, lid, gamma, efix=0, chfix=0):
    return (f"{j:g} {pi:d} {e:.6f} {efix} 1 {P['key']} {int(2*s)} {int(2*l)} {lid} 1 {chfix} {gamma:.10g} "
            f"{P['j1']:g} {P['pi1']} {P['j2']:g} {P['pi2']} 0 {P['m1']} {P['m2']} {P['z1']} {P['z2']} "
            f"{ENTRANCE_SEP[P['key']]} {P['sep']} 0 0 0.0 0 {P['rad']} 0 0 0")


def levels(fixed=False, rad2=None, t0_in_4li=False):
    f = int(fixed)
    p2 = dict(P2, rad=rad2) if rad2 else P2
    L = []
    # level 1, 2- T=1: 4Li (reference), 4H, 4He
    L.append(line(P1, 2, -1, 5.0, 1, 1, 1, 2.0e6, f, f))
    L.append(line(p2, 2, -1, 4.14, 1, 1, 1, 2.0e6))
    L.append(line(P3, 2, -1, 24.1736, 1, 1, 1, 1.0e6))
    L.append(line(P4, 2, -1, 24.1736, 1, 1, 1, 1.0e6))
    L.append(line(P5, 2, -1, 24.1736, 1, 1, 1, 0.0))
    L.append("")
    # level 2, 2- T=0: 4He only
    if t0_in_4li:
        L.append(line(P1, 2, -1, 6.0, 1, 1, 2, 1.0e6))
    L.append(line(P3, 2, -1, 24.5, 1, 1, 2, 0.5e6, f, f))
    L.append(line(P4, 2, -1, 24.5, 1, 1, 2, -0.5e6, f, f))
    L.append(line(P5, 2, -1, 24.5, 1, 1, 2, 0.3e6, f, f))
    L.append("")
    return "<levels>\n" + "\n".join(L) + "\n</levels>"


def block(shift_h=-0.86, shift_he=-0.64, vary=0, mix=None, extra=""):
    b = ["<isospin>",
         "reference 3 4 1",
         f"shift 1 4 2 {shift_h} {vary} -0.86 0",
         f"shift 2 4 3 {shift_he} {vary} -0.64 0",
         "level 1 1",
         "level 2 0"]
    if mix is not None:
        b.append(f"mix 1 5 1 1 {mix[0]} {mix[1]} 0")
    if extra:
        b.append(extra)
    b.append("</isospin>")
    return "\n".join(b)


def model(name, levels_text, iso_text, test=TEST, data=""):
    parts = [CONFIG, levels_text]
    if iso_text:
        parts.append(iso_text)
    parts += ["<segmentsData>\n" + data + ("\n" if data else "") + "</segmentsData>",
              "<segmentsTest>\n" + test + "\n</segmentsTest>",
              "<targetInt>\n</targetInt>"]
    with open(os.path.join(OUT, name), "w") as f:
        f.write("\n".join(parts) + "\n")


def standalone_4h():
    # Placeholders: check.sh writes the formal E and gamma of iso.azr's 4H member in.
    lv = "<levels>\n" + line(P2, 2, -1, 0.0, 1, 1, 1, 0.0) + "\n\n</levels>"
    iso = "<isospin>\nreference 1 4 2\nlevel 1 1\n</isospin>"
    model("standalone_4H.azr", lv, iso, test="1 2 2 2.0 12.0 0.05 0 0 0 0")


def main():
    os.makedirs(OUT, exist_ok=True)
    model("iso.azr", levels(), block())
    model("iso_mix.azr", levels(), block(mix=(0.05, 0)))
    model("iso_fit.azr", levels(fixed=True), block(-0.5 + 0.0, -0.4, vary=1, mix=(0.0, 1)), test=TEST, data=DATA)
    standalone_4h()
    model("bad_undeclared.azr", levels(), "")
    model("bad_radius.azr", levels(rad2=5.0), block())
    model("bad_mix.azr", levels(), block(extra="mix 1 3 1 1 0.01 1 0"))
    model("bad_T.azr", levels(t0_in_4li=True), block())


if __name__ == "__main__":
    main()
