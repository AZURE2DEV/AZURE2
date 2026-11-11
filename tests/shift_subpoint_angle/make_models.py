#!/usr/bin/env python3
"""Models for the shifted sub-point angle test (tests/shift_subpoint_angle).

template.azr is the 16O + p model of the 2026-10-25 audit (p + 16O elastic, 17F g.s. and 0.495 MeV capture pairs),
with the segments and target-effect blocks replaced here.  Every model has a 0.1 keV Gaussian energy convolution
(20 sub-points) on all of its segments:

  a0 : 16O(p,p) lab-frame differential (Harris data, 11 angles), energy shift 0
  a1 : the same with a 1e-7 MeV shift on every segment
  av : the same with shift 0 but VARIED (the fit path is active from the first evaluation)
  g0 : a0 plus a lab-frame differential 16O(p,g0) segment (made-up data, data/pg.dat), shift 0
  g1 : g0 with the 1e-7 MeV shift on the capture segment

Run from this directory: python3 make_models.py   (writes models/)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PP = "1 1 1 0 10 {ang} {ang} 1 1 0 0 {shift} 0 {vary} data/Harris_pp.dat"
PG = "1 1 {gkey} 0 10 {ang} {ang} 1 1 0 0 {shift} 0 0 data/pg.dat"
ANGLES = [37.05, 51.89, 60.68, 86.39, 99.9, 111.29, 122.2, 128.56, 138.36, 147.49, 160.05]
GANGLES = [45.0, 90.0, 135.0]
MP, M16 = 1.007276, 15.99053


def build(name, shift, vary, gamma_shift=None):
    tpl = open(os.path.join(HERE, "template.azr")).read()
    head = tpl.split("<segmentsData>")[0]
    segs = [PP.format(ang=a, shift=shift, vary=vary) for a in ANGLES]
    if gamma_shift is not None:
        segs += [PG.format(gkey=GKEY, ang=a, shift=gamma_shift) for a in GANGLES]
    n = len(segs)
    txt = (head + "<segmentsData>\n" + "\n".join(segs) + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n"
           + f'<targetInt>\n1 "1-{n}" 20 1 0.0001 0 0 "" 0 0 0 0 "" 0\n</targetInt>\n')
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


if __name__ == "__main__":
    # the g.s. gamma pair key from the template (pType 10, Ex 0)
    GKEY = None
    for l in open(os.path.join(HERE, "template.azr")):
        t = l.split()
        if len(t) >= 31 and t[26] == "10" and float(t[16]) == 0.0:
            GKEY = int(t[5]); break
    assert GKEY, "no g.s. gamma pair in template.azr"
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    with open(os.path.join(HERE, "data", "pg.dat"), "w") as f:
        for e in (1.0, 1.5, 2.0, 2.5, 3.0):
            for a in GANGLES:
                f.write(f"{e * (MP + M16) / M16:.6f}  {a}  1.0e-6  1.0e-7\n")
    build("a0", 0, 0)
    build("a1", "1e-07", 0)
    build("av", 0, 1)
    build("g0", 0, 0, gamma_shift=0)
    build("g1", 0, 0, gamma_shift="1e-07")
