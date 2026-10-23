#!/usr/bin/env python3
"""Models for tests/gamma_lab_frame: 12C(a,g)16O(g.s.) external capture only (E1 from l = 1 and E2 from l = 2, ANC
100 fm^-1/2, a = 5.43 fm, one 30 MeV level with a 1e-6 eV alpha width in 1- and in 2+), whose E1-E2 interference
makes the angular distribution strongly asymmetric.

  cm.azr  : one c.m. differential segment (isDiff 4) at E_LAB x ANGLES (data/cm.dat, dummy cross sections)
  lab.azr : the same reaction as a lab differential segment (isDiff 1) reading data/lab.dat, which check.sh writes
            from cm's output with an explicit Lorentz boost of the photon

Run from this directory: python3 make_models.py   (writes models/)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MA, M12, SA = 4.002603, 12.0, 7.16192
E_LAB = [2.0, 3.0, 4.5, 6.0]
ANGLES = [10, 40, 70, 90, 110, 140, 170]

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
</config>
"""
ALPHA = f"0 1 0 1 0 {MA} {M12} 2 6 {SA} {SA} 0 0 0.0 0 5.43 0 0 0"
GAMMA = f"1 1 0 1 0 0 15.9949 0 8 {SA} 0 0 0 0.0 10 0 0 0 5"


def line(J, pi, E, key, s2, l2, lid, width, tail):
    """J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    return f"{J} {pi} {E} 1 1 {key} {s2} {l2} {lid} 1 1 {width} {tail}"


LEVELS = [
    [line(0, 1, 0, 1, 0, 0, 1, 100, ALPHA)],
    [line(1, -1, 30, 1, 0, 2, 2, 1e-6, ALPHA), line(1, -1, 30, 2, 0, 2, 2, 0, GAMMA)],
    [line(2, 1, 30, 1, 0, 4, 3, 1e-6, ALPHA), line(2, 1, 30, 2, 0, 4, 3, 0, GAMMA)],
]


def write(name, isdiff, datafile):
    seg = f"1 1 2 {E_LAB[0] - 0.01} {E_LAB[-1] + 0.01} 0 180 {isdiff} 1 0 0 0 0 0 data/{datafile} 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in LEVELS) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    write("cm", 4, "cm.dat")
    write("lab", 1, "lab.dat")
    with open(os.path.join(HERE, "models", "data", "cm.dat"), "w") as f:
        for e in E_LAB:
            for a in ANGLES:
                f.write(f"{e}  {a}  1e-9  1e-10\n")
