#!/usr/bin/env python3
"""Model for tests/ec_threshold: 16O(n,g)17O(g.s.) across the alpha + 13C threshold of 17O.

Pairs: n + 16O (key 1, entrance), alpha + 13C (key 2; opens at E_cm(n) = 6.35869 - 4.14313 = 2.21556 MeV), photon to the
17O 5/2+ ground state (key 3, E1 external capture).  The g.s. has ANCs in both particle channels; a 5/2- level at
Ex = 9 MeV carries n (l = 3), alpha (l = 2) and an E1 photon width, so the capture includes channel capture through
the alpha + 13C channel, whose Coulomb functions over- and underflow within ~1 keV above its threshold.
The data grid is dense across the threshold (E_cm 2.195-2.235 MeV, 0.25 keV steps, plus points 0.01-0.5 keV around it).

Run from this directory: python3 make_models.py   (writes models/)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MN, M16, SN, SA = 1.00866, 15.9949, 4.14313, 6.35869
THR = SA - SN

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
N16 = f"0.5 1 0 1 0 {MN} {M16} 0 8 {SN} {SN} 0 0 0.0 0 4.4 0 0 0"
A13 = f"0 1 0.5 -1 0 4.0026 13.0033 2 6 {SN} {SA} 0 0 0.0 0 5.4 0 0 0"
GAM = f"1 1 2.5 1 0 0 17.0 0 8 {SN} 0 0 0 0.0 10 0 0 0 1"


def line(J, pi, E, key, s2, l2, lid, width, tail):
    """J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    return f"{J} {pi} {E} 1 1 {key} {s2} {l2} {lid} 1 1 {width} {tail}"


LEVELS = [
    [line(2.5, 1, 0, 1, 1, 4, 1, 0.9, N16), line(2.5, 1, 0, 2, 1, 6, 1, 3.0, A13)],          # g.s. ANCs (fm^-1/2)
    [line(2.5, -1, 9, 1, 1, 6, 2, 3e5, N16), line(2.5, -1, 9, 2, 1, 4, 2, 5e4, A13),
     line(2.5, -1, 9, 3, 1, 2, 2, 1.0, GAM)],                                                  # 5/2-: n, alpha, E1
]

if __name__ == "__main__":
    ecm = [round(2.195 + 0.00025 * i, 7) for i in range(161)]
    ecm += [THR + d for d in (-5e-4, -2e-4, -1e-4, -1e-5, 1e-5, 1e-4, 2e-4, 3e-4, 5e-4)]
    ecm = sorted(set(ecm))
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    with open(os.path.join(HERE, "models", "data", "ng.dat"), "w") as f:
        for e in ecm:
            f.write(f"{e * (MN + M16) / M16:.10f}  0  1e-6  1e-7\n")
    elab0, elab1 = ecm[0] * (MN + M16) / M16, ecm[-1] * (MN + M16) / M16
    seg = f"1 1 3 {elab0 - 1e-4:.6f} {elab1 + 1e-4:.6f} 0 0 0 1 0 0 0 0 0 data/ng.dat 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in LEVELS) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", "ng.azr"), "w").write(txt)
