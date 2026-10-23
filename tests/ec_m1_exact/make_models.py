#!/usr/bin/env python3
"""Build the external-capture (EC) test models: p + 16O -> 17F(0.495 MeV, 1/2+) + gamma, external capture only.

The final state is an s-wave proton bound by 105 keV (ANC 84 fm^-1/2).  Entrance J-groups 1/2-, 3/2- (p-wave, E1) and
1/2+ (s-wave, M1) hold one far level each (30 MeV, proton width 1e-6 eV), so the resonant and resonant-external
terms are negligible (< 1e-12 relative) and the cross section is the hard-sphere external capture alone.  The photon
pair's multipolarity mask selects the term:

  e1.azr : E1 only (mask 1);   m1.azr : M1 only (mask 2)

Each has one angle-integrated segment of 16O(p,g)17F(0.495) at the energies in ECM.  truth.py computes the same cross
sections independently; check.sh runs both models with and without --no-long-wavelength.

Run from this directory: python3 make_models.py   (writes models/)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MP, M16 = 1.007276, 15.99053
ECM = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0]

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
# pair columns: j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask
PROTON = f"0.5 1 0 1 0 {MP} {M16} 1 8 0.599607 0.599607 0 0 0.0 0 5 5.5857 0 0"
def gamma(mask):
    return f"1 1 0.5 1 0.4953 0 17.0021 0 9 0.599607 0 0 0 0.0 10 0 0 0 {mask}"


def line(J, pi, E, key, s2, l2, lid, width, tail):
    """J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    return f"{J} {pi} {E} 1 1 {key} {s2} {l2} {lid} 1 1 {width} {tail}"


def model(name, mask):
    g = gamma(mask)
    levels = [
        [line(0.5, 1, 0.4953, 1, 1, 0, 1, 84, PROTON), line(0.5, 1, 0.4953, 3, 1, 2, 1, 0, g)],   # final state, ANC
        [line(0.5, -1, 30, 1, 1, 2, 2, 1e-6, PROTON), line(0.5, -1, 30, 3, 1, 2, 2, 0, g)],
        [line(1.5, -1, 30, 1, 1, 2, 3, 1e-6, PROTON), line(1.5, -1, 30, 3, 1, 2, 3, 0, g)],
        [line(0.5, 1, 30, 1, 1, 0, 4, 1e-6, PROTON), line(0.5, 1, 30, 3, 1, 2, 4, 0, g)],
    ]
    elab = [e * (MP + M16) / M16 for e in ECM]
    seg = f"1 1 3 {elab[0] - 0.01:.6f} {elab[-1] + 0.01:.6f} 0 0 0 1 0 0 0 0 0 data/ecm.dat 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in levels) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)
    with open(os.path.join(HERE, "models", "data", "ecm.dat"), "w") as f:
        for e in elab:
            f.write(f"{e:.10g}  0  1  0.1\n")


if __name__ == "__main__":
    model("e1", 1)
    model("m1", 2)
