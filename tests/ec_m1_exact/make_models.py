#!/usr/bin/env python3
"""Build the external-capture (EC) test models: p + 16O -> 17F(0.495 MeV, 1/2+) + gamma, external capture only.

The final state is an s-wave proton bound by 105 keV (ANC 84 fm^-1/2).  Entrance J-groups 1/2-, 3/2- (p-wave, E1) and
1/2+ (s-wave, M1) hold one far level each (30 MeV, proton width 1e-6 eV), so the resonant and resonant-external
terms are negligible (< 1e-12 relative) and the cross section is the hard-sphere external capture alone.  The photon
pair's multipolarity mask selects the term:

  e1.azr : E1 only (mask 1);   m1.azr : M1 only (mask 2)

Each has one angle-integrated segment of 16O(p,g)17F(0.495) at the energies in ECM.

  a12c.azr : 12C(a,g)16O(g.s.), external capture only, E1 (l = 1) and E2 (l = 2) (mask 5), ANC 100 fm^-1/2,
             a = 5.43 fm, one 30 MeV level (alpha width 1e-6 eV) in 1- and in 2+; energies ECM_A.  An N = Z pair: the
             long-wavelength E1 charge nearly vanishes (-0.00098), so the exact form and its convection-current term
             are large here.

  p14n_s1.azr, p14n_s3.azr : p + 14N -> 15O(1/2+) or 15O(3/2+), a fictitious s-wave proton state bound by 0.2 MeV
             (ANC 5 fm^-1/2, a = 5 fm), M1 external capture only (mask 2) from the s-wave entrance J-groups 1/2+ and
             3/2+.  Both clusters carry a moment (g_p 5.5857, g(14N) 0.4038) and the target spin is 1, so the M1
             spin term has a g1 S1 and a g2 S2 piece and channel spin 1/2 <-> 3/2 transitions.  With l = 0 the
             orbital term vanishes.  Final channel spin 1/2 (s1) or 3/2 (s3).  Energies ECM.

truth.py computes the same cross sections independently; check.sh runs the models with and without
--no-long-wavelength and --no-ec-convection-current.

Run from this directory: python3 make_models.py   (writes models/)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MP, M16 = 1.007276, 15.99053
ECM = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0]
MA, M12, SA = 4.002603, 12.0, 7.16192
M14, M15, SP15, B15, ANC15, G14N = 13.999234, 15.003066, 7.2971, 0.2, 5.0, 0.4038
ECM_A = [0.6, 1.0, 1.5, 2.0, 3.0]

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


def model_a12c():
    alpha = f"0 1 0 1 0 {MA} {M12} 2 6 {SA} {SA} 0 0 0.0 0 5.43 0 0 0"
    gam = f"1 1 0 1 0 0 15.9949 0 8 {SA} 0 0 0 0.0 10 0 0 0 5"
    levels = [
        [line(0, 1, 0, 1, 0, 0, 1, 100, alpha)],                                                  # 16O g.s., ANC
        [line(1, -1, 30, 1, 0, 2, 2, 1e-6, alpha), line(1, -1, 30, 2, 0, 2, 2, 0, gam)],          # E1 entrance
        [line(2, 1, 30, 1, 0, 4, 3, 1e-6, alpha), line(2, 1, 30, 2, 0, 4, 3, 0, gam)],            # E2 entrance
    ]
    elab = [e * (MA + M12) / M12 for e in ECM_A]
    seg = f"1 1 2 {elab[0] - 0.01:.6f} {elab[-1] + 0.01:.6f} 0 0 0 1 0 0 0 0 0 data/ecm_a.dat 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in levels) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", "a12c.azr"), "w").write(txt)
    with open(os.path.join(HERE, "models", "data", "ecm_a.dat"), "w") as f:
        for e in elab:
            f.write(f"{e:.10g}  0  1  0.1\n")


def model_p14n(name, jf):
    """p + 14N -> 15O(jf+), s-wave final state, M1 only (mask 2); entrance J-groups 1/2+ and 3/2+ (s-wave)."""
    ex = SP15 - B15
    prot = f"0.5 1 1 1 0 {MP} {M14} 1 7 {SP15} {SP15} 0 0 0.0 0 5 5.5857 {G14N} 0"
    gam = f"1 1 {jf} 1 {ex:.4f} 0 {M15} 0 8 {SP15} 0 0 0 0.0 10 0 0 0 2"
    s2f = int(round(2 * jf))
    levels = [[line(jf, 1, f"{ex:.4f}", 1, s2f, 0, 1, ANC15, prot)]]                                  # final state, ANC
    for n, J in enumerate((0.5, 1.5), start=2):
        lv = [line(J, 1, 30, 1, int(round(2 * J)), 0, n, 1e-6, prot), line(J, 1, 30, 3, s2f, 2, n, 0, gam)]
        if J == jf:   # the final state's own J-group: the 30 MeV level shares it
            levels[0].append(line(jf, 1, f"{ex:.4f}", 3, s2f, 2, 1, 0, gam))
        levels.append(lv)
    elab = [e * (MP + M14) / M14 for e in ECM]
    seg = f"1 1 3 {elab[0] - 0.01:.6f} {elab[-1] + 0.01:.6f} 0 0 0 1 0 0 0 0 0 data/ecm_p14n.dat 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in levels) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)
    with open(os.path.join(HERE, "models", "data", "ecm_p14n.dat"), "w") as f:
        for e in elab:
            f.write(f"{e:.10g}  0  1  0.1\n")


if __name__ == "__main__":
    model("e1", 1)
    model("m1", 2)
    model_a12c()
    model_p14n("p14n_s1", 0.5)
    model_p14n("p14n_s3", 1.5)
