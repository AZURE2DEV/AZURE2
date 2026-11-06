#!/usr/bin/env python3
"""Models for ec_e1m1_phase (committed): p + 16O external capture with E1 and M1 into

  s : 17F(0.4953, 1/2+), the 2s_1/2 state of tests/ec_m1_exact (ANC 84 fm^-1/2, bound by 104 keV): E1 from l = 1
      (J = 1/2-, 3/2-), M1 from l = 0 (J = 1/2+, the spin term alone).  For a spin-0 target and an unpolarized beam the
      two cannot interfere in dsigma/dOmega: the E1 amplitude is spin-diagonal and the M1 spin operator is traceless.
  p : a fictitious 1/2- state at E_x = 0.4996 MeV (l_f = 1, bound by 100 keV, ANC 20): E1 from l = 0 (J = 1/2+) and
      l = 2 (J = 3/2+), M1 from l = 1 (J = 1/2-, 3/2-) with the orbital term g_l L and the spin term g_p S_p; the
      orbital M1 amplitude is spin-diagonal and interferes with E1 (odd Legendre orders).

Each with the photon pair's multipolarity mask 1 (E1), 2 (M1) and 3 (both), and one c.m. differential capture
segment (isDiff 4) at the energies and angles of truth_dsdo.py.  Entrance J-groups hold one 30 MeV level with a
1e-6 eV proton width (negligible internal amplitude), as in tests/ec_m1_exact.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MP, M16 = 1.007276, 15.99053
SEP = 0.599607
ENERGIES = (0.3, 1.0, 2.0)
ANGLES = (0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0)
# J, pi, E_x, l_f, ANC (fm^-1/2).  The p-wave state needs a smaller ANC than the 17F* halo: 84 fm^-1/2 on a state bound by
# 100 keV with l_f = 1 would put 130 % of the probability outside 5 fm (AZURE2 then warns "Denominator less than zero while
# transforming" and the cross section no longer scales as ANC^2); 20 fm^-1/2 puts 7 % outside.
FINAL = {"s": (0.5, 1, 0.4953, 0, 84.0), "p": (0.5, -1, 0.4996, 1, 20.0)}
CONFIG = """<config>
true                                                                                                #Perform A-Matrix Calculation
output/                                                                                             #Full Path to Output Directory
checks/                                                                                             #Full Path to Checks Directory
none                                                                                                #Compond Nucleus Check
none                                                                                                #Boundary Condition Check
none                                                                                                #Data Check
none                                                                                                #Lo Matrix Check
none                                                                                                #Legendre Polynomials Check
none                                                                                                #External Capture Check
none                                                                                                #Coulomb Amplitude Check
none                                                                                                #Angular Distribution Check
</config>
"""
# pair columns: j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask
PROTON = f"0.5 1 0 1 0 {MP} {M16} 1 8 {SEP} {SEP} 0 0 0.0 0 5 5.5857 0 0"


def gamma(final, mask):
    J, pi, ex, _, _ = FINAL[final]
    return f"1 1 {J} {pi} {ex} 0 17.0021 0 9 {SEP} 0 0 0 0.0 10 0 0 0 {mask}"


def line(J, pi, E, key, s2, l2, lid, width, tail):
    """J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    return f"{J} {pi} {E} 1 1 {key} {s2} {l2} {lid} 1 1 {width} {tail}"


def model(final, mask):
    J, pi, ex, lf, anc = FINAL[final]
    g = gamma(final, mask)
    levels = [[line(J, pi, ex, 1, 1, 2 * lf, 1, anc, PROTON), line(J, pi, ex, 3, 1, 2, 1, 0, g)]]
    if final == "s":
        entrance = [(0.5, -1, 1), (1.5, -1, 1), (0.5, 1, 0)]          # E1 from l = 1, M1 from l = 0
    else:
        entrance = [(0.5, 1, 0), (1.5, 1, 2), (0.5, -1, 1), (1.5, -1, 1)]   # E1 from l = 0, 2; M1 from l = 1
    for n, (Je, pe, le) in enumerate(entrance, 2):
        levels.append([line(Je, pe, 30, 1, 1, 2 * le, n, 1e-6, PROTON), line(Je, pe, 30, 3, 1, 2, n, 0, g)])
    seg = "1 1 3 0.1 3.0 0 180 4 1 0 0 0 0 0 data/dsdo.dat 0 0"
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in levels) + "\n\n</levels>\n<segmentsData>\n"
           + seg + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    with open(os.path.join(HERE, "models", f"{final}_{mask}.azr"), "w") as f:
        f.write(txt)


if __name__ == "__main__":
    for final in FINAL:
        for mask in (1, 2, 3):
            model(final, mask)
    with open(os.path.join(HERE, "models", "data", "dsdo.dat"), "w") as f:
        for E in ENERGIES:
            for th in ANGLES:
                f.write(f"{E * (MP + M16) / M16:.9f} {th:.1f} 1 0.1\n")
    print("models written")
