#!/usr/bin/env python3
"""Models for tests/ratio_frame (audit item 15, 2026-10-27): a lab-frame RATIO composite is the ratio of the two
terms' LAB cross sections, each term taken to the lab with its OWN pair's Jacobian.

Two reactions, all segments test segments (mode 3; the ratio is read off the lab output, kinflag.dat = 1, and
compared with the ratio of the plain segments' lab values):

  particle      13C + alpha in 17O (the tests/sum_lab_frame model: exits n + 16O(6.05, 0+) key 3 and
                n + 16O(6.13, 3-) key 5, elastic key 2; two levels, no external capture), lab 5.95-6.30 MeV:
                T1 a -> n0+ (0-180 deg), T2 a -> n3- (0-180), T3 a -> a at 60 deg, T4 a -> n3- at 30 deg,
                T5 a -> n3- at 0 deg, T6/T7 the c.m. (isDiff 5) versions of T1/T2, T8/T9 angle-integrated n0+/n3-;
                R1 = T1 / T2 (same angles), R2 = T1 / (a,a)(60), R3 = T1 / (a,n3-)(30), R4 = T1 / (a,n3-)(0),
                R5 = T1 / (a,n0+)(60) (the same pair), R6 = c.m. ratio T6 / T7, R7 = angle-integrated T8 / T9
  particle_rel  the same with relativistic=1
  gamma         12C(a,g)16O external capture (the tests/gamma_lab_frame model: alpha key 1, gamma key 2; E1 + E2
                interference), lab 2-6 MeV: G1 a -> g (10-170 deg), G2 a -> g at 90 deg, G3 a -> a at 60 deg,
                G4 a -> a (10-170); RG1 = G1 / G2 (W(theta)/W(90)), RG2 = G1 / (a,a)(60), RG3 = G4 / (a,g)(90)
  ratio_data    mode 1, relativistic=1: the RATIO data segment n0+ / (a,a)(60 deg) on the particle grid, reading
                data/ratio.dat (written by check.sh from the plain segments' lab values: chi2 must be ~0)
  ratio_free    ratio_data with the seven widths free, for the AZURE_GRAD_CHECK gradient test

Run from this directory: python3 make_models.py  (writes models/).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

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
{extra}</config>
"""

# ---- 13C + alpha (tests/sum_lab_frame) ------------------------------------------------------------------------------
ALPHA = "0  1  0.5  -1  0  4.0026  13.0033  2  6  4.14313  6.35869  0  0  0.0  0  5.4  0  0  0"
N_0P = "0.5  1  0  1  6.0494  1.00866  15.9949  0  8  4.14313  4.14313  0  0  0.0  0  4.4  0  0  0"
N_3M = "0.5  1  3  -1  6.12989  1.00866  15.9949  0  8  4.14313  4.14313  0  0  0.0  0  4.4  0  0  0"


def channel(J, E, key, s, l, width, tail, lid=1, free=False):
    """One <levels> line: J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    fix = 0 if free else 1
    return f"{J}  1  {E}  1  1  {key}  {int(round(2 * s))}  {2 * l}  {lid}  1  {fix}  {width}  {tail}"


def particle_levels(free=False):
    return [
        [channel(2.5, 11.0, 2, 0.5, 3, 3000.0, ALPHA, 1, free), channel(2.5, 11.0, 3, 0.5, 2, 15000.0, N_0P, 1, free),
         channel(2.5, 11.0, 5, 2.5, 1, 20000.0, N_3M, 1, free), channel(2.5, 11.0, 5, 3.5, 1, -10000.0, N_3M, 1, free)],
        [channel(1.5, 11.1, 2, 0.5, 1, 5000.0, ALPHA, 2, free), channel(1.5, 11.1, 3, 0.5, 2, -12000.0, N_0P, 2, free),
         channel(1.5, 11.1, 5, 2.5, 1, -15000.0, N_3M, 2, free)],
    ]


PARTICLE_LEVELS = particle_levels()

# ---- 12C(a,g) (tests/gamma_lab_frame) --------------------------------------------------------------------------------
MA, M12, SA = 4.002603, 12.0, 7.16192
G_ALPHA = f"0 1 0 1 0 {MA} {M12} 2 6 {SA} {SA} 0 0 0.0 0 5.43 0 0 0"
G_GAMMA = f"1 1 0 1 0 0 15.9949 0 8 {SA} 0 0 0 0.0 10 0 0 0 5"


def gline(J, pi, E, key, s2, l2, lid, width, tail):
    return f"{J} {pi} {E} 1 1 {key} {s2} {l2} {lid} 1 1 {width} {tail}"


GAMMA_LEVELS = [
    [gline(0, 1, 0, 1, 0, 0, 1, 100, G_ALPHA)],
    [gline(1, -1, 30, 1, 0, 2, 2, 1e-6, G_ALPHA), gline(1, -1, 30, 2, 0, 2, 2, 0, G_GAMMA)],
    [gline(2, 1, 30, 1, 0, 4, 3, 1e-6, G_ALPHA), gline(2, 1, 30, 2, 0, 4, 3, 0, G_GAMMA)],
]


def write(name, levels, tests, extra="", data=()):
    txt = (CONFIG.format(extra=extra) + "<levels>\n" + "\n\n".join("\n".join(b) for b in levels)
           + "\n\n</levels>\n<segmentsData>\n" + "".join(s + "\n" for s in data) + "</segmentsData>\n<segmentsTest>\n"
           + "".join(s + "\n" for s in tests) + "</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


def main():
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    # test line: active entrance exit minE maxE stepE minA maxA stepA isDiff [isAdvanced op ncomp ent exit angle]
    e = "5.95  6.30  0.05"
    full = f"{e}  0  180  15"
    particle = [f"1  2  3  {full}  1",                      # T1 a -> n0+, lab
                f"1  2  5  {full}  1",                      # T2 a -> n3-, lab
                f"1  2  2  {e}  60  60  15  1",             # T3 a -> a at 60 deg
                f"1  2  5  {e}  30  30  15  1",             # T4 a -> n3- at 30 deg
                f"1  2  5  {e}  0  0  15  1",               # T5 a -> n3- at 0 deg
                f"1  2  3  {full}  5",                      # T6 c.m. n0+
                f"1  2  5  {full}  5",                      # T7 c.m. n3-
                f"1  2  3  {e}  0  0  15  0",               # T8 angle-integrated n0+
                f"1  2  5  {e}  0  0  15  0",               # T9 angle-integrated n3-
                f"1  2  3  {full}  1  1  1  1  2  5  -999",  # R1 n0+ / n3-, same angles
                f"1  2  3  {full}  1  1  1  1  2  2  60",    # R2 n0+ / (a,a)(60)
                f"1  2  3  {full}  1  1  1  1  2  5  30",    # R3 n0+ / n3-(30)
                f"1  2  3  {full}  1  1  1  1  2  5  0",     # R4 n0+ / n3-(0)
                f"1  2  3  {full}  1  1  1  1  2  3  60",    # R5 n0+ / n0+(60), the same pair
                f"1  2  3  {full}  5  1  1  1  2  5  -999",  # R6 c.m. ratio n0+ / n3-
                f"1  2  3  {e}  0  0  15  0  1  1  1  2  5  -999"]  # R7 angle-integrated ratio
    write("particle", PARTICLE_LEVELS, particle)
    write("particle_rel", PARTICLE_LEVELS, particle, "relativistic=1\n")

    ge = "2.0  6.0  1.0"
    gfull = f"{ge}  10  170  40"
    gamma = [f"1  1  2  {gfull}  1",                       # G1 a -> g, lab
             f"1  1  2  {ge}  90  90  40  1",              # G2 a -> g at 90 deg
             f"1  1  1  {ge}  60  60  40  1",              # G3 a -> a at 60 deg
             f"1  1  1  {gfull}  1",                       # G4 a -> a, lab
             f"1  1  2  {gfull}  1  1  1  1  1  2  90",     # RG1 W(theta) / W(90)
             f"1  1  2  {gfull}  1  1  1  1  1  1  60",     # RG2 (a,g) / (a,a)(60)
             f"1  1  1  {gfull}  1  1  1  1  1  2  90"]     # RG3 (a,a) / (a,g)(90)
    write("gamma", GAMMA_LEVELS, gamma)

    # mode-1 models (relativistic=1): a RATIO data segment n0+ / (a,a)(60 deg) on the extrapolation grid;
    # check.sh writes data/ratio.dat from the plain segments' lab values of the particle_rel run
    # data line: active entrance exit minE maxE minA maxA isDiff norm varyNorm normErr shift varyShift ? file [composite]
    dseg = "5.9  6.4  0  180  1  1.0  0  0  0  0  0  data/ratio.dat  1  1  1  2  2  60"
    write("ratio_data", particle_levels(False), [], "relativistic=1\n", [f"1  2  3  {dseg}"])
    write("ratio_free", particle_levels(True), [], "relativistic=1\n", [f"1  2  3  {dseg}"])


if __name__ == "__main__":
    main()
