#!/usr/bin/env python3
"""Models for the Gaussian kernel-window test (tests/convolution_window).

Base: tests/ec_m1_exact/models/e1.azr (p + 16O -> 17F(0.495) + gamma; here with external capture OFF, mask 0, so the
only cross section is the R-matrix one), plus one NARROW 1/2- resonance at E_x = EX_NARROW (E_cm EX_NARROW - 0.5996,
E_lab x 17.0078/16): proton width GP eV, gamma width GG eV.  Far from it the capture cross section is the 30 MeV
levels' tails, ~10^-9 b, so the resonance is >> 10^3 x its surroundings -- the situation in which a Gaussian kernel
cut at 3 sigma and renormalised (the pre-2026-10-26 TargetEffect::convolutionRange = 3) loses real cross section.

  bare.azr : the angle-integrated 16O(p,g) excitation function on a fine lab-energy grid (STEP) over EMIN..EMAX;
             truth.py convolves it with full Gaussians;
  conv.azr : the same at the test points E_res + k sigma_lab, k = KS, with a Gaussian beam spread SIGMA_LAB
             (constant-sigma field) and NSUB sub-points.

Run from this directory: python3 make_models.py   (writes models/)
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..", "ec_m1_exact", "models", "e1.azr")
MP, M16 = 1.007276, 15.99053
F = M16 / (MP + M16)
SP, EX_NARROW, GP, GG = 0.599607, 1.4, 1000.0, 1.0
ERES_LAB = (EX_NARROW - SP) / F
SIGMA_LAB = 0.003
KS = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
EMIN, EMAX, STEP = round(ERES_LAB - 0.05, 4), round(ERES_LAB + 0.05, 4), 0.00005
NSUB = 50
PROTON = f"0.5 1 0 1 0 {MP} {M16} 1 8 {SP} {SP} 0 0 0.0 0 5 5.5857 0 0"
GAMMA = f"1 1 0.5 1 0.4953 0 17.0021 0 9 {SP} 0 0 0 0.0 10 0 0 0 0"      # ecMultMask 0: no external capture
LEVEL = [f"0.5 -1 {EX_NARROW} 1 1 1 1 2 5 1 1 {GP} {PROTON}", f"0.5 -1 {EX_NARROW} 1 1 3 1 2 5 1 1 {GG} {GAMMA}"]


def model(name, segments, target_line):
    txt = open(BASE).read()
    txt = re.sub(r"(\S+ \S+ \S+ 0 17\.0021 0 9 \S+ 0 0 0 0\.0 10 0 0 0) 1$", r"\1 0", txt, flags=re.M)   # mask 1 -> 0
    txt = txt.replace("</levels>", "\n".join(LEVEL) + "\n\n</levels>")
    txt = re.sub(r"<segmentsData>.*?</segmentsData>", "<segmentsData>\n</segmentsData>", txt, flags=re.S)
    txt = re.sub(r"<segmentsTest>.*?</segmentsTest>", "<segmentsTest>\n" + "\n".join(segments) + "\n</segmentsTest>", txt, flags=re.S)
    txt = re.sub(r"<targetInt>.*?</targetInt>", "<targetInt>\n" + target_line + "\n</targetInt>", txt, flags=re.S)
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    model("bare", [f"1 1 3 {EMIN} {EMAX} {STEP} 0 0 0 0"], "")
    pts = [ERES_LAB + k * SIGMA_LAB for k in KS]
    segs = [f"1 1 3 {e:.6f} {e:.6f} 0.001 0 0 0 0" for e in pts]
    model("conv", segs, f'1 "1-{len(segs)}" {NSUB} 1 {SIGMA_LAB} 0 0 "" 0 0 0 0 "" 0 0 0.04 20 50')
    with open(os.path.join(HERE, "models", "points.txt"), "w") as f:
        for k, e in zip(KS, pts):
            f.write(f"{k:5.1f}  {e:.6f}\n")
    print("E_res(lab) =", round(ERES_LAB, 6), "grid", EMIN, EMAX, STEP, "sigma_lab", SIGMA_LAB)
