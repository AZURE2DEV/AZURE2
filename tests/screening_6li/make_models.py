#!/usr/bin/env python3
"""Models and data of tests/screening_6li: electron screening of 6Li(p,a)3He, Cruz et al. 2005 (Phys. Lett. B 624, 181;
EXFOR O1282005-007) on a Li2WO4 insulator, Li metal and a PdLi alloy, against J. Skowronski's 7Be R-matrix evaluation held
fixed.

The levels are that evaluation's best fit (Supplemental 7Be.azr + param.sav, written into <levels> with pyazr save_fit and
checked on re-read), every energy and width fixed; the evaluation's own segments are removed (it has no 6Li(p,a) data below
E_p = 100 keV, so these points are out of sample).  Only the three Cruz sets remain, each with its own <targetInt> line
carrying a screening block.  Written next to this script:

  data/cruz_{metal,insulator,pdli}.dat   E_lab (MeV), 0, sigma (b), dsigma (b) from the EXFOR S-factors
  bare.azr     U_e = 0 on all three, norms 1
  scr.azr      U_e and norms at the benchmark minimum, fixed
  fit.azr      U_e (start 500 eV) and the norms free, no priors

Run:  python3 tests/screening_6li/make_models.py path/to/7Be_fit.azr
(the snapshot is R-matrix/6Li+p/3He+a/AZURE2_evaluation/10-21-26_screening_benchmark/7Be_fit.azr in the deBoer archive.)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

M1, M2 = 1.00782, 6.0151                            # p + 6Li in the evaluation's pair lines
UCONV, FSTRUC = 931.4940880, 1.0 / 137.0359996790   # include/Constants.h

# EXFOR O1282005 (Li metal), 006 (Li2WO4), 007 (PdLi): E (keV, c.m. -- EXFOR tags the column EN, but the paper's
# 30-100 keV proton beam and its S(E_cm) figures make it c.m.; read as lab, the metals come out 26-27 % below the
# published U_e), S and dS (MeV b, statistical).  Published U_e (ADD-RES): 1320 +- 70, 320 +- 110, 3760 +- 260 eV.
SETS = [
    ("metal", [(24.2, 5.06, .14), (28.3, 4.64, .11), (32.5, 4.61, .09), (39.7, 4.02, .04), (48.0, 3.81, .03),
               (56.4, 3.69, .02), (64.8, 3.50, .04), (73.2, 3.41, .05), (81.7, 3.42, .09)]),
    ("insulator", [(28.3, 3.58, .08), (32.5, 3.69, .10), (39.8, 3.52, .10), (48.0, 3.45, .05), (56.4, 3.34, .05),
                   (64.8, 3.32, .06), (73.2, 3.26, .04)]),
    ("pdli", [(31.3, 7.45, .23), (39.6, 6.00, .21), (47.9, 5.13, .15), (56.3, 4.81, .18), (64.6, 3.98, .22),
              (73.1, 3.51, .29), (81.6, 3.61, .22)]),
]
# The benchmark minimum (chi2 with each norm profiled; AZURE2 dev 334f104, engine and pyazr fits agree to 1e-4 eV).
BEST_UE = (1211.1217, 322.6758, 3432.9061)
BEST_NORM = (1.015325, 1.010056, 1.034700)


def eta(ecm):
    return math.sqrt(UCONV / 2.0) * FSTRUC * 3 * math.sqrt(M1 * M2 / (M1 + M2) / ecm)


def write_data():
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    for key, rows in SETS:
        out = []
        for e_kev, s, ds in rows:
            ecm = e_kev * 1e-3
            g = math.exp(-2 * math.pi * eta(ecm)) / ecm
            out.append(f"{ecm * (M1 + M2) / M2:.8e}  0  {s * g:.8e}  {ds * g:.8e}")
        open(os.path.join(HERE, "data", f"cruz_{key}.dat"), "w").write("\n".join(out) + "\n")


def model(snapshot, ues, norms, vary):
    """The snapshot with every level parameter fixed, its segments, test segments, parameter settings and MCMC block
    dropped, and the three Cruz segments with one screening line each."""
    out, block = [], None
    for line in open(snapshot).read().splitlines():
        s = line.strip()
        if s in ("<parameterSettings>", "<mcmc>"):
            block = "drop"
            continue
        if block == "drop":
            if s in ("</parameterSettings>", "</mcmc>"):
                block = None
            continue
        if s in ("<levels>", "<segmentsData>", "<segmentsTest>", "<targetInt>"):
            block = s
            out.append(line)
            if s == "<segmentsData>":
                for (key, _), n in zip(SETS, norms):
                    out.append(f"1  2  1  0  0.2  0  180  0  {n:.6f}  {int(vary)}  0  0  0.05  0  data/cruz_{key}.dat  0  0")
            if s == "<targetInt>":
                for i, ue in enumerate(ues):
                    out.append(f'1  "{i + 1}"  200  0  0  0  0  ""  0  0  0  0  ""  0  0  0.04  20  50  '
                               f"screening  {ue:.4f}  {int(vary)}  0")
            continue
        if s.startswith("</"):
            block = None
        elif block == "<levels>" and s:
            tok = line.split()
            tok[3] = "1"       # level energy fixed
            tok[10] = "1"      # channel width fixed
            line = "  ".join(tok)
        elif block in ("<segmentsData>", "<segmentsTest>", "<targetInt>") and s:
            continue
        elif s.startswith("output/") or s.startswith("checks/"):
            pass
        out.append(line)
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    snap = sys.argv[1]
    write_data()
    for name, ues, norms, vary in (("bare.azr", (0.0,) * 3, (1.0,) * 3, False),
                                   ("scr.azr", BEST_UE, BEST_NORM, False),
                                   ("fit.azr", (500.0,) * 3, (1.0,) * 3, True)):
        open(os.path.join(HERE, name), "w").write(model(snap, ues, norms, vary))
        print("wrote", name)
