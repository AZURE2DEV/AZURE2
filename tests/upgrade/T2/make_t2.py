#!/usr/bin/env python3
"""T2 generator: p + 15N -> 16O, reactions through an apparatus (docs/upgrade/ARCHITECTURE.md §10).

Writes, from the true parameters in TRUTH / NUISANCE below:

  T2_true.azr        the truth: every parameter fixed at its true value (the data are made from it)
  T2.azr             the fit project: the same data, the free parameters at their start values,
                     the nuisances at their nominal values with their priors
  data/*.dat         synthetic data (lab frame, the engine's input convention), noise drawn with SEED
  cov_diff.txt       covariance of the two differential segments (diagonal + one correlated vector)
  prior_cov.txt      a two-parameter prior covariance (an external normalization measurement)
  truth.json         the true values of every fitted parameter, in the units of parameters.out /
                     param.sav, for check_recovery.py

The noiseless model at every data point is computed by the engine itself (today's AZURE2) from
T2_true.azr: one run per segment (or per group of segments that share a parameter), the model read from the
AZUREOut_*.out file and taken back to the lab frame with the data's own lab->c.m. factor
(column 6 over the datum).  So T2 checks that a fit through the apparatus recovers the parameters
the data were made with (identifiability and the fit machinery), not that the engine's physics is
right; the identity checks and tests/reference do that.

    AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/T2/make_t2.py [--asimov OUTDIR]

--asimov OUTDIR writes the noiseless data (and the covariance built from them) to OUTDIR/data
and OUTDIR/cov_diff.txt instead of the committed files, for the Asimov recovery check.
numpy only.  Runs in about 30 s on one core (11 small engine runs).
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261106

# --------------------------------------------------------------------------------------------
# Particle pairs (atomic masses in u, AME2020; separation energies of 16O in MeV).
# Pair 1 carries the Trojan-horse binding energy B(d -> p + n) = 2.224566 MeV and the
# gammaIsRWA flag 0 (columns 32-33): only THM segments read them.
M_P, M_15N, M_A, M_12C = 1.007825, 15.000109, 4.002603, 12.0
SIN = 12.1274                     # S_p(16O): the entrance pair's separation energy
PAIRS = {
    1: dict(j1=0.5, pi1=1, j2=0.5, pi2=-1, e2=0.0, m1=M_P, m2=M_15N, z1=1, z2=7, S=12.1274, R=4.85,
            thm="2.224566 0"),
    2: dict(j1=0.0, pi1=1, j2=0.0, pi2=1, e2=0.0, m1=M_A, m2=M_12C, z1=2, z2=6, S=7.16192, R=5.43, thm=""),
    3: dict(j1=0.0, pi1=1, j2=2.0, pi2=1, e2=4.43891, m1=M_A, m2=M_12C, z1=2, z2=6, S=7.16192, R=5.43, thm=""),
}
LAB = (M_P + M_15N) / M_15N       # E_lab(p) / E_cm


def ep(ecm):
    return ecm * LAB


# --------------------------------------------------------------------------------------------
# The true levels (Brune parameters: observed energy E_x in MeV, partial widths in eV, the sign
# that of the reduced-width amplitude).  The 1- pair and the two 3- levels follow the 16O levels
# of tests/15N_p_a (fitted to Schardt 1952 and La Cognata 2007); the alpha_1 (12C 2+ 4.439 MeV)
# widths are invented for T2.  'free' marks the parameters the fit varies, with their start values.
#   channel: (pair, 2s, 2l, width_eV, start_or_None)
TRUTH = [
    dict(J=1, pi=-1, E=12.4208, Estart=12.4258, ch=[
        (1, 2, 0, 496.0, None), (1, 2, 4, 0.0, None),
        (2, 0, 2, 80900.0, 75000.0),
        (3, 4, 2, 50.0, None)]),
    dict(J=1, pi=-1, E=13.1050, Estart=None, ch=[
        (1, 2, 0, 36050.0, 40000.0), (1, 2, 4, 0.0, None),
        (2, 0, 2, -689600.0, None),
        (3, 4, 2, 20000.0, 15000.0)]),
    dict(J=1, pi=-1, E=17.0, Estart=None, ch=[          # background pole, fixed
        (1, 2, 0, 8305810.0, None), (1, 2, 4, 0.0, None),
        (2, 0, 2, -3369465.0, None),
        (3, 4, 2, 0.0, None)]),
    dict(J=3, pi=-1, E=13.0734, Estart=None, ch=[
        (1, 2, 4, 4648.0, None), (1, 2, 8, 0.0, None),
        (2, 0, 6, 43940.0, None),
        (3, 4, 2, 5000.0, None)]),
    dict(J=3, pi=-1, E=13.2692, Estart=13.2652, ch=[
        (1, 2, 4, 1702.0, None), (1, 2, 8, 0.0, None),
        (2, 0, 6, 12920.0, 14000.0),
        (3, 4, 2, 0.0, None)]),
]

# --------------------------------------------------------------------------------------------
# Nuisances: true value, nominal (start value and prior centre) and prior width.
NUISANCE = {
    "s1_norm": dict(true=1.04, nominal=1.0, prior="5 %"),
    "s1_shift": dict(true=0.0020, nominal=0.0, prior="0.003 MeV"),
    "D_norm": dict(true=0.95, nominal=1.0, prior="8 %"),              # normgroup D (segments 2, 3)
    "G_shift": dict(true=-0.0015, nominal=0.0, prior="0.003 MeV"),    # shiftgroup G (segments 2, 3)
    "G_fn_a0": dict(true=0.0040, nominal=0.0, prior="0.01"),          # shiftfn a0*(x-0.8) on group G
    "s4_norm": dict(true=0.98, nominal=1.0, prior="3 %"),
    "s4_sqrt": dict(true=0.0015, nominal=0.0, prior="0.003 MeV^1/2"),
    "s5_density": dict(true=2.2e17, nominal=2.0e17, prior="10 %"),    # implanted layer, atoms/cm^2
    "s6_norm": dict(true=1.02, nominal=1.0, prior="5 %"),
    "s7_norm": dict(true=0.97, nominal=1.0, prior="prior_cov.txt"),
    "s8_norm": dict(true=1.03, nominal=1.0, prior="prior_cov.txt"),
    "abundance": dict(true=0.00364, nominal=0.0037, prior="0.0002"),  # 15N in natural N (sample natN)
    "s13_norm": dict(true=1.0 / 1.35, nominal=1.0, prior="none"),    # OES points: free scale
}
# The external normalization measurement of segments 7 and 8 (prior_cov.txt): sigmas and correlation.
PRIOR_SIG = (0.03, 0.04)
PRIOR_RHO = 0.5

# --------------------------------------------------------------------------------------------
# Segments.  Each: key, exit, isDiff, angle, lab energies, relative noise, and the apparatus.
#   For THM (isDiff 10) and OES points the grid is in c.m. energy of p + 15N, converted to lab.
SEGMENTS = [
    dict(key=1, name="s01_int_thin", exit=2, isdiff=0, angle=0.0, E=np.linspace(0.20, 1.40, 40), rel=0.04,
         what="(p,a0) angle-integrated, thin TiN target: stopping, Bohr straggling, 1 keV beam spread; "
              "norm (5 % prior), constant shift (3 keV prior)"),
    dict(key=2, name="s02_diff090", exit=2, isdiff=1, angle=90.0, E=np.linspace(0.25, 1.35, 30), rel=0.04,
         what="(p,a0) dsigma/dOmega lab 90 deg, 2 keV resolution, Q attenuation; normgroup D, shiftgroup G "
              "with shiftfn a0*(x-0.8), covariance cov_diff.txt"),
    dict(key=3, name="s03_diff150", exit=2, isdiff=1, angle=150.0, E=np.linspace(0.25, 1.35, 30), rel=0.04,
         what="(p,a0) dsigma/dOmega lab 150 deg, as segment 2 (same groups and covariance)"),
    dict(key=4, name="s04_thick", exit=2, isdiff=0, angle=0.0, E=np.linspace(0.25, 1.35, 30), rel=0.03,
         what="(p,a0) thick-target yield (density 0), no beam spread; norm (3 %), sqrt(E) shift"),
    dict(key=5, name="s05_implanted", exit=2, isdiff=0, angle=0.0, E=np.linspace(1.17, 1.27, 26), rel=0.03,
         what="(p,a0) across the narrow 3- level, 15N implanted in Ta (stopping per 15N atom), fitted "
              "areal density (10 % prior), straggling 0.6 keV^1/2"),
    dict(key=6, name="s06_gas", exit=2, isdiff=0, angle=0.0, E=np.linspace(0.26, 0.42, 25), rel=0.04,
         what="(p,a0) extended 15N2 gas target, free-gas Doppler 300 K, energy-dependent resolution "
              "a0*x (lab); norm (5 %)"),
    dict(key=7, name="s07_profile", exit=2, isdiff=0, angle=0.0, E=None, rel=0.03,
         what="(p,a0) broad skewed beam profile, reconstructed-energy windows (columns 5-6); norm "
              "(prior_cov.txt)"),
    dict(key=8, name="s08_sum", exit=2, isdiff=0, angle=0.0, E=np.linspace(0.60, 1.40, 30), rel=0.04,
         what="SUM (p,a0) + (p,a1) angle-integrated, 2 keV resolution; norm (prior_cov.txt)"),
    dict(key=9, name="s09_ratio", exit=3, isdiff=0, angle=0.0, E=np.linspace(0.70, 1.40, 20), rel=0.05,
         what="RATIO sigma(p,a1)/sigma(p,a0) angle-integrated, no apparatus"),
    dict(key=10, name="s10_natural", exit=2, isdiff=0, angle=0.0, E=np.linspace(0.25, 0.45, 15), rel=0.06,
         what="(p,a0) on a natural-nitrogen target (sample natN): the 15N abundance is a fit parameter "
              "(prior 0.0002)"),
    dict(key=11, name="s11_thm_hoes_a", exit=2, isdiff=10, angle=0.0, Ecm=np.linspace(0.03, 0.60, 30),
         rel=0.07, what="THM HOES 2H(15N,a12C)n, run a: experiment[LC09], 20 keV c.m. resolution"),
    dict(key=12, name="s12_thm_hoes_b", exit=2, isdiff=10, angle=0.0, Ecm=np.linspace(0.04, 0.58, 28),
         rel=0.07, what="THM HOES run b: same experiment (shared profiled norm), same resolution"),
    dict(key=13, name="s13_thm_oes", exit=2, isdiff=0, angle=0.0, Ecm=np.linspace(0.05, 0.50, 16),
         rel=0.10, what="THM on-shell-equivalent sigma (OES), free scale, no resolution fold (TODO, README)"),
]
# Segment 7: one skewed-Gaussian beam (lab MeV) and the reconstructed-energy windows of its points.
PROFILE = dict(xi=0.34, omega=0.035, alpha=1.5, w=1.0, s=0.004)
RUNS = [(1,), (2, 3), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (11, 12), (13,)]
PROFILE_WINDOWS = [(a, a + 0.015) for a in np.arange(0.27, 0.42, 0.015)]

STOP_TIN = '"a0*x^a1" 2 1.6e-20 -0.45'      # p in TiN, MeV cm^2 per N atom (lab x in MeV)
STOP_TA = '"a0*x^a1" 2 6.0e-20 -0.30'       # p in Ta with implanted 15N, per 15N atom
STOP_GAS = '"a0*x^a1" 2 0.9e-20 -0.50'      # p in N2 gas, per N atom
GRID = "10 20"                               # resonance-width multiplier, points per width (sub-point grids)
SIGMA_THM_CM = 0.020                         # THM resolution, c.m. MeV
SIGMA_THM_LAB = SIGMA_THM_CM * LAB           # the <targetInt> line takes lab (converted back by m2/(m1+m2))

CONFIG = """<config>
true                                                                                                #Perform A-Matrix Calculation
output/                                                                                             #Full Path to Output Directory
checks/                                                                                             #Full Path to Checks Directory
none                                                                                                #Compond Nucleus Check
none                                                                                                #Boundary Condition Check
none                                                                                                #Data Check
none                                                                                                #Lo-Matrix and Penetrability Check
none                                                                                                #Legendre Polynomial Check
none                                                                                                #Coulomb Amplitudes Check
none                                                                                                #Reaction Pathway Check
none                                                                                                #Angular Distributions Check
{extra}</config>
<potential>
useHybridPotential=0
useAdaptiveGrid=1
potentialType=0
V0=150
R=3.6
a=0.6
</potential>
"""
TAIL = """<lastRun>
0
"output/param.sav"
"output/intEC.dat"
1 2
0 ""
0.07 10 -1.
</lastRun>
"""


# A level's widths are entered as observed partial widths at the level's energy, but what the fit holds fixed is the
# reduced-width amplitude gamma (MeV^1/2) of the Brune basis (param.sav).  For a level whose energy starts off the
# truth, the fixed widths that give the true gamma at the start energy differ from the true widths; match_fixed()
# finds them with the engine and keeps them here, (level, channel) -> width in eV at the start energy.
FIXED_AT_START = {}


def level_lines(fit):
    """The <levels> block.  fit=True: free parameters at their start values, the others fixed."""
    out = []
    for lid, lev in enumerate(TRUTH, start=1):
        E, fixE = lev["E"], 1
        if fit and lev["Estart"] is not None:
            E, fixE = lev["Estart"], 0
        for ci, (pair, s2, l2, width, start) in enumerate(lev["ch"], start=1):
            p = PAIRS[pair]
            g, chfix = width, 1
            if fit and start is not None:
                g, chfix = start, 0
            elif fit:
                g = FIXED_AT_START.get((lid, ci), width)
            out.append(f"{lev['J']} {lev['pi']} {E:.6f} {fixE} 1 {pair} {s2} {l2} {lid} 1 {chfix} {g:.12g} "
                       f"{p['j1']} {p['pi1']} {p['j2']} {p['pi2']} {p['e2']} {p['m1']} {p['m2']} {p['z1']} "
                       f"{p['z2']} {SIN} {p['S']} 0 0 0.0 0 {p['R']} 0 0 0" + (f" {p['thm']}" if p['thm'] else ""))
        out.append("")
    return "<levels>\n" + "\n".join(out) + "</levels>\n"


def nz(name, fit):
    """A nuisance's value in the truth project (its true value) or the fit project (its nominal)."""
    return NUISANCE[name]["nominal" if fit else "true"]


def segment_lines(fit, only=None, data_dir="data", cov=True):
    """The <segmentsData> block.  only=(k, ...): every other segment inactive (the truth runs)."""
    v = 1 if fit else 0                     # vary flags only in the fit project
    rows = []
    for s in SEGMENTS:
        k = s["key"]
        act = 1 if (only is None or k in only) else 0
        lo, hi = seg_range(s)
        th = (s["angle"], s["angle"]) if s["isdiff"] == 1 else (0, 180)
        f = f"{data_dir}/{s['name']}.dat"
        norm, vn, nerr, sh, sherr, vs, tail = 1.0, 0, 0, 0.0, 0, 0, "0 0"
        if k == 1:
            norm, vn, nerr = nz("s1_norm", fit), v, 5
            sh, sherr, vs = nz("s1_shift", fit), 0.003, v
        elif k in (2, 3):
            if k == 2:
                norm, vn, nerr = nz("D_norm", fit), v, 8
                sh, sherr, vs = nz("G_shift", fit), 0.003, v
                a0 = nz("G_fn_a0", fit)
                tail += f' normgroup D shiftgroup G shiftfn "a0*(x-0.8)" 1 {a0:g} 0.01 {v}'
            else:
                tail += " normgroup D shiftgroup G"
            if fit and cov:
                tail += " covariance cov_diff.txt"
        elif k == 4:
            norm, vn, nerr = nz("s4_norm", fit), v, 3
            tail += f" sqrtshift {nz('s4_sqrt', fit):g} 0.003 {v}"
        elif k == 6:
            norm, vn, nerr = nz("s6_norm", fit), v, 5
        elif k == 7:
            norm, vn = nz("s7_norm", fit), v
        elif k == 8:
            norm, vn = nz("s8_norm", fit), v
            tail = "1 0 1 1 3 -999"                                   # SUM: (1,2) + (1,3)
        elif k == 9:
            tail = "1 1 1 1 2 -999"                                   # RATIO: (1,3) / (1,2)
        elif k == 10:
            tail += " sample natN"
        elif k in (11, 12):
            vn = 1                                                     # THM: profiled free scale
        elif k == 13:
            norm, vn = nz("s13_norm", fit), v
        rows.append(f"{act} 1 {s['exit']} {lo:.6f} {hi:.6f} {th[0]:g} {th[1]:g} {s['isdiff']} {norm:.6g} {vn} "
                    f"{nerr} {sh:.6g} {sherr:g} {vs} {f} {tail}")
    return "<segmentsData>\n" + "\n".join(rows) + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n"


def target_lines(fit):
    dens = nz("s5_density", fit)
    dfit = f" densityfit 1 10% {3.0e17:.6g}" if fit else ""
    q = "1 5 1 0.995 0.98 0.96 0.93"
    p = PROFILE
    rows = [
        f'1 "1" 20 1 0.001 1 5e17 {STOP_TIN} 0 0 0 "" 0 1 0 {GRID}',                         # thin + spread + Bohr
        f'1 "2,3" 20 1 0.002 0 0 "" 0 {q} 0 "" 0',                                     # resolution + Q
        f'1 "4" 40 0 0 1 0 {STOP_TIN} 0 0 0 "" 0 0 0 {GRID}',                                     # thick target
        f'1 "5" 40 1 0.0008 1 {dens:.6g} {STOP_TA} 0 0 0 "" 0 1 0.6 {GRID}{dfit}',            # implanted
        f'1 "6" 40 0 0 1 1e18 {STOP_GAS} 0 0 1 "a0*x" 1 0.004 0 sigmalab doppler 300',  # gas + Doppler
        f'1 "7" 60 0 0 0 0 "" 0 0 0 0 "" 0 0 0 20 50 beamprofile 1 {p["xi"]} {p["omega"]} {p["alpha"]} '
        f'{p["w"]} {p["s"]} 0 0',                                                      # beam profile
        f'1 "8" 20 1 0.002 0 0 "" 0 0 0 0 "" 0',                                       # SUM resolution
        f'1 "11,12" 30 1 {SIGMA_THM_LAB:.6f} 0 0 "" 0 0 0 0 "" 0',                      # THM resolution
    ]
    return "<targetInt>\n" + "\n".join(rows) + "\n</targetInt>\n"


def abundance_lines(fit):
    a = nz("abundance", fit)
    return f"<abundances>\nnatN: 8 16 {a:g} {1 if fit else 0} {0.0002 if fit else 0}\nnatN: sum 0\n</abundances>\n"


THM_BLOCK = """<thm>
# La Cognata et al. 2009: 2H(15N,a12C)n at 60 MeV; two HOES runs share one profiled scale
experiment[LC09] segments=11,12 beam=15N target=d spectator=n Ebeam=60
</thm>
"""


def project(fit, only=None, data_dir="data", cov=True):
    extra = "priorCovariance=prior_cov.txt\n" if (fit and cov) else ""
    return (CONFIG.format(extra=extra) + level_lines(fit) + abundance_lines(fit) + segment_lines(fit, only, data_dir, cov)
            + target_lines(fit) + "<parameterSettings>\n</parameterSettings>\n" + TAIL + THM_BLOCK)


def seg_range(s):
    if s["key"] == 7:
        c = [0.5 * (a + b) for a, b in PROFILE_WINDOWS]
        return min(c) - 1e-6, max(c) + 1e-6
    E = grid(s)
    return E.min() - 1e-6, E.max() + 1e-6


def grid(s):
    if s["key"] == 7:
        return np.array([0.5 * (a + b) for a, b in PROFILE_WINDOWS])
    if "Ecm" in s:
        return ep(s["Ecm"])
    return s["E"]


def write_data(path, s, values, errors):
    E = grid(s)
    with open(path, "w") as fh:
        fh.write(f"# T2 synthetic data, segment {s['key']}: {s['what']}\n")
        fh.write(f"# generated by make_t2.py (seed {SEED}); lab frame: E_lab(MeV) angle(deg) value error"
                 + (" window_lo window_hi (lab MeV)" if s["key"] == 7 else "") + "\n")
        for i, e in enumerate(E):
            row = f"{e:.6f} {s['angle']:g} {values[i]:.8e} {errors[i]:.8e}"
            if s["key"] == 7:
                row += f" {PROFILE_WINDOWS[i][0]:.6f} {PROFILE_WINDOWS[i][1]:.6f}"
            fh.write(row + "\n")


def run_engine(binary, workdir, azr):
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    os.makedirs(os.path.join(workdir, "checks"), exist_ok=True)
    r = subprocess.run([binary, "--no-gui", "--no-readline", azr], input="1\n\n\n7\n", cwd=workdir,
                       capture_output=True, text=True, timeout=600)
    return r.stdout + r.stderr


def read_out(path):
    """Rows of an AZUREOut file: (E_cm, Ex, theta_cm, fit_cm, S_fit, data_cm, err_cm, ...)."""
    rows = []
    for line in open(path):
        f = line.split()
        if len(f) >= 7:
            try:
                rows.append([float(x) for x in f[:9]])
            except ValueError:
                pass
    return np.array(rows)


def model_lab(binary, tmp):
    """The noiseless model of every segment, in the data file's units (lab, divided by the true norm)."""
    os.makedirs(os.path.join(tmp, "data"), exist_ok=True)
    for s in SEGMENTS:                                  # placeholder data: value 1, error 0.1
        n = len(grid(s))
        write_data(os.path.join(tmp, "data", s["name"] + ".dat"), s, np.ones(n), 0.1 * np.ones(n))
    out = {}
    # Segments that share a parameter are computed together: 2-3 (normgroup D, shiftgroup G and its shiftfn,
    # all declared on segment 2's line) and 11-12 (one profiled THM scale).  Alone, a follower would take its
    # own line's norm and shift, and each THM run its own scale.
    for group in RUNS:
        d = os.path.join(tmp, "run" + "_".join(map(str, group)))
        os.makedirs(d)
        shutil.copytree(os.path.join(tmp, "data"), os.path.join(d, "data"))
        open(os.path.join(d, "true.azr"), "w").write(project(fit=False, only=group))
        log = run_engine(binary, d, "true.azr")
        files = [f for f in os.listdir(os.path.join(d, "output")) if re.match(r"AZUREOut_aa=\d+_R=\d+\.out$", f)]
        if len(files) != 1:
            sys.exit(f"segments {group}: expected one AZUREOut file, found {files}\n{log[-3000:]}")
        rows = read_out(os.path.join(d, "output", files[0]))
        sizes = [len(grid(SEGMENTS[k - 1])) for k in group]
        if len(rows) != sum(sizes):
            sys.exit(f"segments {group}: {len(rows)} output rows for {sum(sizes)} points\n{log[-3000:]}")
        i = 0
        for k, n in zip(group, sizes):                 # rows come segment by segment in key order
            out[k] = rows[i:i + n, 3] / rows[i:i + n, 5]   # fit_cm * (datum_lab / datum_cm), datum_lab = 1
            i += n
    return out


def param_par(binary, workdir, azr_text):
    """param.par of a calculation of azr_text: name -> value (Brune basis: E in MeV, gamma in MeV^1/2)."""
    if os.path.isdir(os.path.join(workdir, "output")):
        shutil.rmtree(os.path.join(workdir, "output"))
    open(os.path.join(workdir, "p.azr"), "w").write(azr_text)
    log = run_engine(binary, workdir, "p.azr")
    path = os.path.join(workdir, "output", "param.par")
    if not os.path.isfile(path):
        sys.exit("no param.par\n" + log[-3000:])
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 2 and f[0].startswith(("width_", "energy_")):
            out[f[0]] = float(f[1])
    return out


def match_fixed(binary, tmp):
    """Fixed widths of the levels that start off their true energy, so that their gammas are the true ones."""
    truth = param_par(binary, tmp, project(fit=False))
    for it in range(10):
        start = param_par(binary, tmp, project(fit=True, cov=False))
        worst = 0.0
        for lid, lev in enumerate(TRUTH, start=1):
            if lev["Estart"] is None:
                continue
            for ci, (pair, s2, l2, width, st) in enumerate(lev["ch"], start=1):
                name = f"width_{lid}_{ci}"
                if st is not None or truth[name] == 0.0:
                    continue
                r = truth[name] / start[name]
                worst = max(worst, abs(r - 1.0))
                FIXED_AT_START[(lid, ci)] = FIXED_AT_START.get((lid, ci), width) * r * r
        if worst < 1e-6:            # param.par carries 8 significant digits; the iteration ends near 2e-7
            return it
    sys.exit(f"fixed widths at the start energies did not converge (worst gamma ratio {worst:.3g})")


def covariance(m2, m3, E2, E3):
    """cov_diff.txt terms: 4 % diagonal and one fully correlated energy-dependent 3 % tilt."""
    m = np.concatenate([m2, m3])
    E = np.concatenate([E2, E3])
    diag = (0.04 * m) ** 2
    vec = 0.03 * m * (E - 0.8) / 0.55
    return diag, vec


def write_cov(path, n2, n3, diag, vec, data_dir="data"):
    with open(path, "w") as fh:
        fh.write("# T2: covariance of segments 2 and 3 (make_t2.py): 4 % uncorrelated + a 3 % efficiency tilt\n")
        fh.write("# (E-0.8)/0.55 common to both angles, in the data's units ((b/sr)^2, lab)\n")
        fh.write(f"file {data_dir}/s02_diff090.dat {n2}\nfile {data_dir}/s03_diff150.dat {n3}\n")
        fh.write("diagonal\n" + "\n".join(f"{x:.8e}" for x in diag) + "\n")
        fh.write("vector\n" + "\n".join(f"{x:.8e}" for x in vec) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default=os.environ.get("AZURE2_BINARY", ""), help="AZURE2 (default $AZURE2_BINARY)")
    ap.add_argument("--asimov", metavar="OUTDIR", help="write noiseless data to OUTDIR instead")
    args = ap.parse_args()
    if not args.binary or not os.access(args.binary, os.X_OK):
        sys.exit("set AZURE2_BINARY or pass --binary path/to/AZURE2")
    binary = os.path.abspath(args.binary)
    rng = np.random.default_rng(SEED)
    tmp = tempfile.mkdtemp(prefix="t2gen.")
    try:
        model = model_lab(binary, tmp)
        match_fixed(binary, os.path.join(tmp, "run1"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    asimov = args.asimov is not None
    root = os.path.abspath(args.asimov) if asimov else HERE
    ddir = os.path.join(root, "data")
    os.makedirs(ddir, exist_ok=True)

    # covariance segments first (correlated draw), then the rest in key order: one stream, one seed
    s2, s3 = SEGMENTS[1], SEGMENTS[2]
    m2, m3 = model[2], model[3]
    diag, vec = covariance(m2, m3, grid(s2), grid(s3))
    C = np.diag(diag) + np.outer(vec, vec)
    z = rng.multivariate_normal(np.zeros(len(diag)), C)
    if asimov:
        z[:] = 0.0
    d23 = np.concatenate([m2, m3]) + z
    err = np.sqrt(diag + vec ** 2)
    write_data(os.path.join(ddir, s2["name"] + ".dat"), s2, d23[:len(m2)], err[:len(m2)])
    write_data(os.path.join(ddir, s3["name"] + ".dat"), s3, d23[len(m2):], err[len(m2):])
    write_cov(os.path.join(root, "cov_diff.txt"), len(m2), len(m3), diag, vec)
    for s in SEGMENTS:
        k = s["key"]
        if k in (2, 3):
            continue
        m = model[k]
        e = s["rel"] * np.abs(m)
        noise = rng.normal(0.0, 1.0, len(m)) * e
        write_data(os.path.join(ddir, s["name"] + ".dat"), s, m + (0.0 if asimov else noise), e)

    # the external normalization measurement of segments 7 and 8 (drawn about the truth)
    cov = np.array([[PRIOR_SIG[0] ** 2, PRIOR_RHO * PRIOR_SIG[0] * PRIOR_SIG[1]],
                    [PRIOR_RHO * PRIOR_SIG[0] * PRIOR_SIG[1], PRIOR_SIG[1] ** 2]])
    truth_pm = np.array([NUISANCE["s7_norm"]["true"], NUISANCE["s8_norm"]["true"]])
    mean = truth_pm + (0.0 if asimov else rng.multivariate_normal(np.zeros(2), cov))
    with open(os.path.join(root, "prior_cov.txt"), "w") as fh:
        fh.write("# T2: external normalization of segments 7 and 8 (make_t2.py; mean drawn about the truth)\n")
        fh.write("segment_7_norm segment_8_norm\n")
        fh.write(f"{mean[0]:.6f} {mean[1]:.6f}\n")
        for r in cov:
            fh.write(" ".join(f"{x:.6e}" for x in r) + "\n")

    open(os.path.join(root, "T2.azr"), "w").write(project(fit=True))
    if not asimov:
        open(os.path.join(HERE, "T2_true.azr"), "w").write(project(fit=False))
        truth = dict(seed=SEED, levels=TRUTH, nuisance=NUISANCE, prior_cov_mean=list(mean),
                     fixed_widths_at_start={f"level {k[0]} channel {k[1]}": v for k, v in FIXED_AT_START.items()},
                     profile=PROFILE, sigma_thm_cm=SIGMA_THM_CM)
        json.dump(truth, open(os.path.join(HERE, "truth.json"), "w"), indent=1, default=float)
    print("wrote", root)


if __name__ == "__main__":
    main()
