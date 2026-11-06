"""Models of tests/relativistic (run once; the .azr files are committed).

na_hs   n + 4He, s wave only (includeAllJPi=1, maxOrbitalMomentum=0) and one far 1/2+ level of zero width:
        pure hard-sphere scattering, sigma = 4 pi sin^2(k a) / k^2, a test of k and pi/k^2.
pa_hs   p + 4He, s wave only, the same: differential elastic at 90 and 150 deg, a test of eta and k through
        the Coulomb functions and phases.
he5     5He: d + t and n + 4He (Q = 17.589 MeV) with the 3/2+ level near the d+t threshold; t(d,n)a and its
        inverse 4He(n,d)t on test grids of equal invariant mass, for detailed balance.
Each is written twice: <name>.azr (non-relativistic) and <name>_rel.azr (relativistic=1).
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
U = 931.4940880
M = dict(n=1.00866491595, p=1.00782503207, d=2.01410177812, t=3.0160492779, a=4.00260325413,
         li7=7.016003434, be7=7.016928717)

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
TAIL = """<targetInt>
</targetInt>
<lastRun>
9997
""
""
0 0
0 ""
-1. -1. -1.
</lastRun>
"""

def line(j, pi, ex, fix, s, l, lid, gamma, pair, chfix=0):
    """One <levels> line (NucLine column order)."""
    (j1, pi1, j2, pi2, m1, m2, z1, z2, sep, rad) = pair
    return (f"{j:g} {pi:d} {ex:.6f} {fix:d} 1 {PAIRKEY[pair]} {int(2*s)} {int(2*l)} {lid} 1 {chfix} {gamma:g} "
            f"{j1:g} {pi1} {j2:g} {pi2} 0 {m1} {m2} {z1} {z2} {sep} {sep} 0 0 0.0 0 {rad} 0 0 0")

NA = (0.5, 1, 0.0, 1, M['n'], M['a'], 0, 2, -0.735, 3.0)       # n + 4He: 5He S_n = -0.735 MeV
PA = (0.5, 1, 0.0, 1, M['p'], M['a'], 1, 2, -1.965, 3.0)       # p + 4He: 5Li S_p = -1.965 MeV
DT = (1.0, 1, 0.5, 1, M['d'], M['t'], 1, 1, 16.854, 5.0)       # d + t above n + 4He by Q
PAIRKEY = {}

def write(name, levels, tests, extra=""):
    for rel in (False, True):
        cfg = CONFIG.format(extra=extra + ("relativistic=1\n" if rel else ""))
        txt = cfg + "<levels>\n" + "\n".join(levels) + "\n\n</levels>\n<segmentsData>\n</segmentsData>\n<segmentsTest>\n" + \
            "\n".join(tests) + "\n</segmentsTest>\n" + TAIL
        with open(os.path.join(HERE, "models", name + ("_rel" if rel else "") + ".azr"), "w") as fh:
            fh.write(txt)

# pure hard sphere, s wave: a far 1/2+ level with zero width defines the pair
PAIRKEY.clear(); PAIRKEY[NA] = 1
write("na_hs", [line(0.5, 1, 60.0, 1, 0.5, 0, 1, 0.0, NA, 1)],
      ["1 1 1 0.5 30.0 0.5 0 0 0 0 0"], "includeAllJPi=1\nmaxOrbitalMomentum=0\n")
PAIRKEY.clear(); PAIRKEY[PA] = 1
write("pa_hs", [line(0.5, 1, 60.0, 1, 0.5, 0, 1, 0.0, PA, 1)],
      ["1 1 1 1.0 30.0 0.5 90 150 60 1"], "includeAllJPi=1\nmaxOrbitalMomentum=0\n")
# 5He: pair 1 = n + 4He (l = 2, s = 1/2), pair 2 = d + t (l = 0, s = 3/2); the 3/2+ level at E_x 16.84 MeV.
# Data segments (mode 1 evaluates them): t(d,n)a at E_d and 4He(n,d)t at the E_n of the SAME excitation energy,
# computed here non-relativistically (E_x additive) for he5.azr and from the invariant mass for he5_rel.azr.
def cm_nr(m1, m2, e): return e * m2 / (m1 + m2)
def lab_nr(m1, m2, t): return t * (m1 + m2) / m2
def cm_rel(m1, m2, e):
    a, b = m1 * U, m2 * U; m = a + b
    return 2 * b * e / ((m * m + 2 * b * e) ** 0.5 + m)
def lab_rel(m1, m2, t):
    a, b = m1 * U, m2 * U
    return t * (2 * (a + b) + t) / (2 * b)
ED = [0.02 + 0.01 * i for i in range(200)]
os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
with open(os.path.join(HERE, "models", "data", "dt_fwd.dat"), "w") as fh:
    for e in ED: fh.write(f"{e:.10f} 0.0 1.0 0.1\n")
for tag, cm, lab in (("nr", cm_nr, lab_nr), ("rel", cm_rel, lab_rel)):
    with open(os.path.join(HERE, "models", "data", f"na_rev_{tag}.dat"), "w") as fh:
        for e in ED:
            ex = cm(M['d'], M['t'], e) + DT[8]               # the excitation energy of the forward point
            fh.write(f"{lab(M['n'], M['a'], ex - NA[8]):.12f} 0.0 1.0 0.1\n")
PAIRKEY.clear(); PAIRKEY[NA] = 1; PAIRKEY[DT] = 2
LEV = [line(1.5, 1, 16.84, 0, 0.5, 2, 1, 0.35, NA), line(1.5, 1, 16.84, 0, 1.5, 0, 1, 1.20, DT)]
for rel in (False, True):
    tag = "rel" if rel else "nr"
    cfg = CONFIG.format(extra="relativistic=1\n" if rel else "")
    segs = ["1 2 1 0.0 5.0 0 180 0 1 0 0 0 0 0 data/dt_fwd.dat 0 0",
            f"1 1 2 0.0 40.0 0 180 0 1 0 0 0 0 0 data/na_rev_{tag}.dat 0 0"]
    txt = cfg + "<levels>\n" + "\n".join(LEV) + "\n\n</levels>\n<segmentsData>\n" + "\n".join(segs) + \
        "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n" + TAIL
    with open(os.path.join(HERE, "models", "he5" + ("_rel" if rel else "") + ".azr"), "w") as fh:
        fh.write(txt)
# lab-frame differential data (value 1): the engine's c.m. data value is dOmega_lab/dOmega_cm at its c.m. angle
with open(os.path.join(HERE, "models", "data", "dt_ang.dat"), "w") as fh:
    for e in (0.1, 1.0, 5.0):
        for a in (0.0, 20.0, 60.0, 90.0, 120.0, 160.0, 180.0): fh.write(f"{e:.4f} {a:.1f} 1.0 0.1\n")
with open(os.path.join(HERE, "models", "data", "pa_ang.dat"), "w") as fh:
    for e in (5.0, 20.0):
        for a in (20.0, 60.0, 90.0, 150.0, 170.0): fh.write(f"{e:.4f} {a:.1f} 1.0 0.1\n")
for rel in (False, True):
    cfg = CONFIG.format(extra="relativistic=1\n" if rel else "")
    PAIRKEY.clear(); PAIRKEY[NA] = 1; PAIRKEY[DT] = 2
    segs = ["1 2 1 0.0 10.0 0 180 1 1 0 0 0 0 0 data/dt_ang.dat 0 0"]
    txt = cfg + "<levels>\n" + "\n".join(LEV) + "\n\n</levels>\n<segmentsData>\n" + "\n".join(segs) + \
        "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n" + TAIL
    with open(os.path.join(HERE, "models", "he5_ang" + ("_rel" if rel else "") + ".azr"), "w") as fh:
        fh.write(txt)
    PAIRKEY.clear(); PAIRKEY[PA] = 1
    txt = cfg.replace("</config>", "includeAllJPi=1\nmaxOrbitalMomentum=0\n</config>") + "<levels>\n" + \
        line(0.5, 1, 60.0, 1, 0.5, 0, 1, 0.0, PA, 1) + "\n\n</levels>\n<segmentsData>\n" + \
        "1 1 1 0.0 30.0 0 180 1 1 0 0 0 0 0 data/pa_ang.dat 0 0" + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n" + TAIL
    with open(os.path.join(HERE, "models", "pa_ang" + ("_rel" if rel else "") + ".azr"), "w") as fh:
        fh.write(txt)
# 7Li(p,n)7Be just above its threshold (Q = -1.644 MeV from the masses): two c.m. angles reach each lab angle
# (g > 1) and the lab angle has a maximum, which relativistically is 0.3-0.5 deg below the Marion-Young one.
QLI = (M['p'] + M['li7'] - M['n'] - M['be7']) * U
PLI = (0.5, 1, 1.5, -1, M['p'], M['li7'], 1, 3, 17.2551, 4.0)          # 8Be S_p
NBE = (0.5, 1, 1.5, -1, M['n'], M['be7'], 0, 4, 17.2551 - QLI, 4.0)   # S_n = S_p - Q
with open(os.path.join(HERE, "models", "data", "li7pn_ang.dat"), "w") as fh:
    for e, angles in ((1.89, (0.0, 10.0, 20.0, 29.0)), (1.90, (0.0, 15.0, 30.0, 44.0)), (1.95, (5.0, 40.0, 80.0)),
                      (3.0, (30.0, 90.0, 150.0))):
        for a in angles: fh.write(f"{e:.4f} {a:.1f} 1.0 0.1\n")
with open(os.path.join(HERE, "models", "data", "li7pn_edge.dat"), "w") as fh:
    fh.write("1.9000 44.9 1.0 0.1\n")   # between the relativistic (44.60) and non-relativistic maxima (atomic masses)
for rel in (False, True):
    cfg = CONFIG.format(extra="relativistic=1\n" if rel else "")
    PAIRKEY.clear(); PAIRKEY[PLI] = 1; PAIRKEY[NBE] = 2
    lev = [line(2, -1, 19.5, 1, 2.0, 0, 1, 0.5, PLI), line(2, -1, 19.5, 1, 2.0, 0, 1, 0.5, NBE)]
    for name, dat in (("li7pn_ang", "li7pn_ang.dat"), ("li7pn_edge", "li7pn_edge.dat")):
        txt = cfg + "<levels>\n" + "\n".join(lev) + "\n\n</levels>\n<segmentsData>\n" + \
            f"1 1 2 0.0 10.0 0 180 1 1 0 0 0 0 0 data/{dat} 0 0" + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n" + TAIL
        with open(os.path.join(HERE, "models", name + ("_rel" if rel else "") + ".azr"), "w") as fh:
            fh.write(txt)
print("models written")
