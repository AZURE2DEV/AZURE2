"""Writes the two p+12C elastic models of this check (run once; the .azr files are committed).

with_dummies.azr : real levels plus a fixed zero-width "dummy" level for every J^pi the pair allows with
                   l <= LMAX -- the convention users followed to get every partial wave's hard-sphere
                   scattering. includeAllJPi absent (off).
no_dummies.azr   : the real levels only, includeAllJPi=1 and maxOrbitalMomentum=LMAX in <config>.
The two must give identical cross sections (check.sh).
"""
import itertools, os
LMAX = 4
# p (1/2+) + 12C (0+): separation energy 1.9435 MeV, channel radius 5.0 fm; channel spin 1/2, so l <= 4 gives
# nine J^pi (1/2+- ... 9/2+), six of them without a level below.
PAIR = dict(j1=0.5, pi1=1, j2=0.0, pi2=1, e2=0.0, m1=1.00782, m2=12.0, z1=1, z2=6, sep=1.9435, ptype=0, rad=5.0)
# 13N levels (E_x MeV, reduced width amplitude column): 1/2+ 2.365, 3/2- 3.502, 5/2+ 3.547
REAL = [(0.5, 1, 2.365, 0.60), (1.5, -1, 3.502, 0.25), (2.5, 1, 3.547, 0.20)]

def channels(j, pi):
    out = []
    s_lo, s_hi = abs(PAIR["j1"] - PAIR["j2"]), PAIR["j1"] + PAIR["j2"]
    s = s_lo
    while s <= s_hi + 1e-9:
        for l in range(LMAX + 1):
            if PAIR["pi1"] * PAIR["pi2"] * (-1) ** l != pi: continue
            # triangle rule, and J runs from |s-l| in integer steps (J = 1 does not couple to s = 1/2, l = 1)
            if abs(s - l) - 1e-9 <= j <= s + l + 1e-9 and abs((j - abs(s - l)) - round(j - abs(s - l))) < 1e-9: out.append((s, l))
        s += 1.0
    return out

def line(j, pi, ex, fix, s, l, lid, active, chfix, gamma):
    P = PAIR
    return (f"{j:g} {pi:d} {ex:g} {fix:d} 1 1 {int(2*s)} {int(2*l)} {lid} {active} {chfix} {gamma:g} "
            f"{P['j1']:g} {P['pi1']} {P['j2']:g} {P['pi2']} {P['e2']:g} {P['m1']} {P['m2']} {P['z1']} {P['z2']} "
            f"{P['sep']} {P['sep']} 0 0 0.0 {P['ptype']} {P['rad']} 0 0 0")

def levels(with_dummies):
    out, lid = [], 1
    seen = set()
    for j, pi, ex, g in REAL:
        for s, l in channels(j, pi): out.append(line(j, pi, ex, 0, s, l, lid, 1, 0, g))
        seen.add((j, pi)); lid += 1
        out.append("")
    if with_dummies:
        for j2 in range(0, 2 * (LMAX + 1) + 1):
            for pi in (1, -1):
                j = j2 / 2
                if (j, pi) in seen or not channels(j, pi): continue
                for s, l in channels(j, pi): out.append(line(j, pi, 40.0, 1, s, l, lid, 1, 1, 0.0))
                lid += 1
                out.append("")
    return "\n".join(out)

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
TAIL = """<segmentsData>
</segmentsData>
<segmentsTest>
1 1 1 1.0 3.0 0.5 60 160 20 1
1 1 1 1.0 3.0 0.25 0 0 0 0 0
</segmentsTest>
<targetInt>
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
for name, dummies, extra in (("with_dummies", True, ""), ("no_dummies", False, f"includeAllJPi=1\nmaxOrbitalMomentum={LMAX}\n")):
    with open(os.path.join("models", name + ".azr"), "w") as f:
        f.write(CONFIG.format(extra=extra) + "<levels>\n" + levels(dummies) + "\n</levels>\n" + TAIL)
    print(name, "written")
