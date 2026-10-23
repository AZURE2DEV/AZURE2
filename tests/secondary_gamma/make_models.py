#!/usr/bin/env python3
"""Build the secondary-gamma (UPOS) test models and the expected angular distributions.

13C + alpha in 17O.  Exit pairs: n + 16O(6.13, 3-) (key 5; channel spins 5/2 and 7/2, E3 to the 0+ ground state)
and alpha + 13C(3.68, 3/2-) (key 6; channel spin 3/2, M1/E2 to the 1/2- ground state).  Single-level models are
evaluated exactly at the level energy, where the exit T-matrix elements of one l' are proportional to
sign(Gamma) sqrt|Gamma| with the same energy factor and phase, so the distribution follows from the widths alone
(truth.py).  Every model also has an angle-integrated segment of the same reaction, which normalises the UPOS cross
section: W(theta) = 4 pi (dsigma/dOmega_gamma) / sigma.

  spin2_plus / spin2_minus / spin2_zero : 5/2+, l' = 1 in s' = 5/2 (Gamma 20 keV) and 7/2 (+10, -10, 0 keV)
  spin2_multil                          : as plus, with l' = 3 in both channel spins as well (-5 and +8 keV)
  mixing                                : 5/2+ -> alpha + 13C(3/2-), l' = 1 and 3; M1/E2, delta = -0.094
  twolevel                              : 5/2+ and 3/2+ -> alpha + 13C(3/2-) (one channel spin), J1 != J2 interference;
                                          expected = the pre-fix engine (dev 6c6932d), unchanged by the fix
  grad                                  : 5/2+ and 3/2+ -> n + 16O(3-), free energies and widths, UPOS data at
                                          three energies -- for the AZURE_GRAD_CHECK gradient test

Run from this directory: python3 make_models.py  (needs numpy, scipy and sympy; writes models/ and expected/).
"""
import os
import numpy as np
import truth

HERE = os.path.dirname(os.path.abspath(__file__))
MA, M13 = 4.0026, 13.0033
S_ALPHA = 6.35869
ANGLES = np.arange(0.0, 181.0, 15.0)

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
ALPHA = "0  1  0.5  -1  0  4.0026  13.0033  2  6  4.14313  6.35869  0  0  0.0  0  5.4  0  0  0"
N_3M = "0.5  1  3  -1  6.12989  1.00866  15.9949  0  8  4.14313  4.14313  0  0  0.0  0  4.4  0  0  0"
A_32 = "0  1  1.5  -1  3.6845  4.0026  13.0033  2  6  4.14313  6.35869  0  0  0.0  0  5.4  0  0  0"


def channel(J, E, key, s, l, width, tail, free=False, lid=1):
    """One <levels> line: J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    fix = 0 if free else 1
    return f"{J}  1  {E}  {fix}  1  {key}  {int(round(2 * s))}  {2 * l}  {lid}  1  {fix}  {width}  {tail}"


def elab(ex):
    return (ex - S_ALPHA) * (MA + M13) / M13


def write(name, levels, segments, data):
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    blocks = []
    for lv in levels:
        blocks.append("\n".join(lv))
    txt = (CONFIG + "<levels>\n" + "\n\n".join(blocks) + "\n\n</levels>\n<segmentsData>\n" + "\n".join(segments)
           + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)
    for fname, rows in data.items():
        with open(os.path.join(HERE, "models", "data", fname), "w") as f:
            for r in rows:
                f.write("  ".join(f"{v:.10g}" for v in r) + "\n")


def single_level(name, exit_key, exit_tail, widths, upos, E=11.0):
    """widths: {(s', l'): Gamma_eV}; upos: (L, Ic, delta)."""
    lv = [channel(2.5, E, 2, 0.5, 3, 1000.0, ALPHA)]
    lv += [channel(2.5, E, exit_key, sp, lp, g, exit_tail) for (sp, lp), g in widths.items()]
    e = elab(E)
    L, Ic, delta = upos
    segs = [f"1  2  {exit_key}  {e - 0.5:.6f}  {e + 0.5:.6f}  0  180  4  1  0  0  0  0.02  1  data/{name}_upos.dat  0  1  {L}  {Ic}  {delta}",
            f"1  2  {exit_key}  {e - 0.5:.6f}  {e + 0.5:.6f}  0  0  0  1  0  0  0  0.02  1  data/{name}_int.dat  0  0",
            f"1  2  {exit_key}  {e - 0.5:.6f}  {e + 0.5:.6f}  0  180  4  1  0  0  0  0.02  1  data/{name}_part.dat  0  0"]
    data = {f"{name}_upos.dat": [(e, a, 1.0, 0.1) for a in ANGLES],
            f"{name}_int.dat": [(e, 0.0, 1.0, 0.1)],
            f"{name}_part.dat": [(e, a, 1.0, 0.1) for a in ANGLES]}
    write(name, [lv], segs, data)


def expected(name, T, Ib, IB, IC, upos):
    L, Ic, delta = upos
    th = np.radians(ANGLES)
    P = truth.populations(T, 0, 0.5, Ib, IB)
    w = truth.gamma_distribution_pure(P, IB, IC, L, th) if delta == 0 else \
        truth.gamma_distribution_mixed(P, IB, IC, L, delta, th)
    os.makedirs(os.path.join(HERE, "expected"), exist_ok=True)
    with open(os.path.join(HERE, "expected", name + ".dat"), "w") as f:
        f.write("# theta_gamma(deg)  W = 4 pi dsigma/dOmega_gamma / sigma  (truth.py)\n")
        for a, v in zip(ANGLES, w):
            f.write(f"{a:6.1f}  {v:.12f}\n")


def amp(g):
    return np.sign(g) * np.sqrt(abs(g))


def main():
    e3 = (3, 0, 0)
    for name, g72 in (("spin2_plus", 1.0e4), ("spin2_minus", -1.0e4), ("spin2_zero", 0.0)):
        w = {(2.5, 1): 2.0e4, (3.5, 1): g72}
        single_level(name, 5, N_3M, w, e3)
        expected(name, {(2.5, sp, lp, 0.5, 3): amp(g) for (sp, lp), g in w.items()}, 0.5, 3, 0, e3)
    w = {(2.5, 1): 2.0e4, (3.5, 1): 1.0e4, (2.5, 3): -5.0e3, (3.5, 3): 8.0e3}
    single_level("spin2_multil", 5, N_3M, w, e3)
    expected("spin2_multil", {(2.5, sp, lp, 0.5, 3): amp(g) for (sp, lp), g in w.items()}, 0.5, 3, 0, e3)
    mix = (1, 0.5, -0.094)
    w = {(1.5, 1): 1.0e4, (1.5, 3): -3.0e3}
    single_level("mixing", 6, A_32, w, mix)
    expected("mixing", {(2.5, sp, lp, 0.5, 3): amp(g) for (sp, lp), g in w.items()}, 0, 1.5, 0.5, mix)

    # two levels, one exit channel spin: J1 != J2 interference (expected = pre-fix engine)
    lv1 = [channel(2.5, 11.0, 2, 0.5, 3, 3000.0, ALPHA), channel(2.5, 11.0, 6, 1.5, 1, 2.0e4, A_32),
           channel(2.5, 11.0, 6, 1.5, 3, -4.0e3, A_32)]
    lv2 = [channel(1.5, 11.1, 2, 0.5, 1, 5000.0, ALPHA, lid=2), channel(1.5, 11.1, 6, 1.5, 1, -1.5e4, A_32, lid=2),
           channel(1.5, 11.1, 6, 1.5, 3, 6.0e3, A_32, lid=2)]
    energies = [elab(e) for e in (10.95, 11.0, 11.05, 11.1, 11.15)]
    segs = [f"1  2  6  {energies[0] - 0.1:.6f}  {energies[-1] + 0.1:.6f}  0  180  4  1  0  0  0  0.02  1  data/twolevel_upos.dat  0  1  1  0.5  -0.094",
            f"1  2  6  {energies[0] - 0.1:.6f}  {energies[-1] + 0.1:.6f}  0  180  4  1  0  0  0  0.02  1  data/twolevel_upos2.dat  0  1  1  0.5  0"]
    rows = [(e, a, 1.0, 0.1) for e in energies for a in ANGLES]
    write("twolevel", [lv1, lv2], segs, {"twolevel_upos.dat": rows, "twolevel_upos2.dat": rows})

    # gradient test: two free levels into both channel spins of n + 16O(3-), UPOS data off the model
    lv1 = [channel(2.5, 11.0, 2, 0.5, 3, 3000.0, ALPHA, True), channel(2.5, 11.0, 5, 2.5, 1, 2.0e4, N_3M, True),
           channel(2.5, 11.0, 5, 3.5, 1, -1.0e4, N_3M, True), channel(2.5, 11.0, 5, 2.5, 3, 4.0e3, N_3M, True),
           channel(2.5, 11.0, 5, 3.5, 3, 3.0e3, N_3M, True)]
    lv2 = [channel(1.5, 11.1, 2, 0.5, 1, 5000.0, ALPHA, True, 2), channel(1.5, 11.1, 5, 2.5, 1, -1.5e4, N_3M, True, 2),
           channel(1.5, 11.1, 5, 2.5, 3, 6.0e3, N_3M, True, 2), channel(1.5, 11.1, 5, 3.5, 3, -5.0e3, N_3M, True, 2)]
    energies = [elab(e) for e in (10.97, 11.03, 11.08, 11.13)]
    segs = [f"1  2  5  {energies[0] - 0.1:.6f}  {energies[-1] + 0.1:.6f}  0  180  4  1  0  0  0  0.02  1  data/grad_upos.dat  0  1  3  0  0"]
    rng = np.random.default_rng(7)
    rows = [(e, a, 0.002 * (1 + 0.3 * rng.standard_normal()), 0.0002) for e in energies for a in ANGLES[1:-1:2]]
    write("grad", [lv1, lv2], segs, {"grad_upos.dat": rows})


if __name__ == "__main__":
    main()
