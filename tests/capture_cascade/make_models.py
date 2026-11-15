#!/usr/bin/env python3
"""Build the capture-cascade test models and the expected secondary-gamma angular distributions.

A capture cascade is a UPOS segment on a CAPTURE pair: A(a,gamma1)B, B -> C + gamma2, gamma1 unobserved, gamma2
observed.  The expected distributions W(theta) = 4 pi (dsigma/dOmega) / sigma come from magnetic-substate sums that use
neither Eq. (22) of Brune and deBoer (2020) nor the R_k (../secondary_gamma/truth.py): an unobserved photon summed over
direction and helicity is an incoherent sum over the multipole states (L, mu, electric/magnetic).

Resonant capture, 13C + alpha -> 17O* (channel spin 1/2), evaluated at the level energy, where T ~ g_in g_out / Gamma_in
(photon channels carry no width in the level matrix):
  single_e1   : 7/2+ (l = 3) -E1-> B(5/2-) -E2-> C(1/2-)                       (orders 0, 2, 4)
  multipoles  : 5/2+ (l = 3) -M1 and E2-> B(3/2+) -M1/E2, delta = 0.35-> C(1/2+)  (the two primary multipoles add
                incoherently)
  twoj_plus / twoj_minus : 7/2+ and 5/2+ at one energy, both E1 -> B(5/2-): J1 != J2 interference, both signs.  The
                PRIMARY distribution of the same amplitudes is checked too (expected/<model>_prim.dat): it fixes the
                relative sign convention of the two pathways with code that existed before.
External capture, p + (1/2-) target -> B by E1, no resonance (zero-width levels carry the J^pi groups):
  ec_s        : B = 1- (s = 1, l_f = 0), p waves through J = 0+, 1+, 2+.  Every J alone leaves B oriented; their
                coherent sum must not (an s state with a spectator channel spin): W = 1.
  ec_p        : B = 2+ (s = 1, l_f = 1), d waves only (the model has no l = 0 entrance channel), J = 1-, 2-, 3-:
                the single-particle d -> p substate sum.
Resonance and external capture together:
  ecres_both  : ec_p plus a 2- resonance (l = 2) with an E1 width to B, at three energies across it.  In the substate
                basis the amplitude is c a_ec + z a_res; |c|^2 is the cross section without the resonance (ecres_hs),
                and the two unknowns Re(z c*), |z|^2 follow from the angle-integrated cross section and the PRIMARY
                distribution of ecres_both itself.  The secondary distribution is then a prediction.  This needs engine
                output, so expected/ecres_both.dat is written only when a binary is given (see below); the primary fit
                residual is printed and must be ~1e-11.
  grad        : ecres_both with free level energy, widths and ANC and cascade data off the model (AZURE_GRAD_CHECK).
  bad_linewindow : a cascade segment with "linewindow" (refused).

  python3 make_models.py [path/to/AZURE2]      (numpy, scipy, sympy; writes models/ and expected/)
"""
import glob
import os
import subprocess
import sys
import tempfile
import shutil
import numpy as np
from sympy import S

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "secondary_gamma"))
import truth  # noqa: E402

MA, M13 = 4.0026, 13.0033
S_ALPHA = 6.35869
ANGLES = np.arange(0.0, 181.0, 15.0)
TH = np.radians(ANGLES)
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
P13 = "0.5  1  0.5  -1  0  1.00728  13.00335  1  6  7.55056  7.55056  0  0  0.0  0  4.5  0  0  0"
ERES = 7.55056 + 0.93
ECRES_E = (0.96, 1.0, 1.04)


def gpair17(JB, piB, ExB):
    """Capture pair 17O(B) + gamma, no external capture."""
    return f"1  1  {JB}  {piB}  {ExB}  0  17.0  0  8  4.14313  0  0  0  0.0  10  0  0  0  0"


def gpair14(JB, piB, ExB, ec):
    """Capture pair 14N(B) + gamma; ec = external-capture multipole mask (1 = E1)."""
    return f"1  1  {JB}  {piB}  {ExB}  0  14.00307  0  7  7.55056  0  0  0  0.0  10  0  0  0  {ec}"


def channel(J, pi, E, key, s, l, width, tail, lid=1, free=False):
    """One <levels> line: J pi E fixE aa key 2s 2l levelID active fixWidth width <pair columns>."""
    fix = 0 if free else 1
    return f"{J}  {pi}  {E}  {fix}  1  {key}  {int(round(2 * s))}  {2 * l}  {lid}  1  {fix}  {width}  {tail}"


def elab(ex):
    return (ex - S_ALPHA) * (MA + M13) / M13


def write(name, levels, segments, data):
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    txt = (CONFIG + "<levels>\n" + "\n\n".join("\n".join(lv) for lv in levels) + "\n\n</levels>\n<segmentsData>\n"
           + "\n".join(segments) + "\n</segmentsData>\n<segmentsTest>\n</segmentsTest>\n<targetInt>\n</targetInt>\n")
    with open(os.path.join(HERE, "models", name + ".azr"), "w") as f:
        f.write(txt)
    for fname, rows in data.items():
        with open(os.path.join(HERE, "models", "data", fname), "w") as f:
            for r in rows:
                f.write("  ".join(f"{v:.10g}" for v in r) + "\n")


def segments(name, aa, key, energies, upos, lo, hi, tail=""):
    """Cascade (UPOS) distribution, angle-integrated, primary distribution, angle-integrated UPOS."""
    L, Ic, delta = upos
    u = f"0  1  {L}  {Ic}  {delta}"
    segs = [f"1  {aa}  {key}  {lo:.6f}  {hi:.6f}  0  180  4  1  0  0  0  0.02  1  data/{name}_upos.dat  {u}{tail}",
            f"1  {aa}  {key}  {lo:.6f}  {hi:.6f}  0  0  0  1  0  0  0  0.02  1  data/{name}_int.dat  0  0",
            f"1  {aa}  {key}  {lo:.6f}  {hi:.6f}  0  180  4  1  0  0  0  0.02  1  data/{name}_prim.dat  0  0",
            f"1  {aa}  {key}  {lo:.6f}  {hi:.6f}  0  0  0  1  0  0  0  0.02  1  data/{name}_iupos.dat  {u}"]
    ang = [(e, a, 1.0, 0.1) for e in energies for a in ANGLES]
    tot = [(e, 0.0, 1.0, 0.1) for e in energies]
    return segs, {f"{name}_upos.dat": ang, f"{name}_int.dat": tot, f"{name}_prim.dat": ang, f"{name}_iupos.dat": tot}


def save(name, w):
    os.makedirs(os.path.join(HERE, "expected"), exist_ok=True)
    with open(os.path.join(HERE, "expected", name + ".dat"), "w") as f:
        f.write("# theta(deg)  W = 4 pi dsigma/dOmega / sigma  (make_models.py)\n")
        for row in np.atleast_2d(w):
            for a, v in zip(ANGLES, row):
                f.write(f"{a:6.1f}  {v:.12f}\n")


def primary(T, IB):
    """Observed PRIMARY photon when every pathway has the same multipole: substate sum with |X_L mu|^2."""
    w = np.zeros_like(TH)
    tot = 0.0
    for nu in truth.projections(S(1) / 2):
        for mB in truth.projections(IB):
            amp = {}
            for (J, sp, L, s, l), t in T.items():
                mu = nu - mB
                if abs(mu) > L or abs(nu) > S(J):
                    continue
                amp[(L, mu)] = amp.get((L, mu), 0) + (np.sqrt(2 * l + 1) * truth.cg(l, s, J, 0, nu, nu)
                                                      * truth.cg(L, sp, J, mu, mB, nu) * t)
            for (L, mu), a in amp.items():
                w += abs(a) ** 2 * truth._x2(L, int(mu), TH)
                tot += abs(a) ** 2
    return w * 4 * np.pi / tot


def resonant_expected(name, groups, IB, upos, with_primary):
    """groups: list of {(J, J_B, L, s, l): amplitude}, one per primary multipole (summed incoherently)."""
    L2, Ic, delta = upos
    P = {m: 0.0 for m in truth.projections(IB)}
    for g in groups:
        for m, p in truth.populations(g, 0, 0.5, 0, IB).items():   # the photon as a spin-0 "particle" with l' = L
            P[m] += p
    save(name, truth.gamma_distribution_pure(P, IB, Ic, L2, TH) if delta == 0
         else truth.gamma_distribution_mixed(P, IB, Ic, L2, delta, TH))
    if with_primary:
        save(name + "_prim", primary(groups[0], IB))


def amp(g):
    return np.sign(g) * np.sqrt(abs(g))


def sp_patterns(li, lf, s, JB, J, L):
    """Substate patterns (nu, M_B, mu): single-particle E1 l_i -> l_f, and the pathway (s, l_i, J, L)."""
    keys, aec, ares = [], [], []
    for nu in truth.projections(s):
        for mu in (-1, 0, 1):
            MB = nu - mu
            if abs(MB) > JB or abs(mu) > lf:
                continue
            keys.append((nu, MB, mu))
            aec.append(truth.cg(lf, 1, li, -mu, mu, 0) * truth.cg(lf, s, JB, -mu, nu, MB))
            ares.append(np.sqrt(2 * li + 1) * truth.cg(li, s, J, 0, nu, nu) * truth.cg(L, JB, J, mu, MB, nu)
                        if J is not None and abs(nu) <= J else 0.0)
    aec, ares = np.array(aec, float), np.array(ares, float)
    return keys, aec / np.linalg.norm(aec), ares / (np.linalg.norm(ares) or 1.0)


def weights_to_w(keys, w2, JB, upos):
    L2, Ic, _ = upos
    P = {m: 0.0 for m in truth.projections(JB)}
    wp = np.zeros_like(TH)
    for (nu, MB, mu), w in zip(keys, w2):
        P[S(MB)] += w
        wp += w * truth._x2(1, int(mu), TH)
    return truth.gamma_distribution_pure(P, JB, Ic, L2, TH), wp * 4 * np.pi / np.sum(w2)


def ec_model(name, JB, piB, lf, li, Js, upos, energies, ec=1, res=None, free=False, cascade_tail=""):
    """p + (1/2-) target, B with one final channel (s = 1, l_f), zero-width levels for the J^pi of entrance l_i."""
    B = gpair14(JB, piB, 5.0, ec)
    pi_c = -piB
    lv = [[channel(float(JB), piB, 5.0, 1, 1, lf, 2.0, P13, 1, free)]]
    if res:
        gp, gg = res
        lv.append([channel(2.0, pi_c, ERES, 1, 1, li, gp, P13, 5, free), channel(2.0, pi_c, ERES, 2, JB, 1, gg, B, 5, free)])
    for i, J in enumerate(Js):
        lv.append([channel(float(J), pi_c, 20.0, 1, 1, li, 0.0, P13, i + 2), channel(float(J), pi_c, 20.0, 2, JB, 1, 0.0, B, i + 2)])
    segs, data = segments(name, 1, 2, energies, upos, 0.5, 1.5, cascade_tail)
    return lv, segs, data


def run_engine(binary, name):
    work = tempfile.mkdtemp(prefix="capture_cascade.")
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copy(os.path.join(HERE, "models", name + ".azr"), os.path.join(work, "t.azr"))
    shutil.copytree(os.path.join(HERE, "models", "data"), os.path.join(work, "data"))
    subprocess.run([binary, "--no-gui", "--no-readline", "t.azr"], input="1\n\n\n7\n", text=True, cwd=work,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
    d = np.loadtxt(glob.glob(os.path.join(work, "output", "AZUREOut_*.out"))[0], ndmin=2)
    shutil.rmtree(work)
    return d[d[:, 0] > 0]


def main(binary=None):
    E = 11.0
    e = elab(E)
    # ---- resonant capture ----
    B = gpair17(2.5, -1, 3.0)
    up = (2, 0.5, 0)
    lv = [channel(3.5, 1, E, 2, 0.5, 3, 1000.0, ALPHA), channel(3.5, 1, E, 5, 2.5, 1, 2.0, B)]
    write("single_e1", [lv], *segments("single_e1", 2, 5, [e], up, e - 0.5, e + 0.5))
    resonant_expected("single_e1", [{(3.5, 2.5, 1, 0.5, 3): 1.0}], 2.5, up, True)

    B = gpair17(1.5, 1, 3.0)
    up = (1, 0.5, 0.35)
    gm1, ge2 = 1.5, -0.7
    lv = [channel(2.5, 1, E, 2, 0.5, 3, 1000.0, ALPHA), channel(2.5, 1, E, 5, 1.5, 1, gm1, B),
          channel(2.5, 1, E, 5, 1.5, 2, ge2, B)]
    write("multipoles", [lv], *segments("multipoles", 2, 5, [e], up, e - 0.5, e + 0.5))
    resonant_expected("multipoles", [{(2.5, 1.5, 1, 0.5, 3): amp(gm1)}, {(2.5, 1.5, 2, 0.5, 3): amp(ge2)}], 1.5, up, False)

    B = gpair17(2.5, -1, 3.0)
    up = (2, 0.5, 0)
    for name, g2 in (("twoj_plus", 3.0), ("twoj_minus", -3.0)):
        ga1, gg1, ga2 = 1000.0, 2.0, 600.0
        lv1 = [channel(3.5, 1, E, 2, 0.5, 3, ga1, ALPHA), channel(3.5, 1, E, 5, 2.5, 1, gg1, B)]
        lv2 = [channel(2.5, 1, E, 2, 0.5, 3, ga2, ALPHA, 2), channel(2.5, 1, E, 5, 2.5, 1, g2, B, 2)]
        write(name, [lv1, lv2], *segments(name, 2, 5, [e], up, e - 0.5, e + 0.5))
        t1, t2 = amp(ga1) * amp(gg1) / ga1, amp(ga2) * amp(g2) / ga2
        resonant_expected(name, [{(3.5, 2.5, 1, 0.5, 3): t1, (2.5, 2.5, 1, 0.5, 3): t2}], 2.5, up, True)

    # ---- external capture ----
    up = (1, 0, 0)
    write("ec_s", *ec_model("ec_s", 1, -1, 0, 1, (0, 1, 2), up, [1.0]))
    keys, aec, _ = sp_patterns(1, 0, 1, 1, None, 1)
    ws, wp = weights_to_w(keys, aec ** 2, 1, up)
    save("ec_s", ws)       # = 1 at every angle
    save("ec_s_prim", wp)
    up = (2, 0, 0)
    write("ec_p", *ec_model("ec_p", 2, 1, 1, 2, (1, 2, 3), up, [1.0]))
    keys, aec, ares = sp_patterns(2, 1, 1, 2, 2, 1)
    ws, wp = weights_to_w(keys, aec ** 2, 2, up)
    save("ec_p", ws)
    save("ec_p_prim", wp)

    # ---- resonance + external capture ----
    res = (20000.0, 2.0e-4)
    write("ecres_both", *ec_model("ecres_both", 2, 1, 1, 2, (1, 2, 3), up, ECRES_E, res=res))
    write("ecres_hs", *ec_model("ecres_hs", 2, 1, 1, 2, (1, 2, 3), up, ECRES_E))
    lv, segs, data = ec_model("grad", 2, 1, 1, 2, (1, 2, 3), up, ECRES_E, res=res, free=True)
    rng = np.random.default_rng(11)
    data = {"grad_upos.dat": [(en, a, 1.2e-8 * (1 + 0.3 * rng.standard_normal()), 1.0e-9) for en in ECRES_E for a in ANGLES[1:-1:2]]}
    write("grad", lv, segs[:1], data)
    lv, segs, data = ec_model("bad_linewindow", 2, 1, 1, 2, (1, 2, 3), up, [1.0], cascade_tail="  linewindow 30")
    write("bad_linewindow", lv, segs[:1], {k: v for k, v in data.items() if k.endswith("_upos.dat")})

    if binary:
        n = len(ANGLES)
        b, h = run_engine(binary, "ecres_both"), run_engine(binary, "ecres_hs")
        X = np.array([truth._x2(1, int(mu), TH) for (_, _, mu) in keys])
        out = []
        for i in range(len(ECRES_E)):
            ne = len(ECRES_E)
            sb, sh = b[ne * n + i, 3], h[ne * n + i, 3]
            dprim = b[ne * n + ne + i * n: ne * n + ne + (i + 1) * n, 3]
            A = np.vstack([np.c_[2 * (aec * ares) @ X, (ares ** 2) @ X], [[2 * aec @ ares, 1.0]]])
            y = np.r_[dprim - sh * (aec ** 2) @ X, sb - sh]
            (u, v), *_ = np.linalg.lstsq(A, y, rcond=None)
            print(f"ecres_both E_lab = {ECRES_E[i]}: |z|^2/|c|^2 = {v / sh:.4f}, cos = {u / np.sqrt(v * sh):+.4f}, "
                  f"primary fit residual {np.abs(A @ [u, v] - y).max() / sb:.1e}")
            out.append(weights_to_w(keys, sh * aec ** 2 + 2 * u * aec * ares + v * ares ** 2, 2, up)[0])
        save("ecres_both", np.array(out))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
