#!/usr/bin/env python3
"""T3: 8Be -- identical bosons, neutrons, a beta pair, a hybrid potential channel.  Holds the true parameters.

    python3 generate.py [path/to/AZURE2]      (or AZURE2_BIN / AZURE2_BINARY in the environment)

Writes T3.azr (the truth), T3_start.azr (free parameters moved: the start of the recovery fit), data/*.dat (the
engine's observables at the truth with Gaussian noise, seed SEED) and the <scatteringLengths> constraints (the
truth's own n + 7Be s-wave lengths with noise of the stated sigma).  The level scheme is schematic: it borrows the
8Be states near the thresholds, but the widths and the background poles are chosen for coverage, not taken from an
evaluation.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import upgrade_project as up  # noqa: E402

SEED = 8
M = up.MASS
S_A = up.separation(M["4He"], M["4He"], M["8Be"])   # -0.0918 MeV (8Be is unbound to alpha + alpha)
S_P = up.separation(M["p"], M["7Li"], M["8Be"])     # 17.2545 MeV
S_N = up.separation(M["n"], M["7Be"], M["8Be"])     # 18.8995 MeV
Q_B = (M["8Li"] - M["8Be"]) * up.U_MEV              # 16.0045 MeV, 8Li(beta-)8Be with atomic masses
SL_SIGMA = 0.1                                      # fm, uncertainty of the scattering-length constraints

# The hybrid nuclear potential (Woods-Saxon) bends the alpha + alpha channel functions only.
POTENTIAL = ("useHybridPotential=0\nuseAdaptiveGrid=1\npotentialType=0\nV0=10\nR=3.0\na=0.6\n"
             "pair=1\nuseHybridPotential=1\npotentialType=0\nV0=10\nR=3.0\na=0.6\n")


def model(scattering_lengths=""):
    aa = up.Pair(1, 0.0, 1, 0.0, 1, M["4He"], M["4He"], 2, 2, S_A, R=5.5)
    p7li = up.Pair(2, 0.5, 1, 1.5, -1, M["p"], M["7Li"], 1, 3, S_P, R=4.0, g1=5.5856946893, g2=2.1709)
    n7be = up.Pair(3, 0.5, 1, 1.5, -1, M["n"], M["7Be"], 0, 4, S_N, R=4.0, g1=-3.82608545, g2=-0.933)
    beta = up.Pair(4, 0.5, 1, 2.0, 1, M["e"], M["8Li"], -1, 3, Q_B, ptype=20, R=5.5)
    C = up.Channel
    levels = [
        # 2+ first excited state: alpha + alpha d wave, fed by the 8Li 2+ Gamow-Teller decay (beta l = 1)
        up.Level(2.0, 1, 3.03, [C(aa, 0, 4, 1.5e6, fixed=False), C(beta, 2, 2, 1.0, fixed=False)],
                 fixed=False, name="2+ 3.03"),
        up.Level(4.0, 1, 11.35, [C(aa, 0, 8, 3.5e6)], name="4+ 11.35"),
        # 1+ just above the p + 7Li threshold (unnatural parity: no alpha + alpha channel)
        up.Level(1.0, 1, 17.64, [C(p7li, 2, 2, 11000.0, fixed=False), C(p7li, 4, 2, 1000.0)],
                 fixed=False, name="1+ 17.64"),
        # 2- just above the n + 7Be threshold: s waves in both nucleon channels, s = 2
        up.Level(2.0, -1, 19.00, [C(p7li, 4, 0, 120000.0, fixed=False), C(n7be, 4, 0, 60000.0, fixed=False)],
                 fixed=False, name="2- 19.00"),
        # background poles (fixed); the 2+ one carries the 7Li(p,alpha) and 7Be(n,alpha) strength
        up.Level(0.0, 1, 25.0, [C(aa, 0, 0, 1.0e7)], name="0+ background"),
        up.Level(2.0, 1, 22.0, [C(aa, 0, 4, 4.0e6), C(p7li, 4, 2, 1.0e6), C(n7be, 4, 2, 1.0e6),
                                 C(beta, 2, 2, -0.5)],
                 name="2+ background"),
        up.Level(1.0, -1, 22.0, [C(p7li, 2, 0, 2.0e6), C(n7be, 2, 0, 2.0e6)], name="1- background"),
    ]
    S = up.DataSegment
    aa_e = up.grid(1.0, 8.0, 0.1)
    data = [
        S("aa_cm050", 1, 1, 4, [(e, 50.0) for e in aa_e], ("rel", 0.03), amin=50, amax=50),
        S("aa_cm090", 1, 1, 4, [(e, 90.0) for e in aa_e], ("rel", 0.03), amin=90, amax=90),
        S("aa_phase_d", 1, 1, 2, [(e, 0.0) for e in up.grid(1.0, 8.0, 0.25)], ("abs", 1.0), phase_jl=(2.0, 2),
          keep=lambda t: abs(t) < 80.0),
        # energy column: the lab energy of one alpha from the decay at rest (E_cm(alpha + alpha) = 2 E_alpha)
        S("beta_alpha", 4, 1, 0, [(e, 0.0) for e in up.grid(0.4, 6.0, 0.08)], ("rel", 0.03)),
        S("p7li_pp_lab150", 2, 2, 1, [(e, 150.0) for e in up.grid(0.30, 2.40, 0.025)], ("rel", 0.03),
          amin=150, amax=150),
        S("p7li_pn", 2, 3, 0, [(e, 0.0) for e in up.grid(1.90, 2.40, 0.01)], ("rel", 0.05)),
        S("p7li_pa", 2, 1, 0, [(e, 0.0) for e in up.grid(0.30, 2.40, 0.05)], ("rel", 0.05)),
        S("n7be_total", 3, 3, 0, [(e, 0.0) for e in up.grid(0.01, 0.40, 0.01)], ("rel", 0.03), tail="total 1"),
        S("n7be_np", 3, 2, 0, [(e, 0.0) for e in up.grid(0.01, 0.40, 0.01)], ("rel", 0.05)),
    ]
    T = up.TestSegment
    tests = [
        T(3, 3, 0.01, 0.40, 0.01, extra="total 1"),       # n + 7Be total (optical theorem) ...
        T(3, 3, 0.01, 0.40, 0.01),                        # ... and its parts: elastic,
        T(3, 2, 0.01, 0.40, 0.01),                        # (n,p),
        T(3, 1, 0.01, 0.40, 0.01),                        # (n,alpha)
        T(1, 1, 1.0, 8.0, 0.25, isdiff=2, extra="0 0"),   # alpha + alpha s, d, g phase shifts
        T(1, 1, 1.0, 8.0, 0.25, isdiff=2, extra="2 2"),
        T(1, 1, 1.0, 8.0, 0.25, isdiff=2, extra="4 4"),
        T(1, 1, 4.0, 4.0, 0.0, 10.0, 170.0, 10.0, isdiff=1),  # alpha + alpha angular distribution (symmetric)
        T(4, 1, 0.1, 7.5, 0.1),                           # beta-delayed alpha spectrum (alpha energy, as data)
    ]
    return up.Project("T3", [aa, p7li, n7be, beta], levels, data, tests, potential=POTENTIAL,
                      scattering_lengths=scattering_lengths, rate=(2, 1))


START_SHIFT = {"E": 0.01, "W": 1.15}


def start_model(sl):
    m = model(sl)
    for _, ref in m.free_parameters():
        v = m.get(ref)
        m.set(ref, v + START_SHIFT["E"] if ref[0] == "E" else v * START_SHIFT["W"])
    return m


def truth_scattering_lengths(binary, m):
    """Complex a0 of n + 7Be for s = 1 and 2 from the engine at the truth (output/scatteringLengths.out)."""
    work = tempfile.mkdtemp(prefix="T3_sl_")
    try:
        ones = {s.name: [(e, a, 1.0, 1.0) for e, a in s.points] for s in m.data}
        azr = m.write(work, data=ones)
        up.run(binary, work, azr, up.CALC)
        out = {}
        with open(os.path.join(work, "output", "scatteringLengths.out")) as f:
            txt = f.read()
        for line in txt.splitlines():
            s = line.split()
            # "  Pair  Particles  J^pi  x  Re a0 (fm)  Im a0 (fm)": J = s for l = 0
            if len(s) >= 6 and s[0] == "3" and s[2] in ("1-", "2-"):
                out[int(s[2][0])] = (float(s[4]), float(s[5]))
        if set(out) != {1, 2}:
            raise RuntimeError("scatteringLengths.out: n + 7Be lengths not found:\n" + txt)
        return out
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    binary = up.find_binary(sys.argv[1] if len(sys.argv) > 1 else None)
    m = model()
    truth = up.truth_values(binary, m)
    for seg in m.data:
        if seg.keep is not None:
            keep = [seg.keep(t) for t in truth[seg.name]]
            seg.points = [p for p, k in zip(seg.points, keep) if k]
            truth[seg.name] = truth[seg.name][keep]
    data = up.synthesize(m, truth, SEED)
    a0 = truth_scattering_lengths(binary, m)
    rng = np.random.default_rng(SEED + 1000)
    sl = ("# active pairKey s Re sigma Im sigma; truth %s fm; Gaussian noise sigma %g fm, seed %d\n"
          % ("; ".join(f"s={s}: {a0[s][0]:.6f} {a0[s][1]:+.6f}i" for s in sorted(a0)), SL_SIGMA, SEED + 1000))
    for s in sorted(a0):
        re_, im_ = a0[s][0] + rng.normal(0.0, SL_SIGMA), a0[s][1] + rng.normal(0.0, SL_SIGMA)
        sl += f"1 3 {s} {re_:.6f} {SL_SIGMA} {im_:.6f} {SL_SIGMA}\n"
    m.scattering_lengths = sl
    m.write(HERE)
    start_model(sl).write(HERE, "T3_start")
    for seg in m.data:
        kind, size = seg.noise
        hdr = (f"# T3 synthetic data ({seg.name}), generate.py: engine truth + Gaussian noise "
               f"({'%g %%' % (100 * size) if kind == 'rel' else '%g absolute' % size}), seed {SEED}\n"
               "# E_lab (MeV; for the beta-delayed spectrum the energy of one alpha, decay at rest)  angle (deg)  value  "
               "uncertainty   [lab frame unless the segment is c.m.: alpha + alpha cross sections]\n")
        up.write_data(os.path.join(HERE, "data", seg.name + ".dat"), data[seg.name], seg, hdr)
    with open(os.path.join(HERE, "temperatures.dat"), "w") as f:   # 7Li(p,alpha)alpha rate, <lastRun> 2 -> 1
        f.write("0.01\n0.03\n0.1\n0.3\n1.0\n3.0\n")
    n = sum(len(s.points) for s in m.data)
    print(f"T3: {len(m.data)} segments, {n} points, {len(m.free_parameters())} free parameters, a0 = {a0}")


if __name__ == "__main__":
    main()
