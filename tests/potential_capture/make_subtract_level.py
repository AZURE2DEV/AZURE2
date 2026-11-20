#!/usr/bin/env python3
"""The level that replaces a subtracted box pole, for check.sh section 10 (potential_background.md section 5.11).

The model: models/capture_a40.azr.in at a_c = 4.87 fm (VD, the s1/2 background with capture interior, B_bg = 0).  Its
s-wave box state with 2 nodes (E_p = 7.59 MeV, gamma_p^2 = 1.65 MeV; the engine's numbers, Session.potential_box_states)
is subtracted from the background and put back as a 1/2+ level:

  - particle channel: the formal level (E_p, gamma_p) with B = B_bg, given in Brune's parameters (they do not depend
    on B): E_b solves E_p - E_b - gamma_p^2 (S~(E_b) - B_bg) = 0, the observed width Gamma_p = 2 P~(E_b) gamma_p^2 /
    (1 + gamma_p^2 S~'(E_b)), with the dressed S~, P~ of the subtracted background (Session.channel_functions);
  - photon channels (E1 to 3/2- and to 1/2-): the formal amplitude of section 5.11,
        gamma_gamma = (1/2) (mu a)^1/2 / (hbar c) R_f gamma_p J_p,   J_p = int_0^a phi_p(r) u_f(r) r dr,
    phi_p = u_p/u_p(a) the box state at E_p (scipy DOP853, independent of the engine's Numerov), u_f the normalised
    bound state of make_reference.py, R_f = sqrt(alpha hbar c) Z_eff sqrt(48)/3 sqrt(2 J_f + 1) W(1 1 1/2 1/2; 0 J_f)
    the E1 factor of ECIntegral.cpp (long-wavelength).  For one level Brune's amplitudes are the formal ones
    (E_p - E - gamma^2 (S~(E) - B) = E_b - E - gamma^2 (S~(E) - S~(E_b)) by the definition of E_b): gamma~ = gamma, and
    only an observed width carries 1/(1 + gamma^2 S~'(E_b)) -- which the particle width Gamma_p does.

The engine reads a photon width as the observed total of the internal and the complex external amplitude at E_b and
subtracts the latter (physics/ subtract_external_capture); one run of the engine gives the external amplitude, from
which the script writes the input that gives the internal gamma~_gamma, into models/subtract_level.azr.in.

Run (needs the engine module and AZURE2_BIN): python3 make_subtract_level.py
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0]]
sys.path.insert(0, HERE)
import make_reference as m  # noqa: E402  (the potential, the bound states, six_j)

import pyazr  # noqa: E402

A_C, B_BG, NODES = 4.87, 0.0, 2
S_SEP = 1.4548867432  # the threshold (MeV of E_x): the p3/2 binding energy
BIN = os.environ["AZURE2_BIN"]


def base_model():
    src = open(os.path.join(HERE, "models", "capture_a40.azr.in")).read()
    out, inl = [], False
    for line in src.split("\n"):
        if line.startswith("<levels>"):
            inl = True
        elif line.startswith("</levels>"):
            inl = False
        t = line.split()
        if inl and len(t) >= 31 and t[26] == "0":
            t[27] = str(A_C)
            line = " ".join(t)
        out.append(line)
    return "\n".join(out)


def subtracted(text):
    old = "wave            l 0 J 1/2 depth 67.5 MeV background interior boundary 0.0 capture interior"
    assert old in text
    return text.replace(old, old.replace("boundary 0.0", f"boundary 0.0 subtract {NODES}"))


def with_level(text, Ex, gamma_p_eV, g32_eV, g12_eV):
    """The 1/2+ level of the group (the 40-MeV placeholder) at Ex with its widths."""
    out = []
    for line in text.split("\n"):
        t = line.split()
        if len(t) >= 31 and t[0] == "0.5" and t[1] == "1" and t[2] == "40":
            t[2] = repr(Ex)
            ir = t[5]
            t[11] = repr({"1": gamma_p_eV, "2": g32_eV, "3": g12_eV}[ir])
            line = " ".join(t)
        out.append(line)
    return "\n".join(out)


def run(text, workdir):
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    os.makedirs(os.path.join(workdir, "checks"), exist_ok=True)
    with open(os.path.join(workdir, "m.azr"), "w") as f:
        f.write(text)
    env = dict(os.environ, OMP_NUM_THREADS="2")
    subprocess.run([BIN, "extrapolate", "--no-ignore-externals", "m.azr"], cwd=workdir, env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)
    return open(os.path.join(workdir, "output", "parameters.out")).read()


def g_line(params, pair):
    """(g_int, Re g_ext, Im g_ext) of the photon channel to `pair` of the 1/2+ level (parameters.out)."""
    block = params.split("J = 0.5+")[1].split("J = ")[0]
    for line in block.split("\n"):
        if re.search(rf"R =\s+{pair}\s", line):
            g = float(re.search(r"g_int =\s+(\S+)", line).group(1))
            ext = re.search(r"g_ext =\s*\(([^,]+),([^)]+)\)", line)
            return g, float(ext.group(1)), float(ext.group(2))
    raise SystemExit(f"no photon channel to pair {pair} in parameters.out")


def g_int(params, pair):
    return g_line(params, pair)[0]


def main():
    tmp = tempfile.mkdtemp(prefix="subtract_level_")
    try:
        S0 = subtracted(base_model())
        s0 = os.path.join(tmp, "S0.azr")
        os.makedirs(os.path.join(tmp, "output"))
        open(s0, "w").write(S0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pyazr.Session(s0, extrapolations=True) as s:
                box = s.potential_box_states(1, l=0, J=0.5, radius=A_C, boundary=B_BG, e_min=-60.0, e_max=9.0)
                k = list(box["nodes"]).index(NODES)
                Ep, g2 = float(box["E"][k]), float(box["gamma2"][k])

                def chan(E):
                    cf = s.channel_functions(1, np.array([E]), l=0, J=0.5)
                    return float(cf["S"][0]), float(cf["P"][0])

                # E_b: E_p - E - gamma_p^2 (S~(E) - B) = 0, secant from E_p
                f = lambda E: Ep - E - g2 * (chan(E)[0] - B_BG)
                e0, e1 = Ep, Ep - 0.1
                f0, f1 = f(e0), f(e1)
                for _ in range(60):
                    e2 = e1 - f1 * (e1 - e0) / (f1 - f0)
                    e0, f0, e1, f1 = e1, f1, e2, f(e2)
                    if abs(e1 - e0) < 1e-12:
                        break
                Eb = e1
                h = 1e-4
                dS = (chan(Eb + h)[0] - chan(Eb - h)[0]) / (2 * h)
                P = chan(Eb)[1]
        norm = 1.0 + g2 * dS
        gt = math.sqrt(g2 / norm)
        Gamma_p_eV = 2.0 * P * gt ** 2 * 1e6

        # the photon amplitudes, independent of the engine's radial functions
        mu = m.MU
        o = m.outward(0, 1, Ep, A_C)
        ua = o.y[0, -1]
        phi_p = lambda r: o.sol(r)[0] / ua
        zeff = m.Z1 * m.M2 / (m.M1 + m.M2) - m.Z2 * m.M1 / (m.M1 + m.M2)
        gammas = {}
        for two_jf, pair in ((3, 2), (1, 3)):
            B, C, u_in, _ = m.bound_state(1, two_jf, 0.3, 4.0)
            Jp = m.quad(lambda r: phi_p(r) * u_in(r) * r, 0.0, A_C, epsabs=0, epsrel=1e-12, limit=400)[0]
            Jf = Fraction(two_jf, 2)
            W = (-1) ** 3 * m.six_j(Fraction(1), Fraction(1), Fraction(0), Fraction(1, 2), Fraction(1, 2), Jf)
            R = math.sqrt(m.ALPHA * m.HBARC) * zeff * math.sqrt(48.0) / 3.0 * math.sqrt(two_jf + 1) * W
            gg = 0.5 * math.sqrt(mu * A_C) / m.HBARC * R * math.sqrt(g2) * Jp
            gammas[pair] = gg  # Brune's amplitude of a single level is the formal one (only the observed width has 1/N)
        Ex = Eb + S_SEP
        print(f"# E_p = {Ep!r} MeV, gamma_p^2 = {g2!r} MeV; E_b = {Eb!r} MeV, S~'(E_b) = {dS!r} /MeV, P~ = {P!r}, "
              f"Gamma_p = {Gamma_p_eV!r} eV; Brune photon amplitudes (internal) {gammas}", file=sys.stderr)

        # The photon inputs.  The engine reads a photon width as the observed total and subtracts the complex
        # external amplitude g_ext at E_b (physics/ subtract_external_capture): g_int = sign(g) sqrt(g^2 - Im^2) - Re,
        # g = sign(Gamma) sqrt(|Gamma| / (2 P_gamma)).  One run gives g_ext and 2 P_gamma (from its g_int); the input
        # for the target t is then Gamma = sign(t + Re) 2 P_gamma ((t + Re)^2 + Im^2).
        p1 = run(with_level(S0, Ex, Gamma_p_eV, 100.0, 100.0), os.path.join(tmp, "r1"))
        inputs = {}
        for pair in (2, 3):
            gi, re_, im_ = g_line(p1, pair)
            twoP = 100.0 / ((gi + re_) ** 2 + im_ ** 2)
            t = gammas[pair]
            inputs[pair] = math.copysign(twoP * ((t + re_) ** 2 + im_ ** 2), t + re_)
        final = with_level(S0, Ex, Gamma_p_eV, inputs[2], inputs[3])
        p3 = run(final, os.path.join(tmp, "r3"))
        for pair in (2, 3):
            print(f"# photon channel to pair {pair}: input {inputs[pair]!r} eV, g_int {g_int(p3, pair)!r} "
                  f"(target {gammas[pair]!r})", file=sys.stderr)
        header = ("# tests/potential_capture/models/subtract_level.azr.in: made by make_subtract_level.py "
                  "(potential_background.md section 5.11)\n")
        with open(os.path.join(HERE, "models", "subtract_level.azr.in"), "w") as f:
            f.write(final)
        print(header.strip(), file=sys.stderr)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
