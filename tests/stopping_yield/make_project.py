#!/usr/bin/env python3
"""Writes the synthetic projects of tests/stopping_yield with LunaSRIM's stopping (run once; the outputs are committed).

    python3 tests/stopping_yield/make_project.py        (any Python that imports lunasrim 1.0.0)

models/stopping_yield.azr  the narrow resonance of tests/layer_yield on TiN layers with `stopping evaluation`
models/fit.azr             two scans on two TiN layers, N1, N2 and the one stopping scale of (TiN, p) free
models/data/*.dat          their data files
models/truth.json          every number the projects were made from, and LunaSRIM's numbers the check compares with

The stopping is SRIM-2013's, computed here by LunaSRIM (J. Skowronski, github.com/skowrons94/LunaSRIM; cite SRIM:
J. F. Ziegler et al., NIM B 268 (2010) 1818), never by the engine: per atom of TiN, Bragg's rule over Ti and N with the
phase factor of a solid (lunasrim.stopping.electronic_stopping, phase 0) and SRIM's nuclear stopping with SRIM's atomic
weights; the layers' exit energies by lunasrim.straggling.energy_loss.  Per active nucleus (one N atom per formula
unit, atom fraction 1/2): eps_active = 1e-21 x S_atom / (1/2) MeV cm^2.
"""
import json
import math
import os
import random

from lunasrim import data_files, stopping, straggling

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
DATA = os.path.join(MODELS, "data")

HBARC = 197.3269804      # MeV fm (physics/Constants.h)
UMEV = 931.49410242      # MeV per u

# p + 15N -> 16O -> alpha + 12C: one isolated narrow 1- level (as tests/layer_yield), and the n + 12C hard sphere (a
# neutral projectile for a refusal)
P15N = dict(j1=0.5, pi1=1, j2=0.5, pi2=-1, m1=1.007825, m2=15.000109, z1=1, z2=7, S=12.1274, R=4.85)
A12C = dict(j1=0.0, pi1=1, j2=0.0, pi2=1, m1=4.002603, m2=12.0, z1=2, z2=6, S=7.16192, R=5.43)
N12C = dict(j1=0.5, pi1=1, j2=0.0, pi2=1, m1=1.008665, m2=12.0, z1=0, z2=6, S=4.94636, R=0.1)
LEVEL = dict(J=1.0, pi=-1, Er_cm=0.9, Gp_eV=60.0, Ga_eV=40.0)

TIN = [(22, 1.0), (7, 1.0)]      # (Z, atoms per formula unit); the active element is N (Z = 7)
F_ACTIVE = 0.5                   # its atom fraction
WEIGHT = {z: data_files.load_scoef()[z].natural_mass for z in (22, 7)}   # SRIM's atomic weights


def eps_atom(e_mev):
    """LunaSRIM's stopping of the proton (mass of <levels>) per atom of solid TiN, eV/(1e15 atoms/cm^2)."""
    ek, m = 1000.0 * e_mev, P15N["m1"]
    tot = sum(c for _, c in TIN)
    return sum(c / tot * (stopping.electronic_stopping(1, m, z, ek, 0) + stopping.nuclear_stopping(1, m, z, WEIGHT[z], ek))
               for z, c in TIN)


def eps_active(e_mev):
    """Per active nucleus, lab MeV cm^2."""
    return 1e-21 * eps_atom(e_mev) / F_ACTIVE


def exit_energy(e_mev, n_active, scale=1.0):
    """LunaSRIM's energy (MeV) after n_active N atoms/cm^2 of TiN (2 n_active atoms of the material), stopping s eps."""
    els = [(z, c, WEIGHT[z]) for z, c in TIN]
    r = straggling.energy_loss(1, P15N["m1"], 1000.0 * e_mev, els, 5.4, scale * n_active / F_ACTIVE, unit="atoms/cm2",
                               phase=0, steps=4000)
    return r["energy_out"] / 1000.0


def entrance_for_exit(e_out, n_active, scale=1.0):
    """The beam energy that leaves the layer at e_out (bisection on exit_energy)."""
    lo, hi = e_out, e_out + 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if exit_energy(mid, n_active, scale) < e_out:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def f_cm(p):
    return p["m2"] / (p["m1"] + p["m2"])


def bw_yield(E0, e_out, eps_r, scale=1.0):
    """The Breit-Wigner integral over the layer, widths and k constant, epsilon = s eps_r across the resonance
    (README.md, eq. (S1)); per incident particle."""
    p, L = P15N, LEVEL
    f = f_cm(p)
    mu = p["m1"] * p["m2"] / (p["m1"] + p["m2"]) * UMEV
    k2 = 2.0 * mu * L["Er_cm"] / HBARC ** 2
    omega = (2 * L["J"] + 1) / ((2 * p["j1"] + 1) * (2 * p["j2"] + 1))
    G = (L["Gp_eV"] + L["Ga_eV"]) * 1e-6
    wg = omega * L["Gp_eV"] * L["Ga_eV"] / (L["Gp_eV"] + L["Ga_eV"]) * 1e-6
    pref = 1e-24 * (2.0 * math.pi / k2 * 0.01) * wg / (f * scale * eps_r)
    return pref * (math.atan(2 * (f * E0 - L["Er_cm"]) / G) - math.atan(2 * (f * e_out - L["Er_cm"]) / G))


def level_line(J, pi, ex, fix, ir, two_s, two_l, lid, chfix, gamma, p, sin):
    return (f"{J:g} {pi:d} {ex:.6f} {fix:d} 1 {ir} {two_s} {two_l} {lid} 1 {chfix} {gamma:.10g} "
            f"{p['j1']:g} {p['pi1']} {p['j2']:g} {p['pi2']} 0.0 {p['m1']} {p['m2']} {p['z1']} {p['z2']} "
            f"{sin} {p['S']} 0 0 0.0 0 {p['R']} 0 0 0")


def levels():
    ex = LEVEL["Er_cm"] + P15N["S"]
    return "\n".join([level_line(1, -1, ex, 1, 1, 2, 0, 1, 1, LEVEL["Gp_eV"], P15N, P15N["S"]),
                      level_line(1, -1, ex, 1, 2, 0, 2, 1, 1, LEVEL["Ga_eV"], A12C, P15N["S"]), "",
                      level_line(0.5, 1, 30.0, 1, 3, 1, 0, 2, 1, 0.0, N12C, N12C["S"]), ""])


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
format=3
</config>
"""


def write_data(name, energies, values=None, errors=None):
    with open(os.path.join(DATA, name), "w") as f:
        f.write("# E_lab(MeV) unused value error (make_project.py)\n")
        for i, e in enumerate(energies):
            v = 1.0 if values is None else values[i]
            d = 1.0 if errors is None else errors[i]
            f.write(f"{e:.12g} 0 {v:.10g} {d:.10g}\n")


def dataset(name, target, data, points=40):
    return (f"dataset {name}\n    reaction       15N(p,a)12C\n    quantity       yield\n    unit           1/particle\n"
            f"    file           data/{data}\n    columns        E_lab[MeV] _ value error\n    norm           1 fixed\n"
            f"    target         {target}\n    quadrature     points {points}\n")


def target(name, n_active, stopping_line, fit=""):
    return (f"target {name}\n    kind           layer\n    active         15N\n    composition    Ti 1 N 1\n"
            f"    areal_density  {n_active:g} atoms/cm2 {fit or 'fixed'}\n    stopping       {stopping_line}\n"
            f"    straggling     none\n")


def main():
    os.makedirs(DATA, exist_ok=True)
    f = f_cm(P15N)
    Er_lab = LEVEL["Er_cm"] / f
    eps_r = eps_active(Er_lab)
    truth = dict(levels=dict(p15N=P15N, a12C=A12C, level=LEVEL), Er_lab=Er_lab, eps_active_Er=eps_r,
                 eps_atom_Er=eps_atom(Er_lab), f_active=F_ACTIVE, sets={},
                 source="lunasrim " + __import__("lunasrim").__version__)
    sets = truth["sets"]
    data, targets = "", ""

    # (S1) the plateau of a thick layer: N = 4e18 (Delta ~ 60 keV), E0 above the resonance
    N_p = 4e18
    E_p = [round(Er_lab + d, 10) for d in (0.002, 0.004, 0.008, 0.016)]
    write_data("plateau.dat", E_p)
    out_p = [exit_energy(e, N_p) for e in E_p]
    sets["plateau"] = dict(N=N_p, E0=E_p, E_out=out_p, Y=[bw_yield(e, o, eps_r) for e, o in zip(E_p, out_p)])
    data += dataset("plateau", "TP", "plateau.dat")
    targets += target("TP", N_p, "evaluation")
    # (S3) the same layer with the stopping as a table of the library's own values (check.sh writes data/tin_p.tab)
    data += dataset("plateau_table", "TT", "plateau.dat")
    targets += target("TT", N_p, 'table "data/tin_p.tab" MeV eV/1e15atoms/cm2')

    # (S2) the step of a layer of N = 5e17 (Delta ~ 8 keV): its two edges every 10 eV and its middle
    N_s = 5e17
    E_back = entrance_for_exit(Er_lab, N_s)
    edge = [0.00001 * i for i in range(-40, 41)]
    E_s = sorted([round(Er_lab + d, 10) for d in edge] + [round(E_back + d, 10) for d in edge] +
                 [round(Er_lab + (E_back - Er_lab) * x, 10) for x in (0.3, 0.4, 0.5, 0.6, 0.7)])
    write_data("step.dat", E_s)
    # the Breit-Wigner yield of the scan with LunaSRIM's E_out at every energy: its half-height points are measured
    # as the engine's are (the Lorentzian tails lower the middle of the step, which moves both points by ~0.3 eV)
    Y_s = [bw_yield(e, exit_energy(e, N_s), eps_r) for e in E_s]
    sets["step"] = dict(N=N_s, E_front=Er_lab, E_back=E_back, width=E_back - Er_lab, E0=E_s, Y=Y_s)
    data += dataset("step", "TS", "step.dat")
    targets += target("TS", N_s, "evaluation")

    head = ("<data>\n# Synthetic (make_project.py): a narrow resonance on TiN layers with SRIM-2013's stopping.  The data\n"
            "# values are placeholders; check.py compares the model values with LunaSRIM's numbers in truth.json.\n")
    with open(os.path.join(MODELS, "stopping_yield.azr"), "w") as out:
        out.write(CONFIG + "<levels>\n" + levels() + "</levels>\n" + head + data + "</data>\n<targets>\n" + targets +
                  "</targets>\n")

    # (S4) the fit: two layers of TiN, one stopping scale s_true = 1.03 (prior 1 4 %), 1 % noise
    rng = random.Random(20261111)
    s_true = 1.03
    fit = dict(s_true=s_true, s_prior=[1.0, 0.04], seed=20261111, sets={})
    fdata, ftargets = "", ""
    for name, tname, n_true, n_start in (("scan1", "F1", 6e17, 5e17), ("scan2", "F2", 9e17, 1e18)):
        back = entrance_for_exit(Er_lab, n_true, s_true)
        E = [round(Er_lab - 0.0006 + 0.0003 * i, 10) for i in range(int((back - Er_lab + 0.0012) / 0.0003) + 1)]
        Y = [bw_yield(e, exit_energy(e, n_true, s_true), eps_r, s_true) for e in E]
        err = [0.01 * max(y, 0.02 * max(Y)) for y in Y]
        obs = [y + rng.gauss(0.0, d) for y, d in zip(Y, err)]
        write_data(f"{name}.dat", E, obs, err)
        fit["sets"][name] = dict(target=tname, N_true=n_true, N_start=n_start, E0=E)
        fdata += dataset(name, tname, f"{name}.dat")
        ftargets += target(tname, n_start, "evaluation scale 1 fit prior 1 4%", fit="fit")
    truth["fit"] = fit
    with open(os.path.join(MODELS, "fit.azr"), "w") as out:
        out.write(CONFIG + "<levels>\n" + levels() + "</levels>\n<data>\n# Synthetic scans (make_project.py): N1, N2 "
                  "and the stopping scale of (TiN, p) free\n" + fdata + "</data>\n<targets>\n" + ftargets + "</targets>\n")
    with open(os.path.join(MODELS, "truth.json"), "w") as out:
        json.dump(truth, out, indent=1)
        out.write("\n")


if __name__ == "__main__":
    main()
