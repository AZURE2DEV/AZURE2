#!/usr/bin/env python3
"""Writes the synthetic projects of tests/layer_yield (run once; the outputs are committed).

models/layer_yield.azr  one data set per limit of docs/upgrade/energy_grid.md section 5.3 (README.md, "Data sets")
models/fit.azr          a resonance scan with noise from a fixed seed; N and Gamma_p free
models/data/*.dat       their data files
models/truth.json       every number the projects were made from (levels, masses, targets, energies, the fit's truth)

The data values of layer_yield.azr are placeholders (1 +- 1): check.py compares the engine's model values with closed
forms it computes itself.  The fit's data are the closed form of the narrow resonance (check.py's arctan, written out
here a second time on purpose: the generator and the check must not share code) times 1 + 2 % Gaussian noise.

    python3 tests/layer_yield/make_project.py
"""
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
DATA = os.path.join(MODELS, "data")

HBARC = 197.3269804      # MeV fm (physics/Constants.h)
UMEV = 931.49410242      # MeV per u

# p + 15N -> 16O -> alpha + 12C: one isolated narrow 1- level (proton s wave, alpha p wave)
P15N = dict(j1=0.5, pi1=1, j2=0.5, pi2=-1, m1=1.007825, m2=15.000109, z1=1, z2=7, S=12.1274, R=4.85)
A12C = dict(j1=0.0, pi1=1, j2=0.0, pi2=1, m1=4.002603, m2=12.0, z1=2, z2=6, S=7.16192, R=5.43)
# n + 12C -> 13C: a 1/2+ s-wave channel of radius 0.1 fm and a level at E_x = 30 MeV without width (hard sphere only)
N12C = dict(j1=0.5, pi1=1, j2=0.0, pi2=1, m1=1.008665, m2=12.0, z1=0, z2=6, S=4.94636, R=0.1)
LEVEL = dict(J=1.0, pi=-1, Er_cm=0.9, Gp_eV=60.0, Ga_eV=40.0)

# stopping formulas (lab MeV cm^2 per active nucleus, x = E_lab / MeV)
EPS_CONST = 5e-21                      # the resonance layers: a constant
EPS_POW = dict(a0=4.5e-20, a1=-0.6)    # thin and thick layers: a0 x^a1
EPS_INV = dict(a0=1.2e-24, a1=-1.0)    # the neutron layers: a0 / x


def f_cm(p):
    """E_cm / E_lab of a pair (non-relativistic, the engine's Segment::CMEnergyFactor)."""
    return p["m2"] / (p["m1"] + p["m2"])


def level_line(J, pi, ex, fix, ir, two_s, two_l, lid, chfix, gamma, p, sin):
    return (f"{J:g} {pi:d} {ex:.6f} {fix:d} 1 {ir} {two_s} {two_l} {lid} 1 {chfix} {gamma:.10g} "
            f"{p['j1']:g} {p['pi1']} {p['j2']:g} {p['pi2']} 0.0 {p['m1']} {p['m2']} {p['z1']} {p['z2']} "
            f"{sin} {p['S']} 0 0 0.0 0 {p['R']} 0 0 0")


def levels(gp_fixed=1):
    ex = LEVEL["Er_cm"] + P15N["S"]
    return "\n".join([
        level_line(1, -1, ex, 1, 1, 2, 0, 1, gp_fixed, LEVEL["Gp_eV"], P15N, P15N["S"]),
        level_line(1, -1, ex, 1, 2, 0, 2, 1, 1, LEVEL["Ga_eV"], A12C, P15N["S"]),
        "",
        level_line(0.5, 1, 30.0, 1, 3, 1, 0, 2, 1, 0.0, N12C, N12C["S"]),
        ""])


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


def range_pow(E0, N, a0, a1):
    """Exact energy loss (lab MeV) of a layer of N atoms/cm^2 with epsilon = a0 x^a1: int_{E_out}^{E0} dE/eps = N."""
    b = 1.0 - a1
    return E0 - (E0 ** b - a0 * b * N) ** (1.0 / b)


def write_data(name, energies, values=None, errors=None):
    with open(os.path.join(DATA, name), "w") as f:
        f.write("# E_lab(MeV) unused value error (make_project.py)\n")
        for i, e in enumerate(energies):
            v = 1.0 if values is None else values[i]
            d = 1.0 if errors is None else errors[i]
            f.write(f"{e:.12g} 0 {v:.10g} {d:.10g}\n")


def dataset(name, reaction, quantity, unit, data, target=None, beam=None, extra="", points=20):
    s = f"dataset {name}\n    reaction       {reaction}\n    quantity       {quantity}\n    unit           {unit}\n"
    s += f"    file           data/{data}\n    columns        E_lab[MeV] _ value error\n    norm           1 fixed\n"
    if target:
        s += f"    target         {target}\n"
    if beam:
        s += f"    beam           {beam}\n"
    if target or beam:
        s += f"    quadrature     points {points}\n"
    return s + extra


def bw_layer_yield(E0, Delta):
    """Arctan form of a narrow resonance on a layer of constant lab stopping (README.md, eq. (R3)), per particle."""
    p, L = P15N, LEVEL
    f = f_cm(p)
    mu = p["m1"] * p["m2"] / (p["m1"] + p["m2"]) * UMEV
    k2 = 2.0 * mu * L["Er_cm"] / HBARC ** 2                     # fm^-2
    omega = (2 * L["J"] + 1) / ((2 * p["j1"] + 1) * (2 * p["j2"] + 1))
    G = (L["Gp_eV"] + L["Ga_eV"]) * 1e-6
    wg = omega * L["Gp_eV"] * L["Ga_eV"] / (L["Gp_eV"] + L["Ga_eV"]) * 1e-6   # MeV
    pref = 1e-24 * (2.0 * math.pi / k2 * 0.01) * wg / (f * EPS_CONST)          # pi/k^2 in b
    return pref * (math.atan(2 * (f * E0 - L["Er_cm"]) / G) - math.atan(2 * (f * (E0 - Delta) - L["Er_cm"]) / G))


def main():
    os.makedirs(DATA, exist_ok=True)
    f = f_cm(P15N)
    Er_lab = LEVEL["Er_cm"] / f
    truth = dict(levels=dict(p15N=P15N, a12C=A12C, n12C=N12C, level=LEVEL), eps_const=EPS_CONST, eps_pow=EPS_POW,
                 eps_inv=EPS_INV, Er_lab=Er_lab, sets={})
    data = ""
    sets = truth["sets"]

    # 1. constant cross section: n + 12C hard sphere, epsilon = a0/x, N = 1e20 (no spread; with a 0.3 keV spread;
    #    with an isotope abundance 0.5)
    N_n = 1e20
    E_n = [0.018, 0.02, 0.025, 0.03]
    write_data("const.dat", E_n)
    for name, target, beam in (("const", "LN", None), ("const_spread", "LN", "B03"), ("const_iso", "LNI", None)):
        data += dataset(name, "12C(n,n)12C", "yield", "1/particle", "const.dat", target, beam, points=100 if beam else 10)
        sets[name] = dict(E0=E_n, N=N_n)
    # 2. thin layer: N = 5e15 of epsilon = a0 x^a1 (Delta ~ 0.3 keV), its target-averaged sigma, and bare sigma at
    #    E0 - Delta/2 and +- 5 keV around it
    N_thin = 5e15
    E_thin = [0.60, 0.65, 0.70, 0.75, 0.80]
    write_data("thin.dat", E_thin)
    mid = [e - 0.5 * range_pow(e, N_thin, EPS_POW["a0"], EPS_POW["a1"]) for e in E_thin]
    bare = []
    for m in mid:
        bare += [m - 0.005, m, m + 0.005]
    write_data("thin_bare.dat", bare)
    data += dataset("thin", "15N(p,a)12C", "yield", "1/particle", "thin.dat", "TT")
    data += dataset("thin_sigma", "15N(p,a)12C", "sigma", "b", "thin.dat", "TT")
    data += dataset("thin_bare", "15N(p,a)12C", "sigma", "b", "thin_bare.dat")
    sets["thin"] = dict(E0=E_thin, N=N_thin, bare=bare)
    # 3. thick layer: N = 1.1e19 (E_out 0.30-0.45 MeV) and the infinite target at E0 and at E_out
    N_thick = 1.1e19
    E_thick = [1.00, 1.05, 1.10]
    E_out = [e - range_pow(e, N_thick, EPS_POW["a0"], EPS_POW["a1"]) for e in E_thick]
    write_data("thick.dat", E_thick)
    write_data("thick_inf.dat", E_thick + E_out)
    data += dataset("thick", "15N(p,a)12C", "yield", "1/particle", "thick.dat", "TK", points=60)
    data += dataset("thick_inf", "15N(p,a)12C", "yield", "1/particle", "thick_inf.dat", "TI", points=60)
    sets["thick"] = dict(E0=E_thick, N=N_thick, E_out=E_out)
    # 4. narrow resonance on a layer of constant stopping: N = 1e18 (Delta = 5 keV lab), E0 every 0.2 keV
    N_res = 1e18
    E_res = [round(Er_lab - 0.0015 + 0.0002 * i, 10) for i in range(41)]
    write_data("res.dat", E_res)
    data += dataset("res", "15N(p,a)12C", "yield", "1/particle", "res.dat", "TR", points=40)
    sets["res"] = dict(E0=E_res, N=N_res)
    # 5. per unit charge: 1/uC with charge state 1 and 2 on the same layer, five energies of the scan
    E_q = E_res[10:35:6]
    write_data("res_q.dat", E_q)
    for q in (1, 2):
        data += dataset(f"res_q{q}", "15N(p,a)12C", "yield", "1/uC", "res_q.dat", "TR", points=40,
                        extra=f"    charge_state   {q}\n")
    sets["res_q"] = dict(E0=E_q, N=N_res)
    # 6. the derivative with respect to N: the same layer with N free (prior 10 %), E_out on the resonance; the bare
    #    sigma at those E_out
    Delta = N_res * EPS_CONST
    E_d = [round(Er_lab + Delta + d, 10) for d in (-0.00015, -0.00005, 0.0, 0.00005, 0.00015)]
    write_data("deriv.dat", E_d)
    write_data("deriv_bare.dat", [e - Delta for e in E_d])
    data += dataset("deriv", "15N(p,a)12C", "yield", "1/particle", "deriv.dat", "TD", points=40)
    data += dataset("deriv_bare", "15N(p,a)12C", "sigma", "b", "deriv_bare.dat")
    sets["deriv"] = dict(E0=E_d, N=N_res)
    # 7. the bare sigma across the resonance, to measure how far the R-matrix line shape is from the Breit-Wigner of
    #    eq. (R1) at E_cm - E_r = 0, +-1, +-3, +-10, +-30 Gamma
    G = (LEVEL["Gp_eV"] + LEVEL["Ga_eV"]) * 1e-6
    E_b = [(LEVEL["Er_cm"] + m * G) / f for m in (-30, -10, -3, -1, 0, 1, 3, 10, 30)]
    write_data("res_bare.dat", E_b)
    data += dataset("res_bare", "15N(p,a)12C", "sigma", "b", "res_bare.dat")
    sets["res_bare"] = dict(E0=E_b)

    targets = f"""<targets>
target LN
    kind           layer
    active         12C
    areal_density  {N_n:g} atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0={EPS_INV['a0']:g} a1={EPS_INV['a1']:g} MeV*cm2/atom
    straggling     none
target LNI
    kind           layer
    active         12C
    areal_density  {N_n:g} atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0={EPS_INV['a0']:g} a1={EPS_INV['a1']:g} MeV*cm2/atom
    straggling     none
    isotope        12C 0.5 fixed
target TT
    kind           layer
    active         15N
    areal_density  {N_thin:g} atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0={EPS_POW['a0']:g} a1={EPS_POW['a1']:g} MeV*cm2/atom
    straggling     none
target TK
    kind           layer
    active         15N
    areal_density  {N_thick:g} atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0={EPS_POW['a0']:g} a1={EPS_POW['a1']:g} MeV*cm2/atom
    straggling     none
target TI
    kind           infinite
    active         15N
    stopping       formula "a0*x^a1" a0={EPS_POW['a0']:g} a1={EPS_POW['a1']:g} MeV*cm2/atom
    straggling     none
target TR
    kind           layer
    active         15N
    areal_density  {N_res:g} atoms/cm2 fixed
    stopping       formula "a0" a0={EPS_CONST:g} MeV*cm2/atom
    straggling     none
target TD
    kind           layer
    active         15N
    areal_density  {N_res:g} atoms/cm2 fit prior {N_res:g} 10%
    stopping       formula "a0" a0={EPS_CONST:g} MeV*cm2/atom
    straggling     none
beam B03
    spread         0.3 keV lab
</targets>
"""
    head = ("<data>\n# Synthetic: one data set per limit of the yield of a layer (README.md).  The data values are\n"
            "# placeholders; check.py compares the model values with closed forms.\n")
    with open(os.path.join(MODELS, "layer_yield.azr"), "w") as out:
        out.write(CONFIG + "<levels>\n" + levels() + "</levels>\n" + head + data + "</data>\n" + targets)

    # the fit: a scan of the resonance on a layer of N_true = 1.2e18 (Delta = 6 keV), every 0.25 keV, 2 % noise
    rng = random.Random(20261111)
    N_true = 1.2e18
    E_fit = [round(Er_lab - 0.002 + 0.00025 * i, 10) for i in range(41)]
    Y_true = [bw_layer_yield(e, N_true * EPS_CONST) for e in E_fit]
    err = [0.02 * max(y, 0.02 * max(Y_true)) for y in Y_true]
    Y_obs = [y + rng.gauss(0.0, s) for y, s in zip(Y_true, err)]
    write_data("fit.dat", E_fit, Y_obs, err)
    truth["fit"] = dict(N_true=N_true, N_start=1.0e18, Gp_true=LEVEL["Gp_eV"], Gp_start=45.0, seed=20261111, E0=E_fit)
    fit_levels = levels(gp_fixed=0).replace(f" {LEVEL['Gp_eV']:.10g} ", f" {45.0:.10g} ", 1)
    fit = (CONFIG + "<levels>\n" + fit_levels + "</levels>\n<data>\n# Synthetic resonance scan (make_project.py): "
           "the arctan form at N = 1.2e18, Gamma_p = 60 eV, times 1 + 2 % noise (seed 20261111)\n"
           + dataset("scan", "15N(p,a)12C", "yield", "1/particle", "fit.dat", "TF", points=40) + "</data>\n"
           f"""<targets>
target TF
    kind           layer
    active         15N
    areal_density  1e18 atoms/cm2 fit
    stopping       formula "a0" a0={EPS_CONST:g} MeV*cm2/atom
    straggling     none
</targets>
""")
    with open(os.path.join(MODELS, "fit.azr"), "w") as out:
        out.write(fit)
    with open(os.path.join(MODELS, "truth.json"), "w") as out:
        json.dump(truth, out, indent=1)
        out.write("\n")


if __name__ == "__main__":
    main()
