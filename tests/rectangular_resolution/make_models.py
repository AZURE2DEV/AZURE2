#!/usr/bin/env python3
"""Write models/bare.azr and models/rect.azr: the level scheme of tests/convolution_window (16O(p,g)17F, a narrow
1/2- resonance near 0.84 MeV lab on a smooth background), unfolded on a fine lab grid (bare) and folded with a
rectangular energy resolution of full width 20 keV lab at points E_res + k w (rect); and models/fit.azr: two data
sets folded with one detector whose width is a fit parameter (start 14 keV).

  make_models.py --data <bare AZUREOut_aa=1_R=3.extrap> [<directory> <noise>]
writes models/data/d1.dat, d2.dat (or <directory>/d1.dat, d2.dat): the pseudo-data at the true width W_TRUE_KEV, each
point the window average of the bare cross section (trapezoid on the bare nodes, no engine kernel) times 1 + <noise>
(default 0.01) times a Gaussian deviate (seeded), error 1 %."""
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "..", "convolution_window", "models", "bare.azr")
E_RES = 0.8508  # MeV lab, the peak of the resonance (c.m. 0.80038 MeV)
WIDTH_KEV = 20.0

head = open(SOURCE).read().split("<data>")[0]


def extrapolation(name, energies, folded):
    lines = [f"extrapolation {name}", "    reaction       16O(p,g)17F[0.4953]", "    quantity       sigma",
             "    unit           b", "    energies       list " + " ".join(f"{e:.6f}" for e in energies) + " MeV lab"]
    if folded:
        lines += ["    detector       det_rect", "    quadrature     points 50 width_multiplier 20 per_width 50"]
    return "\n".join(lines) + "\n"


def write(name, body, targets):
    with open(os.path.join(HERE, "models", name), "w") as f:
        f.write(head + "<data>\n</data>\n<extrapolations>\n" + body + "</extrapolations>\n<targets>\n" + targets +
                "</targets>\n")


fine = [0.79 + 0.00005 * i for i in range(2001)]  # 0.79 .. 0.89 MeV lab, 0.05 keV steps
write("bare.azr", extrapolation("fine", fine, False), "")
w = WIDTH_KEV * 1e-3
points = [E_RES + k * w for k in (-1.0, -0.6, -0.3, 0.0, 0.25, 0.5, 1.0)]
body = "".join(extrapolation(f"r{i + 1}", [e], True) + "\n" for i, e in enumerate(points))
write("rect.azr", body, f"detector det_rect\n    energy_resolution rectangular width {WIDTH_KEV:g} keV fixed\n")

# ---- the fit of the width (check.sh part c)
W_TRUE_KEV, W_START_KEV = 20.0, 14.0
SETS = {"d1": [0.805 + 0.003 * i for i in range(25)], "d2": [0.8065 + 0.003 * i for i in range(25)]}


def dataset(name, points):
    # the two data sets differ in their quadrature, so they are two effects that share the detector's one width
    return "\n".join([f"dataset {name}", "    reaction       16O(p,g)17F[0.4953]", "    quantity       sigma",
                      "    unit           b", f"    file           data/{name}.dat", "    columns        E_lab[MeV] _ value error",
                      "    norm           1 fixed", "    detector       det_fit",
                      f"    quadrature     points {points} width_multiplier 20 per_width 50"]) + "\n"


with open(os.path.join(HERE, "models", "fit.azr"), "w") as f:
    f.write(head + "<data>\n# Pseudo-data at a 20 keV rectangular resolution "
            "(make_models.py --data); the width starts at 14 keV and is fitted\n" + dataset("d1", 50) +
            dataset("d2", 60) + "</data>\n<extrapolations>\n</extrapolations>\n<targets>\ndetector det_fit\n"
            f"    energy_resolution rectangular width {W_START_KEV:g} keV fit\n</targets>\n")


def window_average(e, s, lo, hi):  # trapezoid on the nodes inside [lo, hi] plus the two edges, linear between nodes
    def at(x):
        for i in range(len(e) - 1):
            if e[i] <= x <= e[i + 1]:
                return s[i] + (x - e[i]) / (e[i + 1] - e[i]) * (s[i + 1] - s[i])
        raise ValueError(x)
    xs = [lo] + [x for x in e if lo < x < hi] + [hi]
    ys = [at(x) for x in xs]
    return sum(0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1)) / (hi - lo)


if len(sys.argv) in (3, 5) and sys.argv[1] == "--data":
    where = sys.argv[3] if len(sys.argv) == 5 else os.path.join(HERE, "models", "data")
    noise = float(sys.argv[4]) if len(sys.argv) == 5 else 0.01
    rows = [list(map(float, l.split())) for l in open(sys.argv[2]) if len(l.split()) >= 4]
    e, s = [r[0] for r in rows], [r[3] for r in rows]
    f_cm = e[0] / fine[0]  # c.m./lab energy factor of the pair, from the bare run itself
    rng = random.Random(20261111)
    os.makedirs(where, exist_ok=True)
    for name, energies in SETS.items():
        with open(os.path.join(where, name + ".dat"), "w") as out:
            for lab in energies:
                cm, half = f_cm * lab, 0.5 * f_cm * W_TRUE_KEV * 1e-3
                y = window_average(e, s, cm - half, cm + half)
                out.write(f"{lab:.6f} 0 {y * (1 + noise * rng.gauss(0, 1)):.6e} {0.01 * y:.6e}\n")
