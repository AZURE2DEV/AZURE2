#!/usr/bin/env python3
"""Write models/bare.azr and models/rect.azr: the level scheme of tests/convolution_window (16O(p,g)17F, a narrow
1/2- resonance near 0.84 MeV lab on a smooth background), unfolded on a fine lab grid (bare) and folded with a
rectangular energy resolution of full width 20 keV lab at points E_res + k w (rect)."""
import os

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
write("rect.azr", body, f"detector det_rect\n    energy_resolution rectangular width {WIDTH_KEV:g} keV\n")
