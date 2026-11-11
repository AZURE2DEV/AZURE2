#!/usr/bin/env python3
"""Writes the target_density models from tests/13N/13N.azr (see README.md).

12C(p,p) at 84.3 and 114.5 deg lab (Meyer, segments 3 and 4) from 0.3 to 0.9 MeV, every level fixed, every other
segment inactive, a thin target with a constant stopping cross section of 1e-20 MeV cm^2/atom (10 keV at 1e18).

    python3 tests/target_density/make_models.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "13N", "13N.azr")
STOP = '"a0" 1 1e-20'


def model(target, data_dir="data", second=False):
    lines = open(SRC).read().split("\n")
    l0, l1 = lines.index("<levels>"), lines.index("</levels>")
    for i in range(l0 + 1, l1):
        f = lines[i].split()
        if len(f) > 20:                                  # fix every level energy and width
            f[3], f[10] = "1", "1"
            lines[i] = " ".join(f)
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    for k, i in enumerate(range(i0 + 1, i1), start=1):
        f = lines[i].split()
        if not f:
            continue
        f[0] = "1" if k == 3 or (second and k == 4) else "0"
        if k in (3, 4):
            f[3], f[4] = "0.3", "0.9"
            f[10] = "0"                                  # no norm error: the norm stays fixed at 1
            f[14] = f[14].replace("data/", data_dir + "/")
        lines[i] = " ".join(f)
    j0, j1 = lines.index("<targetInt>"), lines.index("</targetInt>")
    lines[j0 + 1:j1] = target
    return "\n".join(lines)


def thin(seg, n, tail=""):
    return f'1 "{seg}" 50 0 0 1 {n} {STOP} 0 0 0 "" 0 0' + (f" {tail}" if tail else "")


models = {
    # the truth: N = 1.1e18 on both segments (data are made from it by check.sh)
    "true.azr": model([thin(3, "1.1e+18"), thin(4, "1.1e+18")], second=True),
    # one segment, density free from 1e18 without a prior (grids to 1.5e18)
    "fit.azr": model([thin(3, "1e+18", "densityfit 1 0")], data_dir="fitdata"),
    # the same with a 2 % prior about 1e18
    "fit_prior.azr": model([thin(3, "1e+18", "densityfit 1 2%")], data_dir="fitdata"),
    # both segments, one target: a density group (the second line's own 2e18 is replaced by the group's)
    "fit_group.azr": model([thin(3, "1e+18", "densityfit 1 0 1.5e18 densitygroup t"),
                            thin(4, "2e+18", "densitygroup t")], data_dir="fitdata", second=True),
    # refused: densityfit on a line without target integration
    "bad.azr": model(['1 "3" 50 1 0.002 0 0 "" 0 0 0 0 "" 0 densityfit 1 5%']),
}
for name, text in models.items():
    open(os.path.join(HERE, name), "w").write(text)
    print("wrote", name)
