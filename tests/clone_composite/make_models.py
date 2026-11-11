#!/usr/bin/env python3
"""Models for tests/clone_composite, built at run time (check.sh) into a work directory.

From tests/target_density (format 3; 12C(p,p) at 84.3 deg lab, data set meyer_84, every level fixed, a thin target
with a fitted areal density), with meyer_84 made a sum of its own elastic term and one identical elastic component
(2 sigma, "12C(p,p)12C + 12C(p,p)12C"), the only active data set, so the fitted density acts on a data set that holds
components:

  true.azr   N = 1.1e18 (check.sh makes the fit data from its output)
  fit.azr    N free from 1e18, data from true.azr

    python3 make_models.py <tests/target_density> <workdir>
"""
import sys

src, work = sys.argv[1], sys.argv[2]
SUM = "12C(p,p)12C + 12C(p,p)12C"


def sum_meyer_84(text, data_dir):
    """Every data set but meyer_84 inactive; meyer_84 active, the sum, its data file in data_dir."""
    out, name, found = [], None, 0
    for line in text.split("\n"):
        f = line.split()
        if line.startswith("dataset "):
            name = f[1]
            out.append(line)
            out.append("    active         " + ("yes" if name == "meyer_84" else "no"))
            continue
        if not f or line.startswith("<"):
            name = None
        if name and f[0] == "active":
            continue  # written after the header
        if name == "meyer_84" and f[0] == "reaction":
            if f[1:] != ["12C(p,p)12C"]:
                sys.exit("make_models: meyer_84 is not a plain elastic data set")
            line = "    reaction       " + SUM
            found += 1
        elif name == "meyer_84" and f[0] == "file":
            line = "    file           " + data_dir + "/" + f[1].split("/")[-1]
        out.append(line)
    if found != 1:
        sys.exit("make_models: data set meyer_84 not found")
    return "\n".join(out)


open(work + "/true.azr", "w").write(sum_meyer_84(open(src + "/true.azr").read(), "data"))
open(work + "/fit.azr", "w").write(sum_meyer_84(open(src + "/fit.azr").read(), "fitdata"))
