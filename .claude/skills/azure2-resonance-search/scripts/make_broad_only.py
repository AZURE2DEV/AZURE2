#!/usr/bin/env python3
"""Broad-only copy of a fit, for the blur method: every level whose total width (summed over the channels of
parameters.out) is below GMAX keV and whose E_x is at least --ex-min gets its energy fixed and every channel's width
set to 0 and fixed, so the free parameters of the fit never touch it.

Reads:  the fitted project and the fit's parameters.out (the program's output: "J = ... E_level = ..." blocks with
        "G = ... keV" lines).
Writes: OUT (a project; write it in the fit's directory, the data paths stay relative).

    python3 make_broad_only.py fit.azr output/parameters.out broad.azr --gmax 50 --ex-min 4.2

--ex-min keeps bound and sub-threshold states (an ANC only, no width) out of it.
Exit status: 0 written, 1 a level of the project has no unique match in parameters.out (listed; it is kept).
"""
import argparse
import re
import sys


def read_par(path):
    """[(J, parity, E_x, total width in keV)] of parameters.out."""
    out, cur = [], None
    for line in open(path):
        m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
        if m:
            cur = [float(m.group(1)[:-1]), m.group(1)[-1], float(m.group(2)), 0.0]
            out.append(cur)
            continue
        m = re.search(r"G\s+=\s+(\S+) keV", line)
        if m and cur is not None:
            cur[3] += abs(float(m.group(1)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("parameters")
    ap.add_argument("out")
    ap.add_argument("--gmax", type=float, default=50.0, help="keV")
    ap.add_argument("--ex-min", type=float, default=0.0, help="MeV")
    a = ap.parse_args(argv)
    from pyazr import Project
    par = read_par(a.parameters)
    p = Project.open(a.project)
    n, unmatched = 0, 0
    for lv in p.levels():
        if not lv.active:
            continue
        sign = "+" if lv.parity > 0 else "-"
        c = [L for L in par if abs(L[2] - lv.energy) < 2e-3 and abs(L[0] - lv.J) < 1e-6 and L[1] == sign]
        if len(c) != 1:
            print(f"no unique match in {a.parameters}: {lv.label} ({len(c)})")
            unmatched += 1
            continue
        if c[0][3] < a.gmax and lv.energy >= a.ex_min:
            p.set_level(lv.label, fixed=True)
            for ch in lv.channels:
                p.set_channel(lv.label, ch.pair, ch.l, ch.s, width=0.0, fixed=True)
            n += 1
    p.save(a.out)
    print(f"{n} levels with Gamma < {a.gmax} keV zeroed and fixed -> {a.out}")
    return 1 if unmatched else 0


if __name__ == "__main__":
    sys.exit(main())
