#!/usr/bin/env python3
"""Compare the phase shifts of AZUREOut_*.extrap (one block per extrapolation, in the record order of the models: l, 2J
= (0,1) (1,1) (1,3) (2,3) (2,5) (3,5) (3,7), energies 0.5 1 2 4 MeV c.m.; columns E_cm, E_x, angle, delta[deg]) with
data/phases.csv (make_phases.py).  Prints the largest difference (deg); exit 1 above TOL (deg).  Standard library only.

    phase_check.py EXTRAP PHASES.csv TOL            one run against the reference
    phase_check.py EXTRAP EXTRAP2 TOL --pair       two runs against each other
    phase_check.py EXTRAP PHASES.csv MIN --differs  exit 1 unless the s wave differs from the reference by >= MIN at some energy
"""
import sys

WAVES = [(0, 1), (1, 1), (1, 3), (2, 3), (2, 5), (3, 5), (3, 7)]


def extrap(path):
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.lstrip().startswith("#")]
    out = {}
    for i, r in enumerate(rows):
        out[(WAVES[i // 4], round(float(r[0]), 6))] = float(r[3])
    return out


def reference(path):
    out = {}
    for ln in open(path).read().split("\n")[1:]:
        if ln.strip():
            l, tj, e, d = ln.split(",")
            out[((int(l), int(tj)), round(float(e), 6))] = float(d)
    return out


def wrap(d):
    while d > 90.0:
        d -= 180.0
    while d <= -90.0:
        d += 180.0
    return d


def main():
    a = extrap(sys.argv[1])
    mode = sys.argv[4] if len(sys.argv) > 4 else ""
    b = extrap(sys.argv[2]) if mode == "--pair" else reference(sys.argv[2])
    tol = float(sys.argv[3])
    if len(a) != 28 or set(a) != set(b):
        print(f"expected 28 phases of the same waves and energies; got {len(a)}")
        return 1
    diffs = {k: abs(wrap(a[k] - b[k])) for k in a}
    if mode == "--differs":
        s = max(v for k, v in diffs.items() if k[0] == (0, 1))
        print(f"largest s-wave difference {s:.3e} deg")
        return 0 if s >= tol else 1
    worst = max(diffs.values())
    print(f"largest difference {worst:.2e} deg")
    return 0 if worst <= tol else 1


if __name__ == "__main__":
    sys.exit(main())
