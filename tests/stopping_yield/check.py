#!/usr/bin/env python3
"""The model values of a menu-1 run of models/stopping_yield.azr against LunaSRIM's numbers in models/truth.json
(README.md, eqs. (S1)-(S3)).  Prints one line per check; exit 1 when one is outside its tolerance.

    python3 check.py models/truth.json <run>/output
"""
import json
import math
import sys


def blocks(path):
    """The rows of an AZUREOut file, one list per data set (blank lines separate them), each row a list of floats."""
    out, cur = [], []
    with open(path) as f:
        for line in f:
            if line.strip():
                cur.append([float(x) for x in line.split()])
            elif cur:
                out.append(cur)
                cur = []
    if cur:
        out.append(cur)
    return out


def finite(r):
    return r if math.isfinite(r) else math.inf


def crossing(points, half, rising):
    """The energy where the yield crosses `half`, by linear interpolation between the two points that bracket it."""
    for (e1, y1), (e2, y2) in zip(points, points[1:]):
        if (y1 < half <= y2) if rising else (y1 > half >= y2):
            return e1 + (half - y1) * (e2 - e1) / (y2 - y1)
    return math.nan


def main():
    truth = json.load(open(sys.argv[1]))
    plateau, table, step = blocks(sys.argv[2] + "/AZUREOut_aa=1_R=2.out")
    failures = []

    def check(name, value, tol, text):
        ok = finite(value) <= tol
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}: {text}; {value:.2e} (tolerance {tol:g})")
        if not ok:
            failures.append(name)

    # (S1) the plateau of a thick layer: (lambda_r^2/2) omega gamma (m_p + m_N)/(m_N eps(E_r)), eps from LunaSRIM
    s = truth["sets"]["plateau"]
    worst = max(finite(abs(row[3] / y - 1.0)) for row, y in zip(plateau, s["Y"]))
    check("plateau", worst, 3e-4, f"thick TiN layer, eps(E_r) = {truth['eps_active_Er']:.6e} MeV cm^2 per 15N (LunaSRIM); "
          f"largest |Y/Y_SRIM - 1| of {len(plateau)} energies")
    # (S3) the same layer with the stopping as a table of the library's own values
    worst = max(finite(abs(t[3] / p[3] - 1.0)) for t, p in zip(table, plateau))
    check("table", worst, 1e-8, "stopping table of the library's values against stopping evaluation, largest |Y_table/Y_eval - 1|")
    # (S2) the step of a layer: the distance between its half-height points is E_in - E_out of LunaSRIM
    t = truth["sets"]["step"]
    def half_points(pts):
        mid = [y for e, y in pts if t["E_front"] + 0.25 * t["width"] < e < t["E_back"] - 0.25 * t["width"]]
        half = 0.5 * sum(mid) / len(mid)
        return (crossing([p for p in pts if p[0] < t["E_front"] + 0.5 * t["width"]], half, True),
                crossing([p for p in pts if p[0] > t["E_front"] + 0.5 * t["width"]], half, False))
    front, back = half_points([(e, row[3]) for e, row in zip(t["E0"], step)])
    bw_front, bw_back = half_points(list(zip(t["E0"], t["Y"])))
    width, bw_width = back - front, bw_back - bw_front
    check("step width", abs(width / bw_width - 1.0), 2e-5,
          f"half-height points {front:.8f} and {back:.8f} MeV, width {1e3 * width:.6f} keV; the Breit-Wigner scan with "
          f"LunaSRIM's E_out {1e3 * bw_width:.6f} keV (E_in - E_out at E_r: {1e3 * t['width']:.6f} keV)")
    if failures:
        print("check.py: failed: " + ", ".join(failures))
        sys.exit(1)


if __name__ == "__main__":
    main()
