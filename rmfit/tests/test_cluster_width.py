#!/usr/bin/env python3
"""cluster_width / auto_gamma_kev (2026-09-24, the evaluator).

The level search never estimated a width from the data: the template screen
sampled Gamma x [0.3..10] around a guessed Gamma (5 keV in one narrow-level
campaign, so 1.5-50 keV) and sat at the top of that grid 95 times out of 120,
while the released fits in the same window came out at 0.6-2 MeV.  The
residual runs already carry the extent of the structure; it is now reported
per cluster (width = chi2-weighted mean run extent, span = full extent) and
``--gamma auto`` centres the template grid on it.  Pure functions, no engine.
Run:  python3 rmfit/tests/test_cluster_width.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.levelscan import cluster_width, auto_gamma_kev, ExCandidate  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def run(e_lo, e_hi, chi2, seg=1, sign=1, diff=True, ne=None, cover=0.2):
    return dict(segment=seg, e_lo=e_lo, e_hi=e_hi, n=5, sign=sign, zmax=4.0, chi2=chi2, centre=0.5 * (e_lo + e_hi),
                diff=diff, ne=(ne if ne is not None else (5 if e_hi > e_lo else 1)), cover=cover)


def main():
    # 1. one differential run: width = span = its extent
    cw = cluster_width([run(11.60, 12.20, 100.0)])
    check("single run: width = extent", abs(cw["width"] - 0.6) < 1e-9, cw)
    check("single run: span = extent", abs(cw["span"] - 0.6) < 1e-9, cw)
    check("single run: no total width", cw["width_total"] == 0.0, cw)

    # 2. the real first-try failure: the n-total's one heavy 560 keV run must NOT own the
    #    width -- it is a total segment and goes to width_total; the differential runs decide
    runs = [run(10.901, 11.460, 255624.0, seg=93, diff=False),
            run(11.023, 11.084, 6597.0, seg=112), run(11.084, 11.130, 4414.0, seg=112),
            run(10.915, 11.097, 4131.0, seg=51), run(11.008, 11.191, 3527.0, seg=111),
            run(11.073, 11.090, 2547.0, seg=97), run(11.022, 11.097, 2284.0, seg=77)]
    cw = cluster_width(runs)
    check("n-total run does not set width", cw["width"] < 0.2, cw)
    check("n-total run sets width_total", abs(cw["width_total"] - 0.559) < 1e-9, cw)
    check("span still covers everything", abs(cw["span"] - 0.559) < 1e-9, cw)
    check("n_runs counted", cw["n_runs"] == 7, cw)

    # 3. chi2-weighted MEDIAN, not mean: two heavy 60 keV runs and one light 600 keV run -> 60 keV
    cw = cluster_width([run(11.60, 11.66, 500.0), run(11.70, 11.76, 500.0, seg=2), run(11.0, 11.6, 10.0, seg=3)])
    check("weighted median resists one wide light run", abs(cw["width"] - 0.06) < 1e-9, cw)
    cw = cluster_width([run(11.60, 11.66, 10.0), run(11.0, 11.6, 500.0, seg=3)])
    check("weighted median follows the heavy run", abs(cw["width"] - 0.6) < 1e-9, cw)

    # 4. angular distributions at ONE energy (many angles, zero extent) do not count
    cw = cluster_width([run(12.474, 12.474, 7018.0, seg=128, ne=1), run(12.323, 12.476, 4398.0, seg=112)])
    check("zero-extent single-energy run ignored", abs(cw["width"] - 0.153) < 1e-9, cw)
    cw = cluster_width([run(12.40, 12.50, 100.0, ne=2)])
    check("run over < 3 distinct energies ignored -> width 0", cw["width"] == 0.0, cw)

    # 4b. a run covering (nearly) its whole segment is a normalization deficit, not a shape:
    #     the Ge1 (a,a1 gamma) run 11.848-13.236 (chi2 31470) at the 12.863 cluster
    cw = cluster_width([run(11.848, 13.236, 31470.0, seg=117, cover=0.98), run(12.323, 12.476, 4398.0, seg=112)])
    check("whole-segment run ignored for the width", abs(cw["width"] - 0.153) < 1e-9, cw)
    check("...but still inside the span", abs(cw["span"] - 1.388) < 1e-9, cw)
    cw = cluster_width([run(11.848, 13.236, 31470.0, seg=117, cover=0.98)])
    check("only a whole-segment run -> width 0", cw["width"] == 0.0, cw)

    # 5b. only total runs qualify -> width falls back to them
    cw = cluster_width([run(11.0, 11.4, 100.0, diff=False)])
    check("fallback to total runs when no differential run qualifies", abs(cw["width"] - 0.4) < 1e-9, cw)

    # 6. a dispersion pair: two adjacent lobes -> span is the whole shape, width one lobe
    cw = cluster_width([run(11.50, 11.80, 500.0, sign=1), run(11.80, 12.10, 500.0, sign=-1)])
    check("dispersion pair: width = one lobe", abs(cw["width"] - 0.3) < 1e-9, cw)
    check("dispersion pair: span = both lobes", abs(cw["span"] - 0.6) < 1e-9, cw)

    # 7. old-style runs (no diff/ne keys) and empty input are harmless
    cw = cluster_width([dict(segment=1, e_lo=11.6, e_hi=11.9, chi2=5.0)])
    check("runs without diff/ne count as differential", abs(cw["width"] - 0.3) < 1e-9, cw)
    check("no runs -> zeros", cluster_width([]) == dict(width=0.0, span=0.0, width_total=0.0, n_runs=0))

    # 5. auto_gamma_kev: clamp into bounds, default when no width
    check("auto: 600 keV run -> 600 keV template", auto_gamma_kev(dict(width=0.6)) == 600.0)
    check("auto: 3 MeV run clamps to 1000 keV", auto_gamma_kev(dict(width=3.0)) == 1000.0)
    check("auto: 2 keV run clamps to 10 keV", auto_gamma_kev(dict(width=0.002)) == 10.0)
    check("auto: custom bounds", auto_gamma_kev(dict(width=3.0), bounds_kev=(50, 2000)) == 2000.0)
    check("auto: --ex candidate (no width) -> default", auto_gamma_kev(dict(ex=11.6, score=float("nan"))) == 10.0)
    check("auto: width 0 -> default", auto_gamma_kev(dict(width=0.0), default_kev=25.0) == 25.0)

    # 6. the candidate record carries both numbers, and old records without them still load
    c = ExCandidate(ex=11.7, score=1.0, npts=5, segments=[1], runs=[], width=0.4, span=0.7, width_total=0.9)
    j = c.to_json()
    check("to_json carries width/span/width_total", j["width"] == 0.4 and j["span"] == 0.7 and j["width_total"] == 0.9, j)
    old = ExCandidate(ex=11.7, score=1.0, npts=5, segments=[1])
    check("width/span default to 0 on old-style candidates", old.width == 0.0 and old.span == 0.0)

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all cluster-width checks passed")


if __name__ == "__main__":
    main()
