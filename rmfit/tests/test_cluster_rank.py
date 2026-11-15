#!/usr/bin/env python3
"""cluster_density and the rank="density" order (2026-09-24, the evaluator).

The narrow campaign never tested a secondary-gamma feature: clusters were ordered
by raw chi2, and a 12-point Ge run cannot outscore a 791-point n-total run.
Density = chi2 per MeV of structure puts a sharp feature ahead of a smooth deficit
of similar chi2.  Pure functions, no engine.
Run:  python3 rmfit/tests/test_cluster_rank.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.levelscan import cluster_density, ExCandidate  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    # the measured clusters of 2026-09-24: n-total-dominated 11.26 (317,802 chi2, 75 keV
    # differential width) vs a gamma-only cluster (say 20,000 chi2 over 30 keV)
    check("density = chi2 / width", abs(cluster_density(20000, 0.03, 0.2) - 20000 / 0.03) < 1e-6)
    check("span used when width is 0", abs(cluster_density(1000, 0.0, 0.1) - 1000 / 0.1) < 1e-6)
    check("floor stops a zero-extent cluster", abs(cluster_density(1000, 0.0, 0.0, floor_kev=20) - 1000 / 0.02) < 1e-6)
    check("floor applies to a tiny width", abs(cluster_density(1000, 0.001, 0.5, floor_kev=20) - 1000 / 0.02) < 1e-6)

    broad = ExCandidate(ex=12.2, score=170000.0, npts=1000, segments=[93], width=0.6, span=1.4)
    sharp = ExCandidate(ex=12.47, score=40000.0, npts=40, segments=[128], width=0.03, span=0.05)
    for c in (broad, sharp):
        c.density = cluster_density(c.score, c.width, c.span)
    by_chi2 = sorted([broad, sharp], key=lambda c: -c.score)
    by_dens = sorted([broad, sharp], key=lambda c: -c.density)
    check("chi2 order puts the broad deficit first", by_chi2[0] is broad)
    check("density order puts the sharp feature first", by_dens[0] is sharp)
    check("to_json carries density", abs(sharp.to_json()["density"] - sharp.density) < 1e-9)
    check("old-style candidate density defaults to 0", ExCandidate(ex=1, score=1, npts=1, segments=[]).density == 0.0)

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all cluster-rank checks passed")


if __name__ == "__main__":
    main()
