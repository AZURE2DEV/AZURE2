#!/usr/bin/env python3
"""Self-test of the azure2-check scripts: chi2_report.py, widths_check.py, radius_scan.py.

Runs each script on a copy of the synthetic project T1 (tests/upgrade/T1: p + 12C -> 13N, elastic scattering, A_y,
a phase shift and capture, generated from known parameters) in a temporary directory, and checks numbers that mean
something:

- chi2_report.py reproduces the Total-Chi-Squared of the program's own calculation (`calculate`) to 1e-6, and each data
  set's chi2 is the sum of its squared pulls (T1 has no data covariance);
- widths_check.py gives gamma^2_W = 3 hbar^2/(2 mu a_c^2) as computed here from the pair's masses and radius
  (CODATA 2018); the Brune norm J of the 1/2+ and 3/2- levels equals 1 - Gamma S'(E)/(2 P(E)) from the engine's
  physics functions (0.2553: the 1/2+ proton width of T1, 33 keV, is three quarters of the most a level there can
  have, 2P/S' = 44.3 keV, so it is the one flag of T1); a fixed 60 MeV width on the 1/2+ background pole is flagged
  theta^2 > 1; on T2 (tests/upgrade/T2: levels with fixed and free widths in several open channels) the norm of
  each level with a free width equals J = 1 - sum_c Gamma_c S'_c/(2 P_c) over all its channels (0.797, 0.745,
  0.972), and Gamma_W of a free channel is 2 P gamma^2_W;
- radius_scan.py, from the perturbed start of T1 at its own radius, refits to the truth's chi2 neighbourhood
  (objective lowered by the refit, the 1/2+ energy back to 2.365 MeV within 2 keV), and a second radius is fitted
  from the first's result.

Needs: pyazr importable with its engine (PYAZR_ENGINE_DIR), AZURE2_BIN (the program) and AZURE2_COULOMB_TABLE as the
pyazr tests have them. Prints PASS/FAIL lines; exit 0 when all pass, 1 on a failure, 77 (skip) without the engine.
Run time about 30 s at one thread.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
T1 = REPO / "tests" / "upgrade" / "T1"
FAILED = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


def run(script, *args, cwd):
    t = time.time()
    r = subprocess.run([sys.executable, str(HERE / script), *map(str, args)], cwd=cwd, capture_output=True, text=True)
    print(f"  ran {script} {' '.join(map(str, args))} in {time.time() - t:.1f} s (exit {r.returncode})")
    if r.returncode not in (0, 1):
        print(r.stdout[-2000:], r.stderr[-2000:])
    return r


def copy_t1(dest):
    shutil.copytree(T1 / "data", dest / "data")
    for f in ("T1.azr", "T1_start.azr"):
        shutil.copy(T1 / f, dest / f)
    for d in ("output", "checks"):
        (dest / d).mkdir()


def test_chi2_report(work):
    azure2 = os.environ.get("AZURE2_BIN")
    if not azure2:
        check("chi2_report: AZURE2_BIN set", False, "the comparison with the program needs it")
        return
    calc = subprocess.run([azure2, "calculate", "--no-gui", "T1.azr"], stdin=subprocess.DEVNULL, cwd=work,
                          capture_output=True, text=True)
    check("program calculates T1", calc.returncode == 0, f"exit {calc.returncode}")
    r = run("chi2_report.py", "T1.azr", "--chisquared-out", "output/chiSquared.out", "--json", "rep.json", cwd=work)
    check("chi2_report: data chi2 = chiSquared.out total (1e-6)", r.returncode == 0, r.stdout.splitlines()[-1])
    rep = json.loads((work / "rep.json").read_text())
    check("chi2_report: every data set's chi2 = sum of its squared pulls",
          not any(x["covariance"] for x in rep["rows"]))
    check("chi2_report: chi2/nu uses the free parameters", rep["n_free"] > 0 and rep["chi2_per_dof"] is not None
          and abs(rep["chi2_per_dof"] - rep["chi2"] / (rep["n_points"] - rep["n_free"])) < 1e-12)


def test_widths(work):
    import pyazr
    from pyazr import physics
    r = run("widths_check.py", "T1.azr", "--json", "w.json", cwd=work)
    w = json.loads((work / "w.json").read_text())
    flags = [l for l in r.stdout.splitlines() if l.startswith("FLAG")]
    check("widths_check: on the T1 truth the only flag is the saturated 1/2+ level", r.returncode == 1
          and len(flags) == 1 and "1/2+ 2.3650" in flags[0], "; ".join(flags))
    pair = next(q for q in pyazr.Project.open(work / "T1.azr").pairs() if q.key == 1)
    hbarc, amu = 197.3269804, 931.49410242  # MeV fm, MeV/u (CODATA 2018)
    mu = pair.light_mass * pair.heavy_mass / (pair.light_mass + pair.heavy_mass) * amu
    g2w = 3 * hbarc ** 2 / (2 * mu * pair.radius ** 2)
    rows = [c for c in w["channels"] if not c["photon"] and c["pair"] == 1 and c["wigner_gamma2"]]
    worst = max(abs(c["wigner_gamma2"] / g2w - 1) for c in rows)
    check("widths_check: gamma^2_W = 3 hbar^2/(2 mu a^2)", worst < 1e-6, f"{g2w:.6f} MeV, worst rel. diff {worst:.1e}")
    # the Brune norm of a one-channel level, J = 1 - Gamma S'(E)/(2 P(E)), from the engine's physics functions
    args = (1, 6, pair.light_mass, pair.heavy_mass, pair.radius)
    for jpi, ex, l, gamma in (("1/2+", 2.365, 0, 33000.0), ("3/2-", 3.502, 1, 55000.0)):
        e = ex - pair.separation_energy
        j_ref = 1 - gamma * 1e-6 * float(physics.shift_derivative(l, *args, e)) / (
            2 * float(physics.penetrability(l, *args, e)))
        lv = next(lv for lv in w["levels"] if lv["jpi"] == jpi and abs(lv["level_energy"] - ex) < 1e-3)
        check(f"widths_check: J of {jpi} {ex} = 1 - Gamma S'/(2P)", abs(lv["norm"] - j_ref) < 1e-5,
              f"{lv['norm']:.6f} vs {j_ref:.6f}")
    p = pyazr.Project.open(work / "T1.azr")
    lv = p.level(("1/2+", 15.0))
    c = next(c for c in lv.channels if c.pair == 1)
    p.set_channel(lv, c.pair, c.l, c.s, width=60e6)
    p.save(work / "T1_wide.azr")
    r = run("widths_check.py", "T1_wide.azr", "--background-above", "10", cwd=work)
    flags = [l for l in r.stdout.splitlines() if l.startswith("FLAG")]
    check("widths_check: a fixed 60 MeV proton width of the 1/2+ pole at 15 MeV is flagged theta^2 > 1",
          any("theta^2" in f for f in flags), "; ".join(flags) or "no flag")


def test_widths_fixed(work):
    """A level with fixed and free widths: theta^2 / theta^2_formal is the level's whole norm J."""
    import pyazr
    from pyazr import physics
    t2 = work / "T2"
    shutil.copytree(REPO / "tests" / "upgrade" / "T2", t2, ignore=shutil.ignore_patterns("output", "checks"))
    for d in ("output", "checks"):
        (t2 / d).mkdir()
    r = run("widths_check.py", "T2.azr", "--json", "w.json", cwd=t2)
    w = json.loads((t2 / "w.json").read_text())
    with pyazr.Session(str(t2 / "T2.azr")) as s:
        pairs = {q.number: q for q in s.pairs}
    worst, worst_width, n = 0.0, 0.0, 0
    for lv in w["levels"]:
        if lv["norm"] is None:
            continue
        j_ref = 1.0
        for c in w["channels"]:
            if c["photon"] or c["jpi"] != lv["jpi"] or abs(c["level_energy"] - lv["level_energy"]) > 1e-9:
                continue
            q = pairs[c["pair"]]
            args = (c["L"], q.Z1, q.Z2, q.M1, q.M2, q.channel_radius, c["level_energy"] - q.sep_energy - q.excitation)
            two_p = 2 * float(physics.penetrability(*args))
            j_ref -= abs(c["value"]) * 1e-6 * float(physics.shift_derivative(*args)) / two_p
            if c["free"] and c["wigner_width"]:
                worst_width = max(worst_width, abs(c["wigner_width"] / (two_p * c["wigner_gamma2"] * 1e6) - 1))
        n += 1
        worst = max(worst, abs(lv["norm"] - j_ref))
    check("widths_check: on T2 the norm of a level with fixed and free widths is J = 1 - sum Gamma S'/(2P)",
          r.returncode in (0, 1) and n == 3 and worst < 1e-5, f"{n} levels, largest difference {worst:.1e}")
    check("widths_check: on T2 Gamma_W of a free channel is 2 P gamma^2_W", worst_width < 1e-6, f"{worst_width:.1e}")


def test_radius_scan(work):
    import pyazr
    r = run("radius_scan.py", "T1_start.azr", "--pair", "1", "--radii", "3.6", "3.9", "--table", "scan.txt",
            cwd=work)
    check("radius_scan: both radii fitted", r.returncode == 0, r.stdout.strip().splitlines()[-1])
    rows = {float(t[0]): [float(v) for v in t[1:4]] for t in
            (l.split() for l in (work / "scan.txt").read_text().splitlines() if not l.startswith("#"))}
    f0, f, chi2 = rows[3.6]
    check("radius_scan: the refit lowers the objective at the project's radius", f < 0.5 * f0,
          f"{f0:.2f} -> {f:.2f}")
    with pyazr.Session(str(work / "T1_start_a1_3.6.azr")) as s:
        x = s.parameters.read(work / "T1_start_a1_3.6.sav")
        e = next(v for p, v in zip(s.parameters.free, x) if p.kind == "energy" and abs(v - 2.365) < 0.05)
        again = float(s.objective()(x))
    check("radius_scan: the 1/2+ energy returns to the truth (2.365 MeV, 2 keV)", abs(e - 2.365) < 2e-3, f"{e:.5f}")
    check("radius_scan: the saved vector gives the reported objective", abs(again - f) < 1e-6 * max(f, 1),
          f"{again:.6f}")
    check("radius_scan: the second radius is refitted from the first", rows[3.9][1] < rows[3.9][0],
          f"{rows[3.9][0]:.2f} -> {rows[3.9][1]:.2f}")


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr  # noqa: F401
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_check_selftest_") as tmp:
        work = Path(tmp)
        copy_t1(work)
        for test in (test_chi2_report, test_widths, test_widths_fixed, test_radius_scan):
            try:
                test(work)
            except Exception as err:
                check(test.__name__, False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-check selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
