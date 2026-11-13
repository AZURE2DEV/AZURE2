#!/usr/bin/env python3
"""Self-test of the azure2-targets scripts: layer.py, target_check.py.

- The energy loss layer.py prints for 2e18 atoms/cm^2 of TiN and 0.4 MeV protons satisfies the range integral
  int_{E_out}^{E0} dE/eps = N, integrated here with scipy's quad on the library's stopping (to 1e-5);
- a quoted loss (5 keV) gives an areal density whose loss is 5 keV again (round trip, 1e-6 keV);
- the first-order eps(E0) N is within 3 % of the range integral for this thin layer and below it (eps rises toward the
  stopping maximum as the protons slow down at 0.4 MeV);
- the target record layer.py prints for a carbon foil is accepted by the engine: added to a copy of T1
  (``Project.add``, p + 12C) and named by a data set, the project validates without an error.

- target_check.py: on that project, with the foil on the 12C(p,g) data set through the 1/2+ resonance
  (Gamma = 33 keV, a loss of 8 keV), the engine's thin-target fold equals the integral computed outside the engine
  at every point (2e-3; measured 1e-4).

Needs pyazr with its engine (PYAZR_ENGINE_DIR). PASS/FAIL lines; exit 0, 1 on a failure, 77 (skip) without the
engine. About 10 s.
"""
import argparse
import re
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
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
    return r


def number(pattern, text):
    m = re.search(pattern, text)
    return float(m.group(1)) if m else float("nan")


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr
        from pyazr import stopping
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    from scipy.integrate import quad
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_targets_selftest_") as tmp:
        work = Path(tmp)
        try:
            n, e0 = 2e18, 0.4
            r = run("layer.py", "--material", "TiN", "--projectile", "p", "--energy", e0, "--areal-density", n,
                    "--active", "15N", cwd=work)
            loss = number(r"range integral\) (\S+) keV", r.stdout) * 1e-3
            first = number(r"first order eps\(E0\) N (\S+)", r.stdout) * 1e-3
            m = stopping.Material("TiN")
            inv = lambda e: 1e21 / float(stopping.stopping(m, "p", e)["total"])
            n_back = quad(inv, e0 - loss, e0, epsabs=0, epsrel=1e-10)[0]
            check("layer: the loss satisfies int dE/eps = N", abs(n_back / n - 1) < 1e-5,
                  f"loss {loss * 1e3:.4f} keV, N back {n_back:.6e}")
            check("layer: first order eps(E0) N is within 3 % and below", 0.97 < first / loss < 1.0,
                  f"ratio {first / loss:.4f}")
            r = run("layer.py", "--material", "TiN", "--projectile", "p", "--energy", e0, "--loss-kev", 5, cwd=work)
            check("layer: a quoted 5 keV loss round-trips", abs(number(r"range integral\) (\S+) keV", r.stdout) - 5)
                  < 1e-6, r.stdout.splitlines()[3])

            r = run("layer.py", "--material", "C", "--projectile", "p", "--energy", e0, "--areal-density", 1e18,
                    "--active", "12C", cwd=work)
            fields = dict(re.findall(r"^    (\w+)\s+(.+)$", r.stdout, re.M))
            shutil.copytree(T1 / "data", work / "data")
            p = pyazr.Project.open(T1 / "T1.azr")
            p.add("target", "tin", **fields)
            name = "cap_g0"  # 12C(p,g)13N, angle-integrated
            p.set(name, target="tin", quadrature="points 20")  # a folded data set states its quadrature
            path = p.save(work / "T1_target.azr")
            errors = [d for d in pyazr.validate(path) if d["severity"] == "error"]
            check("layer: the printed target record validates in a project", not errors,
                  "; ".join(f"{d['code']}: {d['message']}" for d in errors) or f"{len(fields)} fields")
            # the engine's thin-target fold of the capture data set through the 1/2+ resonance against the integral
            p.set(name, quadrature="points 50")
            p.set("tin", straggling="none")
            for d in ("output", "checks"):
                (work / d).mkdir(exist_ok=True)
            path = p.save(work / "T1_target.azr")
            r = run("target_check.py", "T1_target.azr", "--dataset", name, cwd=work)
            ratios = [float(l.split()[3]) for l in r.stdout.splitlines() if re.match(r"^ +\d", l)]
            check("target_check: the engine's fold equals the independent integral (2e-3)", r.returncode == 0
                  and len(ratios) > 10, f"{len(ratios)} points, ratios {min(ratios, default=0):.6f} to "
                  f"{max(ratios, default=0):.6f}")
        except Exception as err:
            check("layer", False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-targets selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
