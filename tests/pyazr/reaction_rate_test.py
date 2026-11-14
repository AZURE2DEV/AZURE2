#!/usr/bin/env python3
"""The reaction rate through pyazr (pyazr/rate.py, the binding's Session.reaction_rate) against the program.

  1. 12C(p,g) (tests/reaction_rate/models/p12C.azr): pyazr.rate on the program's grid (0.01 to 1 GK, step 0.07) and
     on a temperature list gives every rate the program's `rate` command writes, digit for digit (reactionrates.out
     holds six significant digits), in the same order; T1 (12C(p,g0)) and T3 (7Li(p,a)a) give their golden
     reactionrates.out (tests/upgrade/golden);
  2. a parameter file (tests/reaction_rate: the narrow level moved by moved.par) gives the program's rates with
     --parameter-file;
  3. exit -1 (total radiative capture) equals the rate into the one capture pair to the last bit;
  4. the band from draws equals a direct loop over the same drawn vectors (rates of a fresh session, percentiles);
     the band from a covariance equals g^T C g with g by central differences in a direct loop;
  5. the Gamow window against the closed form E0 = 0.1220 (Z1^2 Z2^2 mu T9^2)^(1/3) MeV,
     width = 0.2368 (Z1^2 Z2^2 mu T9^5)^(1/6) MeV (Rolfs and Rodney; three digits); NaN for a neutral entrance;
  6. the CSV table names the project, the parameters and the units; a comparison table is read back;
  7. the example pyazr/examples/reaction_rate.py runs on T1 and prints the golden rate.

Run from anywhere:  python3 tests/pyazr/reaction_rate_test.py  (needs AZURE2_BIN, the program)
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
MODELS = os.path.join(ROOT, "tests", "reaction_rate", "models")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    from pyazr import rate
    from pyazr import _azure2  # noqa: F401
except Exception as err:
    print(f"skip: pyazr or its engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)
BIN = os.environ.get("AZURE2_BIN")
if not BIN or not os.path.isfile(BIN):
    print("skip: AZURE2_BIN (the program) not set")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="pyazr_rate_")


def place(model, name):
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    shutil.copy(os.path.join(MODELS, model + ".azr"), os.path.join(d, "t.azr"))
    return d


def program(d, *args):
    out = subprocess.run([BIN, "rate", *args, "t.azr"], cwd=d, capture_output=True, text=True, timeout=600)
    if out.returncode != 0:
        raise RuntimeError(f"the program's rate failed in {d}: {out.stdout[-400:]}")
    with open(os.path.join(d, "output", "reactionrates.out")) as fh:
        rows = [ln.split() for ln in fh.read().splitlines()[1:] if ln.strip()]
    return [r[0] for r in rows], [r[1] for r in rows]


print("1. 12C(p,g): pyazr.rate against the program's rate command")
d = place("p12C", "grid")
t_prog, r_prog = program(d, "--rate-entrance=1", "--rate-exit=2", "--temperature-min=0.01", "--temperature-max=1",
                         "--temperature-step=0.07")
r = rate.reaction_rate(os.path.join(d, "t.azr"), 1, 2, t_min=0.01, t_max=1, t_step=0.07)
check(f"the grid: {len(r_prog)} temperatures, every rate as the program writes it",
      [f"{v:g}" for v in r.rate] == r_prog and [f"{v:g}" for v in r.T9] == t_prog,
      list(zip([f"{v:g}" for v in r.rate], r_prog)))
check("... with the units and the reaction", r.units["rate"] == "cm^3 mol^-1 s^-1" and r.units["T9"] == "GK"
      and r.reaction == "12C(p,g)13N", (r.units, r.reaction))
d = place("p12C", "file")
temps = [0.3, 0.01, 2.0, 0.1, 0.03]
with open(os.path.join(d, "temperatures.dat"), "w") as fh:
    fh.write("".join(f"{t!r}\n" for t in temps))
t_prog, r_prog = program(d, "--rate-entrance=1", "--rate-exit=2", "--temperature-file=temperatures.dat")
r_list = rate.reaction_rate(os.path.join(d, "t.azr"), 1, 2, temperatures=temps)
check("a temperature list: the program's rates, in the list's order",
      [f"{v:g}" for v in r_list.rate] == r_prog and np.array_equal(r_list.T9, temps),
      list(zip([f"{v:g}" for v in r_list.rate], r_prog)))

for name, (entrance, exit_key) in (("T1", (1, 2)), ("T3", (2, 1))):
    golden = os.path.join(ROOT, "tests", "upgrade", "golden", name, "rate", "reactionrates.out")
    d = os.path.join(work, name)
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", name), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    temps = np.loadtxt(os.path.join(d, "temperatures.dat"), ndmin=1)
    got = rate.reaction_rate(os.path.join(d, f"{name}.azr"), entrance, exit_key, temperatures=temps)
    with open(golden) as fh:
        want = [ln.split() for ln in fh.read().splitlines()[1:] if ln.strip()]
    check(f"{name} ({got.reaction}): the golden reactionrates.out of the program, digit for digit "
          f"({len(want)} temperatures)", [[f"{a:g}", f"{b:g}"] for a, b in zip(got.T9, got.rate)] == want,
          list(zip([f"{v:g}" for v in got.rate], [w[1] for w in want])))

print("\n2. a parameter file")
d = place("moved", "moved")
shutil.copy(os.path.join(ROOT, "tests", "reaction_rate", "expected", "moved.par"), d)
with open(os.path.join(d, "temps.dat"), "w") as fh:
    fh.write("0.1\n0.3\n1\n")
t_prog, r_prog = program(d, "--rate-entrance=1", "--rate-exit=2", "--temperature-file=temps.dat",
                         "--parameter-file=moved.par")
r_par = rate.reaction_rate(os.path.join(d, "t.azr"), 1, 2, temperatures=[0.1, 0.3, 1.0],
                           parameter_file=os.path.join(d, "moved.par"))
check("moved.par: the program's rates with --parameter-file", [f"{v:g}" for v in r_par.rate] == r_prog,
      list(zip([f"{v:g}" for v in r_par.rate], r_prog)))
r_azr = rate.reaction_rate(os.path.join(d, "t.azr"), 1, 2, temperatures=[0.1, 0.3, 1.0])
check("... and not the .azr's (the level moved: 6 % to 80 % apart)",
      np.all(np.abs(r_azr.rate / r_par.rate - 1) > 0.05),
      (r_azr.rate, r_par.rate))
check("... the source names the file", "moved.par" in r_par.source, r_par.source)

print("\n3. exit -1, the total radiative capture")
p12c = os.path.join(work, "grid", "t.azr")
total = rate.reaction_rate(p12c, 1, -1, grid="table")
one = rate.reaction_rate(p12c, 1, 2, grid="table")
check(f"exit -1 = exit 2 on the table grid ({total.T9.size} temperatures), bit for bit",
      np.array_equal(total.rate, one.rate) and total.T9[0] == 0.01 and total.T9[-1] == 10.0)
for bad, words in (((1, 0), "not 0"), ((1, 1), "different from the entrance"), ((1, -2), "exit -1")):
    try:
        rate.RateSession(p12c, *bad)
        check(f"pairs {bad} refused", False)
    except ValueError as e:
        check(f"pairs {bad} refused ({words})", words in str(e), str(e))

print("\n4. the bands against direct loops")
temps = [0.05, 0.2, 1.0]
with rate.RateSession(p12c, 1, 2) as s:
    names = [n for n, f in zip(s.names, s.fixed) if not f][:3]
    pos = [s.names.index(n) for n in names]
    mean = s.values[pos]
    C = np.diag((0.02 * np.abs(mean) + 1e-3) ** 2)
    C[0, 1] = C[1, 0] = 0.3 * np.sqrt(C[0, 0] * C[1, 1])
    draws = np.random.default_rng(7).multivariate_normal(mean, C, size=9)
    band = rate.rate_band(s, draws=draws, names=names, temperatures=temps)
    lin = rate.rate_band(s, covariance=C, names=names, temperatures=temps)
    base = s.values.copy()
with rate.RateSession(p12c, 1, 2) as fresh:
    loop = []
    for row in draws:
        v = base.copy()
        v[pos] = row
        loop.append(fresh.rate(v, temperatures=temps).rate)
    lo, mid, hi = np.percentile(np.array(loop), [16, 50, 84], axis=0)
    from pyazr.bands import step_sizes
    h = step_sizes(base, pos, 1e-4)
    g = []
    for c, hk in zip(pos, h):
        up, down = base.copy(), base.copy()
        up[c] += hk
        down[c] -= hk
        g.append((np.log(fresh.rate(up, temperatures=temps).rate) - np.log(fresh.rate(down, temperatures=temps).rate))
                 / (2 * hk))
    g = np.array(g)
    sigma = np.sqrt(np.einsum("it,ij,jt->t", g, C, g))
    best = fresh.rate(base, temperatures=temps).rate
check("draws: 16/50/84 % of 9 vectors = the direct loop's, bit for bit",
      np.array_equal(band.low, lo) and np.array_equal(band.median, mid) and np.array_equal(band.high, hi)
      and band.draws == 9, (band.low, lo))
check("... the band holds the best rate's neighbourhood (low < high)", np.all(band.low < band.high))
check("covariance: best x exp(-+sqrt(g^T C g)) = the direct loop's",
      np.allclose(lin.low, best * np.exp(-sigma), rtol=1e-14, atol=0)
      and np.allclose(lin.high, best * np.exp(sigma), rtol=1e-14, atol=0) and np.all(sigma > 0), (lin.low, sigma))
check("... a name the rate has not (a norm) is left out", rate.rate_band(
    rate.RateSession(p12c, 1, 2), covariance=np.diag([1e-4, 1e-4]), names=[names[0], "s01.norm"],
    temperatures=[0.2]).varied == [names[0]])

print("\n5. the Gamow window")
g = rate.gamow_window(p12c, 1, [0.1, 1.0, 3.0])
mu = 1.00782503 * 12.0 / (13.00782503)
e0 = 0.1220 * (36 * mu * g["T9"] ** 2) ** (1 / 3)
w = 0.2368 * (36 * mu * g["T9"] ** 5) ** (1 / 6)
check(f"12C+p: E0 {np.round(g['E0'], 4)} MeV, width {np.round(g['width'], 4)} MeV, the closed form to 1e-3",
      np.allclose(g["E0"], e0, rtol=1e-3) and np.allclose(g["width"], w, rtol=1e-3), (g["E0"], e0, g["width"], w))
n = rate.gamow_window(os.path.join(MODELS, "n12C.azr"), 1, [0.1])
check("a neutron entrance: no Gamow peak (NaN), kT given", np.isnan(n["E0"][0]) and abs(n["kT"][0] - 0.008617333262)
      < 1e-15, n)

print("\n6. the table")
ref = os.path.join(work, "ref.dat")
with open(ref, "w") as fh:
    fh.write("# T9 rate unc\nT9 rate\n" + "".join(f"{t} {v * 1.1} {v * 0.1}\n" for t, v in zip(one.T9, one.rate)))
cmp_ = rate.read_rate_table(ref)
check("a three-column table is read (header words skipped)", cmp_["T9"].size == one.T9.size
      and cmp_["uncertainty"] is not None)
out = rate.write_rate_table(os.path.join(work, "rate.csv"), one, comparison=cmp_)
text = open(out).read()
check("the CSV's header names the project, the parameters' source and the units",
      p12c in text and "parameters: the project's values" in text and "cm^3 mol^-1 s^-1" in text
      and "T9,rate,comparison,ratio" in text, text[:400])
ratio = np.array([float(ln.split(",")[3]) for ln in text.splitlines() if ln[:1].isdigit()])
check("... the ratio to the comparison is 1/1.1", np.allclose(ratio, 1 / 1.1, rtol=1e-9), ratio[:3])

print("\n7. the example pyazr/examples/reaction_rate.py")
out = subprocess.run([sys.executable, os.path.join(ROOT, "pyazr", "examples", "reaction_rate.py"),
                      os.path.join(work, "T1", "T1.azr"), "--entrance", "1", "--exit", "2", "--csv",
                      os.path.join(work, "ex.csv"), "--reaclib", os.path.join(work, "ex.reaclib")],
                     capture_output=True, text=True, timeout=600,
                     env=dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")))
row = [ln.split() for ln in out.stdout.splitlines() if ln.strip().startswith("0.1 ")]
check("it runs on T1 and prints the golden rate at T9 = 0.1 (5.59562e-07) with its Gamow peak (0.0847 MeV)",
      out.returncode == 0 and len(row) == 1 and row[0][1] == "5.59562e-07" and row[0][2] == "0.0847"
      and os.path.isfile(os.path.join(work, "ex.csv")), out.stdout[-800:] + out.stderr[-400:])
check("... and writes REACLIB sets (2 sets, within 5 %)",
      os.path.isfile(os.path.join(work, "ex.reaclib")) and "REACLIB: 2 set(s)" in out.stdout,
      [ln for ln in out.stdout.splitlines() if "REACLIB" in ln])

shutil.rmtree(work, ignore_errors=True)
print(f"\n{len(failures)} failure(s)" if failures else "\nall ok")
sys.exit(1 if failures else 0)
