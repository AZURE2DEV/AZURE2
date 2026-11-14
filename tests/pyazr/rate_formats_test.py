#!/usr/bin/env python3
"""The rate libraries' formats (pyazr/rate_formats.py).

  1. REACLIB fit: a rate made of one known set is recovered (largest deviation 1e-9, the set's numbers), a rate of
     two known sets is reproduced within 1 % by at most three sets;
  2. REACLIB file: every line 74 characters; read back, the species, label, flag, Q and the sets to the format's
     seven digits, and the rate of the sets read equals the fit's to 1e-6;
  3. refusals with their three parts: a rate not covering 0.01-10 GK (no extrapolation), a tolerance not reached,
     a rate <= 0, a wrong label or flag;
  4. STARLIB: T1's 12C(p,g) on the 60 temperatures with a covariance band: written and read back, T9, the rate and
     the factor uncertainty to the written digits, the header's chapter, species and Q; f.u. = exp(sigma) of the
     band; refused without a band and off the grid;
  5. the project's species and Q: T1 (p + 12C -> 13N, Q = 1.94349 MeV, chapter 4), T3 (p + 7Li -> 4He + 4He, Q =
     the two thresholds' difference, chapter 5);
  6. T1's rate on the table grid as REACLIB sets: the number of sets and the deviation printed, within 5 %.

Run from anywhere:  python3 tests/pyazr/rate_formats_test.py
"""
import os
import shutil
import sys
import tempfile
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    from pyazr import rate, rate_formats as rf
    from pyazr.errors import ParameterError
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="pyazr_rate_formats_")
T = np.array(rate.TEMPERATURE_GRIDS["table"][1])
SPECIES = {"chapter": 4, "reactants": ["p", "c12"], "products": ["n13"], "Q": 1.94349}


def refused(name, call, words):
    try:
        call()
        check(f"refused: {name}", False, "no error")
    except ParameterError as e:
        check(f"refused: {name} ({words})", words in str(e) and e.what and e.allowed and e.fix, str(e))


print("1. the REACLIB fit")
a1 = np.array([17.1, 0.0, -13.69, -0.23, 4.62, -2.0, -0.67])
explicit = np.exp(a1[0] + a1[1] / T + a1[2] * T ** (-1 / 3) + a1[3] * T ** (1 / 3) + a1[4] * T + a1[5] * T ** (5 / 3)
                  + a1[6] * np.log(T))
check("a set's rate is exp(a0 + a1/T9 + a2 T9^-1/3 + a3 T9^1/3 + a4 T9 + a5 T9^5/3 + a6 ln T9)",
      np.allclose(rf.reaclib_rate([a1], T), explicit, rtol=1e-13, atol=0))
one = SimpleNamespace(T9=T, rate=explicit)
f1 = rf.fit_reaclib(one, tolerance=1e-6)
check(f"one known set: recovered with 1 set, largest deviation {f1.max_deviation:.1e}",
      f1.n_sets == 1 and f1.max_deviation < 1e-9 and np.allclose(f1.sets[0], a1, rtol=1e-6, atol=1e-6), f1.sets)
a2 = np.array([-5.0, -2.5, 0.0, 0.0, 0.0, 0.0, -1.5])  # a narrow resonance's set: exp(-E/kT) T9^-3/2
two = SimpleNamespace(T9=T, rate=rf.reaclib_rate([a1, a2], T))
f2 = rf.fit_reaclib(two, tolerance=0.01, max_sets=3)
print("   ", f2.summary())
check("two known sets: reproduced within 1 % by at most three", f2.max_deviation <= 0.01 and f2.n_sets <= 3)

print("\n2. the REACLIB file")
path = rf.write_reaclib(os.path.join(work, "r.dat"), f2, SPECIES, label="azr2", flag="r")
lines = open(path).read().splitlines()
check("every line 74 characters; a header entry and three lines per set",
      all(len(ln) == 74 for ln in lines) and len(lines) == 3 + 3 * f2.n_sets, [len(ln) for ln in lines])
back = rf.read_reaclib(path)
check("read back: species, label, flag, Q", all(b["species"] == ["p", "c12", "n13"] and b["label"] == "azr2"
                                               and b["flag"] == "r" and b["chapter"] == 4 and b["Q"] == 1.94349
                                               for b in back) and len(back) == f2.n_sets, back[:1])
check("... the sets to seven digits, the rate of the sets read = the fit's to 1e-6",
      all(np.allclose(b["a"], a, rtol=1e-6, atol=1e-12) for b, a in zip(back, f2.sets))
      and np.allclose(rf.reaclib_rate([b["a"] for b in back], T), rf.reaclib_rate(f2.sets, T), rtol=1e-5))

print("\n3. refusals")
refused("a rate from 0.05 to 1 GK", lambda: rf.fit_reaclib(SimpleNamespace(T9=T[(T >= 0.05) & (T <= 1)],
                                                                         rate=one.rate[(T >= 0.05) & (T <= 1)])),
        "nothing is extrapolated")
rng = np.random.default_rng(3)
refused("a tolerance not reached", lambda: rf.fit_reaclib(SimpleNamespace(T9=T, rate=one.rate * np.exp(
    rng.normal(0, 0.05, T.size))), tolerance=1e-4, max_sets=1), "do not reproduce the rate")
refused("a rate of zero", lambda: rf.fit_reaclib(SimpleNamespace(T9=T, rate=np.where(T < 0.02, 0.0, one.rate))),
        "not positive")
refused("a label of five characters", lambda: rf.write_reaclib(os.path.join(work, "x"), f1, SPECIES, label="abcde"),
        "four characters")
refused("a flag", lambda: rf.write_reaclib(os.path.join(work, "x"), f1, SPECIES, flag="x"), "resonance flag")

print("\n4. STARLIB, with the engine")
try:
    import pyazr._azure2  # noqa: F401
    have_engine = True
except Exception as e:  # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False
if have_engine:
    d = os.path.join(work, "T1")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    azr = os.path.join(d, "T1.azr")
    sp = rf.species_of(azr, 1, 2)
    check(f"T1: {sp}", sp["chapter"] == 4 and sp["reactants"] == ["p", "c12"] and sp["products"] == ["n13"]
          and abs(sp["Q"] - 1.9434900024034978) < 1e-12, sp)
    t3 = os.path.join(ROOT, "tests", "upgrade", "T3", "T3.azr")
    s3 = rf.species_of(t3, 2, 1)
    check(f"T3: {s3}", s3["chapter"] == 5 and s3["reactants"] == ["p", "li7"] and s3["products"] == ["he4", "he4"]
          and abs(s3["Q"] - (17.254406959108284 + 0.09183948734448354)) < 1e-9, s3)
    with rate.RateSession(azr, 1, 2) as s:
        names = [n for n, f in zip(s.names, s.fixed) if not f][:3]
        pos = [s.names.index(n) for n in names]
        C = np.diag((0.01 * np.abs(s.values[pos]) + 1e-3) ** 2)
        r = s.rate(grid="starlib")
        band = rate.rate_band(s, covariance=C, names=names, grid="starlib")
        table = s.rate(grid="table")
    path = rf.write_starlib(os.path.join(work, "s.dat"), r, band, sp, label="azr2")
    got = rf.read_starlib(path)
    check("the header: chapter, species, label, Q", got["chapter"] == 4 and got["species"] == ["p", "c12", "n13"]
          and got["label"] == "azr2" and abs(got["Q"] - sp["Q"]) < 1e-5, got["header"])
    check("60 rows: T9, the rate to 7 digits, f.u. = exp(sigma) to 5 digits",
          got["T9"].size == 60 and np.allclose(got["T9"], rf.STARLIB_T9, rtol=1e-3)
          and np.allclose(got["rate"], band.best, rtol=1e-6) and np.allclose(got["fu"], band.high / band.best,
                                                                            rtol=1e-4), got["fu"][:3])
    refused("STARLIB without a band", lambda: rf.write_starlib(os.path.join(work, "x"), r, None, sp), "no band")
    refused("STARLIB off its grid", lambda: rf.write_starlib(os.path.join(work, "x"), table, band, sp),
            "60 temperatures")
    print("\n6. T1 as REACLIB sets")
    fit = rf.fit_reaclib(table, tolerance=0.05)
    print("   ", fit.summary())
    check("within 5 %", fit.max_deviation <= 0.05)

shutil.rmtree(work, ignore_errors=True)
print(f"\n{len(failures)} failure(s)" if failures else "\nall ok")
sys.exit(1 if failures else 0)
