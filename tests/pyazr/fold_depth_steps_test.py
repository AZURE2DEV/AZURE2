#!/usr/bin/env python3
"""The keyword fold_depth_steps of pyazr.Session: the panels of a target's depth integral with a beam spread or
straggling (docs/upgrade/energy_grid.md section 3.4; the program's --fold-depth-steps=K, AZURE_FOLD_DEPTH_STEPS).

Model: tests/upgrade/T2 (data set 5: a fitted layer with straggling and a beam spread; data set 1: a layer with Bohr
straggling and a spread).  Checked:
  (a) None (not given) is the default 8 (Config::kFoldDepthStepsDefault): every data set's chi2 equal to the session
      with fold_depth_steps=8; a value that is not the default, 4, differs from both 8 and 0 where a nested kernel is;
  (b) 0, the node intervals, moves the data sets with a nested kernel and no other, by less than 1e-3 relative;
  (c) the keyword equals the program: a session with 0 writes the AZUREOut files of AZURE2 calculate
      --fold-depth-steps=0 byte for byte, one without it those of the program without the option, and the two differ;
  (d) the keyword over the variable: AZURE_FOLD_DEPTH_STEPS=0 in a child process gives the chi2 of 0, and with
      fold_depth_steps=8 or 4 in that process those of 8 or 4;
  (e) refusals: -1, 1001, 2.5, nan, inf, True and "4" raise ValueError naming the keyword and the range;
  (f) Session.options records the keyword as given (None or the number), so **options rebuilds the same model.

Run:  python3 tests/pyazr/fold_depth_steps_test.py   (AZURE2_BIN for (c))
"""
import math
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np                    # noqa: E402
    import pyazr                          # noqa: E402
    from pyazr import Session             # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

os.environ.pop("AZURE_FOLD_DEPTH_STEPS", None)
work = tempfile.mkdtemp(prefix="pyazr_fold_depth_steps_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T2"), os.path.join(work, "T2"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    folder = os.path.join(work, "T2")
    os.makedirs(os.path.join(folder, "output"), exist_ok=True)
    azr = os.path.join(folder, "T2.azr")

    def chi2s(**kw):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", pyazr.EngineWarning)
            with Session(azr, **kw) as s:
                x0 = np.array(s.parameters.x0, float)
                return np.array([d.chi2 for d in s.evaluate(x0).data]), dict(s.options)

    DEFAULT = 8  # Config::kFoldDepthStepsDefault (the PI's decision of 2026-11-28; energy_grid.md section 3.4)
    given_none, options_none = chi2s()
    default, options_default = chi2s(fold_depth_steps=DEFAULT)
    four, options_four = chi2s(fold_depth_steps=4)
    zero, options_zero = chi2s(fold_depth_steps=0)
    check(f"(a) not given: the chi2 of fold_depth_steps={DEFAULT}, every data set", np.array_equal(given_none, default),
          (given_none, default))
    moved = [i for i in range(len(default)) if zero[i] != default[i]]
    rel = max((abs(zero[i] / default[i] - 1.0) for i in moved), default=0.0)
    print(f"    data sets moved by 0 against {DEFAULT}: {[i + 1 for i in moved]}, largest relative change of chi2 {rel:.3e}")
    check("(a) 4, not the default: differs from the default and from 0 in a data set with a nested kernel",
          any(four[i] != default[i] and four[i] != zero[i] for i in moved), (four, default, zero))
    check("(b) 0 moves some data sets (a nested kernel in T2)", len(moved) > 0, zero - default)
    check("(b) ... and not every one (the others have no spread inside a target)", len(moved) < len(default), moved)
    check("(b) ... by less than 1e-3 relative", rel < 1e-3, rel)

    binary = os.environ.get("AZURE2_BIN")
    if binary and os.path.exists(binary):
        def fresh(tag):
            d = os.path.join(work, tag)
            shutil.copytree(folder, d, ignore=shutil.ignore_patterns("output", "checks"))
            os.makedirs(os.path.join(d, "output"))
            os.makedirs(os.path.join(d, "checks"))
            return d

        def outputs(d):
            out = os.path.join(d, "output")
            return {f: open(os.path.join(out, f)).read() for f in sorted(os.listdir(out)) if f.startswith("AZUREOut")}

        def program(tag, *flags):
            d = fresh(tag)
            subprocess.run([binary, "calculate", "--no-gui", *flags, "T2.azr"], cwd=d, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL, check=True,
                           env={k: v for k, v in os.environ.items() if k != "AZURE_FOLD_DEPTH_STEPS"})
            return outputs(d)

        def session(tag, **kw):
            d = fresh(tag)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", pyazr.EngineWarning)
                with Session(os.path.join(d, "T2.azr"), **kw) as s:
                    s.write_outputs()
            return outputs(d)
        p0, p4 = program("p0", "--fold-depth-steps=0"), program("p4")
        s0, s4 = session("s0", fold_depth_steps=0), session("s4")
        check("(c) the program's AZUREOut files at 0 differ from those by default", p0 and p0 != p4,
              [f for f in p0 if p0[f] == p4.get(f)])
        check("(c) a session with 0 writes the program's AZUREOut files at --fold-depth-steps=0, byte for byte",
              p0 == s0, [f for f in p0 if p0[f] != s0.get(f)])
        check("(c) a session without the keyword writes the program's files without the option", p4 == s4,
              [f for f in p4 if p4[f] != s4.get(f)])
    else:
        print("    (c) skipped: no AZURE2_BIN")

    child = ("import sys, warnings, numpy as np; sys.path.insert(0, %r); import pyazr; "
             "warnings.simplefilter('ignore');"
             "k = None if sys.argv[1] == 'None' else int(sys.argv[1]);"
             "s = pyazr.Session(%r, fold_depth_steps=k); x = np.array(s.parameters.x0, float);"
             "print(repr([d.chi2 for d in s.evaluate(x).data]))") % (ROOT, azr)

    def in_child(k):
        env = dict(os.environ, AZURE_FOLD_DEPTH_STEPS="0")
        out = subprocess.run([sys.executable, "-c", child, k], env=env, capture_output=True, text=True, check=True)
        return np.array(eval(out.stdout.strip().splitlines()[-1]))
    check("(d) AZURE_FOLD_DEPTH_STEPS=0, keyword not given: the chi2 of 0", np.array_equal(in_child("None"), zero))
    check(f"(d) the keyword {DEFAULT} over AZURE_FOLD_DEPTH_STEPS=0: the chi2 of {DEFAULT}",
          np.array_equal(in_child(str(DEFAULT)), default))
    check("(d) the keyword 4 over AZURE_FOLD_DEPTH_STEPS=0: the chi2 of 4", np.array_equal(in_child("4"), four))

    for bad in (-1, 1001, 2.5, math.nan, math.inf, True, "4"):
        try:
            chi2s(fold_depth_steps=bad)
            check(f"(e) fold_depth_steps={bad!r} refused", False, "accepted")
        except ValueError as err:
            check(f"(e) fold_depth_steps={bad!r} refused with ValueError",
                  "fold_depth_steps needs a whole number from 0 to 1000" in str(err), str(err))

    check("(f) Session.options holds fold_depth_steps as given",
          options_none.get("fold_depth_steps") is None and options_zero.get("fold_depth_steps") == 0
          and options_four.get("fold_depth_steps") == 4, (options_none, options_zero))
    rebuilt, _ = chi2s(**{k: v for k, v in options_zero.items()})
    check("(f) **options of the session with 0 rebuilds its chi2", np.array_equal(rebuilt, zero), (rebuilt, zero))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
