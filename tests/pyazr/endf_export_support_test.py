#!/usr/bin/env python3
"""Two repairs the ENDF export asked for (docs/upgrade/endf.md section 13.15), outside pyazr.endf:

(a) pyazr.evaluate._legendre_rows (the Legendre frame of an extrapolation) cuts a frame of equal counts by a reshape:
    on frames of equal counts (0 and 3), of varying counts and on an empty frame its arrays equal, value for value,
    those of the loop it replaces (written out here);
(b) pyazr.Project.new(nuclear_part=) drops the rows of <parameterSettings> that name a parameter of a removed record:
    from 13N, a project with one data set opens without one W_PARAMETER_SETTING_UNKNOWN warning (16 before), and
    the rows of the levels stay.
Exits 77 (skipped) without the compiled engine module.  Run from anywhere.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
    from pyazr.evaluate import _legendre_rows
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def loop_rows(flat):
    """The loop _legendre_rows had: per point a count, then that many coefficients."""
    rows, i = [], 0
    flat = np.asarray(flat, float)
    while i < flat.size:
        count = int(round(flat[i]))
        rows.append(flat[i + 1:i + 1 + count].copy())
        i += 1 + count
    return rows


rng = np.random.default_rng(7)
frames = {"empty": np.zeros(0), "1000 points without coefficients": np.zeros(1000),
          "500 points of 3": np.concatenate([[3.0, *rng.normal(size=3)] for _ in range(500)]),
          "varying counts": np.concatenate([[float(c), *rng.normal(size=c)] for c in rng.integers(0, 5, 300)]),
          "first count 1, then 0": np.array([1.0, 0.5, 0.0, 0.0])}
for name, frame in frames.items():
    a, b = loop_rows(frame), _legendre_rows(frame)
    same = len(a) == len(b) and all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(a, b))
    check(f"_legendre_rows, {name}: the arrays of the loop, value for value", same, (len(a), len(b)))

tmp = tempfile.mkdtemp(prefix="export_support_")
d = os.path.join(tmp, "13N")
shutil.copytree(os.path.join(ROOT, "tests", "13N"), d, ignore=shutil.ignore_patterns("output*", "checks*"))
for sub in ("output", "checks"):
    os.makedirs(os.path.join(d, sub), exist_ok=True)
p = pyazr.Project.new(os.path.join(d, "new.azr"), nuclear_part=os.path.join(d, "13N.azr"))
np.savetxt(os.path.join(d, "probe.dat"), [(0.5, 90.0, 1.0, 1.0), (1.0, 90.0, 1.0, 1.0)], fmt="%.17g")
p.add_dataset("probe", reaction=p.reaction(1, 1), quantity="dsigma_domega", unit="b/sr",
              file=os.path.join(d, "probe.dat"), columns="E_lab[MeV] theta_cm[deg] value error",
              norm=pyazr.Nuisance(1, "fixed"))
p.save()
text = open(os.path.join(d, "new.azr")).read()
settings = text[text.index("<parameterSettings>"):text.index("</parameterSettings>")]
# the engine prints its warnings of the set-up: a session in a process of its own, its output read
import subprocess
run = subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); import pyazr; "
                      "pyazr.Session(sys.argv[2]).close()", ROOT, os.path.join(d, "new.azr")],
                     capture_output=True, text=True, timeout=600)
unknown = [line for line in (run.stdout + run.stderr).splitlines() if "W_PARAMETER_SETTING_UNKNOWN" in line]
check("Project.new of 13N: the session opens (exit 0) without a W_PARAMETER_SETTING_UNKNOWN warning",
      run.returncode == 0 and not unknown, (run.returncode, unknown[:3], run.stderr[-300:]))
check("Project.new of 13N: the rows of the levels kept, none naming a data set", "Level 3 Energy" in settings and
      "baumann_1758." not in settings, settings[:400])
shutil.rmtree(tmp, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
