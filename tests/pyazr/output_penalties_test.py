#!/usr/bin/env python3
"""chiSquared.out written through pyazr carries the penalty totals (audit 2026-10-25, item 30; 2026-10-27).

`EData::WriteOutputFiles` prints the Park, abundance and scattering-length totals from fields on the EData object.
The CLI's chi-squared path stores them there; the API's `CalculateChi2RWA` never did (pyazr adds those terms on the
Python side), so a `write_output_files()` from a session wrote 0 for all three -- "Total-Scattering-Length-Chi-Squared"
was missing from the file -- or, after a session had been through the CLI path, a stale value.  `AZUREAPI::WriteOutputFiles`
now evaluates them at the current parameters before writing.

Check on tests/scattering_length/be7p_fit.azr (two phase-shift segments plus a <scatteringLengths> constraint):
the file's "Total-Scattering-Length-Chi-Squared" equals `penalties()["scattering_length"]` (which is > 0 here).
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/output_penalties_test.py
"""
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "scattering_length", "be7p_fit.azr")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

work = tempfile.mkdtemp(prefix="output_penalties.")
try:
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copytree(os.path.join(os.path.dirname(SRC), "data"), os.path.join(work, "data"))
    shutil.copy(SRC, os.path.join(work, "m.azr"))
    s = azure2(os.path.join(work, "m.azr"), cwd=work)
    want = float(s.penalties()["scattering_length"])
    out = s.write_output_files()
    text = open(os.path.join(out, "chiSquared.out")).read()
    m = re.search(r"Total-Scattering-Length-Chi-Squared:\s*(\S+)", text)
    got = float(m.group(1)) if m else None
    check(f"the model has a scattering-length penalty to report ({want:.6g})", want > 1e-6)
    check("chiSquared.out carries Total-Scattering-Length-Chi-Squared", got is not None,
          "the line is missing: the API wrote the file with the penalty at 0")
    check(f"the file's value equals penalties()['scattering_length'] ({got} vs {want:.6g})",
          got is not None and abs(got - want) <= 1e-5 * max(1.0, abs(want)))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
