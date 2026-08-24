#!/usr/bin/env python3
"""A Python session must be able to write the files the CLI writes.

Everything a run produces -- AZUREOut_*, chiSquared.out -- came only from the
executable. A pyazr session could compute all of it and still had to shell out
to the binary to get files a colleague, or the GUI's plot tab, could open.

The catch this pins: chiSquared.out reports the per-segment chi-squared *stored
on each segment*, and only the chi-squared evaluation sets that. Writing after a
bare forward pass produces a well-formed file full of zeros -- which is worse
than no file, because it looks like an answer.

Compares against the binary's own output rather than a recorded reference, so
the two cannot drift apart.

Needs the compiled engine and an AZURE2 binary; skips cleanly without either.

Run from anywhere:  python3 tests/pyazr/output_files_test.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

binary = None
for cand in sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        binary = cand
        break
if binary is None:
    print("skip: no AZURE2 binary to compare against")
    sys.exit(0)

# hybrid_potential is the quick one: three elastic segments, no external
# capture, so a full run is seconds rather than minutes.
SOURCE = os.path.join(ROOT, "tests", "hybrid_potential")


def run(work):
    """Produce output/ with the binary, and return what it wrote."""
    subprocess.run([binary, "--no-gui", "--no-readline", "hybrid_potential.azr"],
                   cwd=work, input="1\n\n\n7\n", text=True,
                   capture_output=True, timeout=900)
    out = os.path.join(work, "output")
    return {os.path.basename(f): open(f).read()
            for f in sorted(glob.glob(os.path.join(out, "*")))
            if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}


with tempfile.TemporaryDirectory() as tmp:
    cli_dir = os.path.join(tmp, "cli")
    py_dir = os.path.join(tmp, "py")
    for d in (cli_dir, py_dir):
        shutil.copytree(SOURCE, d)
        for junk in ("output", "checks"):
            shutil.rmtree(os.path.join(d, junk), ignore_errors=True)
        os.makedirs(os.path.join(d, "output"))
        os.makedirs(os.path.join(d, "checks"))

    from_cli = run(cli_dir)
    check("the binary wrote something", bool(from_cli), "no output files")

    with azure2(os.path.join(py_dir, "hybrid_potential.azr"), cwd=py_dir) as m:
        outdir = m.write_output_files()
    check("write_output_files returned the output directory",
          os.path.isdir(outdir), outdir)
    from_py = {os.path.basename(f): open(f).read()
               for f in sorted(glob.glob(os.path.join(outdir, "*")))
               if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    # EData::WriteOutputFiles owns the data files. A CLI run additionally
    # leaves param.par (written at initialization) and parameters.out (from
    # CNuc::PrintTransformParams at the end); neither belongs to this call, so
    # compare what it does own and say what it does not.
    owned = {"chiSquared.out"} | {n for n in from_cli if n.startswith("AZUREOut_")}
    check("every data file the CLI wrote is there",
          owned <= set(from_py),
          f"missing {sorted(owned - set(from_py))}")
    for name in sorted(owned & set(from_py)):
        check(f"{name} is byte-identical", from_py[name] == from_cli[name])
    print("        not written by this call, and not expected: "
          + ", ".join(sorted(set(from_cli) - owned)))

    chi = from_py.get("chiSquared.out", "")
    total = [l for l in chi.splitlines() if l.startswith("Total-Chi-Squared:")]
    check("chiSquared.out carries a total", bool(total))
    if total:
        value = float(total[0].split()[1])
        print(f"        {total[0].strip()}")
        check("and it is not zero -- the forward pass alone would leave it so",
              value != 0.0, total[0])

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all output-file checks passed")
