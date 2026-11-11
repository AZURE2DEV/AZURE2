#!/usr/bin/env python3
"""The external-capture node table keys on the relativistic switch (audit 2026-10-25, checklist t3-25; 2026-10-27).

The ECTable of external-capture radial integrals lives for the whole process and is keyed on the pair, channel and
run options.  Up to dev b80b857 `relativistic=1` was not in the key, so a pyazr session that evaluated a model and
then the same model with relativistic kinematics reused nodes computed for the other, and its relativistic objective
differed from that of a fresh process (eta, rho and the flux factor of every node change).  The key carries the switch
now (and the intEC signature too: tests/ec_relativistic_signature).

Checks on tests/13N (p + 12C capture with external capture):
  (1) in one process, objective(non-relativistic) then objective(relativistic) equals the relativistic objective of a
      fresh process to 1e-9 relative;
  (2) the two objectives differ by more than 1e-4 relative (teeth).
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/ec_relativistic_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def objective_of(path, work):
    from pyazr import azure2
    return azure2(path, cwd=work).objective()


if __name__ == "__main__" and len(sys.argv) == 3 and sys.argv[1] == "--fresh":
    # child: the relativistic objective in a process of its own
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    print(repr(objective_of(sys.argv[2], os.path.dirname(sys.argv[2]))))
    sys.exit(0)

try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2, AzrModel
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

work = tempfile.mkdtemp(prefix="ec_relativistic.")
try:
    # one directory per evaluation: each writes its own output/intEC.dat, so no run reads another's file (the
    # amplitudes' text round trip moves the objective by ~1e-5, which would mask what is being tested here)
    def model_dir(name, relativistic):
        d = os.path.join(work, name)
        os.makedirs(os.path.join(d, "output"))
        shutil.copytree(os.path.join(os.path.dirname(SRC), "data"), os.path.join(d, "data"))
        m = AzrModel.from_file(SRC)
        m.relativistic = relativistic
        path = os.path.join(d, name + ".azr")
        m.write(path)
        return path, d

    plain, plain_dir = model_dir("plain", False)
    rel, rel_dir = model_dir("rel", True)
    rel2, rel2_dir = model_dir("rel_fresh", True)
    check("the relativistic variant carries the config key", AzrModel.from_file(rel).relativistic and not AzrModel.from_file(plain).relativistic)

    o_plain = objective_of(plain, plain_dir)
    o_rel_after = objective_of(rel, rel_dir)       # same process: the node table already holds non-relativistic nodes
    out = subprocess.run([sys.executable, os.path.abspath(__file__), "--fresh", rel2], capture_output=True, text=True,
                         env=dict(os.environ, OMP_NUM_THREADS="2"), cwd=ROOT)
    lines = [l for l in out.stdout.splitlines() if l.strip()]
    o_rel_fresh = float(lines[-1]) if lines else float("nan")
    check(f"non-relativistic {o_plain:.6f} and relativistic {o_rel_fresh:.6f} objectives differ by more than 1e-4",
          abs(o_rel_fresh / o_plain - 1.0) > 1e-4, f"{abs(o_rel_fresh / o_plain - 1.0):.2e}")
    check(f"relativistic objective after a non-relativistic one in the same process = fresh process ({o_rel_after:.6f} vs {o_rel_fresh:.6f})",
          abs(o_rel_after / o_rel_fresh - 1.0) < 1e-9,
          f"relative difference {abs(o_rel_after / o_rel_fresh - 1.0):.2e}: the node table was shared across the switch")
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
