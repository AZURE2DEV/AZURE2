#!/usr/bin/env python3
"""An external-capture final state must be a level of the project; pyazr sees the engine's new rule (2026-10-26).

The API always switches external capture on, so it is where the two branches of CNuc::ParseExternalCapture show:

  * models/missing_level.azr (tests/ec_final_state_listed): the only data segment exits into the capture pair and
    its final state is not listed -> the engine prints the ERROR line and Initialize fails, so pyazr.azure2(...)
    raises AZURE2Error.  Older engines ran it with a hidden level (a fixed energy_2 and a FREE width_2_1 in no
    .azr line) that shifted every later energy_N name.
  * models/unused.azr: the same levels, but the only segment is elastic 16O(p,p), so nothing uses the capture
    pair -> a WARNING, the project runs, and the parameter list holds exactly the three listed levels.
  * models/listed.azr runs as before.

Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/ec_final_state_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MODELS = os.path.join(ROOT, "tests", "ec_final_state_listed")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def stage(model, work):
    d = os.path.join(work, model)
    os.makedirs(os.path.join(d, "output")); os.makedirs(os.path.join(d, "checks"))
    shutil.copy(os.path.join(MODELS, "models", model + ".azr"), os.path.join(d, "t.azr"))
    shutil.copytree(os.path.join(MODELS, "data"), os.path.join(d, "data"))
    return os.path.join(d, "t.azr")


def main():
    try:
        import pyazr
        from pyazr import azure2
    except ImportError as e:
        print(f"skipping: {e}")
        return 0
    work = tempfile.mkdtemp(prefix="ec_final_state_")
    try:
        # (1) listed: runs
        s = azure2(stage("listed", work))
        names = [p.name for p in s.parameters.energies]
        check("listed final state: session opens", len(names) > 0)
        check("listed final state: four levels in the parameter list", len(names) == 4, f"{names}")
        # (2) missing level, capture segment present: refused
        raised = None
        try:
            azure2(stage("missing_level", work))
        except Exception as e:  # pyazr.AZURE2Error
            raised = e
        check("missing final state used by a segment: pyazr.azure2 raises",
              raised is not None and type(raised).__name__ == "AZURE2Error", repr(raised))
        # (3) missing level, no segment uses the pair: WARNING, runs, no hidden level
        # the engine writes on the C-level stdout, which redirect_stdout does not see: capture file descriptor 1
        cap = os.path.join(work, "unused_stdout.txt")
        saved = os.dup(1); sys.stdout.flush()
        fd = os.open(cap, os.O_WRONLY | os.O_CREAT | os.O_TRUNC); os.dup2(fd, 1); os.close(fd)
        try:
            s = azure2(stage("unused", work))
        finally:
            sys.stdout.flush(); os.dup2(saved, 1); os.close(saved)
        out = open(cap).read()
        check("unused capture pair: session opens", s is not None)
        check("unused capture pair: WARNING names the pair and the final state",
              "WARNING: external capture into pair 3 needs its final state" in out, out[-400:])
        names = [p.name for p in s.parameters.energies]
        check("unused capture pair: exactly the three listed levels, no hidden one", len(names) == 3, f"{names}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("ec_final_state_test:", "PASS" if not failures else f"FAIL ({len(failures)})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
