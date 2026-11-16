#!/usr/bin/env python3
"""An external-capture final state too close above a particle threshold is refused in pyazr as by the program
(E_FINAL_STATE_UNBOUND); one farther above it is a narrow resonance and is computed.

tests/upgrade/T1 with its 3/2- final state of 12C(p,g)13N[1.0] (photon pair 3, external capture E1, M1, E2) moved,
with its level: to 1.94354 MeV, 0.05 keV above S_p = 1.94349 MeV (eta = 129 > 100), Project.validate() (what the web
page's Check shows) lists the code on the data set cap_g1, and Session(...) raises AZURE2Error with the engine's line;
to 2.5 MeV it validates and evaluates a finite chi^2; at 1.94 MeV, below S_p, it validates and opens.  The program's
side: tests/ec_final_state_unbound/check.sh.

Needs the compiled engine; skips cleanly without it.  Run from anywhere:
  python3 tests/pyazr/ec_final_state_unbound_test.py
"""
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T1 = os.path.join(ROOT, "tests", "upgrade", "T1")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def variant(work, name, energy):
    """T1 with the 3/2- level at 1.0 MeV and photon pair 3's final state (column 17) at `energy` (MeV)."""
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    shutil.copytree(os.path.join(T1, "data"), os.path.join(d, "data"))
    lines, inside = [], False
    with open(os.path.join(T1, "T1.azr")) as f:
        for line in f.read().split("\n"):
            inside = (inside or line.startswith("<levels>")) and not line.startswith("</levels>")
            fields = line.split()
            if inside and len(fields) >= 31:
                if fields[:3] == ["1.5", "-1", "1.0"]:
                    fields[2] = str(energy)
                if fields[5] == "3":
                    fields[16] = str(energy)
                line = " ".join(fields)
            lines.append(line.replace("13N[1.0]", f"13N[{energy}]"))
    path = os.path.join(d, "t.azr")
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def main():
    try:
        import pyazr
    except ImportError as e:
        print(f"skipping: {e}")
        return 77
    work = tempfile.mkdtemp(prefix="ec_final_state_unbound_")
    try:
        close = variant(work, "close", 1.94354)
        diagnostics = pyazr.Project.open(close).validate()
        codes = [(d["code"], d.get("record")) for d in diagnostics]
        check("final state 0.05 keV above S_p: Project.validate() lists E_FINAL_STATE_UNBOUND on cap_g1",
              codes == [("E_FINAL_STATE_UNBOUND", "cap_g1")], f"{codes}")
        raised = None
        try:
            pyazr.Session(close)
        except Exception as e:  # pyazr.AZURE2Error
            raised = e
        check("final state 0.05 keV above S_p: Session raises AZURE2Error with the engine's line",
              raised is not None and type(raised).__name__ == "AZURE2Error" and
              re.search(r"E_FINAL_STATE_UNBOUND.*cannot be evaluated this close to the threshold", str(raised)) is not None,
              repr(raised))
        above = variant(work, "above", 2.5)
        check("final state at 2.5 MeV (a narrow resonance): Project.validate() is empty",
              pyazr.Project.open(above).validate() == [])
        s = pyazr.Session(above)
        chi2 = s.evaluate(s.parameters.x0).chi2
        check("final state at 2.5 MeV: the session evaluates a finite chi^2", chi2 == chi2 and chi2 < 1e99, f"{chi2}")
        s.close() if hasattr(s, "close") else None
        below = variant(work, "below", 1.94)
        check("final state at 1.94 MeV: Project.validate() is empty", pyazr.Project.open(below).validate() == [])
        s = pyazr.Session(below)
        check("final state at 1.94 MeV: the session evaluates a finite chi^2", s.evaluate(s.parameters.x0).chi2 < 1e99)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if failures:
        print(f"ec_final_state_unbound (pyazr): FAIL ({len(failures)})")
        return 1
    print("ec_final_state_unbound (pyazr): PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
