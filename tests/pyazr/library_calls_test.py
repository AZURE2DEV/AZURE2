#!/usr/bin/env python3
"""pyazr's library calls are quiet, write only where the caller says, and refuse with three parts.

  (a) quiet: a session whose set-up makes the engine print (a stale external-capture cache: "WARNING: ... holds 40
      external capture amplitudes ... Recalculating"), reaction_rate (the engine's progress bars) print nothing on
      standard output (each run in a child process, its descriptor 1 captured); the engine's lines reach the logger
      pyazr.engine instead (WARNING lines at level WARNING), and the rate's numbers are those of a run that printed;
  (b) the rate writes nothing beside the project (no output/cross_section_grid.dat, no intEC.extrap) and leaves no
      temporary directory once closed;
  (c) a set-up the engine refuses (T3 with its 2+ levels inactive: the phase-shift data set has no channel) raises
      AZURE2Error whose what/allowed/fix are set, in the message, and whose message quotes the engine's ERROR line;
  (d) a fit vector of the wrong length or with a NaN, given to Session.evaluate and to the Objective's calls, is a
      ParameterVectorError with its three parts.

Run from anywhere:  python3 tests/pyazr/library_calls_test.py
"""
import glob
import logging
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
    import numpy as np  # noqa: E402
    import pyazr  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)


def child(code):
    """Run ``code`` in a child Python with this pyazr; returns (exit code, stdout, stderr)."""
    env = dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""))
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=900, env=env)
    return run.returncode, run.stdout, run.stderr


def files(top):
    return {os.path.relpath(os.path.join(a, f), top) for a, _, fs in os.walk(top) for f in fs}


def parts_in_message(e):
    return all(getattr(e, k) for k in ("what", "allowed", "fix")) and all(
        getattr(e, k) in str(e) for k in ("what", "allowed", "fix"))


warnings.simplefilter("ignore", pyazr.EngineWarning)
tmp = tempfile.mkdtemp(prefix="library_calls_")
try:
    for name in ("T1", "T3"):
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", name), os.path.join(tmp, name),
                        ignore=shutil.ignore_patterns("output", "checks"))
    T1 = os.path.join(tmp, "T1", "T1.azr")
    T3 = os.path.join(tmp, "T3", "T3.azr")

    # (a) quiet
    print("(a) nothing on standard output")
    os.makedirs(os.path.join(tmp, "T1", "output"))
    with pyazr.Session(T1):
        pass                                  # output/intEC.dat for every capture data set
    fewer = pyazr.Project.open(T1)
    fewer.set("cap_g0_lab090", active="no")  # the cache now holds amplitudes this set-up does not need
    fewer.save(os.path.join(tmp, "T1", "fewer.azr"))
    code = ("import logging, warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            "logging.basicConfig(level=logging.INFO, format='LOGGED %(name)s %(levelname)s %(message)s')\n"
            f"pyazr.Session({os.path.join(tmp, 'T1', 'fewer.azr')!r}).close()\n")
    rc, out, err = child(code)
    check("a set-up that makes the engine print: exit 0, nothing on standard output", rc == 0 and out == "",
          (rc, out[-500:], err[-500:]))
    check("its lines reach the logger pyazr.engine, the WARNING at level WARNING",
          "LOGGED pyazr.engine WARNING WARNING: 'output/intEC.dat' holds" in err, err[-800:])

    code = ("import logging, warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            "logging.basicConfig(level=logging.INFO, format='LOGGED %(name)s %(levelname)s %(message)s')\n"
            f"r = pyazr.rate.reaction_rate({T1!r}, 1, 2, temperatures=[0.1, 1.0])\n"
            "import sys; print('RATE', *[repr(float(v)) for v in r.rate], file=sys.stderr)\n")
    before = files(os.path.dirname(T1))
    scratch_before = set(glob.glob(os.path.join(tempfile.gettempdir(), "pyazr_rate_*")))
    rc, out, err = child(code)
    check("reaction_rate: exit 0, nothing on standard output", rc == 0 and out == "", (rc, out[-500:], err[-500:]))
    check("reaction_rate: the progress bar's last frame reaches the logger at INFO",
          "LOGGED pyazr.engine INFO \t[*************************] 100%" in err, err[-800:])
    rate = [line for line in err.splitlines() if line.startswith("RATE ")]
    # the rate of the program's run on this project, written by a run that printed its progress (T1: 1 -> 2)
    check("reaction_rate: the numbers of a run that printed",
          rate == ["RATE 5.595623195798498e-07 530.5291324488501"], rate)

    # (b) the rate writes nothing beside the project
    print("(b) the rate's files")
    check("reaction_rate: nothing new beside the project", files(os.path.dirname(T1)) == before,
          sorted(files(os.path.dirname(T1)) - before))
    check("reaction_rate: no temporary directory left",
          set(glob.glob(os.path.join(tempfile.gettempdir(), "pyazr_rate_*"))) <= scratch_before, "")
    with pyazr.rate.RateSession(T1, 1, 2) as rs:
        scratch = rs._scratch
        check("RateSession: its engine writes into its temporary directory", os.path.isdir(scratch), scratch)
        rs.rate(temperatures=[1.0])
        check("RateSession: the grid file is there, not beside the project",
              os.path.isfile(os.path.join(scratch, "cross_section_grid.dat"))
              and not os.path.exists(os.path.join(tmp, "T1", "output", "cross_section_grid.dat")), os.listdir(scratch))
    check("RateSession.close() removes it", not os.path.exists(scratch), scratch)

    # (c) a refused set-up
    print("(c) a refused set-up")
    # the 2+ levels made inactive past the Project's checks, as a hand edit of the file would (the Project edit
    # refuses this edit: tests/pyazr/project_phase_shift_test.py)
    from pyazr import nuclear_part as npart
    p = pyazr.Project.open(T3)
    for lv in [lv for lv in p.levels() if lv.jpi == "2+"]:
        npart._replace(p._doc, npart._engine_call("level", npart._engine_module().set_level_edits, p._doc,
                                                  npart._level_address(lv.label), energy=None, fixed=None,
                                                  active=False))
    inactive = p.save(os.path.join(tmp, "T3", "inactive.azr"))
    logging.getLogger("pyazr.engine").setLevel(logging.CRITICAL)
    try:
        pyazr.Session(inactive)
        check("T3 without its 2+ levels: refused", False, "accepted")
    except pyazr.AZURE2Error as e:
        check("T3 without its 2+ levels: e.what, e.allowed, e.fix set and in the message", parts_in_message(e),
              (e.what, e.allowed, e.fix, str(e)[:500]))
        check("the message quotes the engine's ERROR line",
              "ERROR: phase-shift segment #3 (J = 2, l = 2) asks for channel spin s = 0" in str(e)
              and "ERROR: phase-shift segment #3" in e.what, str(e)[:800])
        check("the message no longer points at an output above", "check the output above" not in str(e), str(e))

    # (d) the fit vector's refusals
    print("(d) ParameterVectorError")
    with pyazr.Session(T1) as s:
        x = s.parameters.x0
        n = x.size
        calls = {"Session.evaluate": s.evaluate, "Objective(x)": s.objective(),
                 "Objective.value_and_grad": s.objective().value_and_grad,
                 "Objective.residuals": s.objective().residuals, "Objective.log_prob": s.objective().log_prob}
        bad = x.copy()
        bad[1] = np.nan
        for label, call in calls.items():
            for what, vector, words in (("one value too many", np.append(x, 1.0),
                                         [f"got {n + 1}", f"expected the free parameter vector ({n} values"]),
                                        ("a NaN", bad, ["entry 1", s.parameters.free.names[1], "is nan"])):
                try:
                    call(vector)
                    check(f"{label} with {what}: refused", False, "accepted")
                except pyazr.ParameterVectorError as e:
                    check(f"{label} with {what}: three parts, in the message, naming the entry",
                          parts_in_message(e) and label in e.what and all(w in str(e) for w in words),
                          (e.what, e.allowed, e.fix, str(e)))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
