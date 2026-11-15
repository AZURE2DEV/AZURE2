#!/usr/bin/env python3
"""pyazr's library calls are quiet, write only where the caller says, and refuse with three parts.

  (a) quiet: a session whose set-up makes the engine print (a stale external-capture cache: "WARNING: ... holds 40
      external capture amplitudes ... Recalculating"), reaction_rate (the engine's progress bars) print nothing on
      standard output (each run in a child process, its descriptor 1 captured); the engine's lines reach the logger
      pyazr.engine instead (WARNING lines at level WARNING), and the rate's numbers are those of a run that printed;
      an Objective's calls likewise (13N's note at the first call; the engine stood in for by an object that
      writes on descriptor 1, for every call);
  (b) the rate writes nothing beside the project (no output/cross_section_grid.dat, no intEC.extrap) and leaves no
      temporary directory once closed;
  (c) a set-up the engine refuses (T3 with its 2+ levels inactive: the phase-shift data set has no channel) raises
      AZURE2Error whose what/allowed/fix are set, in the message, and whose message quotes the engine's ERROR line;
  (d) a fit vector of the wrong length or with a NaN, given to Session.evaluate and to the Objective's calls, is a
      ParameterVectorError with its three parts;
  (e) trimmed_model, extrapolation_bands and tabulate write nothing beside the project, not even while they run
      (a thread watches the directory): trimmed_model's project is in the work directory it returns (a temporary
      one by default), extrapolation_bands and tabulate remove theirs.

Run from anywhere:  python3 tests/pyazr/library_calls_test.py
"""
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
    env = dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""), TMPDIR=CHILD_TMP)
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=900, env=env)
    return run.returncode, run.stdout, run.stderr


def files(top):
    return {os.path.relpath(os.path.join(a, f), top) for a, _, fs in os.walk(top) for f in fs}


def parts_in_message(e):
    return all(getattr(e, k) for k in ("what", "allowed", "fix")) and all(
        getattr(e, k) in str(e) for k in ("what", "allowed", "fix"))


warnings.simplefilter("ignore", pyazr.EngineWarning)
tmp = tempfile.mkdtemp(prefix="library_calls_")
# the temporary directories of the calls (children and this process) in directories of this test, so that what is
# left there is theirs (other processes of the machine make their own in the system's)
CHILD_TMP = os.path.join(tmp, "child_tmp")
os.makedirs(CHILD_TMP)
tempfile.tempdir = os.path.join(tmp, "own_tmp")
os.makedirs(tempfile.tempdir)
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

    # a session's call that makes the engine print (the parameter records of a project whose limit the engine's
    # transformation cannot reach: its "**WARNING: Denominator less than zero"), and a fit whose progress callback
    # prints: the call is quiet, the callback's text is shown
    limited = pyazr.Project.open(T1)
    limited.set_parameter_setting("1/2+ 2.365", (1, 0, 0.5), limits=(16500.0, 66000.0))
    limited.save(os.path.join(tmp, "T1", "limited.azr"))
    code = ("import warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            f"s = pyazr.Session({os.path.join(tmp, 'T1', 'limited.azr')!r})\n"
            "s.parameters\n"
            "s.fit_with_engine(max_calls=3, progress=lambda r: print('PROGRESS', r['calls'], flush=True))\n")
    rc, out, err = child(code)
    check("Session.parameters and fit_with_engine of a project the engine warns about: not on standard output",
          rc == 0 and "Denominator" not in out and "Denominator less than zero" in err, (rc, out[-400:], err[-400:]))
    check("fit_with_engine: the progress callback's prints are shown", "PROGRESS 1" in out, out[-400:])

    # an Objective's calls are held as a session's are.  No project of the tests makes the engine print during one,
    # so the engine is stood in for by an object that writes a line on descriptor 1 and then asks the real engine
    code = ("import logging, os, warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            "logging.basicConfig(level=logging.INFO, format='LOGGED %(name)s %(levelname)s %(message)s')\n"
            f"s = pyazr.Session({T1!r})\n"
            "f = s.objective()\n"
            "want = (f(f.x0), f.value_and_grad(f.x0)[0], f.residuals(f.x0).sum(), f.log_prob(f.x0))\n"
            "class Loud:\n"
            "    def __init__(self, engine): self.engine = engine\n"
            "    def __getattr__(self, name):\n"
            "        call = getattr(self.engine, name)\n"
            "        def said(*a, **k):\n"
            "            os.write(1, ('WARNING: the engine speaks in ' + name + '\\n').encode())\n"
            "            return call(*a, **k)\n"
            "        return said\n"
            "s._sess = Loud(s._sess)\n"
            "got = (f(f.x0), f.value_and_grad(f.x0)[0], f.residuals(f.x0).sum(), f.log_prob(f.x0))\n"
            "f.describe()\n"
            "print('SAME', got == want, flush=True)\n"
            "s._sess = s._sess.engine\n"
            "s.close()\n")
    rc, out, err = child(code)
    check("Objective: what the engine writes during a call is not on standard output; the numbers are the same",
          rc == 0 and out == "SAME True\n", (rc, out[-400:], err[-400:]))
    spoken = ["calculate_objective_rwa", "calculate_objective_grad_rwa", "objective_rows_rwa", "log_prob_rwa",
              "objective_term_names"]
    check("Objective: the lines of each of its engine calls reach the logger pyazr.engine at level WARNING",
          all(f"LOGGED pyazr.engine WARNING WARNING: the engine speaks in {name}" in err for name in spoken),
          err[-800:])
    # and with the engine itself: 13N's sub-point grids make it write a note at the objective's first call
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), os.path.join(tmp, "13N"),
                    ignore=shutil.ignore_patterns("output", "checks", "expected"))
    code = ("import logging, warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            "logging.basicConfig(level=logging.INFO, format='LOGGED %(name)s %(levelname)s %(message)s')\n"
            f"s = pyazr.Session({os.path.join(tmp, '13N', '13N.azr')!r})\n"
            "f = s.objective()\n"
            "print('VALUE', f(f.x0) == s.evaluate(f.x0).objective, flush=True)\n"
            "s.close()\n")
    rc, out, err = child(code)
    check("Objective of 13N: the engine's note is not on standard output and reaches the logger at INFO",
          rc == 0 and out == "VALUE True\n" and "LOGGED pyazr.engine INFO Note: the sub-point grids" in err,
          (rc, out[-400:], err[-600:]))

    code = ("import logging, warnings, pyazr\n"
            "warnings.simplefilter('ignore')\n"
            "logging.basicConfig(level=logging.INFO, format='LOGGED %(name)s %(levelname)s %(message)s')\n"
            f"r = pyazr.rate.reaction_rate({T1!r}, 1, 2, temperatures=[0.1, 1.0])\n"
            "import sys; print('RATE', *[repr(float(v)) for v in r.rate], file=sys.stderr)\n")
    before = files(os.path.dirname(T1))
    rc, out, err = child(code)
    check("reaction_rate: exit 0, nothing on standard output", rc == 0 and out == "", (rc, out[-500:], err[-500:]))
    check("reaction_rate: the progress bar's last frame reaches the logger at INFO",
          "LOGGED pyazr.engine INFO \t[*************************] 100%" in err, err[-800:])
    rate = [line for line in err.splitlines() if line.startswith("RATE ")]
    # the rate of the program's run on this project, written by a run that printed its progress (T1: 1 -> 2)
    # (to 1e-8: these are computed with the Coulomb table; the exact functions, where there is no table file, give
    # 5.595623194312042e-07 and 530.5291324346432, 3e-10 away)
    want = (5.595623195798498e-07, 530.5291324488501)
    got = [float(v) for v in rate[0].split()[1:]] if len(rate) == 1 else []
    check("reaction_rate: the numbers of a run that printed",
          len(got) == 2 and all(abs(g / w - 1) < 1e-8 for g, w in zip(got, want)), rate)

    # (b) the rate writes nothing beside the project
    print("(b) the rate's files")
    check("reaction_rate: nothing new beside the project", files(os.path.dirname(T1)) == before,
          sorted(files(os.path.dirname(T1)) - before))
    check("reaction_rate: no temporary directory left",
          not os.listdir(CHILD_TMP), os.listdir(CHILD_TMP))
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

    # (e) the band and tabulation calls write nothing beside the project
    print("(e) nothing beside the project")
    import threading
    import time

    def watched(call):
        """(result, files that appeared beside T1 at any time during the call)"""
        before, seen, done = files(os.path.dirname(T1)), set(), []

        def watch():
            while not done:
                seen.update(files(os.path.dirname(T1)) - before)
                time.sleep(0.01)
        thread = threading.Thread(target=watch)
        thread.start()
        try:
            result = call()
        finally:
            done.append(True)
            thread.join()
        return result, sorted(seen | (files(os.path.dirname(T1)) - before))

    left_before = set(os.listdir(tempfile.gettempdir()))
    (path, work), seen = watched(lambda: pyazr.trimmed_model(T1, keys=[1]))
    check("trimmed_model: nothing beside the project; its project in the work directory it returns",
          not seen and os.path.dirname(path) == work and not work.startswith(os.path.dirname(T1)), (seen, path, work))
    shutil.rmtree(work)
    with pyazr.Session(T1) as s:
        x0 = s.parameters.x0
        n_rmatrix = sum(1 for p in s.parameters if not p.fixed and p.kind in ("energy", "width"))
    bands, seen = watched(lambda: pyazr.extrapolation_bands(T1, keys=[1], params=x0,
                                                            covariance=np.eye(n_rmatrix) * 1e-6))
    check("extrapolation_bands: a band, nothing beside the project", list(bands) == [1] and not seen, seen)
    tables, seen = watched(lambda: pyazr.tabulate(T1, [(1, 2)], 0.5, 1.0, rel_tol=0.05, min_points=5, max_rounds=1))
    check("tabulate: a table, nothing beside the project", len(tables[(1, 2)].sigma) > 5 and not seen, seen)
    left = set(os.listdir(tempfile.gettempdir())) - left_before
    check("no temporary directory of theirs left", not left, left)
finally:
    tempfile.tempdir = None
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
