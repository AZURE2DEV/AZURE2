#!/usr/bin/env python3
"""Session.sample_with_engine against the program's command mcmc (docs/upgrade/pyazr3.md section 12).

Each case runs the program (AZURE2_BIN: `AZURE2 mcmc --walkers=W --steps=N --sampler-threads=T --basis=B`) on one
copy of a project and the session's call with the same answers on another copy, and compares what the sampler
writes.  The sampler's seed is fixed in the engine, and both go through one engine function (RunSampler), so the
comparison is exact:
- samples.mcmc, walkers.mcmc and param.mcmc are equal byte for byte (and parameters.out and the AZUREOut_* files of
  the best walker, written every 10 steps);
- the returned chain, log_prob, log_likelihood and log_prior are the numbers of samples.mcmc exactly (read here with
  float(); the file has 17 digits), with the shapes (steps, walkers, parameters) and (steps, walkers), and
  log_prob = log_likelihood + log_prior as the file's.
Cases: T3 (tests/upgrade; 8 free parameters, a hybrid potential, a neutron channel) in the rwa basis at 1 and at 2
sampler threads, and in the physical basis at 2; T1 (13 free parameters, capture) in the rwa basis at 2 threads;
T3 in Park's parametrization (--use-park) in the physical basis; T3 from a parameter file with moved free values,
in both bases, against a session of the project given those values as x0; T3 from a parameter file whose fixed
level values are moved too against a session made with that file; T3 with a Gaussian prior in <mcmc> <parameters>;
a chain of 100 steps continued to 200 against the program's continued chain.
Controls: one progress record per step; a stop request set at step 30 ends the run there (status "stopped", the
files complete to that step) and a later call continues it; an exception in progress is raised; each keyword out of
its range is refused with the command line's words; a row of <mcmc> that names no parameter is refused with the
sampler's words; the session's values and objective are unchanged by a run.

Run: python3 tests/pyazr/sample_with_engine_test.py   (AZURE2_BIN: the program).  About 40 s.
"""
import filecmp
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import pyazr
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)
PROGRAM = os.environ.get("AZURE2_BIN")
if not PROGRAM or not os.path.exists(PROGRAM):
    print("SKIP: AZURE2_BIN (the program) is not set")
    sys.exit(77)
if not pyazr.describe()["engine"].get("mcmc", True):
    print("SKIP: this build has no sampler (USE_MCMC)")
    sys.exit(77)

logging.getLogger("pyazr.engine").setLevel(logging.ERROR)   # the sampler's warnings are in result.log
tmp = tempfile.mkdtemp(prefix="sample_with_engine_")


def copy_project(name, label, edit=None):
    """A copy of tests/upgrade/<name> with empty output directories; the path of its project file."""
    work = os.path.join(tmp, label)
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", name), work,
                    ignore=shutil.ignore_patterns("output", "checks", "expected"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    path = os.path.join(work, name + ".azr")
    if edit:
        with open(path) as f:
            text = f.read()
        with open(path, "w") as f:
            f.write(edit(text))
    return path


def program(azr, walkers, steps, threads, basis, *more):
    """The command on the project: (exit status, standard output)."""
    proc = subprocess.run([PROGRAM, "mcmc", f"--walkers={walkers}", f"--steps={steps}",
                           f"--sampler-threads={threads}", f"--basis={basis}", *more, os.path.basename(azr)],
                          stdin=subprocess.DEVNULL, cwd=os.path.dirname(azr), capture_output=True, text=True,
                          timeout=1200)
    return proc.returncode, proc.stdout


def session(azr, **options):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pyazr.Session(azr, **options)


def file_numbers(path):
    """samples.mcmc as (header, rows of floats)."""
    with open(path) as f:
        lines = f.read().splitlines()
    return lines[0].split(","), np.array([[float(v) for v in line.split(",")] for line in lines[1:]])


def compare(label, prog_azr, sess_azr, r, walkers, steps):
    out_p, out_s = os.path.join(os.path.dirname(prog_azr), "output"), os.path.join(os.path.dirname(sess_azr), "output")
    same = {f: os.path.exists(os.path.join(out_p, f)) and os.path.exists(os.path.join(out_s, f))
            and filecmp.cmp(os.path.join(out_p, f), os.path.join(out_s, f), shallow=False)
            for f in ("samples.mcmc", "walkers.mcmc", "param.mcmc")}
    check(f"{label}: samples.mcmc, walkers.mcmc and param.mcmc equal the program's byte for byte", all(same.values()),
          str(same))
    others = sorted(f for f in os.listdir(out_p) if f.startswith("AZUREOut_") or f == "parameters.out")
    differ = [f for f in others if not (os.path.exists(os.path.join(out_s, f))
                                        and filecmp.cmp(os.path.join(out_p, f), os.path.join(out_s, f), shallow=False))]
    check(f"{label}: the best walker's {len(others)} output files equal the program's byte for byte",
          others and not differ, str(differ))
    header, rows = file_numbers(os.path.join(out_p, "samples.mcmc"))
    n = len(r.names)
    shapes = (r.chain.shape == (steps, walkers, n) and r.log_prob.shape == (steps, walkers)
              and r.log_likelihood.shape == (steps, walkers) and r.log_prior.shape == (steps, walkers))
    check(f"{label}: chain {r.chain.shape}, log_prob {r.log_prob.shape}; status {r.status}; "
          f"{len(header) - 5} parameter columns for {n} names",
          shapes and r.status == "completed" and len(header) - 5 == n and (r.steps, r.walkers) == (steps, walkers)
          and rows.shape == (steps * walkers, 5 + n))
    if not shapes or rows.shape != (steps * walkers, 5 + n):
        return
    order = np.lexsort((rows[:, 1], rows[:, 0]))   # by step, then walker
    rows = rows[order]
    exact = (np.array_equal(r.chain.reshape(-1, n), rows[:, 5:]) and np.array_equal(r.log_prob.ravel(), rows[:, 2])
             and np.array_equal(r.log_likelihood.ravel(), rows[:, 3])
             and np.array_equal(r.log_prior.ravel(), rows[:, 4]))
    check(f"{label}: the returned arrays are the program's samples.mcmc number for number "
          f"({rows.size} numbers, exact equality)", exact)
    moved = len(np.unique(r.chain[:, :, 0])) > walkers
    check(f"{label}: log_prob = log_likelihood + log_prior (1e-12); the walkers moved",
          np.allclose(r.log_prob, r.log_likelihood + r.log_prior, rtol=1e-12, atol=1e-9) and moved)
    check(f"{label}: the result names the files written",
          os.path.samefile(r.samples_file, os.path.join(out_s, "samples.mcmc"))
          and os.path.samefile(r.walkers_file, os.path.join(out_s, "walkers.mcmc"))
          and os.path.samefile(r.output_dir, out_s) and any("Acceptance fraction" in line for line in r.log))


try:
    # ---------------------------------------------------------------- the chain of the program
    for name, walkers, basis, threads, park in (("T3", 20, "rwa", 1, False), ("T3", 20, "rwa", 2, False),
                                                ("T3", 20, "physical", 2, False), ("T1", 30, "rwa", 2, False),
                                                ("T3", 20, "physical", 1, True)):
        label = f"{name}, {basis}, {threads} thread{'s' if threads > 1 else ''}" + (", Park" if park else "")
        print(label)
        prog = copy_project(name, f"{name}_{basis}_{threads}_{park}_program")
        sess = copy_project(name, f"{name}_{basis}_{threads}_{park}_session")
        status, out = program(prog, walkers, 100, threads, basis, "--overwrite-samples", *(["--use-park"] * park))
        check(f"{label}: the program ran", status == 0, out[-300:])
        with session(sess, use_park=park) as s:
            x_before = s.parameters.x0.copy()
            f_before = s.evaluate(x_before).objective
            records = []
            r = s.sample_with_engine(walkers, 100, sampler_threads=threads, basis=basis, overwrite_samples=True,
                                     progress=records.append)
            compare(label, prog, sess, r, walkers, 100)
            check(f"{label}: {len(records)} progress records, steps 1 to 100, with the best walker's log-probability",
                  [q["step"] for q in records] == list(range(1, 101)) and all(q["steps"] == 100 for q in records)
                  and all(abs(q["log_prob"] - q["log_likelihood"] - q["log_prior"]) <= 1e-9 * abs(q["log_prob"])
                          for q in records)
                  and abs(records[-1]["log_prob"] - r.log_prob[-1].max()) <= 1e-12 * abs(records[-1]["log_prob"]),
                  str(records[-1:]))
            check(f"{label}: the session's values and objective are unchanged by the run",
                  np.array_equal(s.parameters.x0, x_before) and s.evaluate(x_before).objective == f_before)
            if name == "T3" and basis == "rwa" and threads == 1:
                one_thread = r
            if name == "T3" and basis == "rwa" and threads == 2:
                check("T3, rwa: the chain at 2 sampler threads is the chain at 1, number for number",
                      np.array_equal(r.chain, one_thread.chain) and np.array_equal(r.log_prob, one_thread.log_prob))
                again = s.sample_with_engine(walkers, 100, sampler_threads=2)
                check("T3, rwa: a call on a complete samples.mcmc returns the file as it is",
                      np.array_equal(again.chain, r.chain) and again.status == "completed")

    # ---------------------------------------------------------------- a parameter file, and x0
    print("T3 from a parameter file (the program) against the session's parameter file and x0")

    def write_start(path, s, free, fixed_factor):
        """A parameter file with every digit (Parameters.write gives 8, as param.sav has them): the free values
        `free`, the fixed level values times `fixed_factor`."""
        s.parameters.write(path, free)
        with open(path) as f:
            tag = f.readline()
        with open(path, "w") as f:
            f.write(tag)
            for q in s.parameters:
                level = q.kind in ("energy", "width")
                value = free[q.free_index] if not q.fixed else q.start * (fixed_factor if level else 1.0)
                f.write(f"{q.name:>28s} {float(value)!r:>26s} {0.0:.7e}\n")

    for basis in ("rwa", "physical"):
        prog, sess = copy_project("T3", f"file_{basis}_program"), copy_project("T3", f"file_{basis}_session")
        with session(sess) as s:
            moved = s.parameters.x0 * (1 + 0.002 * np.cos(np.arange(s.parameters.x0.size)))
            write_start(os.path.join(os.path.dirname(prog), "start.sav"), s, moved, 1.0)
            status, out = program(prog, 20, 100, 1, basis, "--overwrite-samples", "--parameter-file=start.sav")
            check(f"x0, {basis}: the program ran from the parameter file",
                  status == 0 and "Reading User Parameter File" in out, out[-300:])
            r = s.sample_with_engine(20, 100, basis=basis, x0=moved, overwrite_samples=True)
            compare(f"x0 against the program's parameter file, {basis}", prog, sess, r, 20, 100)
            if basis == "rwa":
                check("x0: the chain is not the one from the project's values",
                      not np.array_equal(r.chain, one_thread.chain))
    # fixed level values of the file that are not the project's: a session made with the file has them
    prog, sess = copy_project("T3", "fixed_program"), copy_project("T3", "fixed_session")
    start = os.path.join(os.path.dirname(prog), "start.sav")
    with session(sess) as s:
        moved = s.parameters.x0 * (1 + 0.002 * np.cos(np.arange(s.parameters.x0.size)))
        n_fixed = sum(1 for q in s.parameters if q.fixed and q.kind in ("energy", "width") and q.start)
        write_start(start, s, moved, 1.01)
    status, out = program(prog, 20, 100, 1, "rwa", "--overwrite-samples", "--parameter-file=start.sav")
    check(f"parameter file with {n_fixed} fixed level values moved by 1 %: the program ran",
          status == 0 and n_fixed > 0, out[-300:])
    with session(sess, parameter_file=start) as s:
        r = s.sample_with_engine(20, 100, overwrite_samples=True)
        compare("a session made with the parameter file (fixed values moved)", prog, sess, r, 20, 100)

    # ---------------------------------------------------------------- a prior of <mcmc>
    print("T3 with a Gaussian prior of <mcmc> <parameters>")
    prior = "<mcmc>\n<parameters>\nenergy_1 3.03 3.0305 0.0004 1\n</parameters>\n</mcmc>\n"
    prog = copy_project("T3", "prior_program", lambda text: text + prior)
    sess = copy_project("T3", "prior_session", lambda text: text + prior)
    status, out = program(prog, 20, 100, 2, "rwa", "--overwrite-samples")
    check("prior: the program applied it", status == 0 and "Applied 1 Gaussian prior(s)" in out, out[-300:])
    with session(sess) as s:
        r = s.sample_with_engine(20, 100, sampler_threads=2, overwrite_samples=True)
        compare("prior of <mcmc>", prog, sess, r, 20, 100)
        check("prior: the log-prior is not the one without it, and the sampler's lines say it was applied",
              not np.array_equal(r.log_prior, one_thread.log_prior)
              and any("Applied 1 Gaussian prior(s)" in line for line in r.log))

    # ---------------------------------------------------------------- a continued chain
    print("T3, 100 steps continued to 200")
    prog, sess = copy_project("T3", "continued_program"), copy_project("T3", "continued_session")
    for steps, more in ((100, ("--overwrite-samples",)), (200, ())):
        status, out = program(prog, 20, steps, 1, "rwa", *more)
    check("continued: the program resumed its chain", status == 0 and "Resuming from step 100" in out, out[-300:])
    with session(sess) as s:
        s.sample_with_engine(20, 100, overwrite_samples=True)
        r = s.sample_with_engine(20, 200)
        compare("continued to 200 steps", prog, sess, r, 20, 200)
        check("continued: the first 100 steps are the uninterrupted chain's",
              np.array_equal(r.chain[:100], one_thread.chain))

    # ---------------------------------------------------------------- the controls
    print("controls")
    sess = copy_project("T3", "controls")
    with session(sess) as s:
        stop = threading.Event()
        seen = []

        def at_thirty(record):
            seen.append(record["step"])
            if record["step"] == 30:
                stop.set()
        r = s.sample_with_engine(20, 100, overwrite_samples=True, progress=at_thirty, stop=stop)
        header, rows = file_numbers(r.samples_file)
        check(f"stop: a request set at step 30 ends the run there (status {r.status}, {r.steps} steps, "
              f"{len(rows)} rows, last record {seen[-1]})",
              r.status == "stopped" and r.steps == 30 and r.chain.shape == (30, 20, 8) and len(rows) == 600
              and seen[-1] == 30 and np.array_equal(r.chain, one_thread.chain[:30]))
        check("stop: walkers.mcmc holds the 20 walkers of the last step",
              np.array_equal(np.loadtxt(r.walkers_file, delimiter=","), r.chain[-1]))
        r = s.sample_with_engine(20, 100)
        check(f"stop: the next call continues the chain to 100 steps (status {r.status}, {r.steps} steps)",
              r.status == "completed" and r.chain.shape == (100, 20, 8)
              and np.array_equal(r.chain[:30], one_thread.chain[:30]))
        request = pyazr.session_fit.StopRequest()
        request.set()
        r = s.sample_with_engine(20, 100, overwrite_samples=True, stop=request)
        check(f"stop: a StopRequest set before the run ends it after its first step ({r.steps} step)",
              r.status == "stopped" and r.steps == 1)

        class Broken(Exception):
            pass

        def broken(record):
            raise Broken(f"at step {record['step']}")
        try:
            s.sample_with_engine(20, 100, overwrite_samples=True, progress=broken)
            check("progress: an exception in it is raised", False, "no exception")
        except Broken as e:
            check("progress: an exception in it ends the run and is raised", "step 1" in str(e), str(e))

        refusals = (("walkers=1", dict(walkers=1), "--walkers takes a whole number from 2 to 1000"),
                    ("walkers=1001", dict(walkers=1001), "--walkers takes a whole number from 2 to 1000"),
                    ("steps=99", dict(steps=99), "--steps takes a whole number from 100 to 10000000"),
                    ("sampler_threads=0", dict(sampler_threads=0),
                     "--sampler-threads takes a whole number of at least 1"),
                    ("parameter_spread=0", dict(parameter_spread=0.0),
                     "--parameter-spread takes a number above 0 and up to 50"),
                    ("parameter_spread=51", dict(parameter_spread=51.0), "--parameter-spread takes"),
                    ("energy_spread=2", dict(energy_spread=2.0), "--energy-spread takes a number above 0 and up to 1"),
                    ("basis='widths'", dict(basis="widths"), "'rwa'"),
                    ("walkers=20.5", dict(walkers=20.5), "a whole number"),
                    ("x0 of 3 values", dict(x0=[1.0, 2.0, 3.0]), "8 free parameters"))
        before = open(r.samples_file, "rb").read()
        for label, keywords, words in refusals:
            arguments = dict(dict(walkers=20, steps=100), **keywords)
            try:
                s.sample_with_engine(**arguments)
                check(f"refusal: {label}", False, "no error")
            except ValueError as e:
                check(f"refusal: {label} is a ValueError naming the option's range", words in str(e), str(e))
        try:
            s.sample_with_engine(20, 100, stop=object())
            check("refusal: a stop without is_set()", False, "no error")
        except TypeError as e:
            check("refusal: a stop without is_set() is a TypeError", "is_set" in str(e), str(e))
        check("refusal: a refused call leaves samples.mcmc alone", open(r.samples_file, "rb").read() == before)

    sess = copy_project("T3", "bad_prior", lambda text: text + "<mcmc>\n<parameters>\nenergy_99 1 1 0.1 1\n"
                                                                "</parameters>\n</mcmc>\n")
    with session(sess) as s:
        try:
            s.sample_with_engine(20, 100, overwrite_samples=True)
            check("refusal: a row of <mcmc> that names no parameter", False, "no error")
        except pyazr.EngineError as e:
            check("refusal: a row of <mcmc> that names no parameter is an EngineError with the sampler's words",
                  "names no parameter of this model" in str(e) and "energy_99" in str(e), str(e))
        check("refusal: it wrote no samples.mcmc",
              not os.path.exists(os.path.join(os.path.dirname(sess), "output", "samples.mcmc")))
    with session(copy_project("T3", "extrapolations"), extrapolations=True) as s:
        try:
            s.sample_with_engine(20, 100)
            check("refusal: a session of the extrapolations", False, "no error")
        except pyazr.EngineError as e:
            check("refusal: a session of the extrapolations has no data to sample", "extrapolations=True" in str(e))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
