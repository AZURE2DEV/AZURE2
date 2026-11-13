#!/usr/bin/env python3
"""Two sessions evaluating at the same time on two Python threads do not see each other (no hidden state).

The engine calls of a session release the GIL (``py::call_guard<py::gil_scoped_release>`` on the evaluations:
calculate_rwa, segment_chi2, chi2_and_grad and the rest; the set-up in the Session constructor keeps it), so two
Python threads drive two sessions truly at once: one thread's set-up runs beside the other's evaluation, and two
evaluations run side by side.  Every value a calculation keeps between calls must therefore belong to the session
(its Config: the run's scratch, the external-capture amplitudes, the notices), never to the process.

Sessions, each with its own copy of its project directory:
  A   tests/13N (12C+p: elastic scattering and capture, external capture from the energy tables), default options;
  B   tests/18O_p_a_thm (Trojan-horse data on the shared mesh) with the Coulomb table off;
  N   tests/13N with --no-ec-table (the interpolated external-capture amplitudes, ECAmplitudeCache);
  M   tests/13N with --no-ec-table and the channel radius of its first pair at 5.0 fm (the same model's keys as N);
  T   tests/13N with a second, coarser Coulomb table file (made here with the table tool; `--coulomb-table`);
  R   tests/13N as written: its data, output and checks paths relative to the project's directory;
  S   tests/18O_p_a_thm as written (relative paths), the Coulomb table off.
The paths of A, B, N, M, T are made absolute; R and S keep the project's relative paths, and their sessions also
write their output files (write_output_files), whose bytes are part of what is compared: a session resolves its
relative paths against its own directory (Config::baseDirectory), never against the process's working directory,
which two threads share (pyazr changed it per call before, and a session writing while the other had entered its
own directory wrote into the other's).  Pairs run at once on two threads, 20 times each (each time both threads build
their session and evaluate it, a barrier starting them together): A with B, N with M, A with T, R with S.  In one more
process the sessions run one after another: A T A, T A T (two table files in one process).  Every evaluation must give
the bytes of the same session alone in a fresh process.  What is compared: the cross sections of every data point
(calculate_rwa), the chi-squared per segment and the objective with its analytic gradient (chi2_and_grad), as raw
IEEE bytes, and for R and S the output files written.

Before the sessions owned their state, the binding swapped a process-wide Config pointer and amplitude cache per call
(two threads raced on them), and the external-capture term table and Park's shift-derivative memo were keyed by
whether the Coulomb table is used, not by its file.

One OpenMP thread per session (2 threads in all).  Needs the compiled engine (PYAZR_ENGINE_DIR) and the Coulomb table
tool beside the engine (AZURE2_COULOMB_TABLE_TOOL, else azure2-coulomb-table next to the module's build); without the
tool the T cases are skipped.  Skips cleanly without the engine.
Run:  python3 tests/pyazr/concurrent_sessions_test.py
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
REPEATS = 20
# label: (project folder, session options, pair radius edit)
CASES = {
    "A": ("13N", {}, None),
    "B": ("18O_p_a_thm", {"use_coulomb_table": False}, None),
    "N": ("13N", {"use_ec_table": False}, None),
    "M": ("13N", {"use_ec_table": False}, (1, 5.0)),
    "T": ("13N", {"coulomb_table": "@SECOND"}, None),
    "R": ("13N", {}, None),
    "S": ("18O_p_a_thm", {"use_coulomb_table": False}, None),
}
RELATIVE = {"R", "S"}  # the project's paths kept relative; the session writes its output files too
PAIRS = [("A", "B"), ("N", "M"), ("A", "T"), ("R", "S")]
SEQUENCES = [("A", "T", "A"), ("T", "A", "T")]


def table_tool():
    tool = os.environ.get("AZURE2_COULOMB_TABLE_TOOL")
    if tool:
        return tool if os.access(tool, os.X_OK) else None
    engine = os.environ.get("PYAZR_ENGINE_DIR")
    if engine:
        for cand in (os.path.join(os.path.dirname(os.path.abspath(engine)), "src", "azure2-coulomb-table"),
                     os.path.join(engine, "azure2-coulomb-table")):
            if os.access(cand, os.X_OK):
                return cand
    return None


def stage(tmp, label):
    """A copy of the label's project directory of its own (output files of one session are not another's)."""
    folder, _, radius = CASES[label]
    work = os.path.join(tmp, label)
    azr = os.path.join(work, folder + ".azr")
    if not os.path.isdir(work):
        shutil.copytree(os.path.join(ROOT, "tests", folder), work)
        for junk in ("output", "checks"):
            shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
            os.makedirs(os.path.join(work, junk))
        if radius:
            from pyazr import Project
            proj = Project.open(azr)
            proj.set_pair(radius[0], radius=radius[1])
            proj.save(azr)
        if label in RELATIVE:
            return azr
        # Absolute paths for the other cases, as written before relative paths were the session's own.
        with open(azr) as f:
            text = f.read()
        text = re.sub(r"^(output|checks)/", lambda m: os.path.join(work, m.group(1)) + "/", text, flags=re.M)
        text = re.sub(r"^(\s+file\s+)(data/)", lambda m: m.group(1) + os.path.join(work, "data") + "/", text, flags=re.M)
        with open(azr, "w") as f:
            f.write(text)
    return azr


def options(label, second_table):
    opts = dict(CASES[label][1])
    if opts.get("coulomb_table") == "@SECOND":
        opts["coulomb_table"] = second_table
    return opts


def evaluate(azr, opts, write=False):
    import numpy as np
    from pyazr import azure2
    files = []
    with azure2(azr, **opts) as s:
        x = np.asarray(s.params_rwa, float)
        parts = [np.concatenate([np.asarray(v, float).ravel() for v in s.calculate_rwa(x)]),
                 np.asarray(s.segment_chi2(x), float)]
        chi2, grad = s.chi2_and_grad(x)
        parts += [np.asarray([chi2], float), np.asarray(grad, float).ravel()]
        if write:
            out = s.output_dir
            for name in os.listdir(out):  # this session's files only: written afresh below
                if name.startswith("AZUREOut_") or name == "chiSquared.out":
                    os.remove(os.path.join(out, name))
            s.write_output_files(x)
            files = sorted(n for n in os.listdir(out) if n.startswith("AZUREOut_") or n == "chiSquared.out")
            if not files:
                raise RuntimeError(f"no output file written in {out}")
    h = hashlib.sha256()
    for p in parts:
        h.update(np.ascontiguousarray(p, dtype="<f8").tobytes())
    for name in files:
        h.update(name.encode())
        with open(os.path.join(out, name), "rb") as f:
            h.update(f.read())
    return h.hexdigest()[:16]


def worker(tmp, mode, labels, second_table):
    """mode 'seq': the sessions one after another; mode 'pair': two threads at once, REPEATS times."""
    sys.path.insert(0, ROOT)
    out = []
    if mode == "seq":
        for label in labels:
            out.append([label, evaluate(stage(tmp, label), options(label, second_table), label in RELATIVE)])
    else:
        azrs = {label: stage(tmp, label) for label in labels}
        for _ in range(REPEATS):
            barrier = threading.Barrier(len(labels))
            got, errors = {}, []

            def run(label):
                try:
                    barrier.wait()
                    got[label] = evaluate(azrs[label], options(label, second_table), label in RELATIVE)
                except Exception as err:  # reported below, not swallowed
                    errors.append(f"{label}: {type(err).__name__}: {err}")

            threads = [threading.Thread(target=run, args=(label,)) for label in labels]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            if errors:
                raise RuntimeError("; ".join(errors))
            out += [[label, got[label]] for label in labels]
    print("DIGESTS " + json.dumps(out))


def run(tmp, mode, labels, second_table):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    tmp = tempfile.mkdtemp(dir=tmp)  # fresh copies of the projects per process
    r = subprocess.run([sys.executable, os.path.abspath(__file__), "--worker", tmp, mode, ",".join(labels),
                        second_table or ""], env=env, capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("DIGESTS "):
            return json.loads(line[8:])
    raise RuntimeError(f"worker {mode} {labels} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def main():
    try:
        sys.path.insert(0, ROOT)
        import numpy  # noqa: F401
        from pyazr import azure2  # noqa: F401
    except Exception as err:
        print(f"skip: engine not available ({type(err).__name__}: {err})")
        return 77
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        second_table = None
        tool = table_tool()
        if tool:
            # A coarse table (fewer nodes per decade, l <= 4): other interpolation errors than the default file's, so a
            # value computed with one table and served to a session of the other shows in the bytes.
            second_table = os.path.join(tmp, "coarse_coulomb_table.bin")
            r = subprocess.run([tool, "build", second_table, "--lmax", "4", "--per-decade", "24", "--threads", "1"],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"FAIL  the table tool could not build the second table:\n{r.stdout[-1000:]}\n{r.stderr[-1000:]}")
                return 1
        else:
            print("note: no Coulomb table tool (AZURE2_COULOMB_TABLE_TOOL); the T cases are skipped")
        labels = [k for k in CASES if second_table or k != "T"]
        fresh = {}
        for label in labels:
            fresh[label] = run(tmp, "seq", (label,), second_table)[0][1]
            print(f"  alone  {label} ({CASES[label][0]}, {options(label, second_table) or 'defaults'}): {fresh[label]}")
        if second_table and fresh["T"] == fresh["A"]:
            failures.append("the second Coulomb table does not change 13N's numbers: the T cases would test nothing")
        cases = [("pair", p) for p in PAIRS] + [("seq", s) for s in SEQUENCES]
        for mode, group in cases:
            if not second_table and "T" in group:
                continue
            try:
                got = run(tmp, mode, group, second_table)
            except RuntimeError as err:  # a session that could not be built or evaluated (a crash of the worker)
                failures.append(f"{mode} {'-'.join(group)}: {str(err).splitlines()[0]} ({str(err).splitlines()[-1]})")
                continue
            bad = [(k, label, digest) for k, (label, digest) in enumerate(got) if digest != fresh[label]]
            print(f"  {mode:4s} {'-'.join(group):6s} {len(got)} evaluations, {len(bad)} differ from a fresh process")
            for k, label, digest in bad[:5]:
                failures.append(f"{mode} {'-'.join(group)}: evaluation {k + 1} ({label}) {digest} != {fresh[label]}")
    for f in failures:
        print("FAIL  " + f)
    print("FAILED" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    if len(sys.argv) == 6 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3], sys.argv[4].split(","), sys.argv[5] or None)
        sys.exit(0)
    sys.exit(main())
