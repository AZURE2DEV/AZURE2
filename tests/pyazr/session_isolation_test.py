#!/usr/bin/env python3
"""Two calculations in one process do not see each other (no hidden state).

Projects A and B with different compound nuclei, channel radii and channel
kinds: tests/13N (12C+p, elastic scattering and capture with external capture,
radius 4.6 fm) and tests/18O_p_a_thm (18O+p, Trojan-horse data on the shared
mesh, radius 5.6 fm), plus tests/hybrid_potential (a channel with a nuclear
potential) as a third project.  In one process the sessions run in the orders
A B A and B A B (and A C A, and A R A with R project A at another channel radius), each at 1 and at 3 OpenMP threads; every
evaluation of a project must give the same bytes every time it runs and the
same bytes as that project alone in a fresh process at the same thread count.
What is compared: the cross sections of every data point (calculate_rwa), the
chi-squared per segment, and the objective with its analytic gradient
(chi2_and_grad), as raw IEEE bytes.  The same numbers at 1 and at 3 threads are
checked too (results do not depend on the thread count).

A leak shows as an evaluation of A after B differing from A's first one: some
process-wide object (a static, a thread_local memo, a global pointer) kept a
value of B under a key A uses.

Each sequence runs in its own process (the thread count is fixed when the
engine starts); one process at a time, at most 3 threads.

Needs the compiled engine (PYAZR_ENGINE_DIR); skips cleanly without it.
Run:  python3 tests/pyazr/session_isolation_test.py
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PROJECTS = {"A": "13N", "B": "18O_p_a_thm", "C": "hybrid_potential", "R": "13N"}
# R: project A with the channel radius of its first pair set to 5.0 fm in the file (the same nucleus, the same
# objects' sizes, so a memo keyed by an address or by anything but the radius would hand A's values to R)
RADIUS = {"R": (1, 5.0)}
SEQUENCES = [("A",), ("B",), ("C",), ("R",), ("A", "B", "A"), ("B", "A", "B"), ("A", "C", "A"), ("A", "R", "A"),
             ("R", "A", "R")]
THREADS = (1, 3)


def stage(tmp, name):
    work = os.path.join(tmp, name)
    name = name.split("-")[0]
    azr = os.path.join(work, name + ".azr")
    if not os.path.isdir(work):
        shutil.copytree(os.path.join(ROOT, "tests", name), work)
        for junk in ("output", "checks"):
            shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
            os.makedirs(os.path.join(work, junk))
        label = {v + "-r": k for k, v in PROJECTS.items() if k in RADIUS}.get(os.path.basename(work))
        if label:
            from pyazr import Project
            proj = Project.open(azr)
            pair, radius = RADIUS[label]
            proj.set_pair(pair, radius=radius)
            proj.save(azr)
    return azr


def worker(tmp, sequence):
    """Run the sessions of `sequence` one after another in this process; print one digest per evaluation."""
    sys.path.insert(0, ROOT)
    import numpy as np
    from pyazr import azure2

    out = []
    for label in sequence:
        with azure2(stage(tmp, PROJECTS[label] + ("-r" if label in RADIUS else ""))) as s:
            x = np.asarray(s.params_rwa, float)
            parts = [np.concatenate([np.asarray(v, float).ravel() for v in s.calculate_rwa(x)]),
                     np.asarray(s.segment_chi2(x), float)]
            chi2, grad = s.chi2_and_grad(x)
            parts += [np.asarray([chi2], float), np.asarray(grad, float).ravel()]
            h = hashlib.sha256()
            for p in parts:
                h.update(np.ascontiguousarray(p, dtype="<f8").tobytes())
            out.append([label, h.hexdigest()[:16], int(sum(p.size for p in parts))])
    print("DIGESTS " + json.dumps(out))


def run(tmp, sequence, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    # A fresh copy of the projects per process: what one process wrote (output files, the EC integral file) is
    # not seen by the next; within a process a project's second session finds its first session's files.
    tmp = tempfile.mkdtemp(dir=tmp)
    r = subprocess.run([sys.executable, os.path.abspath(__file__), "--worker", tmp, ",".join(sequence)],
                       env=env, capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("DIGESTS "):
            return json.loads(line[8:])
    raise RuntimeError(f"worker {sequence} at {threads} threads failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def main():
    try:
        sys.path.insert(0, ROOT)
        import numpy  # noqa: F401
        import pyazr  # noqa: F401
        from pyazr import azure2  # noqa: F401
    except Exception as err:
        print(f"skip: engine not available ({type(err).__name__}: {err})")
        return 77
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        fresh = {}
        for threads in THREADS:
            for sequence in SEQUENCES:
                got = run(tmp, sequence, threads)
                print(f"  {threads} thread(s) {'-'.join(sequence):6s} " +
                      "  ".join(f"{label}:{digest}" for label, digest, _ in got))
                if len(sequence) == 1:
                    fresh[(sequence[0], threads)] = got[0][1]
                    continue
                for k, (label, digest, n) in enumerate(got):
                    want = fresh[(label, threads)]
                    if digest != want:
                        failures.append(f"{threads} thread(s), sequence {'-'.join(sequence)}: evaluation {k + 1} "
                                        f"({PROJECTS[label]}) differs from a fresh process ({digest} != {want})")
        for label in PROJECTS:
            if fresh[(label, 1)] != fresh[(label, 3)]:
                failures.append(f"{PROJECTS[label]}: 1 thread {fresh[(label, 1)]} != 3 threads {fresh[(label, 3)]}")
    for f in failures:
        print("FAIL  " + f)
    print("FAILED" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3].split(","))
        sys.exit(0)
    sys.exit(main())
