#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): results do not depend on the number of threads.

T1 and T3 are calculated from data (menu 1), extrapolated (menu 3) and, for the reaction rate, integrated (menu 5)
with OMP_NUM_THREADS = 1 and 4, and fitted (menu 2, MIGRAD on the analytic gradient) from their start models with 1
and 4 threads.  Every output file (cross sections, A_y, phase shifts, S factors, chi^2, scattering lengths, rates)
must agree to 1e-9 relative, and the fitted parameter vectors (param.sav) to 1e-12.  The test also reports which
files are not byte-for-byte identical.

With the c1cc92c7 build every file is byte-for-byte identical except T1's AZUREOut_aa=1_R=3.out (capture to the 3/2-
state), which differs between 1 and 4 threads by up to 3.6e-11 relative (one unit in the 11th printed digit), in the
calculation as well as after the fit, and not in every run (two 4-thread runs differ by 1e-14 absolute); other
files occasionally differ only in the sign of a printed zero.  Today's
engine is therefore thread-independent to ~1e-11, not bit for bit; bit-for-bit independence by construction is the
Phase 6 gate (ROADMAP.md), which will tighten this check to identity.  A data race, a reduction whose order depends on the thread schedule, or a cache filled by
whichever thread arrives first shows up here as a difference in the last digits or worse.  The login node of the
machine this was written on had one core; four threads on one core still interleave, but run it on a multi-core
machine for the full check (CI does).

    python3 tests/pyazr/upgrade_thread_independence_test.py   (binary: $AZURE2_BIN, $AZURE2_BINARY or newest)
"""
import filecmp
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "upgrade", "lib"))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    import upgrade_project as up
    binary = up.find_binary()
except (ImportError, SystemExit) as e:
    print(f"skip: {e}")
    sys.exit(77)   # ctest SKIP_RETURN_CODE (tests/pyazr/CMakeLists.txt)

RATE = {"T1": (1, 2), "T3": (2, 1)}


def numeric_difference(a, b):
    """Largest relative difference of the numbers of two text files; problems if their shapes or words differ."""
    if not os.path.exists(b):
        return float("inf"), ["missing"]
    la, lb = [x.split() for x in open(a)], [x.split() for x in open(b)]
    if len(la) != len(lb):
        return float("inf"), ["line count"]
    worst, problems = 0.0, []
    for x, y in zip(la, lb):
        if len(x) != len(y):
            return float("inf"), ["field count"]
        for u, v in zip(x, y):
            try:
                fu, fv = float(u.rstrip(",")), float(v.rstrip(","))
            except ValueError:
                if u != v:
                    problems.append(u)
                continue
            if max(abs(fu), abs(fv)) > 1e-200:
                worst = max(worst, abs(fu - fv) / max(abs(fu), abs(fv)))
    return worst, problems


SKIP = {"intEC.dat", "intEC.dat.sig", "intEC.extrap", "intEC.extrap.sig", "param.sav", "fitStatus.out"}


def run_all(project, threads, top):
    """Outputs of menus 1, 3, 5 in top/<mode>/output, and of the fit in top/fit/output."""
    for mode, answers, azr_name in (("calc", up.CALC, project), ("extrap", up.EXTRAP, project),
                                    ("rate", up.rate_answers(*RATE[project]), project),
                                    ("fit", up.FIT, project + "_start")):
        work = os.path.join(top, mode)
        os.makedirs(work)
        up.project_copy(project, work)
        flags = ["--use-gradient"] if mode == "fit" else []
        up.run(binary, work, azr_name + ".azr", answers, flags=flags, threads=threads)


top = tempfile.mkdtemp(prefix="upgrade_threads.")
try:
    for project in ("T1", "T3"):
        a, b = os.path.join(top, project, "1"), os.path.join(top, project, "4")
        run_all(project, 1, a)
        run_all(project, 4, b)
        for mode in ("calc", "extrap", "rate", "fit"):
            oa, ob = os.path.join(a, mode, "output"), os.path.join(b, mode, "output")
            files = sorted(f for f in os.listdir(oa) if f not in SKIP)
            worst, bad, notbitwise = 0.0, [], []
            for f in files:
                w, problems = numeric_difference(os.path.join(oa, f), os.path.join(ob, f))
                worst = max(worst, w)
                if problems or w > 1e-9:
                    bad.append(f)
                elif not filecmp.cmp(os.path.join(oa, f), os.path.join(ob, f), shallow=False):
                    notbitwise.append(f)
            check(f"{project} {mode}: {len(files)} output files equal to 1e-9, 1 vs 4 threads", not bad,
                  f"largest relative difference {worst:.1e}" + (f"; differ: {', '.join(bad)}" if bad else ""))
            if notbitwise:
                print(f"  info  {project} {mode}: not byte-for-byte identical: {', '.join(notbitwise)}")
        pa = up.read_param_file(os.path.join(a, "fit", "output", "param.sav"))
        pb = up.read_param_file(os.path.join(b, "fit", "output", "param.sav"))
        worst = max(abs(pa[k][0] - pb[k][0]) / max(abs(pa[k][0]), 1e-300) for k in pa)
        check(f"{project} fit: parameters equal, 1 vs 4 threads", worst <= 1e-12, f"largest relative difference {worst:.1e}")
finally:
    shutil.rmtree(top, ignore_errors=True)

print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
