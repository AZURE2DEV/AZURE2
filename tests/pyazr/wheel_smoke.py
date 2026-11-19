#!/usr/bin/env python3
"""The test of an installed (or unpacked) pyazr wheel: the wheel workflow's test command (pyproject.toml,
[tool.cibuildwheel] test-command) and a part of tests/pyazr/wheel_test.py.

    python tests/pyazr/wheel_smoke.py REPOSITORY

REPOSITORY is the source checkout, read for the test projects only (tests/upgrade/T1-T3); pyazr must come from
somewhere else (the check fails if it is imported from REPOSITORY).  Checks: pyazr.describe() (the engine module
next to the package, the SRIM citation, the mass table); the objective of T1, T2, T3 at the start parameters
(objective()(x0)) against the values stored below, to 1e-7 relative (the platforms differ: see OBJECTIVE_TOLERANCE);
the stopping of protons in carbon; the web page's files found; ``python -m pyazr --help``; the program AZURE2 on T1
(on the PATH, or AZURE2_BIN) against the stored total chi-squared of chiSquared.out (6 digits, as printed).  Writes
only to a scratch folder in the working directory; exits 0 when every check passes.
"""
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile

#: objective()(x0) of the test projects (the engine's exact Coulomb routines, no table), Linux x86-64.  T1 and T3
#: re-recorded on 2026-11-18 for the CODATA 2018 m_p of the M1 capture and m_e of the Fermi function (were
#: 482.27523957308824, 1.6e-9 higher, and 518.713434718418, 1.9e-8 lower).
OBJECTIVE_AT_START = {
    "T1": 482.27523879169445,
    "T2": 10910.13380363804,
    "T3": 518.7134446833464,
}
#: How far a platform's objective may be from the stored Linux x86-64 value (relative).  It was 1e-9 until the
#: wheels' first runs on other platforms (2026-11-16, GitHub run 43859564660), which gave
#:   Linux aarch64 (GCC, manylinux_2_28)   T2 10910.133783195419   1.9e-9 from the stored value
#:   macOS arm64 (Apple clang)             T3 518.7134335819021    2.2e-9
#:   macOS x86-64 (Apple clang)            T3 518.7134336276669    2.1e-9
#: and the other projects within 1e-9.  The engine's stated spread between compilers is near 1e-5 (CLAUDE.md,
#: "Determinism"); this test asks whether the wheel computes with the right engine and data files, for which
#: 1e-7 (50 times the largest difference seen) is ample.  Why the platforms differ at 2e-9 was not traced.
OBJECTIVE_TOLERANCE = 1e-7
#: Total-Chi-Squared of the program's chiSquared.out on T1 (menu 1, calculate with the data).
PROGRAM_T1_TOTAL = 482.275

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)
        if os.environ.get("GITHUB_ACTIONS"):
            # One annotation per failed check: the run's log needs a login, its annotations do not.
            text = f"{name} -- {str(detail)[:600]}".replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
            print(f"::error title=wheel_smoke ({platform.machine()})::{text}", flush=True)


def objective_at_start(path):
    """objective()(x0) of a session on the project at ``path``, with the exact Coulomb routines."""
    import pyazr
    session = pyazr.Project.open(path).session(use_coulomb_table=False)
    try:
        f = session.objective()
        return float(f(f.x0))
    finally:
        session.close()


def main(repo, scratch):
    repo = os.path.abspath(repo)
    import pyazr
    where = os.path.dirname(os.path.abspath(pyazr.__file__))
    check("pyazr is not the repository's", not where.startswith(repo + os.sep), where)
    d = pyazr.describe()
    check("describe(): the engine module is the package's", os.path.dirname(d["engine"]["module"]) == where, d["engine"])
    check("describe(): the SRIM citation", "Ziegler" in d["stopping"]["citation"], d["stopping"])
    check("describe(): the mass table found", os.path.isfile(d["masses"].get("file", "")), d["masses"])
    for name, stored in OBJECTIVE_AT_START.items():
        project = os.path.join(scratch, name)
        shutil.copytree(os.path.join(repo, "tests", "upgrade", name), project)
        os.makedirs(os.path.join(project, "output"), exist_ok=True)
        value = objective_at_start(os.path.join(project, name + ".azr"))
        check(f"{name}: the objective at the start parameters, {value!r}, is the stored {stored!r} "
              f"({OBJECTIVE_TOLERANCE:g})",
              stored is not None and math.isclose(value, stored, rel_tol=OBJECTIVE_TOLERANCE), value)
    s = pyazr.stopping.stopping(pyazr.stopping.Material("C"), "p", [0.5, 1.0])
    check("stopping of protons in carbon: positive, finite", all(v > 0 and math.isfinite(v) for v in s["total"]), s)
    from pyazr.ui import routes
    for name in ("index.html", "vendor/plotly-cartesian-4.1.2.min.js"):
        try:
            full, _ = routes.static_file(name)
            check(f"the web page's {name} is found", os.path.isfile(full) and full.startswith(where), full)
        except Exception as e:
            check(f"the web page's {name} is found", False, e)
    p = subprocess.run([sys.executable, "-m", "pyazr", "--help"], capture_output=True, text=True, cwd=scratch)
    check("python -m pyazr --help: exit 0, the commands listed", p.returncode == 0 and "coulomb-table" in p.stdout,
          p.stderr)
    # The program: AZURE2_BIN, else the console script on the PATH (it runs the package's pyazr/bin/AZURE2), else that
    # file itself (an unpacked wheel has no console script).
    from pyazr import program as packaged
    check("the package holds the program (pyazr/bin/AZURE2)", packaged.path() is not None, where)
    program = os.environ.get("AZURE2_BIN") or shutil.which("AZURE2") or packaged.path()
    check("the program AZURE2 is installed", program is not None, os.environ.get("PATH"))
    if program:
        folder = os.path.join(scratch, "T1-program")
        shutil.copytree(os.path.join(repo, "tests", "upgrade", "T1"), folder)
        for sub in ("output", "checks"):
            shutil.rmtree(os.path.join(folder, sub), ignore_errors=True)
            os.makedirs(os.path.join(folder, sub))
        # The package's file run directly loads its DLLs (a Windows wheel: pyazr.libs) from the PATH, which the
        # console script (pyazr.program.main) sets the same way.
        environment = dict(os.environ)
        libs = packaged._libs_directory()
        if os.name == "nt" and libs:
            environment["PATH"] = libs + os.pathsep + environment.get("PATH", "")
        p = subprocess.run([program, "calculate", "--no-gui", "--no-coulomb-table", "T1.azr"], stdin=subprocess.DEVNULL,
                           capture_output=True, text=True, cwd=folder, timeout=600, env=environment)
        total = None
        try:
            with open(os.path.join(folder, "output", "chiSquared.out")) as f:
                for line in f:
                    if line.startswith("Total-Chi-Squared:"):
                        total = float(line.split()[1])
        except OSError:
            pass
        check(f"AZURE2 on T1: exit 0, Total-Chi-Squared {total} is the stored {PROGRAM_T1_TOTAL}",
              p.returncode == 0 and total is not None and PROGRAM_T1_TOTAL is not None
              and math.isclose(total, PROGRAM_T1_TOTAL, rel_tol=1e-5),
              f"exit {p.returncode} of {program}; " + p.stdout[-400:] + p.stderr[-400:])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    with tempfile.TemporaryDirectory(prefix="wheel_smoke.", dir=os.getcwd()) as scratch:
        main(sys.argv[1], scratch)
    if failures:
        print(f"{len(failures)} failure(s)")
        sys.exit(1)
    print("the wheel works")
