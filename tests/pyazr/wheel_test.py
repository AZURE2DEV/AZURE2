#!/usr/bin/env python3
"""pyazr as a wheel (pyproject.toml, scikit-build-core; docs/upgrade/pyazr3.md section 2.9).

  1. a wheel is built from the repository with pip (``pip wheel ROOT --no-deps --no-build-isolation``, the build
     directory and every scratch file under this test's folder, pip kept off the network), and the repository's
     ``git status`` is the same before and after (the build writes nothing into the source);
  2. its files are the expected list: the package's files (pyazr/datafiles.py: every .py file and every data file
     DATA_FILES names, less NOT_SHIPPED), the engine module, the program AZURE2 in the scripts, the metadata with
     the console script ``pyazr`` and the licence files (LICENSE, the SRIM notice, plotly's); a missing or an extra
     file fails;
  3. the wheel is **unpacked** (nothing is installed) and a fresh interpreter imports it with only that folder on
     its path and the repository out of reach (``pyazr.__file__`` is checked): the objective of T2 at the start
     parameters equals the in-tree build's bit for bit, as do the stopping of protons in carbon; the console script
     ``pyazr`` (the entry point of the metadata) answers; tests/pyazr/wheel_smoke.py (the wheel workflow's test
     command) passes with the wheel's AZURE2 on the PATH;
  4. the wheel's AZURE2 on T1 writes output files identical to those of the in-tree program (AZURE2_BIN).

PYAZR_WHEEL=<file> tests a wheel already built instead of building one.  CMAKE_BUILD_PARALLEL_LEVEL (default 2) sets
the compile jobs.  Scratch: ./wheel_test (ctest: the build's tests/pyazr).  Exits 77 (skipped) without the in-tree
engine module, pip or scikit-build-core.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import pyazr
    from pyazr import _azure2  # noqa: F401
    from pyazr import datafiles
    import scikit_build_core  # noqa: F401
    import pip  # noqa: F401
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

failures = []
OFFLINE = {"PIP_NO_INDEX": "1", "http_proxy": "http://127.0.0.1:9", "https_proxy": "http://127.0.0.1:9",
           "HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9"}


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:1500]}"))
    if not ok:
        failures.append(name)


def git_status():
    p = subprocess.run(["git", "-c", "safe.directory=*", "status", "--porcelain", "--untracked-files=all"],
                       cwd=ROOT, capture_output=True, text=True)
    return p.stdout if p.returncode == 0 else None


def build_wheel(scratch):
    out = os.path.join(scratch, "dist")
    for sub in ("dist", "build", "pip-cache", "tmp"):
        shutil.rmtree(os.path.join(scratch, sub), ignore_errors=True)
        os.makedirs(os.path.join(scratch, sub))
    env = dict(os.environ, PIP_CACHE_DIR=os.path.join(scratch, "pip-cache"), TMPDIR=os.path.join(scratch, "tmp"),
               CMAKE_BUILD_PARALLEL_LEVEL=os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "2"), **OFFLINE)
    env.pop("PYAZR_ENGINE_DIR", None)
    with open(os.path.join(scratch, "pip-wheel.log"), "w") as log:
        p = subprocess.run([sys.executable, "-m", "pip", "wheel", ROOT, "--no-deps", "--no-build-isolation",
                            "-w", out, "-C", "build-dir=" + os.path.join(scratch, "build")],
                           stdout=log, stderr=subprocess.STDOUT, env=env)
    wheels = glob.glob(os.path.join(out, "pyazr-*.whl"))
    check("pip wheel: exit 0, one wheel", p.returncode == 0 and len(wheels) == 1,
          f"exit {p.returncode}; see {os.path.join(scratch, 'pip-wheel.log')}")
    return wheels[0] if len(wheels) == 1 else None


def expected_files(version):
    info = f"pyazr-{version}.dist-info"
    files = {"pyazr/" + rel for rel in datafiles.package_files(os.path.join(ROOT, "pyazr"))}
    files |= {"pyazr/_azure2" + sysconfig.get_config_var("EXT_SUFFIX"), f"pyazr-{version}.data/scripts/AZURE2"}
    files |= {f"{info}/{name}" for name in ("METADATA", "WHEEL", "RECORD", "entry_points.txt")}
    files |= {f"{info}/licenses/{name}" for name in ("LICENSE", "stopping/data/srim/README.md",
                                                     "pyazr/ui/static/vendor/plotly-LICENSE.txt")}
    return files


def child_env(unpacked):
    """An environment in which only the unpacked wheel provides pyazr (no PYAZR_ENGINE_DIR, no repository)."""
    env = {k: v for k, v in os.environ.items() if k not in ("PYAZR_ENGINE_DIR", "PYTHONPATH", "PYTHONHOME")}
    env.update(PYTHONPATH=unpacked, OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    return env


CHILD = r"""
import json, sys
import pyazr
from pyazr import stopping
f = pyazr.Project.open(sys.argv[1]).session(use_coulomb_table=False).fit_objective()
s = stopping.stopping(stopping.Material("C"), "p", [0.25, 0.5, 1.0, 2.0])
print(json.dumps({"file": pyazr.__file__, "version": pyazr.__version__, "objective": float(f(f.x0)).hex(),
                  "stopping": [float(v).hex() for v in s["total"]]}))
"""


def in_tree(t2):
    from pyazr import stopping
    f = pyazr.Project.open(t2).session(use_coulomb_table=False).fit_objective()
    s = stopping.stopping(stopping.Material("C"), "p", [0.25, 0.5, 1.0, 2.0])
    return {"objective": float(f(f.x0)).hex(), "stopping": [float(v).hex() for v in s["total"]]}


def run_program(binary, folder):
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), folder)
    for sub in ("output", "checks"):
        shutil.rmtree(os.path.join(folder, sub), ignore_errors=True)
        os.makedirs(os.path.join(folder, sub))
    p = subprocess.run([binary, "--no-gui", "--no-readline", "--no-coulomb-table", "T1.azr"], input="1\n\n\n7\n",
                       capture_output=True, text=True, cwd=folder, env=dict(os.environ, OMP_NUM_THREADS="1"),
                       timeout=600)
    return p.returncode, os.path.join(folder, "output")


def same_tree(a, b):
    names = sorted(os.listdir(a)), sorted(os.listdir(b))
    if names[0] != names[1]:
        return f"files differ: {names}"
    for n in names[0]:
        with open(os.path.join(a, n), "rb") as fa, open(os.path.join(b, n), "rb") as fb:
            if fa.read() != fb.read():
                return f"{n} differs"
    return ""


def main(scratch):
    before = git_status()
    wheel = os.environ.get("PYAZR_WHEEL") or build_wheel(scratch)
    if not wheel:
        return
    check("the repository's git status is unchanged by the build", git_status() == before,
          f"before:\n{before}\nafter:\n{git_status()}")
    version = pyazr.__version__
    with zipfile.ZipFile(wheel) as z:
        names = {n for n in z.namelist() if not n.endswith("/")}
        expected = expected_files(version)
        check(f"the wheel's files are the expected {len(expected)}", names == expected,
              f"missing: {sorted(expected - names)}; extra: {sorted(names - expected)}")
        entry = z.read(f"pyazr-{version}.dist-info/entry_points.txt").decode() if \
            f"pyazr-{version}.dist-info/entry_points.txt" in names else ""
        check("the console script pyazr = pyazr.cli:main", "pyazr = pyazr.cli:main" in entry, entry)
        tag = z.read(f"pyazr-{version}.dist-info/WHEEL").decode()
        check("the wheel is platform-specific (Root-Is-Purelib: false)", "Root-Is-Purelib: false" in tag, tag)
        unpacked = os.path.join(scratch, "unpacked")
        shutil.rmtree(unpacked, ignore_errors=True)
        z.extractall(unpacked)
    program = os.path.join(unpacked, f"pyazr-{version}.data", "scripts", "AZURE2")
    if os.path.isfile(program):
        os.chmod(program, 0o755)  # zipfile does not restore the mode; pip does
    print(f"  wheel: {os.path.basename(wheel)}, {os.path.getsize(wheel)} bytes, {len(names)} files")

    run_dir = os.path.join(scratch, "run")
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(run_dir)
    t2 = os.path.join(run_dir, "T2")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T2"), t2)
    os.makedirs(os.path.join(t2, "output"), exist_ok=True)
    p = subprocess.run([sys.executable, "-s", "-c", CHILD, os.path.join(t2, "T2.azr")], capture_output=True,
                       text=True, cwd=run_dir, env=child_env(unpacked))
    try:
        got = json.loads(p.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        check("the unpacked wheel runs a session on T2", False, p.stdout[-800:] + p.stderr[-1500:])
        return
    check("the unpacked wheel is the pyazr imported", got["file"].startswith(unpacked + os.sep), got["file"])
    t2_tree = os.path.join(run_dir, "T2-tree")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T2"), t2_tree)
    os.makedirs(os.path.join(t2_tree, "output"), exist_ok=True)
    want = in_tree(os.path.join(t2_tree, "T2.azr"))
    check(f"T2: the objective at x0 equals the in-tree build's bit for bit ({float.fromhex(want['objective'])!r})",
          got["objective"] == want["objective"], (got["objective"], want["objective"]))
    check("stopping of p in C at 0.25-2 MeV equals the in-tree build's bit for bit", got["stopping"] ==
          want["stopping"], (got["stopping"], want["stopping"]))
    p = subprocess.run([sys.executable, "-s", "-c", "import sys; from pyazr.cli import main; sys.exit(main())",
                        "describe", "--json"], capture_output=True, text=True, cwd=run_dir, env=child_env(unpacked))
    try:
        doc = json.loads(p.stdout)
        ok = p.returncode == 0 and doc["ok"] and doc["result"]["engine"]["module"].startswith(unpacked + os.sep)
    except (ValueError, KeyError, TypeError):
        ok = False
    check("the console script's entry point (describe --json) answers from the wheel", ok, p.stdout[-600:] +
          p.stderr[-600:])
    env = child_env(unpacked)
    env["PATH"] = os.path.dirname(program) + os.pathsep + env.get("PATH", "")
    p = subprocess.run([sys.executable, "-s", os.path.join(HERE, "wheel_smoke.py"), ROOT], capture_output=True,
                       text=True, cwd=run_dir, env=env)
    print("\n".join("    smoke:" + line for line in p.stdout.splitlines()))
    check("tests/pyazr/wheel_smoke.py (the workflow's test command) passes on the unpacked wheel", p.returncode == 0,
          p.stderr[-800:])

    binary = os.environ.get("AZURE2_BIN")
    if binary and os.path.isfile(program):
        code_w, out_w = run_program(program, os.path.join(run_dir, "T1-wheel"))
        code_t, out_t = run_program(binary, os.path.join(run_dir, "T1-tree"))
        difference = same_tree(out_w, out_t)
        check("AZURE2 of the wheel on T1: exit 0, output files identical to the in-tree program's",
              code_w == 0 and code_t == 0 and not difference, (code_w, code_t, difference))
    else:
        check("the wheel's AZURE2 and AZURE2_BIN are there", False, (program, binary))


if __name__ == "__main__":
    scratch = os.path.join(os.getcwd(), "wheel_test")
    os.makedirs(scratch, exist_ok=True)
    main(scratch)
    if failures:
        print(f"{len(failures)} failure(s)")
        sys.exit(1)
    print("all checks passed")
