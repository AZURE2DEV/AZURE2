#!/usr/bin/env python3
"""Where rmfit imports pyazr from (rmfit/__init__.py AZURE2_ROOT), with rmfit inside an AZURE2 checkout or beside it.

  (a) rmfit at the top of an AZURE2 tree, no $AZURE2_ROOT: AZURE2_ROOT is that tree and pyazr comes from it;
  (b) the same rmfit with $AZURE2_ROOT naming ANOTHER tree that is later on PYTHONPATH (a campaign pinned to a frozen
      engine, the order the job templates write): pyazr comes from $AZURE2_ROOT, not from the tree rmfit sits in;
  (c) rmfit beside AZURE2 (a directory that is not an AZURE2 tree, as in a development archive), $AZURE2_ROOT set:
      pyazr from it; without $AZURE2_ROOT: ~/AZURE2;
  (d) the job scripts' import guard (rmfit first, then pyazr) passes in case (b).
Engine-free (fake pyazr packages), seconds.  Run:  python3 rmfit/tests/test_azure2_root.py"""
import os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))      # .../rmfit
FAILS = []


def check(name, ok, info=""):
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"  [{info}]"))
    if not ok:
        FAILS.append(name)


def fake_tree(path, name, azure2=True):
    os.makedirs(os.path.join(path, "pyazr"), exist_ok=True)
    open(os.path.join(path, "pyazr", "__init__.py"), "w").write(f"NAME = {name!r}\n")
    if azure2:
        os.makedirs(os.path.join(path, "src"), exist_ok=True)
        open(os.path.join(path, "src", "AZURE2.cpp"), "w").write("// marker\n")


def probe(pythonpath, azure2_root=None, code=None):
    env = {k: v for k, v in os.environ.items() if k not in ("AZURE2_ROOT", "PYTHONPATH")}
    env["PYTHONPATH"] = pythonpath
    if azure2_root:
        env["AZURE2_ROOT"] = azure2_root
    code = code or "import rmfit, pyazr; print(rmfit.__file__); print(rmfit.AZURE2_ROOT); print(pyazr.NAME)"
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=tempfile.gettempdir())
    return r.returncode, r.stdout.split("\n"), r.stderr


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_azure2_root_")
    try:
        A, B, arch = os.path.join(tmp, "A"), os.path.join(tmp, "B"), os.path.join(tmp, "archive")
        fake_tree(A, "A"); fake_tree(B, "B"); os.makedirs(arch)
        shutil.copytree(HERE, os.path.join(A, "rmfit"), ignore=shutil.ignore_patterns("tests", "attic", ".git", "__pycache__"))
        shutil.copytree(HERE, os.path.join(arch, "rmfit"), ignore=shutil.ignore_patterns("tests", "attic", ".git", "__pycache__"))

        rc, out, err = probe(A)
        check("(a) rmfit inside an AZURE2 tree: AZURE2_ROOT is that tree, pyazr from it",
              rc == 0 and out[0].startswith(A) and out[1] == A and out[2] == "A", err or out)
        rc, out, err = probe(os.pathsep.join((A, B)), azure2_root=B)
        check("(b) $AZURE2_ROOT names another tree later on the path: pyazr from $AZURE2_ROOT",
              rc == 0 and out[0].startswith(A) and out[1] == B and out[2] == "B", err or out)
        rc, out, err = probe(os.pathsep.join((arch, B)), azure2_root=B)
        check("(c) rmfit beside AZURE2, $AZURE2_ROOT set: pyazr from it", rc == 0 and out[1] == B and out[2] == "B", err or out)
        rc, out, err = probe(arch, code="import os, rmfit; print(rmfit.AZURE2_ROOT == os.path.expanduser('~/AZURE2'))")
        check("(c) rmfit beside AZURE2, no $AZURE2_ROOT: ~/AZURE2", rc == 0 and out[0] == "True", err or out)
        guard = ("import sqlite3, rmfit, pyazr; assert pyazr.__file__.startswith(rmfit.AZURE2_ROOT), pyazr.__file__; "
                 "print('import guard ok', pyazr.__file__)")
        rc, out, err = probe(os.pathsep.join((A, B)), azure2_root=B, code=guard)
        check("(d) the job scripts' import guard passes in case (b)", rc == 0 and "import guard ok" in out[0], err or out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all AZURE2_ROOT checks passed")


if __name__ == "__main__":
    main()
