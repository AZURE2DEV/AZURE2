#!/usr/bin/env python3
"""tools/runall.sh and tools/run_one.sh never leave the results of an earlier run under a tag beside new ones.

A fake engine (a shell script that writes output/value.out with $FAKE_VALUE and checks/c.out) runs a tests tree of two
folders, a/ and b/, in the Brune basis:
  1. a fresh tag: both folders' results, exit 0;
  2. the same tag again: refused (exit 2), the results untouched;
  3. with REPLACE_TAG=1 after b/ was renamed c/: a's value is the new one, .res-b is gone, .res-c is there;
  4. run_one.sh onto a tag whose output/ folder of a file exists and is not empty: replaced (before, mv refused the
     folder silently and the old output stayed beside a new run.log);
  5. a result that cannot be moved (a file where the result folder goes): run_one.sh prints it and exits 1, and
     runall.sh ends with "FAILED" and exit 1, not "runs done".
"""
import os
import shutil
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "tools")
ENGINE = """#!/bin/bash
cat > /dev/null
[ -n "${FAKE_BLOCK:-}" ] && mkdir -p "$(dirname "$FAKE_BLOCK")" && touch "$FAKE_BLOCK"
echo "$FAKE_VALUE" > output/value.out
echo checks > checks/c.out
echo "fake engine $*"
"""

failures = []


def check(name, ok, detail=""):
    print("  %s  %s%s" % ("ok  " if ok else "FAIL", name, "" if ok or not detail else "  -- " + detail))
    if not ok:
        failures.append(name)


def main():
    with tempfile.TemporaryDirectory() as work:
        tests, engine = os.path.join(work, "tests"), os.path.join(work, "engine.sh")
        for d in ("a", "b"):
            os.makedirs(os.path.join(tests, d))
            open(os.path.join(tests, d, d + ".azr"), "w").close()
        with open(engine, "w") as f:
            f.write(ENGINE)
        os.chmod(engine, 0o755)
        res = os.path.join(work, "cmp", "brune", "t")

        # The caller's REPLACE_TAG is not this test's (a run of all projects exports it), and no tool reads the
        # caller's standard input (the fake engine's `cat` would wait on an open socket).
        base = {k: v for k, v in os.environ.items() if k != "REPLACE_TAG"}

        def runall(value, extra=None):
            env = dict(base, SRC_TESTS=tests, WORK=work, JOBS="1", AZURE2_COULOMB_TABLE=engine, FAKE_VALUE=value)
            env.update(extra or {})
            return subprocess.run(["bash", os.path.join(TOOLS, "runall.sh"), engine, "t", ""], env=env,
                                  stdin=subprocess.DEVNULL, capture_output=True, text=True)

        def value(folder, azr):
            path = os.path.join(res, ".res-" + folder, azr, "output", "value.out")
            return open(path).read().strip() if os.path.exists(path) else None

        r = runall("1")
        check("1. a fresh tag runs", r.returncode == 0 and value("a", "a.azr") == "1" and value("b", "b.azr") == "1",
              r.stdout + r.stderr)
        r = runall("2")
        check("2. an existing tag is refused", r.returncode == 2 and "exists" in r.stderr and value("a", "a.azr") == "1",
              "exit %d: %s" % (r.returncode, r.stderr))
        shutil.move(os.path.join(tests, "b"), os.path.join(tests, "c"))
        os.rename(os.path.join(tests, "c", "b.azr"), os.path.join(tests, "c", "c.azr"))
        r = runall("3", {"REPLACE_TAG": "1"})
        check("3. REPLACE_TAG=1 replaces the tag: new values, the renamed folder's old results gone",
              r.returncode == 0 and value("a", "a.azr") == "3" and value("c", "c.azr") == "3" and
              not os.path.exists(os.path.join(res, ".res-b")), "exit %d, a %s, .res-b %s" %
              (r.returncode, value("a", "a.azr"), os.path.exists(os.path.join(res, ".res-b"))))
        env = dict(base, AZURE2_COULOMB_TABLE=engine, FAKE_VALUE="4")
        r = subprocess.run(["bash", os.path.join(TOOLS, "run_one.sh"), tests, os.path.join(work, "cmp"), "t", engine, "",
                            "a"], env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True)
        check("4. run_one.sh replaces an existing, non-empty output folder", r.returncode == 0 and
              value("a", "a.azr") == "4", "exit %d, value %s: %s" % (r.returncode, value("a", "a.azr"), r.stderr))
        block = os.path.join(work, "cmp", "brune", "u", ".res-a")  # a file where the results of a/ go
        env = dict(base, SRC_TESTS=tests, WORK=work, JOBS="1", AZURE2_COULOMB_TABLE=engine, FAKE_VALUE="5",
                   FAKE_BLOCK=block)
        r = subprocess.run(["bash", os.path.join(TOOLS, "runall.sh"), engine, "u", ""], env=env,
                           stdin=subprocess.DEVNULL, capture_output=True, text=True)
        check("5. a result that cannot be moved: printed, runall.sh exits 1 with FAILED",
              r.returncode == 1 and "could not be moved" in r.stderr and "FAILED" in r.stderr and
              "runs done u\n" not in r.stdout, "exit %d: %s %s" % (r.returncode, r.stdout, r.stderr))
    print("FAIL" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
