#!/usr/bin/env python3
"""tools/classify.py: a file that cannot be read is never identical, not even to another that cannot be read.

Two runs of one project are laid out as tools/runall.sh writes them (CMP/<basis>/<tag>/.res-<dir>/<file>.azr/output/),
in a temporary directory: a file equal in both runs, a file that differs, and a file that is a dangling symbolic link
in both runs (unreadable).  The equal file is counted identical; the other two are listed and not counted.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSIFY = os.path.join(HERE, "..", "..", "tools", "classify.py")


def main():
    with tempfile.TemporaryDirectory() as cmp:
        for tag, value in (("ref", "1.0"), ("new", "2.0")):
            out = os.path.join(cmp, "brune", tag, "t", ".res-t", "t.azr", "output")
            os.makedirs(out)
            with open(os.path.join(out, "same.out"), "w") as f:
                f.write("1 2 3\n")
            with open(os.path.join(out, "moved.out"), "w") as f:
                f.write(f"1 {value}\n")
            os.symlink(os.path.join(cmp, "missing"), os.path.join(out, "gone.out"))
        r = subprocess.run([sys.executable, CLASSIFY, cmp, "ref", "new", "brune"], capture_output=True, text=True)
        print(r.stdout, end="")
        failures = []
        if r.returncode != 0:
            failures.append(f"exit {r.returncode}: {r.stderr}")
        if "## brune: 3 files, 1 identical" not in r.stdout:
            failures.append("expected 3 files, 1 identical")
        if not any(line.endswith("gone.out\tunreadable-both") for line in r.stdout.splitlines()):
            failures.append("the unreadable file is not listed as unreadable-both")
        if not any("moved.out\tabsrel=" in line for line in r.stdout.splitlines()):
            failures.append("the differing file is not listed")
    for f in failures:
        print("FAIL:", f)
    print("PASS" if not failures else f"{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
