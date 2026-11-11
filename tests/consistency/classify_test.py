#!/usr/bin/env python3
"""tools/classify.py: a file that cannot be read is never identical, not even to another that cannot be read; the
class lab-frame-rounding accepts a last-digit change only in a listed project and nothing beyond the limit.

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
    failures += lab_frame_class()
    for f in failures:
        print("FAIL:", f)
    print("PASS" if not failures else f"{len(failures)} failure(s)")
    return 1 if failures else 0


def lab_frame_class():
    """A listed project (T2) with one file one printed digit apart (accepted) and one 3e-10 apart at 10 digits (beyond
    the limit, listed); an unlisted project with a last-digit change (listed); without the option nothing is accepted."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("classify", CLASSIFY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    failures = []
    if ".res-upgrade__T2/T2.azr" not in mod.LAB_FRAME_PROJECTS:
        return ["T2 is not in LAB_FRAME_PROJECTS"]
    with tempfile.TemporaryDirectory() as cmp:
        files = {("upgrade", ".res-upgrade__T2", "T2.azr", "flip.out"): ("2.1931573490e+00", "2.1931573489e+00"),
                 ("upgrade", ".res-upgrade__T2", "T2.azr", "far.out"): ("2.1931573490e+00", "2.1931573487e+00"),
                 ("other", ".res-other", "other.azr", "flip.out"): ("2.1931573490e+00", "2.1931573489e+00")}
        for tag, which in (("ref", 0), ("new", 1)):
            for (top, res, azr, name), values in files.items():
                out = os.path.join(cmp, "brune", tag, res, azr, "output")
                os.makedirs(out, exist_ok=True)
                with open(os.path.join(out, name), "w") as f:
                    f.write(f"1.0 {values[which]}\n")
        r = subprocess.run([sys.executable, CLASSIFY, cmp, "ref", "new", "brune", "--class=lab-frame-rounding"],
                           capture_output=True, text=True)
        print(r.stdout, end="")
        lines = r.stdout.splitlines()
        if "## brune: 3 files, 0 identical, 1 lab-frame rounding" not in r.stdout:
            failures.append("class: expected 3 files, 0 identical, 1 lab-frame rounding")
        if not any("T2.azr/output/flip.out\tlab-frame-rounding" in l for l in lines):
            failures.append("class: the last-digit change of T2 is not in the class")
        if not any("T2.azr/output/far.out\tabsrel=" in l for l in lines):
            failures.append("class: the change beyond the limit is not listed as differing")
        if not any("other.azr/output/flip.out\tabsrel=" in l for l in lines):
            failures.append("class: the change of an unlisted project is not listed as differing")
        r = subprocess.run([sys.executable, CLASSIFY, cmp, "ref", "new", "brune"], capture_output=True, text=True)
        if "lab-frame" in r.stdout or "## brune: 3 files, 0 identical" not in r.stdout:
            failures.append("without --class=lab-frame-rounding a change was accepted")
    return failures


if __name__ == "__main__":
    sys.exit(main())
