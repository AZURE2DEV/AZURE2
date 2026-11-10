#!/usr/bin/env python3
"""The expression tree of every data set and extrapolation of every project of tests/ and examples/
(docs/upgrade/data_model.md section 3, step S10), as `AZURE2 --print-predictions` prints it, against the committed
table tests/upgrade/predictions.txt.

The tree is what the model pass evaluates at each point (src/Predictions.cpp, experiment/Prediction.h): a change in
the builder, in the set-up that feeds it (components, pairs, samples, transforms) or in a project shows here line by
line.  Each project is set up in its own directory (data files are read relative to the working directory).

  python3 tests/upgrade/predictions_test.py --binary build/src/AZURE2           # check
  python3 tests/upgrade/predictions_test.py --binary build/src/AZURE2 --write   # regenerate the table
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TABLE = os.path.join(HERE, "predictions.txt")


def projects():
    found = []
    for top in ("tests", "examples"):
        for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, top)):
            dirnames[:] = sorted(d for d in dirnames if d not in ("output", "checks"))
            for name in sorted(filenames):
                if name.endswith(".azr"):
                    found.append(os.path.relpath(os.path.join(dirpath, name), ROOT))
    return sorted(found)


def table(binary):
    out = []
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for rel in projects():
        path = os.path.join(ROOT, rel)
        r = subprocess.run([binary, "--print-predictions", os.path.basename(path)], cwd=os.path.dirname(path),
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env, timeout=600)
        out.append(f"== {rel}")
        out.extend(r.stdout.decode("utf-8", "replace").splitlines())
        if r.returncode != 0:
            out.append(f"exit {r.returncode}")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--write", action="store_true", help="rewrite tests/upgrade/predictions.txt")
    a = ap.parse_args()
    text = table(os.path.abspath(a.binary))
    if a.write:
        with open(TABLE, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"wrote {TABLE}: {text.count(chr(10) + '== ') + 1} projects")
        return 0
    with open(TABLE, encoding="utf-8") as f:
        golden = f.read()
    if text == golden:
        print(f"predictions: {golden.count('== ')} projects, every tree as committed")
        return 0
    new, old = text.splitlines(), golden.splitlines()
    shown = 0
    for k in range(max(len(new), len(old))):
        a_, b_ = (old[k] if k < len(old) else "<none>"), (new[k] if k < len(new) else "<none>")
        if a_ != b_:
            print(f"line {k + 1}:\n  committed: {a_}\n  now:       {b_}")
            shown += 1
            if shown >= 20:
                break
    print("FAIL: the printed trees differ from tests/upgrade/predictions.txt (regenerate with --write only for an "
          "intended change, and say why in the commit)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
