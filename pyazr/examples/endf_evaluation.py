#!/usr/bin/env python3
"""An ENDF-6 evaluation next to a project: read it, see what it holds, draw it over the data, check the file.

The four things an evaluator does with ``pyazr.endf`` (user guide: ENDF-6 evaluations), on any ENDF-6 file and
project; with no arguments, on the hand-written p + 12C file of the tests (``tests/pyazr/data/endf/p-C012_hand.endf``,
made-up coefficients, not physics) and T1 (p + 12C):

1. ``pyazr.endf.read(path)`` and ``print(material.summary())``;
2. ``material.curves(data, session=s)`` for every data set: the file's values at its points in its quantity;
3. ``pyazr.endf.check(path)``: the module's own format check, then CHECKR, FIZCON and PSYCHE when they are found
   (``$AZURE2_ENDF_UTILITIES``);
4. ``pyazr.endf.standardize(path)`` when STANEF is found: what it would change, written nowhere.

Usage:
    python3 endf_evaluation.py [FILE.endf PROJECT.azr]
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

import pyazr
import pyazr.endf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
path, project = (sys.argv[1:3] if len(sys.argv) >= 3 else
                 (os.path.join(ROOT, "tests", "pyazr", "data", "endf", "p-C012_hand.endf"),
                  os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")))

# 1. what the file holds
material = pyazr.endf.read(path)
print(material.summary())

# 2. the evaluation over each data set of the project
with pyazr.Session(project) as s:
    for data in s.evaluate(s.parameters.x0).data:
        try:
            curve = material.curves(data, session=s)
        except pyazr.FormatRefused as e:
            print(f"{data.name}: not in the file ({e.what})")
            continue
        covered = curve.coverage
        ratio = np.median(curve.values[covered] / data.model[covered]) if covered.any() else float("nan")
        print(f"{data.name}: {curve.quantity} from {curve.source}, {covered.sum()} of {data.n} points covered, "
              f"median file/model {ratio:.4g}")

# 3. the file's format: the module's own check, then the checking codes found
print(pyazr.endf.check(path))

# 4. what STANEF would change (a candidate in a scratch directory; nothing replaced)
if any(c.name == "stanef" for c in pyazr.endf.checkers()):
    print(pyazr.endf.standardize(path).summary())
else:
    print("STANEF not found: pyazr.endf.checkers() lists the codes found")
