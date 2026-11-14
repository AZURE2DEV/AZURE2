#!/usr/bin/env python3
"""An evaluation's resonance parameters as a start model, and a model written as an evaluation (pyazr.endf).

1. The ENDF/B-VIII.1 3He + 4He extract of the tests (MF=2/MT=151 in LRF=7): its record printed, a new project made
   of it (Brune basis, every level fixed; separation energy from AME2020), the conversion report, and the file's own
   MF=3/MT=600 against the engine's reconstruction from that project.
2. The project T1 of the tests written as an ENDF-6 file (MF=1, MF=2 in the formal basis, MF=3 from the engine,
   LRP=2), the module's own format check of it, and its MF=3 read back.

Usage:
    python3 endf_start_model.py [WORKDIR]      (default: a new temporary directory)
"""
import os
import shutil
import sys
import tempfile
import warnings

import numpy as np

import pyazr
import pyazr.endf as endf

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
EXTRACT = os.path.join(REPO, "tests", "pyazr", "data", "endf", "h-002_He_004_extract.endf")
T1 = os.path.join(REPO, "tests", "upgrade", "T1")


def main(workdir):
    warnings.simplefilter("ignore", pyazr.EngineWarning)
    # [recipe]
    # 1. Read the resonance parameters (refusals name what, what is allowed, how to fix).
    m = endf.read(EXTRACT)
    record = m.resonances()
    print(record)
    # 2. A new project from them: pairs, levels, channels; the report says what was converted and what was not.
    for d in ("output", "checks"):
        os.makedirs(os.path.join(workdir, d), exist_ok=True)
    project, report = endf.to_project(record, os.path.join(workdir, "he3_he4.azr"), material=m)
    print(report)
    # 3. The engine's reconstruction against the file's own MF=3/MT=600 (an extrapolation at its energies).
    xs = m.cross_section(600)
    np.savetxt(os.path.join(workdir, "e_lab.dat"), np.asarray(xs.table.x) / 1e6, fmt="%.10f")
    project.add_extrapolation("hp", reaction=project.reaction(1, 2), quantity="sigma", unit="b",
                              energies="file e_lab.dat MeV lab")
    with project.session(extrapolations=True) as s:
        model = np.asarray(s.evaluate(s.parameters.x0).data[0].model, float)
    above = np.asarray(xs.table.x) > 7.2e6
    rel = np.abs(model[above] - np.asarray(xs.table.y)[above]) / np.asarray(xs.table.y)[above]
    print(f"MF=3/MT=600 above 7.2 MeV: largest relative difference {rel.max():.2e}")
    # 4. A model written as an evaluation: MF=1, MF=2 (formal basis), MF=3 from the engine; LRP=2.
    t1 = os.path.join(workdir, "T1")
    shutil.copytree(T1, t1, ignore=shutil.ignore_patterns("output*", "checks*"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(t1, d), exist_ok=True)
    with pyazr.Session(os.path.join(t1, "T1.azr")) as s:
        written = endf.write(s, os.path.join(workdir, "T1.endf"), header={"AUTH": "example", "ALAB": "EXAMPLE"})
    print(written)
    print(endf.check_format(written.path))
    print(endf.read(written.path).summary())
    # [/recipe]


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="endf_start_model_"))
