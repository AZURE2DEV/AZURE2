#!/usr/bin/env python3
"""Write a format-3 project from Python: schema -> records -> diagnostics -> describe -> calculate.

The recipe of docs/source/user_guide/project_from_python.rst (and of the azure2-setup skill): a thick-target yield
curve of 15N(p,a)12C in counts per uC (charge state 1) on the levels of T2, the worked example of
docs/upgrade/project_format3.md section 9.1 as tests/units_scales/models/yield_story.azr holds it.  Every record is
checked by the engine as it is added; a refusal is a pyazr.FormatError with the engine's codes.

Usage:
    python3 project_from_python.py [WORKDIR]      (default: a new temporary directory; the project is WORKDIR/yield.azr)
"""
import os
import shutil
import sys
import tempfile

import pyazr
from pyazr import Project, Nuisance, FormatError

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
STORY = os.path.join(REPO, "tests", "units_scales", "models")


def main(workdir):
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    if not os.path.exists(os.path.join(workdir, "data")):
        shutil.copytree(os.path.join(STORY, "data"), os.path.join(workdir, "data"))

    # [recipe]
    # 1. The schema: what a data set of quantity "yield" takes (the engine's table, never spelled in pyazr).
    yield_row = pyazr.schema()["quantities"]["yield"]
    print("yield units:", ", ".join(yield_row["units"]))
    # help(Project.add_dataset) lists every field with its grammar, condition and units.

    # 2. A new project on the nuclear part (<config>, <levels>, ...) of a format-3 file; data files are read
    #    relative to the project's directory.
    p = Project.new(os.path.join(workdir, "yield.azr"), nuclear_part=os.path.join(STORY, "yield_story.azr"))

    # 3. Records, the target before the data set that names it.  A value is the file's text, a number, a
    #    (value, unit) tuple or a pyazr.Nuisance.
    p.add_target("N15T", kind="infinite", active="15N", composition="N 1",
                 stopping='formula "a0*x^a1" a0=4.5e-20 a1=-0.6 MeV*cm2/atom', straggling="none")
    p.add_beam("pbeam", spread=(0.3, "keV", "lab"))
    curve = dict(reaction="15N(p,a)12C", quantity="yield", unit="counts/uC", efficiency=Nuisance(0.0182, "fixed"),
                 file="data/yield_curve.dat", columns="E_lab[keV] _ value error",
                 norm=Nuisance(1, "fit", prior=(1, "3%")), target="N15T", beam="pbeam", quadrature="points 60")

    # 4. The diagnostics: a record the engine refuses raises FormatError, the project unchanged.
    try:
        p.add_dataset("Y_curve", **curve)
    except FormatError as e:
        print("refused:", e.codes, [d["field"] for d in e.diagnostics])   # ['E_FIELD_REQUIRED'] ['charge_state']
    p.add_dataset("Y_curve", charge_state=1, **curve)
    assert p.validate() == []

    # 5. What the engine makes of it, before computing.
    print(p.describe("Y_curve"))

    # 6. Calculate: the project is saved and opened as a pyazr.azure2 session.
    with p.session() as m:
        chi2 = float(sum(m.segment_chi2()))
    print(f"chi2 = {chi2:.6g} on {len(m.energies[0])} points; the project is {p.path}")
    # [/recipe]
    return p, chi2


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    main(sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="project_from_python."))
