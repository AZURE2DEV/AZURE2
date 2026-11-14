"""Turn an R-matrix level (a resonance / pole) on and off.

Two ways, both shown below:

  * RUNTIME  -- zero a level's reduced widths in the parameter vector and
    recalculate.  No file is written; one AZURE2 instance can toggle any level
    instantly.  Use this to *probe* what a level does (compare cross sections
    with it on, off, or alone).

  * FILE     -- Project.set_level(level, active=False) writes a new project
    in which the level's rows are inactive: the engine does not read them.  Use
    this to *persist* a "level removed" model and refit it (the rows stay in
    the file; set_level(level, active=True) brings the level back).

A level is switched off by decoupling it from every channel (all its reduced
widths = 0); its energy then no longer matters.  "Reactivating" is just not
zeroing it -- i.e. running the original parameters again.
"""

import numpy as np

from pyazr import Session, Project


azr = Session("13N.azr")
P = azr.parameters                                # the table of the fit parameters
best = P.x0                                       # the fitted free parameters (the fit vector)

# What are the levels?  Each energy row names its level (level_label, as Project.levels() gives it).
print("physical levels (resonances + background poles):")
for e in P.of_kind("energy"):
    n = len(P.of_level(e.level_label).of_kind("width").free)
    print(f"  {e.level_label}   ({n} free widths)")

# Pick one to study, by its label.
target = next(e for e in P.of_kind("energy") if e.level_energy)   # e.g. the first level above the ground state
# (the ground state is the final state of the external capture: it cannot be switched off in the file)
label = target.level_label


# -- RUNTIME: off / only / full, no file rewrite ----------------------------
widths = P.of_kind("width").free
off = P.put(best, P.of_level(label).of_kind("width").free, 0.0)                     # this level OFF
only = P.put(best, [w for w in widths if w.level_label != label], 0.0)               # ONLY this level (+ background)

xs_full = [d.model for d in azr.evaluate(best).data]                 # one cross section per segment
xs_off = [d.model for d in azr.evaluate(off).data]
xs_only = [d.model for d in azr.evaluate(only).data]

s = 0                                             # inspect the first segment
print(f"\nsegment {s}, first energy point:")
print(f"  full model : {xs_full[s][0]:.5g}")
print(f"  {label} OFF : {xs_off[s][0]:.5g}")
print(f"  {label} ONLY: {xs_only[s][0]:.5g}")
# interference of the target with the rest = (full - off) - (only - background)

# Several levels at once: P.put(best, [w for w in widths if w.level_label not in keep], 0.0).


# -- FILE: persist a "level removed" model to refit -------------------------
model = Project.open("13N.azr")
model.set_level(label, active=False)
print(f"\ndeactivated in file: {label}")
path = model.save("13N_no_level.azr")             # original untouched
azr2 = Session(path)                               # refit / analyze this reduced model
print(f"wrote {path}: {len(azr2.parameters.x0)} free parameters "
      f"(was {len(best)})")

# Reactivating = just use the original model / parameters again (`best`),
# or set_level(..., active=True) on the project.
