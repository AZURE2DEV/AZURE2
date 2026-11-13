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
best = np.asarray(azr.params_rwa, float)          # the fitted free parameters

# What are the levels?  Each is a LevelKey you can hand to the toggles below.
print("physical levels (resonances + background poles):")
for k in azr.physical_levels():
    n = len(azr.level_free_width_indices(level=k))
    print(f"  {k}   ({n} free widths)")

# Pick one to study -- by (jpi, energy), by a LevelKey, or by index.
target = next(k for k in azr.physical_levels() if k.energy)   # e.g. the first level above the ground state
# (the ground state is the final state of the external capture: it cannot be switched off in the file)
# target = azr._match_level(jpi="1/2-", energy=0.0)   # ...or select it explicitly


# -- RUNTIME: off / only / full, no file rewrite ----------------------------
off = azr.without_level(best, level=target)       # this level OFF
only = azr.only_level(best, level=target)         # ONLY this level (+ background)

xs_full = [d.model for d in azr.evaluate(best).data]                 # one cross section per segment
xs_off = [d.model for d in azr.evaluate(off).data]
xs_only = [d.model for d in azr.evaluate(only).data]

s = 0                                             # inspect the first segment
print(f"\nsegment {s}, first energy point:")
print(f"  full model : {xs_full[s][0]:.5g}")
print(f"  {target} OFF : {xs_off[s][0]:.5g}")
print(f"  {target} ONLY: {xs_only[s][0]:.5g}")
# interference of the target with the rest = (full - off) - (only - background)

# You can toggle several at once by intersecting index sets, e.g.:
#   keep = set(azr.level_free_width_indices(level=a)) \
#        | set(azr.level_free_width_indices(level=b))
#   v = best.copy(); ...zero every free width whose index is not in keep...


# -- FILE: persist a "level removed" model to refit -------------------------
model = Project.open("13N.azr")
model.set_level((target.jpi, target.energy), active=False)
print(f"\ndeactivated in file: {target.jpi} at {target.energy} MeV")
path = model.save("13N_no_level.azr")             # original untouched
azr2 = Session(path)                               # refit / analyze this reduced model
print(f"wrote {path}: {len(azr2.params)} free parameters "
      f"(was {len(azr.params)})")

# Reactivating = just use the original model / parameters again (`best`),
# or set_level(..., active=True) on the project.
