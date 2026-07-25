"""Turn an R-matrix level (a resonance / pole) on and off.

Two ways, both shown below:

  * RUNTIME  -- zero a level's reduced widths in the parameter vector and
    recalculate.  No file is written; one AZURE2 instance can toggle any level
    instantly.  Use this to *probe* what a level does (compare cross sections
    with it on, off, or alone).

  * FILE     -- AzrModel.deactivate_level() writes a new .azr in which the
    level's widths are zeroed and fixed.  Use this to *persist* a "level
    removed" model and refit it (the level stays in its J-group's shared
    channel set, so the file is still valid).

A level is switched off by decoupling it from every channel (all its reduced
widths = 0); its energy then no longer matters.  "Reactivating" is just not
zeroing it -- i.e. running the original parameters again.
"""

import numpy as np

from pyazr import azure2, AzrModel


azr = azure2("13N.azr")
best = np.asarray(azr.params_rwa, float)          # the fitted free parameters

# What are the levels?  Each is a LevelKey you can hand to the toggles below.
print("physical levels (resonances + background poles):")
for k in azr.physical_levels():
    n = len(azr.level_free_width_indices(level=k))
    print(f"  {k}   ({n} free widths)")

# Pick one to study -- by (jpi, energy), by a LevelKey, or by index.
target = azr.physical_levels()[0]                 # e.g. the first level
# target = azr._match_level(jpi="1/2-", energy=0.0)   # ...or select it explicitly


# -- RUNTIME: off / only / full, no file rewrite ----------------------------
off = azr.without_level(best, level=target)       # this level OFF
only = azr.only_level(best, level=target)         # ONLY this level (+ background)

xs_full = azr.calculate_rwa(best)                 # one cross section per segment
xs_off = azr.calculate_rwa(off)
xs_only = azr.calculate_rwa(only)

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
model = AzrModel.from_file("13N.azr")
removed = model.deactivate_level(jpi=target.jpi, energy=target.energy)
print(f"\ndeactivated in file: {[str(r) for r in removed]}")
path = model.write("13N_no_level.azr")            # original untouched
azr2 = azure2(path)                               # refit / analyze this reduced model
print(f"wrote {path}: {len(azr2.params)} free parameters "
      f"(was {len(azr.params)})")

# Reactivating = just use the original model / parameters again (`best`),
# or edit the AzrModel channels' gamma back to non-zero before writing.
