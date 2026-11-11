"""Add and remove levels from Python and write a new project (original untouched)."""
from pyazr import azure2, Project

model = Project.open("13N.azr")

model.remove_level(("1/2+", 20))                   # drop a background pole

model.add_level(1.5, +1, 4.1,                      # add a 3/2+ resonance (E_x in MeV)
                channels=[dict(pair=1, l=2, s=0.5, width=1000.0),
                          dict(pair=2, l=1, s=0.5, width=0.1)])

path = model.save("13N_edited.azr")
for level in model.levels():
    print(level.label, [(c.pair, c.l, c.s, c.width) for c in level.channels])

azr = azure2(path)                                 # run the edited scheme
print(f"\n{len(azr.level_scheme.levels)} levels, {len(azr.params)} free parameters")
