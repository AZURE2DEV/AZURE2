"""Add and remove levels from Python and write a new .azr (original untouched)."""
from pyazr import azure2, AzrModel

model = AzrModel.from_file("13N.azr")

model.remove_level(jpi="1/2+", energy=20)          # drop a background pole

model.add_level(J=1.5, parity=1, energy=4.1,       # add a 3/2+ resonance
                channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0),
                          dict(pair=2, L=1, S=0.5, gamma=0.1)])

path = model.write("13N_edited.azr")
print(model)

azr = azure2(path)                                 # run the edited scheme
print(f"\n{len(azr.level_scheme.levels)} levels, {len(azr.params)} free parameters")
