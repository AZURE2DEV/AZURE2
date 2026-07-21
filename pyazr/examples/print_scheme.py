"""Read a model's level scheme and dataset provenance."""
from pyazr import azure2

azr = azure2("13N.azr")

print(azr.level_scheme)          # pairs, J-groups, levels, channels, widths, Wigner
print()
print(azr.datasets.table())      # data file, reaction, observable, systematics

# everything is queryable, e.g.
scheme = azr.level_scheme
print("\nresonances below 5 MeV:")
for lv in scheme.resonances(energy_max=5):
    print(f"  {lv.jpi:>4}  E = {lv.energy:.3f} MeV  ({len(lv.channels)} channels)")
