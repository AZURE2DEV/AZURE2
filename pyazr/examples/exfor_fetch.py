"""Fetch nuclear data from EXFOR/NDS and drop it into an AZURE2 model.

Demonstrates the ``pyazr.nds`` module end to end:

1. search EXFOR for a reaction (``search_exfor``)
2. fetch a dataset and convert it to an AZURE2 data file (``ExforData.to_azr``)
3. attach it to an ``AzrModel`` as a new data segment (``add_data_segment``)
4. pull the level scheme for the compound nucleus from ENSDF (``fetch_levels``)
5. resolve the paper behind the dataset (``reference`` + ``resolve_doi``)

Run from the repository root::

    python pyazr/examples/exfor_fetch.py

The EXFOR server is on the public internet; a network connection is required.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from pyazr import AzrModel
from pyazr import nds as N

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "_exfor_demo")
DATA = os.path.join(WORK, "data")

# 1. Search EXFOR for 13C(p,g)14N total cross sections (S-factor datasets).
print("== search_exfor(C-13, p,g, SIG) ==")
hits = N.search_exfor(target="C-13", reaction="p,g", quantity="SIG", limit=5)
for h in hits:
    print(" ", h)
if not hits:
    raise SystemExit("no EXFOR hits; check connectivity / parameters.")

# 2. Fetch one dataset (O2599004: Skowronski et al. 2023, LUNA) and convert
#    it into an AZURE2 data file in the AZURE2 layout:
#      energy(MeV, lab)  angle(deg)  cross-section(b)  error
os.makedirs(DATA, exist_ok=True)
dataset = N.fetch_exfor("O2599004")
print("\n== fetch_exfor('O2599004') ==")
print("  reaction:", dataset.reaction)
print("  columns:", [(c.name, c.unit) for c in dataset.columns])

kw = dataset.to_azr(DATA, entrance=1, exit=2, observable="total-capture")
print("\n== to_azr(...) segment kwargs ==")
for k, v in kw.items():
    print(f"  {k} = {v}")
print("  data file head (lab E_MeV, angle, sigma_b, err_b):")
with open(kw["data_file"]) as f:
    for line in f.readlines()[:3]:
        print("   ", line.rstrip())

# 3. Attach the new segment to a model.  The 13N.azr compound nucleus is
#    13C+p (pair 1) -> capture (pair 2), which matches the EXFOR entrance.
mdl = AzrModel.from_file(os.path.join(HERE, "..", "..", "tests", "13N", "13N.azr"))
mdl.add_data_segment(**kw)
out_azr = os.path.join(WORK, "13N_nds.azr")
mdl.write(out_azr)
print("\n== wrote model with the new segment: %s ==" % out_azr)

# 4. Level scheme of the compound nucleus (14N) from ENSDF via LiveChart.
print("\n== fetch_levels('14n') ==")
for lv in N.fetch_levels("14n")[:6]:
    print("  %8.3f keV  %-5s  %s" % (lv.energy, lv.jp, lv.half_life))

# 5. The paper behind the dataset.
print("\n== reference('O2599004') ==")
ref = N.reference("O2599004")
print("  entry  :", ref.entry)
print("  title  :", ref.title)
print("  cite   :", ref.reference)
print("  authors:", ref.authors[:80], "...")
doi = N.resolve_doi(ref)
print("  DOI    :", doi, " -> https://doi.org/%s" % doi)

print("\nNOTE: data segments changed the grids; delete output/intEC.dat and")
print("output/intEC.extrap before running this model (or set a fresh output dir).")
