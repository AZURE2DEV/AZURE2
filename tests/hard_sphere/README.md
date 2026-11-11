# tests/hard_sphere

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 738de836:tests/hard_sphere/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Writes the two p+12C elastic models of this check (run once; the .azr files are committed).

with_dummies.azr : real levels plus a fixed zero-width "dummy" level for every J^pi the pair allows with
                   l <= LMAX -- the convention users followed to get every partial wave's hard-sphere
                   scattering. includeAllJPi absent (off).
no_dummies.azr   : the real levels only, includeAllJPi=1 and maxOrbitalMomentum=LMAX in <config>.
The two must give identical cross sections (check.sh).
```
