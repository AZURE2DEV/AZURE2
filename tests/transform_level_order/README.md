# transform_level_order

The standard (non-Brune) parameter transformations hand the eigenvalues to the
levels of a J^pi group by energy, not by their position in the project file
(code audit item 13, fixed 2026-10-27).

`CompoundSystem::TransformIn` turns the physical level parameters into formal ones by a
generalized eigenproblem per J^pi group. GSL returns the eigenvalues in
ascending order and AZURE2 assigned them to the levels in FILE order;
`CompoundSystem::TransformOut` followed the same positional index. A group listed out of
energy order (a 40 MeV background pole first, say) therefore had its energies,
widths and fix flags permuted: the first-listed level carried the lowest
formal energy, `parameters.out` listed the group in ascending order rather than
the file's, and a fit that fixed the pole's energy varied the low level
instead. The cross sections were never affected -- the set of levels is the
same whatever the labels -- which is why nothing warned. Only runs through
the transformation (`--no-brune`, or the Brune box unchecked in the GUI;
pyazr `use_brune=False`) are concerned; the Brune parametrization (the
default everywhere) does not solve the eigenproblem.

Now each level receives the eigenvalue of the same rank as its physical
energy, and `TransformOut` follows the rank of the level's formal energy. For a
file in energy order both mappings are the identity, so every existing
result is unchanged.

## Models (`make_models.py` (deleted; see "How the models were made"))

From `tests/capture_ay_coupling/models/twolevel.azr` (n + 15N, one 1+ group,
every width fixed):

- `asc.azr`: levels 3.491 MeV (energy free) and 3.791 MeV (energy fixed) in energy order;
- `desc.azr`: the same two levels with the 3.791 level listed first;
- `desc_fit.azr`: `desc` with a total-capture data segment on `data/cap.dat`.

## Checks (`check.sh path/to/AZURE2`, all with `--no-brune`)

- (a) `asc` and `desc` give the same total capture cross section (51 energies, 1e-10);
- (b) `desc`: `param.par`'s `energy_1` is the first-listed level's formal energy
  and `energy_2` the 3.491 level's, and they equal `asc`'s two formal energies
  exchanged (3.7921077 and 3.491, to the 10⁻⁷ MeV printed): the boundary
  condition is set at the lowest level whatever the listing order (HANDOFF
  decision 10), so the lowest level's formal energy is 3.491 exactly; `asc`
  unchanged;
- (c) `desc`: `parameters.out` lists the levels in file order;
- (d) a fit of `desc_fit` against data made from `desc` with the 3.491 level's
  formal energy lowered by 0.02 MeV (through a parameter file, so the truth
  lies in the fit's own one-parameter space on every build) moves the free
  level to the truth (3.4712) and keeps the fixed one at 3.791, chi2 ~ 1e-8.

Builds from 2026-10-27 to 2026-11-07 set the boundary condition at the
first-listed level: (b) `energy_1` = 3.791, `energy_2` = 3.4888, not `asc`'s
levels; (d) took its truth at 3.4688.

Old build: (b) `energy_1` = 3.4888, (c) 3.491 listed first, (d) the fixed flag
sat on the 3.491 level, the fit moved the 3.791 level to 3.7912 and chi2
stayed at 203702.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 8c238222:tests/transform_level_order/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for tests/transform_level_order (audit item 13, 2026-10-27).

Derived from tests/capture_ay_coupling/models/twolevel.azr (n + 15N, one 1+ group with two levels,
every width fixed):

  asc.azr       levels in energy order: 3.491 (energy free), 3.791 (energy fixed)
  desc.azr      the same two levels with the 3.791 level listed FIRST
  desc_fit.azr  desc with a total-capture data segment on data/cap.dat (written by check.sh from desc
                with the 3.491 level's formal energy lowered by 0.02 MeV)

Run from anywhere:  python3 tests/transform_level_order/make_models.py
```
