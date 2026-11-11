# transform_level_order

The standard (non-Brune) parameter transformations hand the eigenvalues to the
levels of a J^pi group by energy, not by their position in the project file
(code audit item 13, fixed 2026-10-27).

`CNuc::TransformIn` turns the physical level parameters into formal ones by a
generalized eigenproblem per J^pi group. GSL returns the eigenvalues in
ascending order and AZURE2 assigned them to the levels in FILE order;
`CNuc::TransformOut` followed the same positional index. A group listed out of
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

## Models (`make_models.py`)

From `tests/capture_ay_coupling/models/twolevel.azr` (n + 15N, one 1+ group,
every width fixed):

- `asc.azr`: levels 3.491 MeV (energy free) and 3.791 MeV (energy fixed) in energy order;
- `desc.azr`: the same two levels with the 3.791 level listed first;
- `desc_fit.azr`: `desc` with a total-capture data segment on `data/cap.dat`.

## Checks (`check.sh path/to/AZURE2`, all with `--no-brune`)

- (a) `asc` and `desc` give the same total capture cross section (51 energies, 1e-10);
- (b) `desc`: `param.par`'s `energy_1` is the first-listed level's formal energy
  (3.791 exactly, since the boundary condition is set at the first-listed
  level) and `energy_2` the 3.491 level's (3.4888); `asc` unchanged;
- (c) `desc`: `parameters.out` lists the levels in file order;
- (d) a fit of `desc_fit` against data made from `desc` with the 3.491 level's
  formal energy lowered by 0.02 MeV (through a parameter file, so the truth
  lies in the fit's own one-parameter space on every build) moves the free
  level to the truth (3.4712) and keeps the fixed one at 3.791, chi2 ~ 1e-8.

Old build: (b) `energy_1` = 3.4888, (c) 3.491 listed first, (d) the fixed flag
sat on the 3.491 level, the fit moved the 3.791 level to 3.7912 and chi2
stayed at 203702.
