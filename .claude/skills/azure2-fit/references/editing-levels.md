# Editing the model and snapshotting a fit

The engine reads its model from the `.azr`, so an edited level scheme is applied by writing a new `.azr` and opening a
fresh session on it. `pyazr.AzrModel` parses the `<levels>` block and re-emits everything else verbatim; the original
file is never modified. This page covers adding, removing and deactivating levels and channels, other `.azr` edits,
and `save_fit`, which writes a fit back to a `.azr` plus its `param.sav`. Building a project from nothing (new pairs,
first level scheme) is `azure2-setup`.

## Adding, removing, deactivating

```python
from pyazr import AzrModel, azure2

mdl = AzrModel.from_file("7Be.azr")
print(mdl)                                        # J^pi / energies / channels / fixed flags
mdl.find(jpi="5/2+", energy=10.253, tol=2e-2)     # -> [AzrLevel]

mdl.remove_level(jpi="1/2+", energy=20)           # drop a background pole entirely
mdl.add_level(J=1.5, parity=+1, energy=8.6,       # add a 3/2+ resonance, E_x in MeV
              channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0, fixed=False),
                        dict(pair=2, L=1, S=0.5, gamma=0.1)],   # never gamma=0 on a channel meant to be fitted
              level_fixed=False)                  # level_fixed=False -> the energy is a fit parameter
mdl.deactivate_level(jpi="7/2-", energy=4.572)    # keep in file, zero and fix every width
path = mdl.write("_test.azr")

with azure2(path, cwd=HERE) as m:
    ...
```

Rules:

- **All levels of one J^π share one channel set** (an R-matrix J-group requirement). If a level of that (J, parity)
  exists, `add_level` clones the group's channel structure and applies your `gamma`/`fixed` only to matching
  (pair, L, S); unmentioned channels are added inert (`gamma=0`, fixed). A spec matching none raises and lists the
  allowed channels. Only for a new J^π do your channels define the group.
- `remove_channel(jpi, pair, L, S)` removes that channel from every level of the group; removing it from one level
  makes the engine reject the file.
- `add_level`/`add_channel` reference only a pair that already exists in the file; new particle pairs are added in the
  GUI or by building the project (`azure2-setup`).
- `gamma` in the file is the width input (Γ in eV for an open channel, an ANC in fm^-1/2 for a closed one; see
  `docs/source/user_guide/levels_channels.rst`), not a reduced-width amplitude. Seed a new channel with a small nonzero
  value and set `fixed=False`: a width of exactly 0 is fixed whatever its flag says (SKILL.md, "Free and fixed").
- Levels are renumbered (`levelID`) on every edit; every later `energy_<N>`/`width_<N>_<c>` name moves
  (`references/parameters-and-sav.md`).
- Round level energies to at most 8 decimals when writing a `.azr` programmatically. A level line written by
  `AzrModel.add_level` with a 16-significant-digit energy (e.g. `3.436212548111519`) made the engine misread the
  model without a message (χ² × 4000 on ¹³N); the same energy with 8 decimals, and the 18-character energies that
  `apply_fit` writes, read correctly, so the trigger is the token layout of `add_level`'s lines. rmfit rounds to 1 eV.
- `observable="total-capture"` needs every significant γ-cascade transition set up beforehand as its own
  (particle, γ) pair; the engine sums over them for a segment with this observable. Other observables include
  `analyzing-power` (code 7). Fetching data sets from EXFOR to feed `add_data_segment`: `azure2-setup`.
- **Adding or removing data segments invalidates the external-capture cache.** After `add_data_segment`,
  `remove_data_segments`, `clear_data_segments` or `set_extrapolations`, delete `output/intEC.dat` and
  `output/intEC.extrap` before the next run, or give the edited model its own output directory
  (`mdl.set_output_dir(...)`); otherwise the engine reuses integrals computed for the old grids. In a live session
  `azr.recalculate_external_capture()` forces a recompute (`pyazr/examples/edit_model.py`).
- **Removing or inserting a `<segmentsData>` line renames every later segment parameter**: remap every per-segment
  family or, better, deactivate (`references/parameters-and-sav.md`, `references/drop-dataset.md`).

## `save_fit`: snapshot a fit

```python
azr, sav = m.save_fit("7Be_fit.azr")          # current parameters
azr, sav = m.save_fit("7Be_fit.azr", x_best)  # or an explicit free vector
```

`m.save_fit(path, x=None)` writes the `.azr` and the companion `param.sav`, then verifies: it reopens what it wrote and
compares the physical values transformed back, and the reduced-width amplitudes the file reads for the particle and
β channels of each level, against the fit's (up to a common sign of a level, which the standard parametrization may
flip). If either disagrees it removes both files and raises. `path` is explicit; nothing is written in place. It
writes every level parameter at its current observed value (`transform_all_rwa(full, free_only=False)`,
`apply_fit(..., include_fixed=True)`), so the re-read is exact; `load_sav` of the companion `.sav` gives the same.

- **Fixed parameters are fixed in reduced-width space, not in what the file shows.** A fixed channel's observed width
  or ANC moves whenever the free widths of the same level move (Brune's factor 1 + Σγ² dS/dE, the external width), and
  a fixed width's observed value moves with the level energy through the penetrability P(E). A `<levels>` block with
  the file's old observed values for fixed channels and fitted values for free ones is therefore not the fit (one
  model re-read 21 units off, 98,642 → 98,663, with 19 fixed amplitudes drifting by up to 14 %).
- builds before 2026-10-12 (expires: 2027-02-01): `save_fit` compared only physical values and kept the file's old
  observed values for fixed channels; the engine returned |B| for a β-decay feeding amplitude, so a negative one was
  written back positive (one β-delayed α χ² went 750 → 29,724) and the verify passed. `tests/pyazr/beta_sign_test.py`
  tells in a minute whether a pyazr module has the defect.
- The external-capture integrals depend only on the EC final-level energies (fixed in practice), pair data and the data
  grid, so `recalculate_external_capture()` changes a fit's χ² by about 0.1 at most; it uses the level energies from
  the file, not the session's. `ignore_externals=True` is `--ignore-externals`: it drops the external amplitude only
  where the internal width is zero, it does not switch external capture off.
- **Target-integration sub-point grids** are re-placed on every evaluation from the current level energies and widths
  (`--no-grid-replacement`, pyazr `replace_grids=False`, restores a grid fixed at start-up; segments that need
  external-capture amplitudes keep their grids; `docs/source/user_guide/experimental_effects.rst`). With replacement a
  fresh load and a running session agree (one model: 25,776.076357 from either file).
  builds before 2026-10-19 (expires: 2027-02-01) placed the grid once from the file's levels, so a fit with narrow
  levels under `<targetInt>` minimized against the start model's grid; re-evaluate such a fit's end point in a current
  session before comparing.
  builds before 2026-10-22 (expires: 2027-02-01) also integrated the tail of a level much narrower than the coarse step
  badly (χ² oscillating with a 10 eV period in E_level, analytic gradient of the wrong sign); there, a wide fine region
  (Width Mult. 200, Pts/Width 5) is the work-around. If χ² oscillates or gradient checks fail in one narrow level's
  energy under `<targetInt>`, suspect the grid before the physics.
- **The verify tolerance (`rtol=1e-4`) can hide fragility.** The rwa ↔ physical conversion for an open channel is
  reproducible between sessions to 10⁻⁷–10⁻⁶ relative. Near a destructive-interference minimum (a broad, high-lying
  background pole feeding a low-energy capture channel) that noise can move a data set's χ² by orders of magnitude: a
  5 × 10⁻⁷ relative change to one 887 keV-wide background level's widths turned one segment's χ² of 883 into 358,000.
  Check the total (or `m.segment_chi2`) objective of a reopened snapshot before using it as a refit's start.

## What a snapshot still cannot carry

- **Fitted normalizations and shifts are only in the `.sav`.** `save_fit` leaves `<segmentsData>` alone, so a
  calculation with a blank parameter file is not the fit, and baking norms in is not the fix
  (`references/normalizations-and-shifts.md`).
- **Verify a bake of `<levels>` with a blank parameter file**, not with `param.sav`: mode 1 with an external file
  loads that file's values by name and overrides `<levels>`, so it can only confirm the `.sav`. A blank-file run must
  equal a pyazr evaluation at the fitted R-matrix block with every norm and shift at its nominal value.
- **The check dumps are keywords, not filenames.** `<config>` accepts only `none`, `screen` or `file`
  (`Config::ReadConfigFile`); anything else leaves the check off and `checks/` stays empty. Write `file`.

## `apply_fit` and the example script

`AzrModel.apply_fit` is the lower-level half: it takes `pairs=` (the engine's pair numbering, `m.pairs`), matches
levels on the engine's own (jgroup, level) via `engine_level_keys()`, and refuses any parameter it cannot place. It does
not verify; `save_fit` does. `pyazr/examples/save_fit_to_azr.py` does the whole job: it loads a fit from `param.sav` or
an `.npz`, optionally changes a channel radius (`--radius 1=4.70`, dropping the stale `intEC` caches), writes the `.azr`
plus a companion `param.sav` with the norms, and fails if the result does not round-trip.
