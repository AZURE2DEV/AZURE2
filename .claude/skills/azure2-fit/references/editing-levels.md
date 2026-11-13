# Editing the model and snapshotting a fit

The engine reads its model from the `.azr`, so an edited level scheme is applied by writing a new `.azr` and opening a
fresh session on it. `pyazr.Project` (format 3 only; an old file is refused with `AZURE2 --convert OLD NEW`)
edits the `<levels>` rows, the records, `<config>` and `<thm>` through the engine's grammars; the lines an edit does
not touch keep their bytes and nothing is written until `save`. This page covers adding, removing and deactivating
levels and channels, other edits, and `save_fit`, which writes a fit back to a `.azr` plus its `param.sav`. Building a
project from nothing (new pairs, first level scheme) is `azure2-setup`.

## Adding, removing, deactivating

```python
from pyazr import Project, Session

p = Project.open("7Be.azr")
for lv in p.levels():                             # label, J^pi, E_x, channels (pair, l, s, width), fixed, active
    print(lv.label, [(c.pair, c.l, c.s, c.width) for c in lv.channels])
lv = p.level(("5/2+", 10.253))                    # by J^pi and E_x (within 1 keV), or by its label

p.remove_level(("1/2+", 20))                      # drop a background pole entirely
p.add_level(1.5, +1, 8.6,                         # add a 3/2+ resonance, E_x in MeV
            [dict(pair=1, l=2, s=0.5, width=1000.0, fixed=False),
             dict(pair=2, l=1, s=0.5, width=0.1)],   # never width=0 on a channel meant to be fitted
            fixed=False)                          # fixed=False -> the energy is a fit parameter
p.set_level(("7/2-", 4.572), active=False)        # keep in the file, not read by the engine
path = p.save("_test.azr")

with Session(path, cwd=HERE) as m:
    ...
```

Rules:

- **All levels of one J^π share one channel set** (an R-matrix J-group requirement). If a level of that (J, parity)
  exists, `add_level` clones the group's rows and applies your `width`/`fixed` only to matching (pair, l, s);
  unmentioned channels are added inert (width 0, fixed). A channel outside the set is refused, with the set listed.
  Only for a new J^π do your channels define the group.
- `add_level` references only a pair that already exists in the file; new particle pairs are added by writing
  their level rows or by building the project (`azure2-setup`). Removing one channel from a group is not offered.
- The width column is the width input (Γ in eV for an open channel, an ANC in fm^-1/2 for a closed one; see
  `docs/source/user_guide/levels_channels.rst`), not a reduced-width amplitude. Seed a new channel with a small nonzero
  value and set `fixed=False`: a width of exactly 0 is fixed whatever its flag says (SKILL.md, "Free and fixed").
- A new level gets the next free level ID; the engine's (group, level) numbering of the others is unchanged unless the
  level falls before them in its J^π group (`references/parameters-and-sav.md`).
- Numbers are written as the shortest decimal that reads back as the same double, each right-aligned in the width of
  the token it replaces.
- A capture data set summed over the γ cascade (`reaction 12C(p,g)13N` with every capture pair) needs every
  significant transition set up beforehand as its own (particle, γ) pair. Fetching data sets from EXFOR:
  `azure2-setup`.
- **Adding or removing data sets invalidates the external-capture cache.** After `add_dataset`, `remove`, or an
  extrapolation change, delete `output/intEC.dat` and `output/intEC.extrap` before the next run, or give the edited
  model its own output directory (`p.set_config(output_dir="out2/")`); otherwise the engine may reuse integrals
  computed for the old grids. In a live session `s.recalculate_external_capture()` forces a recompute
  (`pyazr/examples/edit_model.py`).
- Data-set parameters are named by record (`meyer_84.norm`), so adding or removing a data set renames nothing else;
  deactivating one (`p.set(name, active="no")`) keeps its rows of `<parameterSettings>` for when it comes back
  (`references/drop-dataset.md`).

## `save_fit`: snapshot a fit

```python
azr, sav = s.save_fit("7Be_fit.azr")          # current parameters
azr, sav = s.save_fit("7Be_fit.azr", x_best)  # or an explicit free vector
```

`s.save_fit(path, x=None)` writes the `.azr` and the companion `param.sav`, then verifies: it reopens what it wrote and
compares the physical values transformed back, and the reduced-width amplitudes the file reads for the particle and
β channels of each level, against the fit's (up to a common sign of a level, which the standard parametrization may
flip). If either disagrees it removes both files and raises. `path` is explicit; nothing is written in place. It
writes every level parameter at its current observed value (`transform_all_rwa(full, free_only=False)`,
`place_level_values`, fixed parameters included), so the re-read is exact; `load_sav` of the companion `.sav` gives the
same.

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
  where the width is given as exactly zero in the file, it does not switch external capture off.
- **Target-integration sub-point grids** are re-placed on every evaluation from the current level energies and widths
  (`--no-grid-replacement`, pyazr `replace_grids=False`, restores a grid fixed at start-up; segments that need
  external-capture amplitudes keep their grids; `docs/source/user_guide/experimental_effects.rst`). With replacement a
  fresh load and a running session agree (one model: 25,776.076357 from either file).
  builds before 2026-10-19 (expires: 2027-02-01) placed the grid once from the file's levels, so a fit with narrow
  levels under a target integration minimized against the start model's grid; re-evaluate such a fit's end point in a
  current session before comparing.
  builds before 2026-10-22 (expires: 2027-02-01) also integrated the tail of a level much narrower than the coarse step
  badly (χ² oscillating with a 10 eV period in E_level, analytic gradient of the wrong sign); there, a wide fine region
  (Width Mult. 200, Pts/Width 5) is the work-around. If χ² oscillates or gradient checks fail in one narrow level's
  energy under a target integration, suspect the grid before the physics.
- **The verify tolerance (`rtol=1e-4`) can hide fragility.** The rwa ↔ physical conversion for an open channel is
  reproducible between sessions to 10⁻⁷–10⁻⁶ relative. Near a destructive-interference minimum (a broad, high-lying
  background pole feeding a low-energy capture channel) that noise can move a data set's χ² by orders of magnitude: a
  5 × 10⁻⁷ relative change to one 887 keV-wide background level's widths turned one segment's χ² of 883 into 358,000.
  Check the objective (`s.objective()(x)`, or `s.evaluate(x)` per data set) of a reopened snapshot before using it as a
  refit's start.

## What a snapshot still cannot carry

- **Fitted normalizations and shifts** go into their records' nuisances with `save_fit(norms="fitted")` (the
  default; each prior keeps its own centre); with `norms="nominal"` they are only in the `.sav`, and a calculation
  with a blank parameter file is then not the fit (`references/normalizations-and-shifts.md`).
- **Verify a bake of `<levels>` with a blank parameter file**, not with `param.sav`: mode 1 with an external file
  loads that file's values by name and overrides `<levels>`, so it can only confirm the `.sav`. A blank-file run must
  equal a pyazr evaluation at the fitted R-matrix block with every norm and shift at its nominal value.
- **The check dumps are keywords, not filenames.** `<config>` accepts only `none`, `screen` or `file`
  (`Config::ReadConfigFile`); anything else leaves the check off and `checks/` stays empty. Write `file`.

## The placement and the example script

`pyazr.snapshot.place_level_values(project, parameters, pairs, want)` is the lower-level half: it maps the engine's
pair numbers to the file's keys through `pairs=` (`s.pairs`), places each value by the engine's own (group, level) and
the channel's (pair key, l, 2s), and refuses any free parameter it cannot place. It does not verify; `save_fit` does.
`pyazr/examples/save_fit_to_azr.py` does the whole job: it loads a fit from `param.sav` or an `.npz`, optionally changes
a channel radius (`--radius 1=4.70`, dropping the stale `intEC` caches), writes the `.azr` plus a companion `param.sav`
with the norms, and fails if the result does not round-trip.
