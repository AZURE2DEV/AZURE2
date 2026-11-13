# Parameter names, parameter files and remapping

A fit is the `.azr` plus its parameter file (`param.sav`, or `param.par`/`param.fit`). The engine and pyazr match a
parameter file to the model by parameter name only, and they say nothing about a name that matches the wrong
parameter or matches nothing. This page states how the names are built, how to read and load a parameter file, and how
to carry values across a structural edit. Behaviour reference: `docs/source/reference/output_files.rst` (parameter
files) and `docs/source/user_guide/levels_channels.rst` (how level lines are grouped).

## How the names are assigned

- **Level parameters** (`CNuc::FillMnParams`): `energy_<N>` counts levels in J-group order (J-groups in order of first
  appearance in `<levels>`, then the group's levels in file order), not in raw file order and not by `levelID`.
  `width_<N>_<c>` is channel `c` of the J-group's channel set, in the order of that level's channel lines. Only active
  levels are counted.
- **Segment parameters**: `segment_<key>_norm`, `segment_<key>_energy_shift`, `segment_<key>_energy_shift_sqrt`
  (and the shift-function terms), where `<key>` is the 1-based position of the `<segmentsData>` line counting every
  line, active or not (`Segment.key`, matching `Parameter.segment_key`; `WorkingCopy::FillMnParams`).
- **Consequence**: deactivating, removing or adding a level renames every later level parameter, and removing or
  inserting a `<segmentsData>` line renames every later segment parameter. A `param.sav` written before the edit then
  puts every later value on the wrong parameter, without an error.
- Old builds used `j=<J>_la=<L>_energy` and `j=<J>_la=<L>_ch=<C>_rwa`; convert such files first
  (`references/legacy-params.md`).

## Reading a parameter file

- **Read by name, never by row position.** `param.par`/`param.sav` carry a parametrization tag (formal 0, Brune 1,
  Park 2). Three forms exist in archives side by side: untagged; a bare first row `parametrization  1.0  0.0` (written
  by builds of one day; a reader indexing rows is off by one, and in one rate study this took the ground state for the
  Hoyle state and gave a rate exactly 2× high); and a comment line `#parametrization 1 0` (skipped by numpy, read by
  the engine). The engine's own reader and rmfit read by name and were never affected. A name-based reader that
  accepts all three:
  ```python
  vals = {t[0]: float(t[1]) for t in (l.split() for l in open(sav))
          if len(t) >= 2 and not t[0].startswith("#") and t[0] != "parametrization"}
  ```
- `pyazr.bands.best_fit_params` raises "has N parameters but the model has M" on an old-format `.sav` in data mode.

## Loading a parameter file in pyazr

- **Use `x = m.load_sav(path)`**, not a slice of `full_rwa_from_sav`. The CLI's parameter-file reader also sets fixed
  level energies and widths. pyazr's free-vector calls (`evaluate(x)`, `objective()(x)`,
  `transform_rwa`, `save_fit`, ...) take every fixed value from the session, that is from the `.azr`. So
  `x = m.full_rwa_from_sav(sav)[~fixed]` evaluates a different model whenever the `.sav`'s fixed values differ from the
  `.azr`'s (a fit restarted from another fit's `.sav`, or an `.azr` whose fixed widths were never updated). In one case
  10 fixed widths off by up to 1.6 % gave χ² = 9925.03 against the CLI's 9868.04; with `load_sav` every `AZUREOut_*`
  line and `chiSquared.out` matched the CLI exactly.
- `load_sav` pushes the whole file into the engine (`sess.set_params_all_rwa`), sets `m.params_rwa` and returns the
  free vector; `update_rwa_params_from_sav` and `bands.best_fit_params` use it. It does not change the free/fixed
  split (the CLI frees a parameter fixed at 0 that the file sets nonzero) and refuses a file tagged with another basis.
- `extrap_mode()` and `data_mode()` re-initialise the engine from the `.azr`: call `load_sav` again after switching.
- builds before 2026-10-12 (expires: 2027-02-01) have no `load_sav`; there, `m.sess.update_segments_all_rwa(full)`
  (one forward pass) has the same effect.

## Which physical level is `energy_<N>`?

- pyazr's `m.parameters` names need not equal the names in a `param.sav` written from another file: in one model
  pyazr's `energy_39` was the 14.72 MeV level while the `param.sav`'s `energy_39` was the 16.77 MeV one (pyazr showed
  the `.azr` seed, `param.sav` the fit). Identify a `param.sav` name by (a) `parameters.out`, which lists
  `J = ... E_level = ...` blocks in the same order, or (b) a one-line CLI probe: change that entry in a scratch copy,
  run mode 1, and see which `E_level` or width line changed.
- `m.parameters` is the engine's records, one row per parameter (`print(m.parameters.table())`: name, kind, unit,
  start, free, bounds, what it belongs to). A `**WARNING: Denominator less than zero while transforming` printed
  while it is built concerns the transformation of the start values, not the table: cross-check with a plain
  `evaluate(x).data[i].energies` or the CLI before concluding the parameters are corrupted.
- A `LevelKey` prints as `5/2-#2@6.588MeV`; `(jgroup, level)` is its identity (level numbering restarts inside every
  J-group).
- **A `Parameter`'s `pair` is the engine's pair number**, which counts particle pairs in the order `<levels>` first
  mentions them, not the file's pair key. In one model engine pair 1 is file key 2 and file key 1 is engine pair 6.
  When calling `pyazr.snapshot.place_level_values` directly, pass `pairs=m.pairs`, or it cannot translate.

## Remapping after a level edit

- A plain copy of the old `param.sav` is unsafe even as a starting point. After a 3⁺ level was added to one model, the
  new level's energy was set to 40.0 MeV by a name that now belonged to a background pole (the start χ² was
  2.7 × 10⁸ instead of about 7.7 × 10⁵). After one level was switched off in another, 279 names became 271 and
  `energy_31` named the next level. Inserting a background pole into a 367-parameter model matched 327 names by name
  and started the fit at χ² = 3 × 10⁹.
- Safe options: (a) bake the fit into `<levels>` (`m.save_fit`) before the structural edit and seed only the
  `segment_*` lines afterwards; or (b) remap by physical identity (kind, level energy, J^π, pair key, l, radiation
  type; segment key for normalizations), let genuinely new parameters keep the value the new `.azr` carries, and check
  that every mapped pair has identical values before trusting the map. To list names on both models:
  `[p.name for p in m.parameters if p.kind in ("energy", "width")]`.

## Segment keys: removing or inserting a `<segmentsData>` line

`AZUREParams::ReadUserParameters` matches an external parameter file by exact name across all parameters, fixed or
free, and calls `SetValue()` on any match, with no warning on a stale one. Deleting an earlier line shifts every later
key down by one. Two failures seen in one archive fit, both loaded and fitted without a message:

1. a remap script renamed `segment_N_norm` past the cut but not `segment_N_energy_shift`, so one segment's fixed
   energy calibration was applied to its neighbour;
2. a remap regex anchored at the line start missed every `segment_N_norm`, because those lines are indented in
   `param.sav` while `segment_N_energy_shift` lines are not.

Remap every per-segment family by the same shifted key and verify it directly (diff the rows on both sides of the
edit). Better: deactivate instead of deleting (`references/drop-dataset.md`).

## Two levels of one J^π closer than 1 keV

Lines of `<levels>` belong to one level when they carry the same level ID (ninth column) within 1 keV, or the
identical energy; two levels of one J^π closer than 1 keV with different IDs stay two levels and print
`WARNING: level ID ... is within 1 keV` (`docs/source/user_guide/levels_channels.rst`).
builds before 2026-10-24 (expires: 2027-02-01) merged any same-J^π lines within 1 keV into one level whatever their
IDs: the second level's channels became the first's, its parameters disappeared from the fit, and a model written out
and read back had a different level count. A fit that drives two same-J^π levels within 1 keV on such a build cannot
round-trip through a `.azr`; rmfit's bake guard reports it as lost keys.

## Limits and priors written by the Qt GUI

The `<parameterSettings>` of the Fitting tab of released builds' Qt GUI names `Level N` as `energy_N` and `Level N Channel M` as `width_N_M` in the
engine's numbering, and marks the block `# level numbering: engine` (`docs/source/user_guide/fitting.rst`).
builds before 2026-10-21 (expires: 2027-02-01) numbered GUI settings by model row (inactive levels counted, widths by
global channel row) while the engine read engine numbers, so after a GUI save of a file not sorted by J, parity and
energy, or with an inactive level ahead of active ones, limits and priors landed on other levels (a prior centred on
the 3.52 MeV level acted on the 5.0 MeV one, χ² +21,904). A current build reads an unmarked block in whichever
numbering makes more value columns equal the parameters' current values (tie: engine). To audit an old file: is
`Level N Channel M` (M = global channel-line row) the engine's `width_N_M`? The GUI keeps settings only for the
parameters it lists as free, so script-made bounds on fixed parameters vanish on a GUI save.
