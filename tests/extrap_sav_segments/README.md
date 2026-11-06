# extrap_sav_segments

A parameter file's `segment_*` entries are not applied in a run without data.

An extrapolation (mode 3, `<segmentsTest>`) and a reaction-rate run number their test grids 1..N and
`EData::FillMnParams` emits `segment_<k>_energy_shift`, `_energy_shift_sqrt`, `_energy_shift_fn<j>` and
the normalization for each of them, so a `param.sav` written by a fit matched DATA segment k's fitted
shift to TEST segment k by name (`AZUREParams::ApplyUserParameterFile`), and
`EData::FillEnergyShiftsFromParams` applied it (a nonzero constant shift counts as active): the
requested grid moved by the shift, converted to c.m. -- the energy-shift twin of the `<targetInt>` trap
(audit 2026-10-25, item 12).  `ReadUserParameters` now passes whether the run has data; without data the
`segment_*` lines are skipped and the run log says how many.  pyazr's `load_sav` does the same in
extrapolation mode (it pushes the whole vector into the engine, fixed shifts included, so it moved the
grid too).

`check.sh` (model: tests/sum_lab_frame/models/extrap.azr): (a) mode 3 with a `.sav` whose
`segment_1_energy_shift` is 0.01 gives the same grid as without it (values to 1e-6: the file carries
the level parameters with 7 digits), and the note is in the log; (b) with data (mode 1 on a one-point
segment) the same file still shifts the point, as it must; (c) pyazr `load_sav` in extrapolation mode
leaves `calculate_energies` unchanged.  dev 81abfd4: (a) the first grid energy moved from 4.5496 to
4.5572 MeV c.m. (0.01 lab x 13/17).
