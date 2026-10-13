# shift_groups — one constant energy shift shared by several segments

A segment line may end with `shiftgroup <label>`.  Every segment carrying the
same label uses ONE constant energy shift: the parameter of the group's first
segment in file order (its value, uncertainty, vary flag and penalty).  The
other members' own shift value, uncertainty and vary flag are ignored, their
`segment_<key>_energy_shift` is fixed in the fit and written out with the
shared value.  The sqrt(E) term (`sqrtshift`) is not shared.  Introduced
2026-10 for data sets split into per-angle segments that share one beam-energy
calibration (12C+a: the Frawley 15N+p excitation functions, 37 segments).

Model: tests/identical_pp_res via tests/energy_shift_sqrt (p+p, E_cm 0.4-5
MeV, its own values as data; chi2 = 24 per unshifted segment).

`shift_groups.azr`: `xs_40.dat` with a = 0.004 MeV (fixed) and `xs_90.dat`
with its own a = 0, both in group G, so segment 2 is shifted by 0.004 too.
`expected/chiSquared.out` pins both.

`check.sh`:

1. the grouped pair gives exactly the chi2 of the same two segments each given
   a = 0.004 directly (1e-10), segment 2 is really shifted, and the log notes
   the group;
2. `xs_40_off6.dat` and `xs_90_off6.dat` carry labels L = E + 0.006 and the
   model's exact yields; with every level fixed and the group's shift free
   (the follower's own line says "free from 0.01", ignored), MIGRAD must
   return a = -0.006 (|error| < 1e-5), chi2 = 0, the follower written with the
   same value in `param.sav`, and both rows in `shifts.out` -- once with
   Minuit's numerical gradient and once with `--use-gradient` (the group
   leader's shift takes the full-dataset finite difference there, since it
   moves the follower's points too);
3. `xs_40_off6.dat` with `xs_90_off2.dat` (labels 2 keV high): one shift for
   both, strictly between -0.006 and -0.002.
