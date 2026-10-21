# ec_shift — external capture at shifted energies

`check.sh` (run by `run_tests.sh`) checks the external-capture energy tables
(`ECTable`, on by default; `--no-ec-table` turns them off) on the 13N reference
model (`../13N`), with target integration off and only the capture segments 1
(Artemov), 6 and 8 (Skowronski, LUNA HPGe and FELS) and 10 (Kettner, 55°) active:

1. **Shift = moved data.** The segments with an energy shift of +30 and −30 keV
   (lab) give the same cross sections as their data files with every lab energy
   moved by the same amount and no shift (limit 1e-7; the result is 0 to the
   printed digits).  Both runs evaluate the capture amplitudes at the same c.m.
   energies from the same tables.
2. **The shift has an effect** (more than 1e-3), so (1) proves something.
3. **Table vs converged integrals.**  `AZURE2_EC_TABLE_CHECK=1` computes every
   tabulated integral directly as well, converged to 1e-10, and the largest
   relative difference must stay below 1e-7 (it is 1.6e-10 at the default
   50 nodes per decade).
4. **Table vs `--no-ec-table`** without a shift, within 1e-3 (1.1e-6): the
   direct integrals are converged only to 1e-4 relative, which is where the two
   differ.

The `--no-ec-table` difference between shifted and moved runs is printed for
reference: before the tables a shifted point interpolated linearly between the
amplitudes of its segment's own points and held them flat beyond the segment's
ends — 9 % off for +30 keV and a factor 57 for −30 keV, where the lowest LUNA
points move below the segment's lowest energy and the steeply falling amplitude
was held constant.
