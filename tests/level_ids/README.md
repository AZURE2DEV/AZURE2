# level_ids

Level identity in the `<levels>` block: lines of one J^pi group belong to one level when they carry the same level
ID (column 9) within 1 keV, or the identical energy (`CNuc::Fill`).  Before 2026-10-24 any lines within 1 keV were
one level whatever their IDs (`JGroup::IsLevel`, tolerance 1e-3 MeV), so a narrow level next to a broad one of the
same J^pi was merged into it: the later lines' widths replaced the earlier level's and its parameters vanished from
the fit (17O+a/10-22-26_an3_jpiscan: 3/2- 8.2900 MeV, ID 91, swallowed by 3/2- 8.2893 MeV, ID 22; 32 instead of 35
free parameters).

`check.sh` edits the 13N reference model (3/2- ID 5 at 3.5032 MeV) at run time:

| model  | edit                                                         | expected                                   |
|--------|--------------------------------------------------------------|--------------------------------------------|
| zero   | ID 13 at 3.5039 MeV, zero widths, lines last                 | chiSquared.out = reference                 |
| narrow | ID 13 at 3.5039 MeV, 50 eV / -0.01 eV, lines last            | same as `first`, different from reference  |
| first  | the same lines first in `<levels>`                           | same as `narrow`                           |
| same   | ID 13 at exactly 3.5032 MeV with ID 5's widths               | one level (warning), chiSquared.out = ref. |

The build before the fix fails 9 checks (zero: total chi2 1.92e6 instead of 336480; narrow 1.92e6 vs first 315588).
