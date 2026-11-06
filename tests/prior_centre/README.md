# prior_centre — explicit centres of the normalization and shift priors

A `<segmentsData>` norm (shift) field is the start value of a run and, in a
classic file, also the centre of the segment's normalization (energy-shift)
prior.  A `<parameterSettings>` row

    segment_N_norm prior_centre c
    segment_N_energy_shift prior_centre c

(`N` the segment key) moves only the centre (`EData::ReadPriorCentres`, which
sets `ESegment::SetNominalNorm` / `SetNominalEnergyShift`; every penalty reads
`GetNominalNorm` / `GetNominalEnergyShift`).  A file can then store fitted
norms and keep its priors on the experimental values.  Older readers skip
three-token rows.

`check.sh` runs `tests/15N_p_a` in Calculate mode (segment 3: norm
0.9950877877494786, 15 % error):

| case | rows | expected |
|---|---|---|
| none | — | classic total, norm term 0 |
| same | centre = the field | stdout and `chiSquared.out` byte-identical to none |
| one | `segment_3_norm prior_centre 1` | + ((n - 1)/0.15)^2 = 0.00107244 |
| shift | shift free, error 0.01 MeV, `segment_3_energy_shift prior_centre 0.005` | + 0.25 |
| stale | `segment_9_norm prior_centre 1` | warning, classic total |
| bad | `... prior_centre abc`, norm centre 0 | refused |

    ./tests/prior_centre/check.sh path/to/AZURE2
