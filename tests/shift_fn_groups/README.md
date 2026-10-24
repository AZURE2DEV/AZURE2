# shift_fn_groups — one shift function shared by segments with their own constants

A segment line may end with

    shiftfngroup <label>

Every segment carrying the label applies the shift function (`shiftfn`) declared
on the FIRST of them (file order), with that segment's parameters
(`segment_<first>_energy_shift_fn<j>`, penalised once), while each keeps its own
constant energy shift: dE(x) = a_k + f(x; a0, ...). `shiftgroup` shares the
constant as well; use `shiftfngroup` when the per-segment constants are real
(e.g. separate figures digitized with separate offsets) but the energy-dependent
calibration is one (11B+a 10-22-26: one centred sqrt(E) term per experiment).
A segment may not be in both kinds of group; only the first member declares the
function.

Model and data: `tests/shift_fn` (p+p model, E_cm 0.2-4.8 MeV). The check script
also writes `xs_90_lin_mislabelled_p3.dat`, the 90 deg mislabelled data with every
label moved by +3 keV, so that the 40 and 90 deg segments share a0 = 0.006 in
`a0*(x-1)` but have constants -0.004 and -0.004 - 0.003*(1 + 0.006) = -0.007018.

    ./tests/shift_fn_groups/check.sh path/to/AZURE2
