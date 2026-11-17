# spin_observables_dd — the quantity spin_observable on identical spin-1 bosons (d + d)

Nineteen spin observables of a d + d toy (identical bosons: only l + s even; levels 0+, 2+ with
s- and l-mixing, 1-, 2-, 0-, 3-): vector and spherical-tensor analysing powers of beam and target
(`b:y`, `t:y`, `b:iT11`, `b:T20/21/22`, `t:T20/21/22`), Cartesian tensor powers (`b:yy`, `b:xx`,
`b:zz`, `b:xz`), spin correlations including tensor x tensor (`b:y*t:y`, `b:x*t:x`, `b:z*t:z`,
`b:yy*t:yy`), outgoing polarization (`e1:y`) and transfer (`b:y>e1:y`).

**Truth.** As in `spin_observables_p3He`, the data are an independent calculation (here with the
identical-particle symmetrization: nuclear pathways x (1 + (-1)^(l'+s')), Coulomb
f_C(theta) + (-1)^s f_C(pi - theta)) plus an alternating +-0.01 offset with errors 0.01, so
chi2/N = 1.0000 per data set records agreement (1e-7). Beam and target observables of identical
particles come out equal, as they must.

**Origin.** dev's `tests/spin_observables_dd`, written in format 3 by `AZURE2 --convert` (0 notes,
0 questions); the expected file is this engine's and equals dev's.

**Run.** `tests/run_tests.sh` compares `output/chiSquared.out` with `expected/`; the analytic
gradient of this project is checked by `tests/spin_observables_reaction/check.sh`.
