# spin_observables_p3He — the quantity spin_observable on spin-1/2 x spin-1/2 (p + 3He)

Fifteen spin observables of a p + 3He toy with l-mixing (2-, 1+) and singlet-triplet mixing (1-, 2+):
beam and target analysing powers (`b:y`, `t:y`, `b:iT11`, `t:iT11`), spin correlations (`b:y*t:y`,
`b:x*t:x`, `b:z*t:z`, `b:x*t:z`, `b:z*t:x`, and the parity-forbidden `b:x*t:y`), outgoing
polarizations (`e1:y`, `e2:y`) and transfer coefficients (`b:y>e1:y`, `t:y>e2:y`, `b:y>e2:y`).

**Truth.** The data are the values of an independent calculation (mpmath Coulomb functions, sympy
Clebsch-Gordan coefficients, its own S-matrix and spin operators; in the owner's fit archive,
`Aeql4/10-18-26_spin_impl/check/indep_spin.py`) plus an alternating offset of +-0.01, with errors
0.01. Every point therefore contributes chi2 = 1 when the engine agrees with that calculation, and
`expected/chiSquared.out` records chi2/N = 1.0000 for every data set (the engine agrees to 1e-6).
A sign or convention error moves a data set's chi2 by orders of magnitude.

Grid: E_lab 2, 6, 10 MeV, theta_cm 30-150 deg in 20-deg steps (21 points per data set).

**Origin.** dev's `tests/spin_observables_p3He` (isDiff 10 lines), written in format 3 by
`AZURE2 --convert` (0 notes, 0 questions); the expected file is this engine's and equals dev's
(1e-5 of a chi2 printed with six digits).

**Run.** `tests/run_tests.sh` compares `output/chiSquared.out` with `expected/`; the analytic
gradient of this project is checked by `tests/spin_observables_reaction/check.sh`.
