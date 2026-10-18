# spin_observables_p3He — isDiff 9 on spin-1/2 x spin-1/2 (p + 3He)

Fifteen spin observables of a p + 3He toy with l-mixing (2-, 1+) and singlet-triplet mixing (1-, 2+):
beam and target analyzing powers (`b:y`, `t:y`, `b:iT11`, `t:iT11`), spin correlations (`b:y*t:y`,
`b:x*t:x`, `b:z*t:z`, `b:x*t:z`, `b:z*t:x`, and the parity-forbidden `b:x*t:y`), outgoing
polarizations (`e1:y`, `e2:y`) and transfer coefficients (`b:y>e1:y`, `t:y>e2:y`, `b:y>e2:y`).

The data are the values of an **independent** calculation (mpmath Coulomb functions, sympy
Clebsch-Gordan coefficients, its own S-matrix and spin operators; fit archive
`Aeql4/10-18-26_spin_impl/check/indep_spin.py`) plus an alternating offset of +-0.01, with errors 0.01.
Every point therefore contributes chi2 = 1 when the engine agrees with that calculation, and the
recorded chi2/N is 1.0000 for every segment (the engine agrees to 1e-6). A sign or convention error
moves a segment's chi2 by orders of magnitude.

Grid: E_lab 2, 6, 10 MeV, theta_cm 30-150 deg in 20-deg steps (21 points per segment). Built by
`check/make_tests.py`.
