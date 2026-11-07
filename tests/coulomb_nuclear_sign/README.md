# coulomb_nuclear_sign

The sign of the Coulomb–nuclear interference term in elastic scattering, checked against an independent
partial-wave sum (`docs/source/theory/conventions.rst`, "Elastic scattering").

AZURE2 forms the elastic cross section as |C_α|² (π/k²) + Σ B_L P_L / k² + an interference term
(i/√π)(π/k²) Σ g_J C_α(θ) T*_cc P_l(cos θ) (`GenMatrixFunc.cpp`), with C_α = (4π)^(-1/2) η csc²(θ/2)
exp(−2iη ln sin(θ/2)) and T = e^{2iω} − U. Whether that combination equals the physical
|f_C + f_N|² depends on three sign conventions at once (the sign of C_α relative to f_C, the i, and T = e^{2iω} − U
rather than U − e^{2iω}); none of them is visible in an angle-integrated cross section or in a calculation without
Coulomb field.

**Model** (`models/hard_sphere_pc.azr`): `tests/hard_sphere/models/no_dummies.azr` (p + ¹²C, a = 5.0 fm,
`includeAllJPi=1`, l ≤ 4) with the three level widths set to zero, so every partial wave is pure hard-sphere
scattering, U_l = exp(2i(ω_l − φ_l)), independent of J and of the parametrization.

**Truth** (`truth.py`, mpmath; run once, output committed as `expected/dsdo.dat`): f = f_C + (1/2ik) Σ (2l+1)
e^{2iω_l}(e^{−2iφ_l} − 1) P_l with f_C = −η/(2k sin²(θ/2)) exp(−iη ln sin²(θ/2)), on the c.m. grid AZURE2 writes for
the test segment (E_lab = 1–3 MeV, θ_lab = 60–160°; 30 points).

**Checks** (`check.sh path/to/AZURE2`): A-matrix and R-matrix routes, `--use-brune` and `--no-brune`; every point
within 1e-6 relative of the truth. On `feature/upgrade` f8b7bf22 all four agree to 8.4e-11. With the interference
term reversed the truth moves by 7–74 % on this grid.
