# ec_m1_exact — external capture against an independent calculation

Two systems: p + 16O (E1 in three forms, and M1) and alpha + 12C (E1 + E2 in three forms). The three electric
forms are the long-wavelength operator, the extended-Siegert charge term (`--no-long-wavelength
--no-ec-convection-current`), and charge plus convection current (`--no-long-wavelength`, the default exact form
since 2026-10-23; derivation in `docs/source/theory/external_capture_operator.rst`).

`a12c.azr` is 12C(a,g)16O(g.s.) with E1 (l = 1) and E2 (l = 2) external capture, ANC 100 fm^-1/2 and a = 5.43 fm.
alpha + 12C is an N = Z pair, so the long-wavelength E1 charge is only -0.00098; there the charge term is 3.5-4.7
times the LWA E1 amplitude, and the convection current adds another +10 %.

`p + 16O -> 17F(0.495 MeV, 1/2+) + gamma`, hard-sphere external capture only. The final state is an s-wave proton
bound by 104 keV with ANC 84 fm^-1/2. The entrance J-groups 1/2-, 3/2- (E1) and 1/2+ (M1) each hold a single 30 MeV
level with a 1e-6 eV proton width, so the resonant terms are negligible. The photon pair's multipolarity mask
selects E1 (`e1.azr`) or M1 (`m1.azr`).

`truth.py` computes the angle-integrated cross sections without any AZURE2 code. It uses the golden rule with an
energy-normalized continuum, sigma = (pi/k^2) sum_J g_J 2 pi Gamma_J. The bound state is u_f = C W(2 kappa r), and
the continuum is mpmath Coulomb functions with the hard-sphere phase. It uses the textbook one-body E1 and M1
operators, with reduced matrix elements in (l s)J coupling. The finite-k_gamma E1 operator is the charge term of the
extended Siegert theorem (Friar and Fallieros, PRC 29, 1645 (1984)). For point clusters this is
`3/k [Z1 g(k r m2/M) - Z2 g(k r m1/M)]`, with `g(y) = int_0^y j_1(t)/t dt` in closed form; the same operator appears
in Canton and Levchuk, arXiv:0805.2667, Eqs. 33 and 39-41. The script writes `expected/sigma.dat` and
`expected/sigma_a12c.dat`. It takes about 20 minutes, so the files are committed.

`check.sh` runs the models in each form and checks:

| | check | tolerance |
|---|---|---|
| (a) | p+16O E1, long wavelength | 2e-5 |
| (b) | p+16O E1, `--no-long-wavelength` (charge + convection current) | 2e-5 |
| (c) | p+16O E1, `--no-long-wavelength --no-ec-convection-current` (charge term) | 2e-5 |
| (d) | M1 against the truth | 2e-5 |
| (e) | M1 with `--no-long-wavelength` against the same truth | 2e-5 |
| (f), (g), (h) | alpha+12C E1 + E2 in the three forms | 2e-5 |
| (i) | an `intEC.dat` written by a long-wavelength run is not reused by a `--no-long-wavelength` run | identical to fresh |
| (j) | an `intEC.dat` written with the convection current is not reused without it | identical to fresh |

AZURE2 agrees with the truth to about 1e-8.

**What it caught (fixed 2026-10-23, dev 41a468b).** Before the fix, (d) failed by 6.25 %: `src/ECIntegral.cpp` had
`sqrt(16 / 3)`, which integer arithmetic evaluates as sqrt(5). The intended value is the L = 1 multipole factor
sqrt(8(2L+1)(L+1)/L)/(2L+1)!!. (e) failed by a factor 0.294, because the E1 charge function (x0.560 for p + 16O)
multiplied the M1 integrand. (i) failed because the intEC signature did not include the long-wavelength switch.
Both M1 defects date from the 2020 initial commit. The E1 terms are unchanged by the fix.

**A trap in the input.** The ANC-to-amplitude conversion uses the photon pair's final-state energy (`e2` on the
pair line). The radial integral uses the level energy. The 16O+p model these numbers come from had 0.49533 MeV
on the pair and 0.4953 MeV on the level; that 0.03 keV difference moves the external capture by 0.16 %. Keep the
two equal.
