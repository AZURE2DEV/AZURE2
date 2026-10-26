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
| (k), (l) | p+14N M1, both cluster moments, final channel spin 1/2 and 3/2 | 2e-5 |
| (m) | an `intEC.dat` from before the M1 spin-term fix (`models/stale_v2`, case (k) by dev 823369d) is not reused | identical to fresh |
| (n), (o) | p+14N M1 from p-waves into a p-wave state: orbital + spin terms; orbital term alone (g = 0) | 2e-5 |
| (p) | an `intEC.dat` from before the M1 orbital-term fix (`models/stale_v3`, case (n) by dev c23a2d3) is not reused | identical to fresh |

AZURE2 agrees with the truth to about 1e-8.

**What it caught (fixed 2026-10-23, dev 41a468b).** Before the fix, (d) failed by 6.25 %: `src/ECIntegral.cpp` had
`sqrt(16 / 3)`, which integer arithmetic evaluates as sqrt(5). The intended value is the L = 1 multipole factor
sqrt(8(2L+1)(L+1)/L)/(2L+1)!!. (e) failed by a factor 0.294, because the E1 charge function (x0.560 for p + 16O)
multiplied the M1 integrand. (i) failed because the intEC signature did not include the long-wavelength switch.
Both M1 defects date from the 2020 initial commit. The E1 terms are unchanged by the fix.

`p14n_s1.azr` and `p14n_s3.azr` are p + 14N -> 15O, M1 external capture from the s-wave entrance J-groups 1/2+ and
3/2+ into a fictitious s-wave state bound by 0.2 MeV (ANC 5 fm^-1/2, a = 5 fm). Both clusters carry a moment
(g_p = 5.5857, g(14N) = 0.4038) and the target has spin 1, so the M1 spin term has a g1 S1 and a g2 S2 piece and
connects channel spins 1/2 and 3/2. With l = 0 the orbital term vanishes. The truth builds the spin matrix element
from magnetic substates (`python3 truth.py p14n`, about three minutes, writes `expected/sigma_p14n.dat`).

**What it caught (fixed 2026-10-25).** Before the fix, (k) failed by 47 % and (l) by 17 %. The `tau` term in
`src/ECIntegral.cpp` lacked a factor (-1)^(2 j2) on the g1 S1 piece and (-1)^(2 j1) on the g2 S2 piece, so a cluster's
moment term had the wrong sign whenever the other cluster had half-integer spin. When both clusters have
half-integer spin (p+13C, p+15N, p+7Li) the whole spin term had the wrong sign against the orbital term. This test
cannot see that case, because with l = 0 there is no orbital term. Spin-0 targets (p+16O, case d) were right. The
spin and orbital amplitudes now equal an m-scheme calculation times one common constant, 1/sqrt(2 J_i + 1), in all
5,277 combinations of spins 0 to 2, l = 0 to 2 and both g-factors (R-matrix archive `AZURE2_code_audit/10-25-26/fix01_m1_spin`). The
first push of the fix (dev 0465351) still reused an `intEC.dat` written before it, because the signature version
was unchanged, and gave (k) 47 % off again; (m) checks the bump to signature v3.

`p14n_p.azr` and `p14n_p0.azr` are p + 14N -> 15O(3/2-), M1 external capture from the p-wave entrance J-groups 1/2-,
3/2- (channel spins 1/2 and 3/2) and 5/2- (s = 3/2) into a fictitious p-wave state bound by 0.2 MeV (ANC 5 fm^-1/2
in the s = 3/2 channel, 0 in s = 1/2, a = 5 fm). With l = 1 the orbital term g_l L is present (g_l = m_p mu
(Z1/m1^2 + Z2/m2^2) = 0.96) and interferes with the spin term on the s = 3/2 pathways; the s = 1/2 entrance
channels change the channel spin, which only the spin term can do. `p14n_p0` sets both g-factors to 0. The truth
builds the full operator g_l L + g_p S_p + g_14N S_14N from magnetic substates (`python3 truth.py p14n_p`, about
five minutes, writes `expected/sigma_p14n_p.dat`).

**What it caught (fixed 2026-10-26).** Before the fix, (n) was 39 % high and (o) exactly 2x: `src/ECIntegral.cpp`
added the orbital term on every M1 pathway, including those with s_f != s_i, where L cannot act. The spurious term
needs a final state that couples to both channel spins at the same l (J_f = 3/2 here; with J_f = 5/2 and s_i = 1/2 it
vanishes by a triangle rule, which is why a first version of this model saw nothing), so l >= 1 on both sides and
two cluster spins: p+14N, p+13C, p+15N, p+10B, p+6Li, p+23Na M1 capture into p- or d-wave states, never s-wave
capture (the orbital term has sqrt(l(l+1))). (p) checks the signature bump to v4.

**A trap in the input.** The ANC-to-amplitude conversion uses the photon pair's final-state energy (`e2` on the
pair line). The radial integral uses the level energy. The 16O+p model these numbers come from had 0.49533 MeV
on the pair and 0.4953 MeV on the level; that 0.03 keV difference moves the external capture by 0.16 %. Keep the
two equal.
