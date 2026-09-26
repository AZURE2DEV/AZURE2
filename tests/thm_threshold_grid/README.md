# thm_threshold_grid — a THM curve through the entrance threshold

`check.sh` evaluates the THM (HOES) excitation functions of `tests/6Li_d`
(6Li+d) and `tests/7Li_p_a` (p+7Li) on a 1 keV grid over E_cm = 0–0.3 MeV,
with their 30 keV (sigma) Gaussian folding, and passes when every value is
finite and positive and neighbours differ by less than 2 %.

The folding window of a THM segment is ±5 sigma and may extend below the
entrance threshold (EData.cpp), so for every point below 0.15 MeV the
sub-point grid crosses E = 0. Before the fix, a sub-point landing at
round-off distance from threshold (E ≈ −1e-16 MeV) got a NaN shift function
— the plain-double Whittaker function underflows there, and for
|E| < ~1e-5 MeV it returned S = 0 — and the folded cross section of each
affected point came out at 1e16–1e20 against neighbours of ~0.005 (on 6Li_d
every multiple of 9 keV below 0.11 MeV). The negative-energy shift function
is now the analytic Whittaker log-derivative, continued smoothly to the exact
zero-energy limit S_l(0) = −l − (x/2) K_{2l}(x)/K_{2l+1}(x); and a positive
channel energy with Sommerfeld parameter above 100 (P < 1e-273) is treated as
at threshold, with S continued as 2 S(0) − S(−E), since the Coulomb-wave
routines return NaN or garbage there. See `ShftFunc` and `IsCoulombThreshold`
in `src/EPoint.cpp`.

The models are the two pinned THM projects unchanged except for the data file;
this checks behaviour, not a chi-squared.
