# 7Li_p_a — the THM validity test of 7Li(p,alpha), pinned

## Data provenance — read this first

**This project is a regression pin, not a physics benchmark.** Its chi2 says
the engine has not moved; it says nothing about how well AZURE2's THM
formalism describes these data, because the data are not what the model
assumes:

* **What the points are.** EXFOR O1653002 (Tumino et al., EPJA 27 s01, 243
  (2006)) is the 7Li(p,alpha) cross section extracted from the *3He* breakup
  3He(7Li,alpha alpha)2H at 33 MeV, |p_s| <= 30 MeV/c, theta_cm = 50-70 deg,
  E_cm resolution 80-120 keV (FWHM). In the published excitation function the
  "penetrability effects" (s + p waves, Coulomb + centrifugal) were "included
  before the comparison with direct data", and it was normalized to direct data (Engstler 1992, Mani 1964, Cassagnou 1962) at
  E_cm = 2-3 MeV (pp. 244-247 of the paper). It is a penetrability-corrected,
  direct-normalized, on-shell-equivalent cross section, **not a raw
  half-off-shell (HOES) excitation function** — yet the project fits it with
  the HOES observable.
* **The carrier.** The Trojan horse was 3He = p (x) d: the transferred particle
  is the proton and the spectator the deuteron, with **B(3He -> p + d) =
  5.4935 MeV** (AME2020: 7288.971 + 13135.723 - 14931.219 keV). The project
  carries **B = 2.2246 MeV**, the deuteron binding of a d(7Li,alpha alpha)n
  measurement (Lattuada et al. 2001 / Lamia et al. 2012), which is not where
  these points come from.
* **What is pinned.** The recorded chi2 (`expected/`) keeps B = 2.2246 MeV, as
  it always has. The same model with the 3He binding is pinned alongside, in
  `tests/thm_options/check.sh` (chi2 = 1827.89 on 2026-09-26); `run_tests.sh`
  takes the first `.azr` of a directory, so it cannot live here as a second
  file.

Source of this analysis: the THM refits of 2026-09-25/26 (outside the repo,
`thm_refits/REFITS.md` §1 and §9, `literature_values.md` §1a), where treating
the same points as an on-shell cross section with a free scale gave THM
chi2/N = 3.5 against 12.5-12.8 under HOES (with either binding), with direct
norms returning to ~1. A faithful HOES treatment would need the uncorrected
THM yields, which are not in EXFOR.

## The project

The 7Li(p,alpha)4He excitation function extracted by the Trojan Horse method
(from 3He(7Li,alpha alpha)2H breakup, not the d(7Li,alpha alpha)n breakup this
file used to name — see *Data provenance* above), evaluated here as if it were
the half-off-shell one:

A. Tumino et al., Eur. Phys. J. A 27 (S1) 243 (2006) — the THM *validity
test* of this reaction — EXFOR **O1653002**, 66 points,
E_cm = 0.076-6.87 MeV, arbitrary units, theta_cm = 50-70 deg.

Model: the validated Paneru et al. 8Be evaluation (PRC 111, 064609,
Table IV; the same base as `6Li_d` and `7Li_p_ay`), entrance pair 5
(7Li+p), exit the identical-boson alpha+alpha pair. The model takes the
Trojan horse to be the deuteron (it was 3He — see *Data provenance*), so the
entrance-pair channel lines carry B = 2.2246 MeV in the optional 33-field
column — the complementary case to `6Li_d`, whose
carrier is 6Li at 1.4735 MeV. A 30 keV (sigma) Gaussian resolution is
folded through the engine's sub-point convolution on the THM segment.

EXFOR carries no experimental uncertainties for this entry (only a
digitizing error), so the data files hold an assigned 10% (0.01 absolute
floor); the c.m. energies were converted to the lab convention the data
files use. The segment normalization is the THM arbitrary scale, solved
analytically against this model and written into the file.

The recorded chi2 = 4952.9 over 66 points is a **pin, not a fit**: the
Paneru parameters come from on-shell data alone, and the half-off-shell
observable weights the same poles differently. The two-peak structure
(E_cm = 2.6 and 5.1 MeV) sits at the right energies; the misfit
concentrates below 0.3 MeV, where the HOES observable rises toward
threshold against falling data, and in the relative peak heights — the
behaviour a joint THM+direct refit corrects (see the analysis notes).

Note (thm branch): the THM segment's arbitrary overall scale is now profiled out analytically by the engine (a free THM norm is set to its chi2-optimal value each evaluation), so the pinned chi2 is 2740.48 at the optimal scale -- superseding the earlier hand-solved-norm value (4952.9). This is the correct minimal chi2 for the arbitrary THM scale.

Note (Sep 2026): entrance partial waves are now summed incoherently in l (the l cross terms vanish once the exit direction is integrated and the spin projections summed; see src/THMMatrixFunc.cpp). This case has two entrance l in one channel spin, so the pin moves from 2740.48 to 2180.69 with the same parameters. A `<thm>` block with `entranceL=coherent` restores the old sum and reproduces 2740.48.

Note (2026-09-25): the default THM entrance vertex changed from `vertex=perlevel` (per-level S_c(E_lambda) under Brune, as mrmpy) to `vertex=constant` (the formal-R-matrix vertex with the channel boundary B_c = S_c(E_1), applied after the level sum; see docs/source/theory/thm_implementation.rst). Same parameters; the pin moves from 2180.69 to 2137.83 (optimal norm 0.00143722 -> 0.00189034). A `<thm>` block with `vertex=perlevel` reproduces the previous pin, 2180.69 (checked on 2026-09-25; tests/thm_options pins it).

Note (2026-09-25, later): `vertex=constant` first took B_c = S_c at the first level of each J group in file order, so the result depended on level order (this file sorted as the GUI writes it gave a different chi2). It now uses S_c at the lowest-energy level of each J group, independent of order; the pin above is that value.
