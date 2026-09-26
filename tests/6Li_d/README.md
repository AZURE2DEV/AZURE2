# 6Li_d — a second Trojan-Horse regression, different carrier

## Data provenance — read this first

**This project is a regression pin, not a physics benchmark.** Its chi2 says
the engine has not moved; it says nothing about how well AZURE2's THM
formalism describes these data, because the data are not what the model
assumes:

* **What the points are.** EXFOR D0649002 (Pizzone et al., PRC 83, 045801
  (2011)) is the 6Li(d,alpha) cross section extracted from the *3He* breakup
  6Li(3He,alpha alpha)H at 17.5 MeV (|p_s| < 35 MeV/c). The THM yields were
  divided by the kinematic factor lambda_3/lambda_2 |phi|^2, then "corrected
  for the penetrability factor" with an l = 0 term and an l = 2 term (the 2+
  near 25 MeV), each with its own normalization constant, fixed against direct
  data (McClenahan 1975, Mani 1964, Elwyn 1977, Engstler 1992) over
  E_cm = 0.4-5 MeV (§V-VI of the paper). They are penetrability-corrected,
  direct-normalized, on-shell-equivalent cross sections — **not a raw
  half-off-shell (HOES) excitation function** — and the project fits them with
  the HOES observable. (Treated as an on-shell cross section their free scale
  comes out at 1.11: they are essentially absolute.)
* **The carrier.** The Trojan horse was 3He = d (x) p: the transferred
  particle is the deuteron and the spectator the proton, with
  **B(3He -> d + p) = 5.4935 MeV** (AME2020: 7288.971 + 13135.723 -
  14931.219 keV). The project carries **B = 1.4735 MeV**, the alpha + d
  binding of 6Li — the carrier of the *6Li*-breakup measurement
  6Li(6Li,alpha alpha)4He (Spitaleri et al. 2001), which is not where these
  points come from. The section below, written when the file was made,
  describes the model as built, not the experiment.
* **What is pinned.** The recorded chi2 (`expected/`) keeps B = 1.4735 MeV,
  as it always has. The same model with the 3He binding is pinned alongside,
  in `tests/thm_options/check.sh` (chi2 = 483.474 on 2026-09-26);
  `run_tests.sh` takes the first `.azr` of a directory, so it cannot live here
  as a second file.

Source of this analysis: the THM refits of 2026-09-25/26 (outside the repo,
`thm_refits/REFITS.md` §1 and §9, `literature_values.md` §1b). For these points
chi2 does not discriminate HOES from on-shell (0.98 vs 1.32 per point with 23
free parameters), so the provenance, not the fit, is what settles it.

## The project

The THM excitation function of 6Li(d,alpha)4He, on the same 8Be compound as
`7Li_p_ay` — the validated Paneru et al. model (PRC 111, 064609, Table IV
parameters) — entered through the other door: entrance pair 1 (6Li+d), exit
the identical-boson alpha+alpha pair.

Data: R. G. Pizzone et al., PRC 83 (2011) 045801 (EXFOR D0649002), the
THM-derived two-body cross section, 62 points, E_cm = 0.43-4.96 MeV. EXFOR
carries no uncertainties for the entry, so the files hold an assigned 10%.

## What this covers that tests/17O does not

* **A different carrier binding energy.** 17O's Trojan horse is the deuteron
  (B = 2.225 MeV); this model takes the transferred particle to be the
  deuteron bound in 6Li = alpha (x) d with **B = 1.4735 MeV** (the experiment
  used 3He, B = 5.4935 MeV — see *Data provenance*), written on the
  entrance-pair channel lines in field 32, the first of the two optional THM
  fields (field 33 is the gammaIsRWA flag). A
  hard-coded deuteron binding would pass 17O and fail here.
* **The energy-resolution convolution on a THM segment.** The `<targetInt>`
  line folds a 30 keV Gaussian (sigma) through the engine's sub-point
  machinery, as the 17O case does with 21 keV — the convolution range for THM
  segments follows `TargetEffect::thmConvolutionRange`.
* **An identical-particle exit** under the THM observable.

## The recorded chi-squared

The segment's normalization is the THM arbitrary scale, solved analytically
against this model (engine convention: the norm multiplies the *data*) and
written into the file; chi2 = 2270.61 over 62 points is a **pin, not a fit**
— the Paneru parameters were fit to on-shell data only. The number is large
for a physical reason worth keeping visible: the data show a broad structure
near E_cm = 3 MeV (Ex = 25.3 MeV, the 2+ 25.72 region) that the on-shell
model barely produces through the half-off-shell observable. A refit that
frees the 2+ 22.98 / 2+ 25.72 / 0+ 27.49 widths reduces it — see the
analysis in the evaluation notes — but the pinned reference deliberately
stays at the published parameters.

Note (thm branch): the THM segment's arbitrary overall scale is now profiled out analytically by the engine (a free THM norm is set to its chi2-optimal value each evaluation), so the pinned chi2 is 1662.03 at the optimal scale -- superseding the earlier hand-solved-norm value (6135.89). This is the correct minimal chi2 for the arbitrary THM scale.

Note (Sep 2026): entrance partial waves are now summed incoherently in l (the l cross terms vanish once the exit direction is integrated and the spin projections summed; see src/THMMatrixFunc.cpp). This case has two entrance l in one channel spin, so the pin moves from 1662.03 to 1076.06 with the same parameters. A `<thm>` block with `entranceL=coherent` restores the old sum and reproduces 1662.03.

Note (2026-09-25): the default THM entrance vertex changed from `vertex=perlevel` (per-level S_c(E_lambda) under Brune, as mrmpy) to `vertex=constant` (the formal-R-matrix vertex with the channel boundary B_c = S_c(E_1), applied after the level sum; see docs/source/theory/thm_implementation.rst). Same parameters; the pin moves from 1076.06 to 682.099 (optimal norm 0.00735648 -> 0.0319074). A `<thm>` block with `vertex=perlevel` reproduces the previous pin, 1076.06 (checked on 2026-09-25).

Note (2026-09-25, later): `vertex=constant` first took B_c = S_c at the first level of each J group in file order, so the result depended on level order (this file sorted as the GUI writes it gave a different chi2). It now uses S_c at the lowest-energy level of each J group, independent of order; the pin above is that value.

Re-pinned 2026-09-26: the shift-function fix b6cc41b removed ~1e-9 noise from S(E) that ShftFunc::EnergyDerivative amplified into dS/dE, which enters the Brune transformation of sub-threshold levels. 682.099 -> 682.133 (5.0e-5).
