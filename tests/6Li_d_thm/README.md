# 6Li_d_thm — a second Trojan-Horse regression, different carrier

The THM excitation function of 6Li(d,alpha)4He, on the same 8Be compound as
`7Li_p_ay` — the validated Paneru et al. model (PRC 111, 064609, Table IV
parameters) — entered through the other door: entrance pair 1 (6Li+d), exit
the identical-boson alpha+alpha pair.

Data: R. G. Pizzone et al., PRC 83 (2011) 045801 (EXFOR D0649002), the
THM-derived two-body cross section, 62 points, E_cm = 0.43-4.96 MeV. EXFOR
carries no uncertainties for the entry, so the files hold an assigned 10%.

## What this covers that tests/17O does not

* **A different carrier binding energy.** 17O's Trojan horse is the deuteron
  (B = 2.225 MeV); here the transferred particle is the deuteron itself,
  bound in 6Li = alpha (x) d with **B = 1.4735 MeV**, written on the
  entrance-pair channel lines through the optional 33-field column. A
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
