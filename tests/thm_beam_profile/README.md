# thm_beam_profile — the beam-profile kernel on THM points below threshold

`check.sh` folds `examples/o18_lacognata2008` (the 18O(p,a)15N THM spectrum of
La Cognata et al. 2008: 20, 90 and 144 keV levels, 31 points from
E_cm = -45 keV) twice:

- **G**: the Gaussian fold, `<targetInt>` convolution with sigma = 18 keV (lab);
- **B**: the `beamprofile` kernel with a flat beam (one skewed-Gaussian
  component, omega = 100 MeV, no skew), the same 18 keV resolution
  (`tpcSigma`) and a 2 eV energy window around every point (data columns 5-6).
  Its kernel, W(E) = [erf((b - E)/s sqrt 2) - erf((a - E)/s sqrt 2)]/2 times a
  constant, is then the Gaussian to 1e-9.

B must equal G at every point to 1e-4; it does to about 1e-6.

What it guards (fixed 2026-09-28): the beam-profile window of a point was
clipped at E_cm = +1 keV, where the Gaussian fold of a THM point goes down to
its floor 1 keV above the compound threshold (TargetEffect::minIntegrationEnergy
minus the entrance separation energy). A point at -45 keV was folded from the
sub-points above +1 keV only and came out 66 times too large (78 times in this
setup), -15 keV 5.4 times, still +0.5 % at +45 keV. The window also cut the
resolution tails at 4 s, where the THM Gaussian fold keeps 5 sigma: 81 keV
above the 289 eV level at 144 keV the fold was 2 % low. Above ~50 keV the two
kernels had agreed to 1e-4.

The model is the example's, read in place; nothing here depends on its
parameters, since both folds use the same ones. About 10 s.
