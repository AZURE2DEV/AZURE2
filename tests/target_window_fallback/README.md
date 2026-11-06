# target_window_fallback — the effective energy outside the windows of a target integration

Since aeb54b7 (2026-08-31) a `<targetInt>` line can be restricted to lab-energy windows (`"lo-hi,..." blend [tol]`).
A point outside every window was left a plain point at its data energy. The data of a target integration carry the
BEAM energy E0, so that point's model was sigma(E0) instead of the target average, off by the centroid shift (about
half the energy loss). On real data (13C+a ND2021r, 16 lines with windows around the narrow levels,
13C+a/10-30-26_an0_structure/ti/g5n.azr) the windowed model came out worse than no target integration at all:
chi2 421,432 against 400,616 integrated everywhere. With the fallback it is 401,435.

Now (2026-11-05) such a point is sigma at the depth-averaged energy of its own target,
`E_bar = int E dE/eps / int dE/eps` over `[E0 - Delta, E0]` (`TargetEffect::MeanEnergyLoss`): one sub-point at E_bar,
or three (E_bar at 0.01 N0, N0 and the layout density) under `densityfit`. Window edges blend with sigma(E_bar) inside
the linear kernel, so the analytic gradient covers them. The automatic tolerance falls back to sigma(E_bar) and judges
the curvature only. Windows on a thick target, and windows or a tolerance on a line with `doppler`, are refused.

Model: `tests/thick_target_resonance` (p + 14N, external capture off), a thin target with epsilon = a0/x^0.8 and
N = 1e18 atoms/cm^2 (Delta = 18 keV at 360 keV), 57 points over 240-520 keV, the 278 keV (lab) resonance inside the
window 0.27-0.32 MeV, and the bare cross section on a 0.1 keV grid as the truth table. For the power law E_bar has a
closed form. The checks, (a)-(h), are listed at the top of `check.sh`. The ones worth knowing:

- (a) outside the window, every point equals sigma(E_bar) of the truth table (9e-6; the table's interpolation).
- (b) |full integration - sigma(E_bar)| matches |sigma''(E_bar)| Delta^2/24 to within 1 % at every point away from
  the resonance. It is up to 12.6 % at 240 keV, where the Coulomb barrier bends sigma by e every 16 keV across a
  25 keV target. That is the price of the fallback, and the reason windows belong wherever it is not small.

Builds before 2026-11-05 fail 7 checks: (a), (a'), (d), (e), (f) and both refusals of (g).
