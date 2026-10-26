# convolution_window — a Gaussian kernel must reach 5 sigma

`IntegrateTargetEffectImpl` integrates every Gaussian kernel (beam spread, straggling, the pure energy
convolution) over `[E - R sigma, E + R sigma]` and divides by the kernel's integral over that window, with
`R = TargetEffect::convolutionRange`. Until 2026-10-26 `R` was 3: the 0.27 % of the weight beyond 3 sigma was
dropped and the rest rescaled. Next to a resonance a thousand times its surroundings that tail carries real cross
section, so a point 2-3 sigma from the peak came out 2-4 % low and one more than 3 sigma away lost the peak
entirely (17O+a ODeSA targets, reported 2026-10-26; audit fix05). `R` is 5 now (lost weight 5.7e-7).

`make_models.py` builds, from `tests/ec_m1_exact/models/e1.azr` with external capture off, a 16O(p,g) model with
a 1 keV-wide 1/2- resonance at E_lab 0.8508 MeV on a ~1e-7 b background (peak 3.5e-3 b: 10^4 contrast):

| model | what |
|---|---|
| bare | the angle-integrated excitation function on a 0.05 keV lab grid, 0.80-0.90 MeV |
| conv | the points E_res + k sigma_lab, k = +-1.5, 2, 2.5, 3, 3.5, 4, with a 3 keV (lab) Gaussian beam spread, 50 sub-points |

`truth.py` convolves the bare output with full Gaussians (sigma_cm = sigma_lab F) over the whole grid and writes
`expected/truth.dat`; the bare output does not depend on the kernel, so the truth is independent of the window.
`check.sh`: (a) every convolved point within 2e-3 of the truth; (b) the bare resonance unchanged (peak, position,
background against `expected/bare_peak`).

**What it catches.** With `R = 3` the points at |k| = 3, 3.5 and 4 are 5-12 % off (the peak sits outside the
window) and those at |k| = 2.5 by 3-4 %. With `R = 5` the worst point is 1.2e-5 off.
