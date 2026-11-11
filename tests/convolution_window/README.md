# convolution_window — a Gaussian kernel must reach 5 sigma

`IntegrateTargetEffectImpl` integrates every Gaussian kernel (beam spread, straggling, the pure energy
convolution) over `[E - R sigma, E + R sigma]` and divides by the kernel's integral over that window, with
`R = TargetEffect::convolutionRange`. Until 2026-10-26 `R` was 3: the 0.27 % of the weight beyond 3 sigma was
dropped and the rest rescaled. Next to a resonance a thousand times its surroundings that tail carries real cross
section, so a point 2-3 sigma from the peak came out 2-4 % low and one more than 3 sigma away lost the peak
entirely (17O+a ODeSA targets, reported 2026-10-26; audit fix05). `R` is 5 now (lost weight 5.7e-7).

`make_models.py` (deleted; see "How the models were made") builds, from `tests/ec_m1_exact/models/e1.azr` with external capture off, a 16O(p,g) model with
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

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 7d24c881:tests/convolution_window/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for the Gaussian kernel-window test (tests/convolution_window).

Base: tests/ec_m1_exact/models/e1.azr (p + 16O -> 17F(0.495) + gamma; here with external capture OFF, mask 0, so the
only cross section is the R-matrix one), plus one NARROW 1/2- resonance at E_x = EX_NARROW (E_cm EX_NARROW - 0.5996,
E_lab x 17.0078/16): proton width GP eV, gamma width GG eV.  Far from it the capture cross section is the 30 MeV
levels' tails, ~10^-9 b, so the resonance is >> 10^3 x its surroundings -- the situation in which a Gaussian kernel
cut at 3 sigma and renormalised (the pre-2026-10-26 TargetEffect::convolutionRange = 3) loses real cross section.

  bare.azr : the angle-integrated 16O(p,g) excitation function on a fine lab-energy grid (STEP) over EMIN..EMAX;
             truth.py convolves it with full Gaussians;
  conv.azr : the same at the test points E_res + k sigma_lab, k = KS, with a Gaussian beam spread SIGMA_LAB
             (constant-sigma field) and NSUB sub-points.

Run from this directory: python3 make_models.py   (writes models/)
```
