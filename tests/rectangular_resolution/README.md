# rectangular_resolution

The rectangular energy resolution, the table kernel `rectangular` (`experiment/src/KernelRectangular.cpp`): a data set
whose detector has `energy_resolution rectangular width <w>` is folded to the average of its cross section over
[E0 - f w/2, E0 + f w/2] (c.m.; w lab, f = E_cm/E_lab of the pair).

`check.sh` runs `models/bare.azr` (no fold, a 0.05 keV lab grid across the narrow resonance of
`tests/convolution_window`) and `models/rect.azr` (seven points around the resonance, w = 20 keV) and compares each folded
point with the trapezoidal window average of the bare output (1e-3; measured 9e-5 at worst, from the adaptive grid of
50 nodes). `make_models.py` writes both models. The closed forms (sigma linear: sigma(E0); quadratic:
sigma(E0) + sigma'' w^2/24; a step) and dY/dw against central differences are in `experiment/tests/kernels_test.cpp`.

## The width as a fit parameter

`models/fit.azr` has two data sets, `d1` and `d2` (25 points each across the resonance; their quadratures differ, so
they are two effects), that name one detector `det_fit` with `energy_resolution rectangular width 14 keV fit`: one
fit parameter `det_fit.energy_resolution.width` (kind `apparatus`, `docs/upgrade/apparatus_parameters.md`).  Its
pseudo-data `models/data/d1.dat`, `d2.dat` are written by `make_models.py --data <bare output>`: each point the window
average of the bare cross section at w = 20 keV (Python, no engine kernel), times 1 + 1 % Gaussian noise (seed
20261111), error 1 %.  `check.sh` fits (one thread):
- (c) the committed data: one width parameter, w = 20.026 ± 0.013 keV (pull +2.0);
- (d) noise-free data made from the bare run of the check: w within 0.05 σ of 20 keV (measured 19.999998 keV), the
  χ² of both data sets below 0.1 (each is folded with the fitted width);
- (e) `width 21 keV fit prior 21 0.05`: the objective minus the data's χ² is ((w − 21 keV)/0.05 keV)², and
  `chiSquared.out` prints it as `Total-Kernel-Parameter-Chi-Squared`;
- (f) the refusals: `grid_max` below the value, a width 0, no mode, `profiled`,
  the window beside a target.
- (g) the triangular resolution (`experiment/src/KernelTriangular.cpp`, the drill of the recipe "A fit parameter of
  an apparatus effect"): `models/fit.azr` with `triangular width 14 keV fit` on noise-free data made in the check as
  triangle-weighted averages of the bare output at w = 20 keV: w = 19.9995 ± 0.0166 keV.
