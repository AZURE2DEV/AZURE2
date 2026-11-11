# rectangular_resolution

The rectangular energy resolution, the table kernel `rectangular` (`experiment/src/KernelRectangular.cpp`): a data set
whose detector has `energy_resolution rectangular width <w>` is folded to the average of its cross section over
[E0 - f w/2, E0 + f w/2] (c.m.; w lab, f = E_cm/E_lab of the pair).

`check.sh` runs `models/bare.azr` (no fold, a 0.05 keV lab grid across the narrow resonance of
`tests/convolution_window`) and `models/rect.azr` (seven points around the resonance, w = 20 keV) and compares each folded
point with the trapezoidal window average of the bare output (1e-3; measured 9e-5 at worst, from the adaptive grid of
50 nodes). `make_models.py` writes both models. The closed forms (sigma linear: sigma(E0); quadratic:
sigma(E0) + sigma'' w^2/24; a step) and dY/dw against central differences are in `experiment/tests/kernels_test.cpp`.
