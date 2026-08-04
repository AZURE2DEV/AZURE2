# 12C(p,p)12C analyzing power — Baumann et al. (1992)

Exercises the vector analyzing power (observable code 7) against measured data.

## Where the data come from

R. Baumann, G. Keil, N. Kniest, E. Pfaff, M. Preiss, M. Skill and G. Clausnitzer,
*The analyzing power for elastic 12C(p,p)12C scattering below 2.1 MeV*,
Nucl. Phys. **A542** (1992) 53.

The paper publishes no table of A_y. Its results are contour plots (figs. 3 and
4) plus six angular distributions in fig. 1, at E_p = 1.618, 1.658, 1.708,
1.738, 1.758 and 1.779 MeV. `data/baumann_ay.dat` is fig. 1 **digitised**: the
page was rendered at 600 dpi, the axes calibrated on the printed tick labels,
and the drawn curve traced column by column. The curve is the paper's
phase-shift fit, which the plotted points lie on, so these are best understood
as the published analysis sampled at 10-degree intervals rather than as the raw
measurements.

The quoted uncertainty of 0.04 is therefore **digitisation error, not the
experiment's**. Baumann quote statistical errors below 0.004 -- ten times
smaller. 0.04 is roughly the drawn line width (12 px against 466 px per unit of
A_y) plus the axis calibration residual. Angles are centre-of-mass, as printed
on the figure axis, which is what observable 7 expects; energies are laboratory.

Points where the automatic trace could not follow the curve through a steep
section were dropped rather than interpolated, which is why some energies carry
fewer than the full 13 angles. 71 points survive.

## What the model is

`12C_pp_ay.azr` is the `tests/13N` evaluation with its own data segments
switched off and the A_y set added as segment 11. Every level parameter is held
fixed except the two resonances this energy range is sensitive to -- the 3/2-
at E_x = 3.503 MeV and the 5/2+ at 3.545 MeV -- for which the excitation energy
and the proton width are free. Those are the same two states Baumann extract
(their table 1: p3/2 at 3.499 MeV with 57 keV, d5/2 at 3.546 MeV with 50 keV).

The file is stored at the fitted minimum, so the regression runs a plain
calculation.

## Fit result

Fitting those four parameters to the 71 points takes chi-squared from 236.8 to
85.2, i.e. 1.20 per point and 1.27 per degree of freedom:

| parameter | fitted | Baumann table 1 | start |
|---|---|---|---|
| 3/2- E_x     | 3.4954 MeV | 3.499 MeV | 3.5032 MeV |
| 3/2- Gamma_p | 50.5 keV   | 57 keV    | 55.2 keV   |
| 5/2+ E_x     | 3.5430 MeV | 3.546 MeV | 3.5453 MeV |
| 5/2+ Gamma_p | 53.6 keV   | 50 keV    | 49.0 keV   |

Excitation energies agree to 3-4 keV and widths to about 10%, which is what
digitised data can support: the 0.04 uncertainty is an order of magnitude
looser than the measurement, so the fit cannot be expected to recover the
published parameters more tightly than this.

To rerun the fit:

    printf '2\n\nn\n\n' | AZURE2 --no-gui --no-readline 12C_pp_ay.azr

## Reproducing the digitisation

`digitize_fig1.py` regenerates `data/baumann_ay.dat` from the PDF. It needs
`pdftoppm`, numpy, scipy and Pillow, and the paper, which is not redistributed
here.
