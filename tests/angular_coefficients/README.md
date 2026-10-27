# angular_coefficients

Legendre coefficients of an angular distribution as data (segment data type `isDiff` 9, 2026-10-26).

Some experiments publish only the coefficients of a Legendre fit to their angular distributions,
dsigma/dOmega_cm = A_0 [1 + sum_k a_k P_k(cos theta_cm)] (11B(a,p)14C: Wang et al. PRC 43 (1991) 883,
Dayras et al. NPA 261 (1976) 365).  A data segment with `isDiff` 9 fits them: the data file rows are
`E_lab  k  a_k  error`, the order k in the angle column, and the fit value of each point is the model's
coefficient of that order.  The coefficients are the ones the `isDiff` 3 test segments have always
computed (`GenMatrixFunc::CalculateCrossSection`: each interference term of order k divided by the
angle-integrated cross section, so a_0 = 1); the trailing keyword `absolute 1` makes the data the absolute
A_k in b/sr instead (A_0 = sigma_int / 4 pi).

- The points are angular-distribution points; `EPoint::IsCoefficient` / `GetCoefficientOrder` select the
  order, `GetAngDistA0` carries A_0.  Energy shifts move them like any point.
- Target integration and resolution folding act on the absolute A_k, which are linear in the cross
  section, and the ratio is formed afterwards (`EPoint::IntegrateTargetEffectForObservable`, the A_y
  pattern); absolute A_k take the plain path.
- A normalised coefficient has no scale, so a free normalization is refused ("absolute 1" for A_k);
  an elastic pair (the expansion leaves the Coulomb amplitude out), total capture, composites, samples,
  "total" and UPOS are refused, and `isDiff` 9 in `<segmentsTest>` is refused (use 3).
- Analytic gradient (`AMatrixFunc::PointTBar`): a_k = N_k / sigma_int, the differential machinery with a
  Kronecker delta in place of P_l(cos theta) for N_k and the angle-integrated cotangent scaled by
  -a_k / sigma_int; under a target effect the A_y structure with A_0 in place of the cross section
  (`AZUREGrad.cpp`, `GradTargetEffectAdjoint`).

## Models (`make_models.py`, the two-level 13C + alpha model of tests/sum_lab_frame, exit n + 16O(3-))

`extrap` (test segments: the coefficients, the c.m. differential cross section, the angle-integrated one),
`data` (a normalised and an absolute coefficient segment; data written by `verify.py` from `extrap`),
`shift` (a fixed 0.05 MeV shift), `conv` / `conv_diff` (a 30 keV beam spread on the coefficient segments
and on the differential test segment), `fit` (free widths, data a few per cent off), `bad_*` (refusals).

## Checks (`check.sh path/to/AZURE2`)

- (a) the `isDiff` 3 coefficients equal the least-squares Legendre projection (P_0..P_8) of the c.m.
  differential cross section on 0-180 deg, a_k = c_k / c_0, to 1e-6 (the convention and sign);
  (b) c_0 = sigma_int / 4 pi;
- (c) both data segments reproduce their data (chi2 5e-14);
- (d) a 0.05 MeV shift moves the coefficients to E + 0.05 exactly;
- (e) under the beam spread the folded coefficients equal the projection of the folded distribution
  (normalised 2e-9, absolute 3e-9);
- (f) `AZURE_GRAD_CHECK` with `--no-grid-replacement` (the finite differences otherwise include the
  motion of the convolved segment's adaptive grid with the widths, ~1e-5): adjoint and Jacobian-path
  gradients = finite differences to 1.5e-8 over the 7 widths, with normalised, absolute and folded
  segments;
- (g) the refusals.
