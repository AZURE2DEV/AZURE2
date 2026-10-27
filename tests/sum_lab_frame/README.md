# sum_lab_frame

A lab-frame differential SUM composite is the sum of its terms' **lab** cross sections.

A lab-frame differential segment (data `isDiff` 1, test `isDiff` 1) is compared in the c.m. frame of
its base pair: the data are converted once, σ_cm = J_base σ_lab (`ESegment::Fill`,
`EPoint::ConvertCrossSection`).  A SUM component to another exit pair is evaluated in its own c.m.
frame, at its own c.m. angle, with its own Jacobian J_term.  The lab sum Σ s_i σ_cm,i / J_i therefore
reads Σ s_i (J_base / J_i) σ_cm,i in the base frame.  Before 2026-10-26 `CalculateTheoreticalCrossSection`
added the c.m. terms as they were (audit item 7); `ESegment::TermFrameFactor` supplies the ratio now,
in the forward calculation and in both gradient paths (`AZUREGrad.cpp`), and the lab output
(`kinflag.dat` = 1) of the segment is the lab sum.  Each point stores its own conversion factor
(`GetCrossSectionKinFactor`, 1 when no conversion was done, refreshed with the energy shift), so a
c.m. SUM, an angle-integrated SUM, a component of the same pair or gamma terms of one entrance pair
get a factor of 1 and are unchanged.

## Models (`make_models.py`)

13C + α in 17O, two levels (5/2+ at E_x 11.0, 3/2+ at 11.1 MeV), exits n + 16O(6.05, 0+) (key 3) and
n + 16O(6.13, 3−) (key 5): the same masses, different Q values, so the two Jacobians differ at the
same lab angle (2.4 % at most on this grid).  No external capture; a run takes a second.

- `extrap`: test segments T1 (α → n0+), T2 (α → n3−), T3 = SUM(T1 + T2), T4 = SUM of the same pair
  twice, T5–T7 the c.m. versions, T8 = T1 + 0.5 T2 with the scaling token; lab 5.95–6.30 MeV, 0–180°.
- `plain_data` / `sum_data`: the two plain lab data segments and the SUM data segment, each with a
  fixed 0.05 MeV energy shift; `check.sh` writes `data/sum.dat` from the plain segments' fitted lab
  values.
- `sum_free`: `sum_data` with the seven widths free, for the gradient check.

## Checks (`check.sh path/to/AZURE2`)

- (a) lab output: T3 = T1 + T2, T8 = T1 + 0.5 T2, T4 = 2 T1 (1e-9, 104 points); the c.m. SUM T7 =
  T5 + T6; and the frames really differ (the lab SUM in c.m. output is not the plain c.m. sum).
- (b) the SUM data segment under the shift against the shifted plain segments' lab sum: χ² ~ 0.
- (c) `AZURE_GRAD_CHECK` on `sum_free` (data a few % off the model), mapped-components path and
  `AZURE_NO_COMPONENT_MAPPING=1`: adjoint and Jacobian-path gradients = finite differences.

dev 9e34f1c (before the fix): (a) lab SUM off by 2.5 %, the scaled SUM by 1.9 %, the lab SUM equal to
the plain c.m. sum; (b) χ² 47.0 on 36 points at 1 % errors; (c) passes (the gradient matched the wrong
model).  Fixed build: (a) 5e-11, (b) 4e-17, (c) 4e-9.
