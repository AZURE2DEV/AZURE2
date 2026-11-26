# lab_frame_fold -- a point given at a lab angle, folded in the lab

Decision D5 (`docs/upgrade/energy_grid.md` section 4.3): a folded point whose data are given at a lab angle is folded
with the c.m. angle theta_cm(E_n; theta_lab) and the Jacobian J(E_n; theta_lab) = dOmega_lab/dOmega_cm of every node
E_n of its fold, m_lab = sum_n w_n sigma_cm(E_n, theta_cm,n) / J_n (`include/RowSet.h`); before, every node took the
point's own c.m. angle and the comparison its own J.

`models/lab_fold.azr` is the level scheme of `tests/relativistic/models/li7pn_ang.azr` (7Li(p,n)7Be through one far
2- level, a smooth cross section; non-relativistic kinematics) with three data sets at E_lab 2.15-2.45 MeV, 30 deg
lab, a few hundred keV above the threshold (1.881 MeV), where theta_cm runs from 41.7 to 38.6 deg and J from 0.553 to
0.636 across the window:

| data set | what |
|---|---|
| `folded` | E_lab = 2.300 MeV, `detector det_wide` (rectangular resolution, 300 keV lab), 200 quadrature points |
| `fine` | the same lab angle unfolded at 301 energies 1 keV apart (value 1: the c.m. data column is J) |
| `frozen` | the same energies at the folded point's c.m. angle 39.684407814 deg as c.m. data (the fold before D5) |

`check.sh` compares the folded point with Simpson's rule on `fine` (the engine's unfolded points are the reference;
no second formula of the cross section): 2.2e-5 apart (the fold's quadrature: 3.6e-4 at 50 points, 1.4e-6 at 800),
while the frozen angle is 2.2e-2 off.  The build before D5 (milestone 3 of the observables branch) gives 2.2e-2 and
fails.  The data values are 1 with error 0.1 (only the model is compared).

Case (d): `folded` with a fitted energy offset of prior width 0.3 MeV (a variant made in the check).  The fold's planned
nodes then reach below E_c.m. = 1.678 MeV, where 30 deg lab is beyond the kinematic maximum: no neutron reaches the
angle there, so such a node holds a zero row (since 2026-11-25, `experiment::LabBranchesAt`; NaN before, and a fold
that read it was a failed point).  At offset 0 the point is finite (the build of milestone 4 aborted at set-up); at
offset -0.32 MeV the fold reads those energies, which add nothing, and the point is finite and smaller.
