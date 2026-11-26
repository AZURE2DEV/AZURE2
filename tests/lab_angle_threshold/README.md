# lab_angle_threshold — a yield at a lab angle across the threshold of its reaction

`check.sh` (run by `run_tests.sh`) checks the fold of lab-frame data over energies where the reaction has no particle
at the lab angle, or two c.m. angles for it (`experiment::LabBranchesAt`, `Unit::RowSetOver`; since 2026-11-25).
A thick-target yield of an endothermic reaction at a lab angle always reads such energies: below the threshold there
is no ejectile, and just above it the ejectile is slower in the c.m. than the c.m. itself, so it goes into a forward
cone, each lab angle inside the cone reached by a forward and a backward c.m. angle. Before, the point failed.

Model: `../channel_space` (¹²C + p; the inelastic exit ¹²C[0.8]: threshold 0.8672 MeV lab, two-valued up to 0.8733 MeV;
capture), an infinite carbon target.

1. **The yields at 16 lab angles integrate to the angle-integrated yield** (Gauss–Legendre in cos ψ over the sphere),
   at beam energies 2.0, 0.95 and 0.88 MeV; the last is 7 keV above the two-valued range, so the fold crosses all of
   it. Limit 1e-5; measured 2 to 3e-11. The angle-integrated yield has no angle in it, so the check does not use the
   kinematics it tests; it fails if the energies without a particle are not zero or if only one of the two c.m.
   angles is taken. No point fails.
2. **A sum with a term that has no particle at the point.** At 0.8 MeV, below the inelastic threshold, the yield of
   capture + inelastic at 60° is the capture's (to rounding), and the sum with the closed term first is computed too
   (its row stays in the lab frame, `FrameRow::noParticle`).
3. **A point that is not folded is refused as before**, with the conversion's message ("below the threshold of the
   exit pair …").

Limits. With the beam inside the two-valued range the lab cross section has an integrable peak at the edge of the cone
at the top of the fold, and the default quadrature is coarse there: at 0.8725 MeV the integral over the cone is 11 %
below the angle-integrated yield with 400 points, 1.5 % with 4000 and 0.6 % with 40000 (measured 2026-11-25; the yield
there is 1e-43 per particle). Relativistic kinematics, elastic scattering and analysing powers keep the old rule (one
angle; a fold that reads an energy beyond the maximum fails).
