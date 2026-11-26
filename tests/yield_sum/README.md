# yield_sum — the yield of a sum of reactions

`check.sh` (run by `run_tests.sh`) checks that quantity `yield` takes a sum of reactions: a thick-target yield counted
over several exit groups (neutrons to every state of the residual nucleus) is one data set, `r1 + r2 + …`. Until
2026-11-25 format 3 refused it ("takes one reaction, not a sum"), and an old project with such a line (a SUM
segment on a `<targetInt>` line of density 0: ¹⁹F(α,n)²²Na on BaF₂ and CaF₂, nine groups) did not convert.

Model: `../channel_space` (¹²C + p with the inelastic exit ¹²C[0.8] and capture) with one data set on a carbon target.

1. **Angle-integrated, infinite target.** The yield of the sum is the sum of the two terms' yields (limit 1e-9; 5e-11),
   and the order of the terms does not matter.
2. **At 60° lab, on a layer.** The numbers of a data set at a lab angle are written in the c.m. frame of its first
   term, so the sum is a + ρ b with ρ the ratio of the two terms' solid-angle factors; ρ is taken from the reversed sum
   (b + a/ρ), and the sum then equals a + ρ b (limit 1e-9; 3e-11, ρ = 0.905 to 0.923).
3. **Both terms contribute**, and a ratio as a yield is still refused.

The lab-angle case is on a layer because a thick-target yield at a lab angle is not computed for a reaction with a
threshold, with or without a sum: the fold reaches energies where the lab angle is beyond the term's kinematic
maximum and the point fails with a warning (`src/PointFold.cpp`). That is the next thing the ¹⁹F(α,n) data sets at
an angle need.
