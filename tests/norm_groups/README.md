# norm_groups — one normalization shared by several segments

A segment line may end with `normgroup <label>`.  Every segment carrying the
same label uses ONE normalization: that of the group's first segment in file
order (its value, vary flag, uncertainty and penalty).  The other members' own
norm, uncertainty and vary flag are ignored.  A member whose own line varies
its norm keeps its `segment_<key>_norm` parameter, fixed and written with the
shared value, so the parameter numbering does not depend on the groups; it
carries no norm penalty.  Introduced 2026-10-14 for data that share one
absolute scale (11B+a: the Dayras et al. 1976 (a,p) excitation function and
the angular distributions it was integrated from).

Model: the p+p model of tests/shift_groups.  The data files are that model's
exact cross sections (the `*_off6.dat` yields at their true energies) times
0.8 (`xs_40_s08`, `xs_90_s08`) or 0.9 (`xs_90_s09`), with 1% errors, so a
norm of 1.25 (1.111) is the chi2 = 0 point.

`norm_groups.azr`: leader `xs_40_s08` at norm 1.2 (fixed), follower
`xs_90_s09` with its own norm 1, both in group D: the follower is evaluated at
1.2.  `expected/chiSquared.out` pins both.

`check.sh`:

1. the grouped pair gives exactly the chi2 of the same two segments each given
   norm 1.2 directly (1e-10); the follower's own "free, 10%" adds no penalty;
   the log notes the group;
   a fixed leader with a follower whose own line says "free" keeps the group
   fixed: no norm parameter is fitted, the follower stays at the leader's value;
2. both files scaled by 0.8: with every level fixed and the leader's norm
   free, MIGRAD, `--use-gradient` and `--use-lm` return n = 1.25, chi2 = 0,
   and the follower written with the same value in `param.sav`;
3. 0.8 and 0.9: one norm strictly between 1.111 and 1.25;
4. both lines say "1 +- 10%": Total-Norm-Chi-Squared is ((n-1)/0.1)^2 once,
   the follower's row shows 0;
5. `AZURE_GRAD_CHECK=exit`: the leader's norm gradient (adjoint and Jacobian
   paths) agrees with a finite difference of the whole chi2, and only two
   parameters are free (the follower's norm is fixed).

tests/natural_target adds a group over a natural SUM and a plain segment
(gradient check), tests/pyazr/norm_group_test.py the pyazr side.
