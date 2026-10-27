# upos_pair_index

The UPOS secondary-gamma distribution takes the residual nucleus's spin from the exit pair.

`GenMatrixFunc::CalculateBareCrossSection` and its gradient mirror in `AMatrixFunc::PointTBar` formed the
secondary-gamma coefficient `R_k(j2f, L, L', Ic)` with `j2f = GetPair(ir)->GetJ(2)`, where `ir` is the DECAY INDEX
within the entrance pair's list (the function's own comment says so: an earlier multi-compound-nucleus bug had the
same cause).  Decays are appended in pair order and only for pairs that share J-groups with the entrance pair, so
the index equals the pair number only while every earlier pair has a decay from the entrance pair.  A pair of
another compound nucleus listed first, a beta-delayed pair, or a pair whose J-groups are out of the R-matrix breaks
that, and `j2f` becomes the spin of the wrong nucleus (audit item 10, fixed 2026-10-26: `GetPair(exitPairNum)`, as
`CNuc::CalcAngularDists` already did).

`check.sh` prepends a p + 12C pair (key 1, one 1/2- level) to `tests/secondary_gamma/models/spin2_plus.azr`
(alpha + 13C -> n + 16O(3-) with the E3 secondary gamma) and compares the three output blocks of the UPOS/particle
file with the unmodified model: the UPOS distribution (a), the angle-integrated and particle distributions (b) must be
identical, and the extra nucleus must be in the model (c).  dev 17d500b (before the fix): (a) the UPOS values of the
two-nucleus model are NaN (`j2f` = 1/2, the alpha pair's 13C, with L = 3 and I_c = 0 has no R_k); fixed build: identical.
