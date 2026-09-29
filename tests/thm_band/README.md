# thm_band -- the CLI uncertainty band of THM (HOES) points

`check.sh` frees the six R-matrix parameters of `tests/18O_p_a_thm` (two THM
segments with profiled norms), puts a hand-made correlated 6x6
`output/covariance.dat` in place and runs a plain calculation with
`--covariance-band`. The band of every point must equal sqrt(g Sigma g^T),
g taken here by central differences of the CLI's own output files for
q = m n*(p0)/n*(p) (model against the data scaled by the profiled norm; see
docs/source/theory/thm_implementation.rst, "Normalization, gradients and
uncertainty bands"). Tolerance 2e-3 relative; measured 8e-7.

The binary before this test existed used the T-matrix adjoint row for HOES
points; its band was 600-1600 times the finite-difference value. A band of the
model m alone (no profiled-scale term) is off by up to ~50 %, so the check
also pins that term. 13 runs of ~0.15 s.

A second pass puts both segments in one THM experiment with a linear
background (`experiment[A] segments=1,2 background=linear`): the written curve
is the model plus b(E), norm and background are profiled together, and the
band is that of q = (m + b) n*(p0)/n*(p); measured 7e-7. 26 runs in all.
