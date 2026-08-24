# 13N_capture_ay — the capture analyzing power

The same `12C + p` compound nucleus as `13N`, reduced to two segments on one
`(E, theta)` grid — a centre-of-mass differential capture cross section and the
vector analyzing power of the same reaction — so that the two are computed from
the same T-matrix elements and can be compared point by point.

The capture analyzing power follows R. G. Seyler and H. R. Weller,
*Angular distribution theory for particle-capture-gamma reactions*,
Phys. Rev. C **20** (1979) 453, Eqs. (20) and (21) in the channel-spin
representation:

    A_y(theta) = sum_k b_k P_k^1(cos theta) / sum_k a_k P_k(cos theta)

The reference chi-squared is this implementation's own output, so on its own it
locks the numbers against drift rather than proving them right. Three
independent checks do that:

1. **Against the paper.** `seyler_weller_check.py` reproduces every coefficient
   of the worked example on p. 458 — `a_0` through `a_4` and `b_1` through
   `b_4` — to the two decimals printed. That fixes the one convention the paper
   leaves implicit: its `X(l s b; l' s' b'; k 1 k)` is the bare 9-j symbol, not
   the `[prod (2j+1)]^{1/2}`-normalized Fano X.
2. **Against AZURE2 itself.** The denominator `sum_k a_k P_k` must reproduce the
   differential capture cross section AZURE2 builds by the entirely separate
   Blatt-Biedenharn route in `CNuc::CalcAngularDists`, up to the
   angle-independent factor `400 pi / (geom * I1I2)`. Set
   `AZURE2_CAPPOL_DEBUG=1` and the ratio is printed per point; it agrees to
   machine precision, which is what pins the pathway enumeration, the coupling
   order and every hat factor.
3. **Against central differences.** The analytic adjoint in
   `AMatrixFunc::PointAdjoint` agrees with numerical derivatives of `A_y` to
   ~1e-6 of the column scale over all points and all R-matrix parameters,
   external capture included.

Two structural properties are visible in `output/AZUREOut_aa=1_R=2.out`:
`A_y` is exactly zero at 0 and 180 degrees (every `P_k^1` vanishes there), and it
stays inside [-1, 1].

Target effects are switched off in this project on purpose. `A_y` averaged over
a target is cross-section weighted, `<A_y> = int A_y sigma dE / int sigma dE`,
which is a different observable from the thin-target one.

## Segments

| # | data | observable |
|---|---|---|
| 1 | `capture_xs.dat` | differential capture cross section, c.m. |
| 2 | `capture_ay.dat` | **capture analyzing power** |

The data columns are placeholders (`A_y = 0`, uncertainty 1), so segment 2's
chi-squared is just `sum A_y^2` over the grid — which moves if anything in the
calculation moves.
