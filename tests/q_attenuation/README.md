# q_attenuation — attenuation coefficients on a particle exit

A `<targetInt>` line can carry finite-geometry attenuation coefficients Q_K (Rose, Phys. Rev. 91, 610 (1953)). They
belong to the detector, so they act on the **lab** angular distribution. For a particle exit,
`EPoint::CalcLegendreP` now does the following:

1. It expands W_lab = W* dΩ*/dΩ_lab in lab Legendre polynomials. The coefficients are
   b_K = Σ_L M_KL c_L, with M_KL = (2K+1)/2 ∫ P_K(x_lab(x*)) P_L(x*) dx* (96-point Gauss–Legendre).
2. It multiplies b_K by Q_K.
3. It takes the result back to the point's c.m. angle:
   factor_L = P_L(x*) + J Σ_{K<n_Q} (Q_K − 1) M_KL P_K(x_lab), where J = dΩ_lab/dΩ*.

The lab angle comes from g = v_cm/v₃* (with `relativistic=1`, the relativistic boost). Photon exits, UPOS, A_y,
charged-particle elastic scattering and g ≥ 1 multiply P_L by Q_L as entered.

`make_models.py` builds the models from `tests/secondary_gamma/models/grad.azr`. That model is
¹³C(α,n)¹⁶O*(6.13) close to threshold, so g is 0.33–0.37. The output is a c.m. differential test segment at
E_lab 5.9–6.3 MeV and θ* = 0–180° in 5° steps. `check.sh` checks three things:

- (a) 40 coefficients, all equal to 1 (`q1.azr`), give output byte-identical to no coefficients (`plain.azr`).
- (b) The 40 coefficients of a uniform 15° cone (`q.azr`) agree with a brute-force average of the lab distribution
  over that cone, 2-d quadrature in `verify.py`: 8e-10. The reference uses no Legendre algebra for the attenuation.
  Applying Q in the c.m. frame instead would be off by 3.8 %.
- (c) Only Q₀–Q₄ of the cone (`q5.azr`) agree with the same formula, with b_K projected in the **lab** variable:
  1e-9.

**What it caught (fixed 2026-10-25).** The particle-exit branch multiplied P_L by B = U(γ)·Q. U was an O(γ²)
lab→c.m. coefficient matrix, and γ = m₁m₃/(m₂m₄)·E/(E+Q) is about g², not g. B is not the attenuation of the lab
distribution. The dev build before the fix fails all three checks:

- (a) the output changes even with every Q_K = 1;
- (b) it is off by 5.1 %, and its 40 coefficients read past the end of a 20-entry array;
- (c) it is off by 8.5 %.

A non-particle entrance pair also set every factor to 0. Archive fits affected: see
`R-matrix/Brune_vs_Park_Claude_eval/10-25-26_qcoeff_scan/readme` (old ¹³C(α,n) and ¹⁵N(p,α) fits; the 2026
campaigns use photon/UPOS exits only).
