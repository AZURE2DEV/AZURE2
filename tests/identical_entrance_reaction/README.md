# identical_entrance_reaction — (1+δ₁₂) for an identical entrance pair

An identical pair in the entrance channel multiplies every cross section out
of it by 1+δ₁₂ = 2:

    σ(12→34) = (π/k²) Σ_J (2J+1)/[(2i₁+1)(2i₂+1)] (1+δ₁₂) Σ |T^J|²

summed over the channels exchange symmetry allows (ℓ+s even). For spin-0
bosons this is the familiar (π/k²) Σ_ℓ (2ℓ+1)[1+(−1)^ℓ] T_ℓ of 12C+12C
fusion. It is also the convention of the reciprocity theorem
w₁₂k₁₂²σ(12→34)/(1+δ₁₂) = w₃₄k₃₄²σ(34→12)/(1+δ₃₄), with
w = (2i₁+1)(2i₂+1). An identical **exit** pair takes no factor: the cross
section counts reactions, not outgoing particles.

Before this test existed, AZURE2 applied the 2 to elastic scattering only
(commit 7c34992). 12C(12C,α)20Ne, d(d,p)t and α+α → p+7Li came out a
factor 2 low, both angle-integrated and differential. The analytic gradient
also lacked the factor 2 that 7c34992 had put on the elastic angle-integrated
cross section.

`check.sh` runs three one-level projects (`projects/`) at the level energy
E_R. With the boundary condition at S_c(E_R), which is both the default and
`--use-brune`, the cross section there is exactly

    σ(E_R) = (π/k²) g_J (1+δ₁₂) 4 Γ_a Γ_b / Γ²

It compares AZURE2 against that expression, computed in awk with AZURE2's
own constants:

| project | reaction | expected |
|---|---|---|
| c12c12 | 12C+12C → α₀+20Ne, 2⁺ | ×2; dσ/dΩ = σ·5P₂(cos θ)²/4π (so ×2 as well) |
| c12c12 | 12C+12C elastic | ×2, with the hard-sphere phase φ pinned from mpmath |
| dd_pt | d+d → p+t, 2⁺, g = 5/9 | ×2; dσ/dΩ isotropic = σ/4π |
| be8 | p+7Li → α+α | no factor (identical exit) |
| be8 | α+α → p+7Li | ×2, and reciprocity with the forward value |

Tolerance is 10⁻⁶. The two exceptions are the elastic value at 10⁻⁴, because
AZURE2's Coulomb functions set φ ≈ 4.5×10⁻¹² to only about 10⁻⁶, and
reciprocity at 10⁻⁵, because the two runs sit 2.6×10⁻⁷ MeV apart. The
pre-fix binary fails 16 of the 19 checks. The three it passes are the two
elastic values and the identical-exit value.

Also verified when the change went in (not scripted here):

- 12C+12C → α+20Ne from 1.8 to 2.2 MeV, against the full one-level formula
  with P_c(E), S_c(E) and the level shift (mpmath Coulomb functions):
  AZURE2/analytic = 1.0000000–1.0000003, default and Brune.
- 34 single-level projects, one per level of the Tumino et al. (Nature 557,
  687 (2018)) analysis, against the independent Python OES
  (π/k²)(2J+1)(1+δ)ΓₓΓ_c/D: the per-level on-peak ratio went from 0.5002 to
  1.0004. The THM (HOES) curves are unchanged, since their scale is
  arbitrary.
- d+d reaction rate (menu 5): exactly ×2. The rate still excludes the
  1/(1+δ₁₂) of the rate equation, as REACLIB and NACRE do, and AZURE2 now
  says so when it computes such a rate.
- `AZURE_GRAD_CHECK`: the analytic gradient agrees with finite differences
  as closely for these pairs as it does for a mass-perturbed,
  distinguishable copy of the same project.
