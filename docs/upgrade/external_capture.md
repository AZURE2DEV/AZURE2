# External capture: design (split `CNuc`, step 4)

External capture is the part of a radiative-capture amplitude that comes from outside the channel radius a_c (fm):
the entrance scattering wave, the Whittaker tail of the final bound state f, and the electromagnetic multipole
operator between them. It enters the engine in two places: the capture amplitude T^ext to f at each data energy, and
the external part γ^ext_λγ (MeV^½, complex) of the photon amplitude of every R-matrix level λ, which `TransformIn`
subtracts from an observed width and the observed widths add back. This note fixes the part `physics/external_capture`
(ARCHITECTURE §3, `cnuc_split.md` §2.7 and §3 step 4): its inputs and outputs in physical units, the one formula for
the bound-state normalization N_f that the engine now evaluates in six places, and the definition of the
`--ignore-externals` mask, which today depends on the basis and jumps when a fitted amplitude crosses 10⁻⁸ MeV^½.

## 1. The quantities

Final state f: Brune energy Ẽ_f (MeV, E_x of the compound nucleus) and Brune amplitudes γ̃_fc (MeV^½) in its particle
channels c = (α, l_f, s) with radius a_c (fm) and reduced mass μ_c (u). Its bound-state normalization is

    N_f = n₀ + Σ_c (2 μ_c a_c / ħ²) γ̃_fc² I_fc,     I_fc = ∫_{a_c}^∞ [W_{l_f}(2κ_c r) / W_{l_f}(2κ_c a_c)]² dr      (1)

with n₀ = 1 (Brune) or Park's J_f (Park basis, where the amplitude is γ̃_f J_f^½), κ_c the bound-state wave number at
Ẽ_f, and I_fc (fm) the integral of `NFIntegral`, computed once at setup at the input energy. The asymptotic
normalization in channel c, divided by the unit-ANC Whittaker function, is

    s_fc = γ̃_fc N_f^{−½} K_c,     K_c = (2 μ_c a_c / ħ²)^½ / W_{l_f}(2κ_c a_c)                             (2)

(K_c is `ALevel::ECConversionFactor`, fixed at setup). The capture amplitude from entrance channel c to f by the
multipole XL (X = E, M) at the c.m. energy E is

    T^ext_{c→f}(E) = s_fc_f · u_{c c_f XL}(E)                                                                 (3)

where the unit amplitude u (complex) carries the entrance phases, (E_γ/ħc)^{L+½} and the radial integral of
`conventions.rst` Eq. (conv-ext); for channel capture it is multiplied by 2i P_c^½ times the U-matrix bilinear of the
entrance block (section 3). The external part of the photon amplitude of level λ (block b, photon channel γ = XL to f,
energy E_λ) is

    γ^ext_λγ = − Σ_{c ∈ b, particle} Σ_{c_f} s_fc_f γ_λc  I_{c c_f XL}(E_λ, Ẽ_f)                              (4)

over the pairs (c, c_f) of the same partition α allowed by the multipole: for E1, E2 and M1 the orbital rule
|l − L| ≤ l_f ≤ l + L, |s − l_f| ≤ J_f ≤ s + l_f, s = s_f; for M1 also the spin rule |s − 1| ≤ s_f ≤ s + 1,
|l − s_f| ≤ J_f ≤ l + s_f, l = l_f. Only E1, E2 and M1 with the final state's multipolarity flag contribute.
I_{c c_f XL} is `ECIntegral` with the channel-capture form (true): it depends on E_λ, so it is evaluated where E_λ is.

## 2. The part and its interface

`physics/external_capture` (`azure2/physics/ExternalCapture.h`), functions of plain arrays:

| function | inputs | output |
|---|---|---|
| `bound_state_normalization` | n₀; per channel a_c (fm), μ_c (u), γ̃_fc (MeV^½), I_fc (fm) | N_f, and N_f^{−½} (`inverse_sqrt_normalization`) |
| `external_reduced_width` | multipole (X, L), J of λ's block and J_f; per particle channel of λ: (α, l, 2s, γ_λc); per final channel: (α, l_f, 2s_f, s_fc_f); the radial integral as a function of (c, c_f) | γ^ext_λγ, Eq. (4) |
| `couples` | X, L, l, s, l_f, s_f, J_f | whether (c, c_f) enters Eq. (4) |
| `capture_mask` | the Brune-from-formal transform T (n × n) of a block, the masked levels | the projector Π of section 3 |

The part reads no engine state. The radial integral is injected (a `std::function`), because `ECIntegral` needs the
Whittaker functions, which leave GSL only in step 3b (i); the engine supplies it through the adapter
`src/CNucExternalCapture.cpp` (`ExternalReducedWidth`, `SqrtNF`), which reads a level, its group and the final state
and calls the functions above. `CNuc::CalcExternalWidth` and the six copies of N_f^{−½} (`FillCompoundFromParams`,
`FillCompoundFromParamsPhysical`, `TransformIn`, `TransformToBrune`, the Park-norm update of the gradient,
`ECFinalStatesFromFormal`) are replaced by calls of the adapter and deleted.

**Unit amplitudes.** The unit amplitudes u of Eq. (3) are already tabulated once per stage: `EPoint::CalculateECAmplitudes`
fills them per data point and pathway before the fit (with `ECTable` and the `intEC.dat` file as an opt-in grid in E),
from the exact-argument radial Coulomb scratch of step 3b (`channel_table.md` §4) and the per-stage channel table for
P_c^½ and e^{−iφ}; `AMatrixFunc::CalculateTMatrix` and `RMatrixFunc` multiply them by s_fc_f at each evaluation. That
structure is the one ARCHITECTURE §3 asks for and is kept; this step gives the scale s one function. The per-point
storage moves with the observables (step 5 and Phase 2), the radial integral with the Whittaker functions (step 3b (i)).

## 3. The `--ignore-externals` mask

**Today.** With `IGNORE_ZERO_WIDTHS` (`--ignore-externals`, the default of the CLI and pyazr) a photon channel γ of a
level λ whose amplitude has |γ_λγ| < 10⁻⁸ MeV^½ loses its external part: γ^ext_λγ is not subtracted (`TransformIn`,
`subtract_external_capture`) nor added (`ObservedWidths`), and λ is removed from both vectors of the channel-capture
bilinear Σ_{λμ} γ_λc A_λμ γ_μc′ (`AMatrixFunc::GetChannelCaptureBilinear` and two adjoint copies). The test reads the
amplitude in use: the formal one in the formal basis, which mixes the Brune levels, so the bases differ (9.7 × 10⁻² on
the `beam_profile_kernel` data); and it is applied at every evaluation, so a fit or profile that moves γ_λγ across
10⁻⁸ switches the level's external capture and the χ² jumps (AI-R branch `ai-r/fix-ignore-externals-step`, 5b83249a:
d(p,γ)³He χ² 458.0 against 125.2 across the threshold; `tests/13N` 1/2⁺ E1 by a factor 3.4).

**Definition.** The mask is a property of the input, not of the parameters: m_λγ = 1 when the observed photon width
of λ to f is given as exactly zero (a 0 in `<levels>`, a channel the level does not list, or a 0 in the physical vector
handed to `TransformIn`), else 0. The input is the observed (Brune or Park) form in every basis, so m is defined on the
Brune levels; it is recorded when the input is read (`ALevel::ZeroInputWidth`) and read-only afterwards. A masked
channel means "no γ decay of λ to f": γ^int = 0 at the input and no external part. It stays masked if a fit frees its
amplitude (the internal part then follows the fit, the external part stays off). Every other channel keeps its
external part at every amplitude, zero included, so every observable is continuous in every fit parameter. This is
the definition of 5b83249a, which is ported (not merged); it differs from today's only for an input that is not zero
but below 10⁻⁸ MeV^½ (kept now) and in fits that cross 10⁻⁸.

**Channel capture in each basis.** In the Brune and Park bases the levels of the A-matrix are the Brune levels, and
the mask removes masked levels from both vectors: v = P γ_c with P = diag(1 − m_λ). In the formal basis the levels are
formal. The two bases are related by the transform T of `brune_from_formal`: the Brune amplitude of level λ is
γ̃_λc = Σ_μ T_λμ γ_μc (T_λ the unit eigenvector of λ's converged problem), and A = Tᵀ Ã T, so that γ_cᵀ A γ_c′ =
γ̃_cᵀ Ã γ̃_c′ (Brune, PRC 66, 044611 (2002), Eqs. (20)–(30)). The masked Brune bilinear γ̃_cᵀ P Ã P γ̃_c′ is therefore,
in the formal basis,

    (Π γ_c)ᵀ A (Π γ_c′),     Π = T⁻¹ P T = 1 − Σ_{m masked} x_m T_m,  T x_m = e_m                            (5)

a projector (Π² = Π) that is diagonal only when T is. The engine computes Π at each evaluation for each block that
carries channel capture and has a masked level and more than one level (`brune_from_formal` with the transform, to
10⁻¹² MeV), like the final-state amplitudes of step 2, and the analytic gradient falls back to finite differences
there (the derivative of T is a Phase 2 item). With `--no-transform` the levels given are the only ones and P is used.

## 4. Open item 2 of step 2

A final state that is the lowest level of its block has formal parameters equal to its Brune ones only at the energy
where B_c was taken. If a formal-basis fit frees its energy, the two drift apart at second order. The block is then
added to the blocks whose final states are recomputed (`SetUpECFinalBlocks`): a lowest-level final state with a free
energy is treated like one above the lowest level.

## 5. Expected output changes and tests

Against the base (5055bec0), 179 `.azr` × three bases, menu 1: (a) Brune and Park: identical except a photon input
that is not zero but has |γ| < 10⁻⁸ MeV^½ (none expected); (b) formal: the default `--ignore-externals` channel capture
of projects with a masked level in a block of several levels moves to the Brune values (`beam_profile_kernel`); the
formal-basis observed widths of a photon channel whose transformed amplitude is < 10⁻⁸ but whose input is not zero
gain their external part. Anything else stops the step.

Tests: (1) `physics_external_capture_unit`: N_f against a hand value; the selection rules against a table of allowed
and forbidden (l, s, l_f, s_f, J_f) for E1, E2, M1; Eq. (4) with a stub integral against a hand sum; Π for a
diagonal T equals P, Π² = Π, Π T⁻¹ e_m = 0 for a masked m, and on a random two-channel, three-level block the masked
bilinear of Eq. (5) with A = Tᵀ Ã T equals the Brune one to 10⁻¹². A dropped transpose or a wrong T⁻¹ fails it at
order 1. (2) `tests/beam_profile_kernel/check.sh` runs the formal and Brune bases also with the default flags; the
captures must agree to 10⁻⁹ (they differ by 9.7 × 10⁻² before). (3) A continuity check on `tests/13N`: a masked
photon channel's file width set to 10⁻²⁰ eV keeps its external part (unit test of the mask's definition through the
CLI is left to pyazr's `ignore_externals` test of 5b83249a, ported if it fits the time).

## 7. The adjoint and the formal-basis masks (addendum, 2026-11-08)

**Why.** As first built, every formal-basis block with channel capture and several levels made the residual Jacobian
and the model gradients fall back to finite differences, masked level or not, while the χ² gradient
(`AccumulateEGammaGradient`, used by Minuit and pyazr `chi2_and_grad`) had no fallback at all and returned an analytic
gradient without the dependence of Π on the parameters: on `beam_profile_kernel` (formal basis) one component came out
0.103 against −0.394 from central differences, and 1.5 × 10⁻² relative on `tests/13N` with a masked 1/2⁺ level that has
a proton width. The step-2 final-state blocks (`HasFormalECFinalBlocks`) had the same gap on the χ² path. And
`pyazr_zero_width_level` lost its analytic Jacobian in the formal basis for a level whose widths are all zero, where Π
is exactly P.

**Rule.** (a) A masked level whose particle amplitudes are all zero at the current parameters is decoupled: no
element of the level shift matrix links it to another level, so T has e_m as its row and column and Π = P exactly.
When every masked level of a slot is decoupled, no projector is computed and P is used, as in the Brune basis.
(b) The adjoint covers a slot when it has no projector and no masked level has a free particle amplitude; a free
amplitude at zero gives dΠ = [P, dT] ≠ 0 although Π = P at the point. Otherwise, and whenever a formal-basis
final-state block exists, every adjoint entry point (χ² gradient, residual Jacobian, model gradients) falls back to
finite differences through one test (`FormalPartsOutsideAdjoint`, `src/AZUREGrad.cpp`), which asks
`CNuc::FormalCaptureMasksOutsideAdjoint` with the free-parameter map. A block with no masked level keeps the analytic
path (Π = 1).

**Expected.** Outputs of menu 1 unchanged against the branch head except where a decoupled masked level's projector
was used (P instead of a numerically equal Π: differences at the 10⁻¹² level at most). Analytic gradients: unchanged
in the Brune and Park bases; in the formal basis equal to central differences on the probes above (the fallback), and
`pyazr_zero_width_level` passes again (analytic, Π = P).
