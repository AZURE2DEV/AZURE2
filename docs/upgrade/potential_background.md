# The interior of a potential as the background of an R-matrix channel

**Status (2026-11-20): research, side branch `feature/upgrade-9-potential-background`; not for merge without the
PI's word.** Sections 1-6 are the design; section 7 says what is built, section 8 gives the results, section 9 what
is not done.

A channel of a pair with a `potential` record takes its exterior functions (P, S, φ) from the radial equation with a
nuclear potential acting between the channel radius a_c and the matching radius r_m (`hybrid_potential.rst`). This
note adds the interior half: the R-function of the same potential, obtained by integrating the radial equation from
the origin to a_c, is used as the background of that channel's diagonal R-matrix element in place of free background
poles; the fitted levels sit on top of it. With both halves from one potential and no levels, the channel's collision
matrix is that of the potential alone, whatever a_c is (section 2.3); with levels, a_c moves strength between the
levels and the potential. The note states the formalism in AZURE2's conventions (sections 1-2), the rule against
counting a state twice (section 3), what deep cluster potentials with Pauli-forbidden states do to it (section 4),
capture (section 5) and the scope and parameters (section 6).

## 1. The interior R-function of a potential

For a pair α with reduced mass μ (MeV), charge numbers Z₁, Z₂, orbital angular momentum l and (with a spin-orbit term
or a J-dependent depth) total angular momentum J, at the c.m. channel energy E (MeV), the radial equation is that of
`hybrid_potential.rst`,

    u''(r) = (2μ/ħ²) [V_lJ(r) − E] u(r),   V_lJ = V_C + V_N,lJ + ħ² l(l+1)/(2μr²) + v_so ⟨2 l·i⟩.          (1)

Let u_E be the regular solution, u_E(r) ∝ r^{l+1} as r → 0 (V_C of a uniformly charged sphere is finite at the
origin, so the centrifugal term dominates there). Its logarithmic derivative at a_c (fm),

    ℓ(E) = a_c u_E'(a_c) / u_E(a_c)      (dimensionless),                                              (2)

defines, for a boundary-condition constant B (dimensionless), the interior R-function of the potential (Lane and
Thomas, Rev. Mod. Phys. 30 (1958) 257, section IV, single channel; Descouvemont and Baye, Rep. Prog. Phys. 73 (2010)
036301, Eq. (2.13), where the same quantity is computed with a Lagrange mesh, "the calculable R-matrix"):

    R_pot(E; a_c, B) = u_E(a_c) / [a_c u_E'(a_c) − B u_E(a_c)] = 1 / [ℓ(E) − B].                           (3)

R_pot is dimensionless, in the convention of AZURE2's R-matrix R_cc = Σ_λ γ_λc²/(E_λ − E) with γ in MeV^½
(`conventions.rst`). It is a meromorphic function of E with simple poles at the eigenvalues E_p of the potential in the
box [0, a_c] with the condition a_c u'(a_c) = B u(a_c), and residues

    R_pot ≈ γ_p² / (E_p − E),   γ_p² = (ħ²/2μa_c) u_p(a_c)² / ∫₀^{a_c} u_p² dr  (MeV),                      (4)

so that R_pot = Σ_p γ_p²/(E_p − E) exactly (Lane and Thomas, section IV; the box eigenfunctions are complete on
[0, a_c]). Equation (4) follows from the energy derivative of (2), dℓ/dE = −(2μ a_c/ħ²) ∫₀^{a_c} u_E² dr / u_E(a_c)²
(the Wronskian identity), which is negative: R_pot decreases between its poles like any R-function. The deep bound
states of a cluster potential, Pauli-forbidden ones included, are among the E_p (section 4).

Numerics: the Numerov integration of `IntegrateInward` (`physics/HybridPotential.h`) run outward, on the same grid (the
step closest to `step` that puts a_c on a node), started from u(0) = 0 and u(h) = h^{l+1}; the derivative at a_c by the
same five-point difference. One implementation of V(r) (`RadialPotential`) serves both directions.

## 2. Where the background enters the engine's R-matrix

### 2.1 The level matrix with a background

AZURE2 computes, per J^π block, the level matrix (`physics/CollisionMatrix.h`)

    (A⁻¹)_λμ = (E_λ − E) δ_λμ − Σ_c γ_λc γ_μc L⁰_c,   L⁰_c = S_c − B_c + i P_c,
    U_cc' = Ω_c Ω_c' [δ_cc' + 2i P_c^½ P_c'^½ Σ_λμ γ_λc A_λμ γ_μc'],   Ω_c = e^{i(ω_c − φ_c)}.                 (5)

A background in channel c is a term R⁰_c(E) added to the diagonal element: R_cc = Σ_λ γ_λc²/(E_λ − E) + R⁰_c(E).
An energy-dependent R⁰_c is not a pole and has no place in the level sum. It is absorbed instead into the channel's
exterior functions. For one channel with R = R_L + R⁰, where R_L is the level part,

    1 − R L⁰ = (1 − R⁰ L⁰)(1 − R_L L̃⁰),   L̃⁰ ≡ L⁰ / (1 − R⁰ L⁰),                                           (6)

and the collision matrix U = Ω²(1 − R L⁰*)/(1 − R L⁰) becomes

    U = U⁰ + Ω̃² 2i P̃ R_L / (1 − R_L L̃⁰),   U⁰ = e^{2iω} (1 − R⁰L⁰*) / (1 − R⁰L⁰) · e^{−2iφ},                 (7)

with

    P̃ = P / |1 − R⁰ L⁰|²,   φ̃ = φ + arg(1 − R⁰ L⁰),   S̃ − B = Re L̃⁰,   Im L̃⁰ = P̃,   Ω̃ = e^{i(ω − φ̃)},     (8)

so that U⁰ = Ω̃² and U = Ω̃²[1 + 2i P̃ R_L/(1 − R_L L̃⁰)]. The multichannel case is the same with L⁰ diagonal: a
background on the diagonal of one channel factors out of that channel's column. Equation (8) is the whole effect:
**the background is a change of the channel functions (S, P, φ) → (S̃, P̃, φ̃) of the channel it acts on**; the level
matrix, the collision matrix, the T-matrix and every observable keep the form (5). In AZURE2 that is one place, the
channel functions of the pair at a_c.

The same statement on the waves: with H = G + iF the exterior outgoing solution at a_c and ρ = k a_c,

    H̃ = H − R⁰ (ρH' − B H),   ρH̃' = ρH' − B R⁰ (ρH' − B H),                                                  (9)

a real linear map of (u, ρu') with determinant 1, applied to F and G alike; it keeps the Wronskian, so P̃ = ρ/|H̃|²,
φ̃ = arg H̃ and S̃ = Re ρH̃'/H̃ follow from H̃ by the formulas the engine already uses
(`physics/ChannelFunctions.h`). The map acts on the values at a_c only; the waves at r > a_c (the external-capture
integrands) are unchanged (section 5).

### 2.2 The boundary condition of the background

R⁰ = 1/(ℓ − B) depends on B, and the split R = R_L + R⁰ is not invariant under a change of B: a Möbius transformation
of the sum is not the sum of the transforms. B is therefore a number of the background, B_bg, stated in the record,
not the channel's boundary condition of the levels. Once B_bg is fixed, L̃⁰ + B_bg ≡ ℒ(E) is a function of E alone,
and the levels on top see ℒ exactly as they see the exterior log-derivative without a background: the usual
B-invariance of R-matrix theory holds for the level part with ℒ in place of L. In particular

- **Brune's parameters** (Brune, Phys. Rev. C 66 (2002) 044611), the engine's default: the level matrix of
  `fill_level_matrix` reads S_c(E), S_c(E_λ), the chord between two levels and P_c(E); with the background these are
  S̃ and P̃. Brune's derivation uses only that S is a real function of E independent of B, which ℒ is. The on-resonance
  (observed) energy of a level is where (E_λ − E) − γ̃²(S̃(E) − S̃(E_λ)) vanishes, and its observed width
  Γ = 2P̃γ̃²/(1 + γ̃² dS̃/dE): the observed parameters are those of the resonance *on top of the background*, which is what
  a fit is compared with.
- **Formal parameters** (`--no-brune`): the channel's B_c and the formal (E_λ, γ_λ) relative to ℒ; B_bg ≠ B_c is
  allowed.
- **Park's parameters**: S̃'' as S''.
- **Conversion of observed to formal parameters**: the engine's `level_parameters` uses S, S', P at the level energy;
  it uses S̃, S̃', P̃ unchanged. A width read as observed from a fit without background and put into a project with
  background means a different reduced width: the observed width is that of the level on top of the background.

### 2.3 The identity

For a channel with no levels (R_L = 0), U = U⁰ = e^{2iω}(ℓ − L*)/(ℓ − L) e^{−2iφ}, independent of B_bg. If the
exterior functions come from the same potential (the hybrid pair), ℓ − L is the mismatch at a_c between the regular
interior solution and the exterior outgoing one, and U⁰ = e^{2i(ω + δ_pot)} with δ_pot the nuclear phase shift of the
potential (relative to the point-Coulomb phase), *whatever a_c is*. This is the identity the unit test checks
(section 7.2). It holds only when the exterior is the same potential: with a background from a potential and Coulomb
exterior functions, U⁰ is the phase of a potential cut off at a_c.

## 3. Counting a state twice

If a fitted level of the J^π block is also an eigenstate of the potential in that partial wave (a resonance or a bound
state the potential produces), R_pot already contains it as a pole near the level (Eq. (4)), and level plus background
would count it twice. Two rules, each stated in the record:

1. **Subtract the potential's pole.** R⁰ = R_pot − Σ_{p∈P} γ_p²/(E_p − E), with E_p the eigenvalue of the box problem
   with condition B_bg and γ_p² its residue (4), for the eigenstates p named by their node number n (0, 1, ...;
   n = the number of nodes of u_p inside a_c). R⁰ is then regular at E_p; the fitted level replaces the potential's
   pole. With the level at (E_p, γ_p) the potential alone is restored exactly. Consequence: E_p and γ_p depend on a_c
   and B_bg, so the "level on top" is a correction to the potential's state at that radius; its observed (Brune)
   parameters are the physical ones.
2. **Choose a potential with no state there**: a depth for which the partial wave has no eigenstate near the fitted
   level (for instance a potential fitted to the non-resonant part of the phase shift). No subtraction; the level
   carries the whole resonance.

Rule 1 is what the record does when it names `subtract`; rule 2 is the user's choice of depth. A pole that is not
subtracted and lies close to a fitted level makes two poles of nearly equal energy in the channel: the fit then trades
widths between them and their observed parameters are not separately determined. The engine does not detect this.

## 4. Deep cluster potentials and forbidden states

Local potentials for two clusters that respect the Pauli principle bind the forbidden relative states (Buck,
Friedrich and Wheatley, Nucl. Phys. A 275 (1977) 246; Wildermuth's counting N = 2n + l < N_min). For ³He + α, N_min =
3: the s wave has two forbidden states (0s, 1s), the p and d waves one each (0p, 0d), the f wave none; the ⁷Be ground
and first excited states are the 1p states (one node). These are deep bound states of the potential (tens of MeV), and
they are eigenstates of the box too: R_pot has poles near them, at E_p ≈ their binding energies, far below the data.
For R_pot they are distant levels like any other, and they are needed: they fix the number of nodes of the scattering
solution inside a_c and hence ℓ(E) at the energies of the data (Levinson's theorem; the s-wave phase starts at 2π).
They are not subtracted (they are no fitted level) and they are not physical states of ⁷Be: nothing in the engine
reads them as levels. The physical bound states (1p) are poles of R_pot at energies near −1.59 and −1.16 MeV (for a
B_bg near S there); a project that fits them as levels must subtract them (rule 1) or leave the p wave without
background. In this branch the background is refused in a partial wave that the project uses below threshold (section
6), so the p wave of ⁷Be is not given one.

## 5. Capture

The engine's external capture (`ECIntegral.cpp`) has two parts for an entrance channel: the "hard-sphere" term, the
exterior wave cos φ F − sin φ G times the bound state from a_c to infinity, and channel capture, the outgoing wave
√P e^{−iφ} H⁺ fed through the level matrix. Both read φ and P at a_c from the channel functions and F, G at r > a_c
from the pair. With the background (8), φ → φ̃ makes the first term the exterior part of the potential's scattering
wave (cos φ̃ F − sin φ̃ G is that wave beyond a_c, up to normalization), and P̃, φ̃ in the second make the level part
consistent with (7). So the *exterior* capture is consistent with the background with no further change.

What is missing is the *interior* part of the direct capture of the potential's scattering wave, ∫₀^{a_c} u_f r^L u_E dr,
with u_f the bound state inside a_c. Without the background it is carried by the γ-widths of background poles; with
the background those poles are gone. A consistent direct-capture background would add that integral with the
potential's interior solution u_E normalized as in (9) and the potential's bound state normalized to the fitted ANC;
it needs the final state's interior wave function, which an R-matrix level does not have. **Not done in this branch**:
the background acts on the elastic channel functions only, and a capture project with a background has no interior
direct capture from those partial waves. This is the main limit for S(0) (section 9).

## 6. Scope, parameters, validity

- **Channels**: particle channels of a hybrid pair (a pair with a `potential` record), above threshold. The background
  is per partial wave (l, and J where the potential depends on J), on the diagonal; with one spinning particle the
  channel spin is fixed by (l, J) and the background is diagonal in the channels (as the spin-orbit term,
  `spin_orbit_potential.md` §2.2). A partial wave with a background is refused if the project evaluates it below
  threshold (a level of the block below the pair's threshold, a closed channel, an external-capture final state): the
  dressed Whittaker functions of section 4 are not built.
- **What is fixed and what is fitted**: the potential (depth per partial wave, range, R_C) and B_bg are fixed numbers of
  the record; the levels on top are fit parameters as usual. The engine has no potential parameter in the fit vector
  (decided 2026-11-10, `spin_orbit_potential.md` §6 Q7): a depth is fitted by a scan with refits (a profile), which is
  what section 8 does. The gradient of the levels is unchanged (the background changes channel functions, which the
  adjoint reads as constants of the parameters).
- **J-dependent depths**: published cluster potentials for ³He + α use a depth per (l, J) (Dubovichenko, arXiv:1006.4944,
  Table 3.1). The record therefore takes a depth per partial wave, which acts in both halves (interior and exterior),
  so that one potential gives both.
- **Validity**: the identity (2.3) is exact for the pure potential problem. With levels, R_L + R⁰ is a model: the
  potential describes the non-resonant interior, the levels what the potential lacks. The a_c dependence that remains
  comes from the levels and from the missing interior capture; it is measured in section 8.

## 7. What is built

### 7.1 Code

- `physics/HybridPotential.h`: `InteriorSolution` (outward Numerov, the value and derivative at a_c),
  `InteriorRFunction` (Eq. (3)), `InteriorBoxStates` (eigenvalues E_p, residues γ_p², node numbers, Eq. (4)),
  `BackgroundR` (R⁰ with subtraction), `DressWaves` (Eq. (9)).
- (engine and record: see the status line and HANDOFF)

### 7.2 Tests

(filled in with the numbers)

## 8. The radius experiment

(filled in)

## 9. Not done

(filled in)
