# physics/transfer: optical potentials, bound states, distorted waves and DWBA transfer

The Trojan-horse vertices of `THM.md` §3 need distorted waves and a DWBA transfer amplitude that are checked against an
independent reaction code. This note describes the leaf library `physics/transfer` (namespace
`azure2::physics::transfer`, CMake target `azure2_physics_transfer`): an optical potential as a value, the bound state
of two clusters, distorted waves with spin 0, ½ and 1, the transfer A(a, s)F of a cluster x (a = s + x, F = A + x) to
a bound state of F in zero range and in finite range (post form), and the plane-wave limit of the THM entrance
vertex (part 1, §2–§10); the scattering states and continuum bins of x + A, the transfer to them, the narrow-bin limit,
the interior part against the level function, and the prior form with the remnants (part 2, §14–§18). Every part is
compared with a closed form or with FRESCO 3.5 (§9, §18). Nothing of the engine uses the library yet; the wiring is
§12.

Units: MeV, fm, mb/sr; c.m. frame of the pair named; masses in u (×931.49410242 MeV, `physics/Constants.h`).
Kinematics are non-relativistic, as FRESCO's default. The equations below are derived in this note from the
Schrödinger equation and the partial-wave expansions; the only conventions taken from elsewhere are FRESCO's, from its
namelist manual (`fresco-namelist-manual.tex` of FRESCO 3.5: &POT shapes, the spin-orbit factor, the coupling order of
overlaps, KIND 5 and 7), which the library follows so that the benchmarks compare like with like. No equation of
Satchler's *Direct Nuclear Reactions* or of Thompson, Comput. Phys. Rep. 7 (1988) 167, is quoted.

## 1. Files

| File | What |
|---|---|
| `physics/include/azure2/physics/transfer/OpticalPotential.h`, `src/transfer/OpticalPotential.cpp` | Eq. (2.1), radius conventions, the listing `TermsAt` |
| `.../RadialNumerics.h`, `src/transfer/RadialNumerics.cpp` | Numerov outward and inward, Lagrange interpolation, Simpson, Gauss-Legendre, Y_lm(θ, 0), Coulomb functions for matching (§7) |
| `.../BoundState.h`, `src/transfer/BoundState.cpp` | Eq. (3.1), depth search, ANC, D0, analytic tails |
| `.../DistortedWaves.h`, `src/transfer/DistortedWaves.cpp` | Eqs. (4.1)–(4.4) |
| `.../Continuum.h`, `src/transfer/Continuum.cpp` | Eqs. (14.1)–(14.3): scattering states, bins, FRESCO's bin listing |
| `.../Overlap.h`, `src/transfer/Overlap.cpp` | the overlap kinds (bound state, bin, scattering state): one table, one radial form for the kernel |
| `.../Interaction.h`, `src/transfer/Interaction.cpp` | Eq. (17.1): the interaction as a sum of terms (post, prior, remnants), one table |
| `.../Transfer.h`, `src/transfer/Transfer.cpp`, `src/transfer/TransferKernel.{h,cpp}` | Eqs. (5.1)–(6.7); the kernel (6.4) with any interaction of (17.1), per pair of partial waves (l, j) |
| `.../TransferAmplitudes.h`, `src/transfer/TransferAmplitudes.cpp` | Eqs. (19.5)–(19.8): the angular and spin sums, with or without spin-orbit, cross section and iT11 |
| `.../ThmVertex.h`, `src/transfer/ThmVertex.cpp` | Eqs. (8.1)–(8.3) |
| `physics/tests/transfer_*_test.cpp`, `transfer_fresco_common.h`, `physics/tests/data/fresco/`, `physics/tests/data/make_fresco_benchmarks.py` | §9, §18 |

## 2. Optical potential

At the distance r (fm) of the two centres, for the projectile spin s, orbital momentum l and total j = l + s:

```
U(r) = V_C(r) − V f(r; R_V, a_V) − i W f(r; R_W, a_W) − 4 i W_D f(r; R_D, a_D) [1 − f(r; R_D, a_D)]
       + (V_so + i W_so) (ħ/m_π c)² (1/r) df(r; R_so, a_so)/dr ⟨2 l·s⟩,                                   (2.1)
f(r; R, a) = 1/(1 + e^{(r−R)/a}),   (ħ/m_π c)² = 2.000 fm²,   ⟨2 l·s⟩ = j(j+1) − l(l+1) − s(s+1),
```

V, W, W_D > 0 attract and absorb. The surface term is the derivative form, 4f(1 − f) = −4a df/dr, 1 at r = R_D
(FRESCO TYPE=2 SHAPE=0). The spin-orbit term is FRESCO's TYPE=3 SHAPE=0 times its factor 2 l·s; for s = ½,
2 l·s = l·σ, so a strength quoted for l·s without the factor 2 is twice V_so here. Its radial factor is
`SpinOrbitRadialTerm` of `physics/HybridPotential.h` (the same convention), and f is that file's Woods–Saxon form
(`NuclearTerm` of a unit-depth well); V_C is `RadialPotential`'s uniformly charged sphere, Z1Z2e²(3 − r²/R_C²)/(2R_C)
inside R_C and Z1Z2e²/r outside (R_C = 0: point charge), e² = αħc. Each potential states one radius convention for
all its radius parameters r0: R = r0 (`Absolute`), R = r0 A_t^{1/3} (`TargetMass`), R = r0 (A_p^{1/3} + A_t^{1/3})
(`BothMasses`) — FRESCO's TYPE=0 p1 = A_t, p2 = A_p. Explicit numbers only; no global parameter set in this part.
`TermsAt(p, Z1Z2, r)` returns the four terms as FRESCO lists them (fort.34: the Coulomb term per unit Z1Z2 when
called with Z1Z2 = 1; the spin-orbit term before ⟨2 l·s⟩); `Describe` lists the radii.

## 3. Bound state

The radial function u(r) (fm^-1/2) of x relative to the core at binding energy B > 0:

```
u'' + [(2μ/ħ²)(−B − V(r)) − l(l+1)/r²] u = 0,   V(r) = −V0 f(r) + V_C(r) + v_so(r) ⟨2 l·s⟩,                 (3.1)
```

with n interior nodes (FRESCO's NN = n + 1), V0 searched (FRESCO's ISC = 1: only the central depth scales), u
normalised (∫₀^{r_max} u² dr = 1) and positive next to the origin. Beyond the nuclear range u(r) = C W_{−η,l+½}(2κr),
κ = √(2μB)/ħ, η = Z1Z2e²μ/(ħ²κ), W normalised to e^{−κr}(2κr)^{−η} (unit ANC; `physics/coulomb/Whittaker.h`): C
(fm^-1/2) is the single-particle ANC, of sign (−1)ⁿ, FRESCO's convention. Other quantities:

```
r_rms = (∫ r² u² dr)^{1/2},   D0 = ∫ d³r V(r) φ(r) = √(4π) ∫ r V(r) u(r) dr   (l = 0, φ = (u/r) Y_00).        (3.2)
```

The rms radius depends on r_max for a weakly bound state (the deuteron tail beyond 25 fm moves it by 1.5×10⁻³ fm);
the benchmarks use FRESCO's rmatch. Method: Numerov outward from the origin and inward from r_max (starting from W at
the last two nodes); the depth is bracketed by the Sturm count of the outward solution, then found where the
Casoratian z_out(i)z_in(i+1) − z_out(i+1)z_in(i) of z = (1 + h²q/12)u at the matching node (the radius of the volume
term by default), divided by the two norms, vanishes (continuous in V0, zero exactly at an eigenvalue; Illinois
regula falsi, 1e-14 relative). The bracket is widened downward until the Casoratian changes sign, because just above
the eigen-depth the extra zero of the outward solution can lie beyond r_max. Error O(h⁴) for smooth wells, O(h²) for
a step on a node (the mean of the two sides is used there; §9). `SolveBoundStateOfForm` takes any form f and fixed
part; `AnalyticTail` gives the two tails of the engine's THM distortion code, C W(2κr) and C e^{−κr}, zero below
r_min (the engine's φ = u/r without Y_00).

## 4. Distorted waves

```
u_lj'' + [k² − (2μ/ħ²) U_lj(r) − l(l+1)/r²] u_lj = 0,   u(0) = 0,                                           (4.1)
u_lj(r) = (i/2)[H⁻_l(η, kr) − S_lj H⁺_l(η, kr)] = F_l + T_lj H⁺_l,   T = (S − 1)/2i   (r ≥ r_m),            (4.2)
```

k = √(2μE)/ħ, η = Z1Z2e²μ/(ħ²k), H^± = G ± iF; S_lj is the nuclear S-matrix element (the Coulomb factor e^{2iσ_l} is
outside), as FRESCO prints it. The distorted wave of a spinless projectile is
χ⁽⁺⁾(k, r) = (4π/kr) Σ_lm i^l e^{iσ_l} u_l(r) Y*_lm(k̂) Y_lm(r̂), which is e^{ik·r} when U = 0; and
χ⁽⁻⁾*(k, r) = χ⁽⁺⁾(−k, r) = (4π/kr) Σ i^{−l} e^{iσ_l} u_l Y*_lm(k̂) Y_lm(r̂) with the same u (U, not U*). S comes from
two-point matching at r_m and at the node nearest r_m − 1 fm:
S = (H⁻₁u₂ − H⁻₂u₁)/(H⁺₁u₂ − H⁺₂u₁); u is then scaled to Eq. (4.2). For a projectile of spin s (target spin ignored,
as FRESCO without a spin-spin term) the elastic amplitude, derived from the expansion of the incident wave in
|l s j m⟩ and the asymptotic form (4.2), is

```
f_m'm(θ) = f_C(θ) δ_m'm + Σ_lj √(4π(2l+1))/(2ik) e^{2iσ_l} (S_lj − 1) ⟨l 0 s m|j m⟩⟨l m−m' s m'|j m⟩ Y_{l,m−m'}(θ, 0),   (4.3)
f_C = −η/(2k sin²(θ/2)) exp[−iη ln sin²(θ/2) + 2iσ_0],   σ/σ_R = Σ_mm' |f_m'm|² / ((2s+1)|f_C|²),
σ_R(reaction) = (π/k²) Σ_lj (2j+1)/(2s+1) (1 − |S_lj|²).                                                  (4.4)
```

For s = 0 (4.3) is the familiar Σ(2l+1)/(2ik) e^{2iσ_l}(S_l − 1)P_l(cos θ). Without a spin-orbit term one wave per l
is solved and the sums over j run inside (4.3)–(4.4).

## 5. Coordinates of the transfer

With r = r_xA (x relative to A) and ρ = r_sx (s relative to x), α = m_A/m_F, β = m_s/m_a:

```
r_i ≡ r_aA = r + βρ,   r_f ≡ r_sF = αr + ρ,                                                                    (5.1)
r = (r_i − β r_f)/(1 − αβ),   ρ = (r_f − α r_i)/(1 − αβ),   d³r d³ρ = d³r_i d³r_f/(1 − αβ)³.                  (5.2)
```

(5.1) follows from the centres of mass r_a = (m_s r_s + m_x r_x)/m_a and r_F = (m_x r_x + m_A r_A)/m_F; the Jacobian of
(r, ρ) → (r_i, r_f) is the determinant 1 − αβ of [[1, β], [α, 1]], cubed. Masses with m_a = m_s + m_x, m_F = m_A + m_x
(the benchmarks use mass numbers so that this holds exactly).

## 6. Transfer amplitude and cross section

Post form, interaction the binding potential V_sx(ρ) of x in a (Coulomb included, as FRESCO KIND=7 IP1=0 with
KRPOT = KBPOT; `PostInteraction::NuclearPart` drops its Coulomb term, FRESCO with a KRPOT potential that cancels it),
no remnant (FRESCO IP2 = 0):

```
T = ⟨χ⁽⁻⁾_f(r_f) I_F | V_sx(ρ) | I_a χ⁽⁺⁾_i(r_i)⟩,                                                             (6.1)
```

overlap functions coupled as FRESCO couples them, ((l, s_x) j, J_core; J_composite):
I_a = A_a Σ ⟨l_a m_a s_x σ|j_a μ_a⟩⟨j_a μ_a J_s M_s|J_a M_a⟩ φ_{l_a m_a}(ρ) χ_σ, likewise I_F with (l_F, j_F, J_A, J_F),
φ_lm = (u/r) Y_lm, A the signed spectroscopic amplitudes (isospin Clebsch–Gordan factors and √N of identical
nucleons are in A, as for FRESCO; the library adds none). The spin of x is contracted between the two overlaps. With
spin-independent distortions T reduces to orbital amplitudes

```
O_{m_F m_a}(θ) = ∫ d³r d³ρ χ⁽⁻⁾*_f(r_f) φ*_{l_F m_F}(r) V_sx(ρ) φ_{l_a m_a}(ρ) χ⁽⁺⁾_i(r_i).                    (6.2)
```

With k_i along z, k_f at (θ, 0), the expansions of §4, (5.2), φ*_{l_F m_F}(r̂) φ_{l_a m_a}(ρ̂) =
(−1)^{m_F} Σ_LM ⟨l_F −m_F l_a m_a|L M⟩ [Y_{l_F}(r̂) ⊗ Y_{l_a}(ρ̂)]_LM, and the expansion of the rank-L function
G_LM(r_i, r_f) = [u_F(r)/r][V(ρ)u_a(ρ)/ρ][Y_{l_F}(r̂) ⊗ Y_{l_a}(ρ̂)]_LM in bipolar harmonics of (r̂_f, r̂_i),

```
O_{m_F m_a} = (4π)²/(k_i k_f (1−αβ)³) (−1)^{m_F} Σ_{L l_i l_f} ⟨l_F −m_F l_a m_a|L M⟩ i^{l_i−l_f} e^{i(σ_{l_i}+σ_{l_f})}
              √((2l_i+1)/4π) Y_{l_f M}(θ, 0) ⟨l_f M l_i 0|L M⟩ I^L_{l_f l_i},   M = m_a − m_F,                  (6.3)
I^L_{l_f l_i} = ∫∫ dr_i dr_f r_i r_f u_{l_f}(r_f) g^L_{l_f l_i}(r_f, r_i) u_{l_i}(r_i),
g^L_{l_f l_i} = (8π²/(2L+1)) √((2l_i+1)/4π) Σ_M ⟨l_f M l_i 0|L M⟩ ∫_{−1}^{1} dx Y_{l_f M}(x) G_LM(r_i ẑ, r_f r̂_f(x)),   (6.4)
```

x = cos(r̂_i, r̂_f), r̂_f in the xz plane; r̂ then has azimuth π (a factor (−1)^m on Y_{l_F m}) and ρ̂ azimuth 0. The
kernel reduction uses the rotational invariance of Σ_M [Y_{l_f}(r̂_f) ⊗ Y_{l_i}(r̂_i)]*_LM G_LM (∫dΩ_i dΩ_f = 8π²∫dx).
Parity restricts l_i + l_f + l_a + l_F to even values. The spin sum, done explicitly over the Clebsch–Gordan
coefficients (no Racah algebra to get wrong):

```
Σ|T|² = A_a² A_F² Σ_{M_a M_s M_A M_F} | Σ_{σ m_a m_F} ⟨l_a m_a s_x σ|j_a μ_a⟩⟨j_a μ_a J_s M_s|J_a M_a⟩
        ⟨l_F m_F s_x σ|j_F μ_F⟩⟨j_F μ_F J_A M_A|J_F M_F⟩ O_{m_F m_a} |²,                                         (6.5)
dσ/dΩ = (μ_i μ_f/(2πħ²)²)(k_f/k_i) Σ|T|² / ((2J_a+1)(2J_A+1))   (×10 for mb/sr),                              (6.6)
```

μ_i, μ_f the reduced masses of a + A and s + F. Zero range: V_sx φ_00 = D0 δ³(ρ) (l_a = 0; D0 of (3.2) or given), so
r_i = r, r_f = αr (with recoil: FRESCO's INH = 2; FRESCO's default INH = 0 drops α and differs by 25–55 %), and in
general for the integral of the engine's THM weight

```
O_m(θ) = ∫ d³r χ⁽⁻⁾*_f(c_f r) (u(r)/r) Y*_lm(r̂) χ⁽⁺⁾_i(c_i r)
       = (4π)²/(k_i k_f) Σ i^{l_i−l_f} e^{i(σ_i+σ_f)} √((2l_i+1)/4π) Y_{l_f m}(θ,0) √((2l_f+1)(2l_i+1)/(4π(2l+1)))
         ⟨l_f 0 l_i 0|l 0⟩⟨l_f m l_i 0|l m⟩ (c_i c_f)^{−1} ∫ dr u_{l_f}(c_f r) u(r) u_{l_i}(c_i r)/r,            (6.7)
```

post-form ZR: c_i = 1, c_f = α; the engine's M(E) (`src/ThmDistortion.h`): c_i = β, c_f = 1, l = 0, M = √(4π) O_0.

## 7. Numerics

- One outward Numerov (`NumerovOutward`, complex, from the series r^{l+1}(1 + a₁r + a₂r²) with the point-Coulomb
  term a₁ = ηk/(l+1)) and one inward (`NumerovInward`) serve the bound states and the distorted waves. The engine's
  `IntegrateInward` (`physics/HybridPotential.cpp`) is not reused: it integrates a real equation inward from r_m and
  returns only u(a) and u'(a); generalising it to complex coefficients, outward integration and whole waves would
  change the instructions of a routine whose outputs (the hybrid channels of T3) must stay bit-identical. What is
  shared is the potential: the Woods–Saxon form, the spin-orbit factor and the Coulomb sphere are HybridPotential's.
- Matching uses `CoulombForMatching`: COULFG at and beyond the turning point ρ_tp = η + √(η² + l(l+1)), COULCC inside
  (§10 for why not COULCC everywhere).
- Finite range: the (r_i, r_f) grid of step h_K (a multiple of the wave steps) up to r_max, the band
  |r_f − αr_i| ≤ (1 − αβ)ρ_max, the x integral by Gauss–Legendre on the part of [−1, 1] where ρ ≤ ρ_max, the bound
  states by six-point Lagrange interpolation, the radial double integral by the trapezoid rule (the integrand vanishes
  at r = 0 and at r_max). Defaults h_K = 0.1 fm, 48 x nodes, ρ_max = 10 fm; 0.7 s for each FRESCO case below.
- Angular momenta: `physics/AngularCoupling.h` stops at j = 40, so l ≤ 40 (spinless) and l ≤ 39 (spin 1) in (4.3) and
  in (6.3). Enough for the cases here (FRESCO's |S − 1| beyond is < 1e-6); a heavier system needs the limit raised.

## 8. The THM vertex in the plane-wave limit

For a level λ of x + A with level function u (fm^-1/2) inside the channel radius a ((T + V)u = E_λ u, boundary value
B = a u'(a)/u(a)) and the wave w of x relative to A (w = r j_l(pr) without distortion, ε = ħ²p²/2μ):

```
W = ∫₀ᵃ w(r) V(r) u(r) dr,                                                                                  (8.1)
W = (ħ²/2μ) u(a) [B w(a)/a − w'(a)] + (E_λ − ε) ∫₀ᵃ w u dr     (w solves the free equation at ε),          (8.2)
M_l(p) = (B − 1) j_l(pa) − pa j_l'(pa)     (the bracket of (8.2) for w = r j_l(pr)),                         (8.3)
```

(8.2) by replacing Vu with (E_λ − T)u and integrating T by parts twice; (8.3) is `THM.md` §1's M_l, and (8.2) with
(8.3) its full form factor (La Cognata et al., ApJ 708 (2010) 796, Eqs. 20–21, as written in `THM.md` and in the
commit 4882ad0 of the THM branch). `EntranceVertex` takes any wave and its derivative at a, so that the distorted
wave of the next part enters at the same place; `PlaneWaveEntranceVertex` fills in r j_l(pr).

## 9. Validation

Tests (`ctest -R physics_transfer`): `transfer_bound_state`, `transfer_distorted_waves`, `transfer_thm_vertex`,
`transfer_fresco`. Each was run once with an injected error and failed (HANDOFF entry of 2026-11-11).

| Test | Compared with | Asked | Reached (worst, where) |
|---|---|---|---|
| Square well, l = 0, 0 and 1 node (μ = 1 u, R = 3 fm, B = 2 MeV) | root of K cot KR = −κ in the test | 1e-8 in V0 | 6e-9 at h = 0.5 × 10⁻³ fm; error ratio h/(h/2) = 4.0 (O(h²), step on a node) |
| ∫u² dr | trapezoid with end correction | 1 | 1e-9 |
| Tail over the last 5 fm | C e^{−κr} (neutral s), C W(2κr) (p + ¹⁵N p½) | 1e-6 | passes (both) |
| No potential: S, u | 1, F_l(η, kr) (COULFG/COULCC), η = 0 and 0.49, l ≤ 8 | 1e-10 | 2.4e-11 (S), 1.2e-11 (u, rel. to max F) at h = 0.004 fm; at h = 0.001 fm round-off gives 6e-9 |
| Neutral complex square well, l = 0 | e^{2iδ}, K cot KR = k cot(kR + δ) | — | < 1e-6 at h = 2 × 10⁻⁴ fm |
| Real potential (WS + SO + Coulomb), s = ½, 1 | |S| = 1 | 1e-10 | passes |
| M_l(p), l = 0–2, p = 0.1, 0.5, 1.3 fm⁻¹ | j_l in sines and cosines | 1e-6 | passes |
| Full form factor (square-well level) | Lommel integral; (8.1) = (8.2) | 1e-6 | passes |
| FRESCO: constants, E_cm, k, η | its printout | its digits | equal |
| FRESCO: potentials (fort.34), r = 1–10 fm | `TermsAt` | 1e-6 | all cases |
| FRESCO: n-p (d), p-¹⁵N (p½, SO, Coulomb), α-d (2S), α-¹²C (0⁺, 2 nodes): depth, rms, D0, ANC | its listing; ANC from its u/W plateau at 12 fm | its digits | depth 55.30268/55.302 (4 decimals of VR), rms 4.0101/4.010, D0 −130.31/−130.3, ANC 0.905604/0.905604; p-¹⁵N 57.2656/57.265, 2.7440/2.744, 9.634064/9.634060; α-d 77.79654/77.798, 4.0729/4.073, 68.31/68.3, −2.530401/−2.530400; α-¹²C 55.4324/55.43, 3.2503/3.250, 1047.559/1047.56 |
| FRESCO elastic p + ¹⁵N 5 MeV (SO) | S, σ/σ_R, σ_R | 1e-3 | S 8.0e-7; σ/σ_R 2.5e-5 (95°); σ_R 2e-6 |
| FRESCO elastic d + ¹⁵N 10 MeV (spin 1, SO) | same | 1e-3 | S 8.4e-6; σ/σ_R 4.6e-5 (115°) |
| FRESCO elastic ⁶Li + ¹²C 20 MeV | same | 1e-3 | S 1.1e-7; σ/σ_R 3.1e-5 (150°) |
| FRESCO ¹⁵N(d,n)¹⁶O, finite range | dσ/dΩ, 0–90°, > 1 % of peak | 1 % | 3.3e-4 (75°); FRESCO ref − refined 1.7e-4 |
| FRESCO ¹⁵N(d,n)¹⁶O, zero range | same | 1 % | 7.2e-6 (45°) |
| FRESCO ¹²C(⁶Li,d)¹⁶O, finite range (nuclear V_αd) | same | 1 % | 1.8e-3 (35°); FRESCO ref − refined 2.8e-4, its h dependence below |
| FRESCO ¹²C(⁶Li,d)¹⁶O, zero range | same | 1 % | 4.4e-5 (75°) |

FRESCO's convergence (inputs and all runs in `physics/tests/data/fresco/<case>/`): every case was run at the
reference settings and at a refined set (smaller hcm, rintp, hnl or larger rmatch); the test checks that the two agree
to 2e-5 in S, 1e-4 in σ/σ_R and 0.2 % in dσ/dΩ. (d,n) finite range: hcm 0.025 → 0.02 fm, rintp 0.05 → 0.04, hnl
0.025 → 0.02, rnl = 10 fm (FRESCO recommends > 3.4), nnu = 72: change ≤ 1.7e-4; at FRESCO's defaults (nnu, rnl 4)
the cross sections are 0.5 % low. (⁶Li,d) finite range: FRESCO's result moves with hcm linearly (0.05/0.025/0.0125 fm:
0.244189/0.243965/0.243839 mb/sr at 0°, 0.013802/0.013740/0.013722 at 45°), so its remaining step error is about
1e-3; ours moves by < 1e-5 between h = 0.05/0.1 and 0.025/0.05 fm and between ρ_max = 10 and 15 fm. That is the
1.8e-3 left at 35°.

The whole V_αd as interaction (its Coulomb tail 2e²/ρ × u_a(ρ), u_a ~ e^{−0.31ρ}) makes the post-form kernel long
ranged. Ours converges slowly in ρ_max (0°: 0.18650, 0.18367, 0.18422, 0.18423 mb/sr for ρ_max = 10, 15, 20, 25 fm).
FRESCO, with epc = 1e-5, asks for rnl > 19 fm with centre −3.9 fm, and then does not converge: 0.18444, 0.18390,
0.16491, 0.12319 mb/sr at 0° for rnl = 16, 22, 26, 30 fm (`li6d_finite_range_coulomb`, and larger rmatch makes it
worse). At rnl = 22 fm it is 0.2 % (0°), 0.6 % (45°, 90°) from ours. This is a limit of FRESCO's non-local kernel with
a long-ranged vertex, not a difference of physics, so the benchmark uses the nuclear part (KRPOT); the case is kept as
data and not tested.

## 10. Defects and limits found

- **COULCC sign errors (physics/coulomb, outside this branch):** `coulomb::coulomb(l, eta, rho, Method::coulcc)`
  returns F and G with the wrong sign at isolated points beyond the turning point: l = 1, η = 0.49300738446815717,
  ρ = 8.5485146564819541 (F = +0.3622592730 against −0.3622592730 of COULSM and COULFG); a scan of 2.6 × 10⁵ points
  (l ≤ 12, ρ_tp < ρ < ρ_tp + 34) found one such point at η = 0.1 (l = 0, ρ = 7.5943) and one at η = 0.493; none at
  η = 0, 1.55, 3. The physics_coulomb test does not hit them. The engine uses COULSM by default; any user of COULCC
  beyond the turning point is exposed. Not fixed here.
- Not implemented: global optical potentials, relativistic
  kinematics, target spin in the distortion (spin-spin), tensor forces, l > 40, complex (optical or coupled-channel)
  bins. Part 2 added the continuum, the prior form and the remnants (§14–§17); their limits are in §18. Part 3 added
  spin-orbit in the distorted waves of a transfer and the projectile's iT11 (§19); its validation is §22. The ³He Trojan horse is §20; the ¹⁴N
  Trojan horse is not yet benchmarked. The prior form and the remnant terms (§17) with spin-orbit waves go through
  the same sums (the interaction terms are scalar), but are not compared with FRESCO.

## 11. The continuum part

THM needs the transfer to x + A at E > 0: in (6.2) the bound overlap of F is replaced by the x + A scattering function
at E or a bin of them, and the vertex of §8 is the interior part of that amplitude. This is part 2, §14–§18.

## 12. Wiring (later branch)

What the engine's THM observable (Phase 5, σ_HOES on the shared energy mesh) will call, once per fit stage:

- the distortions: `SolveDistortedWaves(OpticalPotential, ScatteringRequest)` for a + A and s + F at each mesh energy,
  the potentials from a global-set table (Koning–Delaroche, An–Cai, GDP08, Avrigeanu: one row each, moved from
  `ThmOptical`, a test per row; never a second copy of the numbers);
- the Trojan horse: `SolveBoundState(well, BoundStateRequest)` → `RadialOverlapOf(BoundState)`;
- per level λ and channel: the level function u_λ (u_E at E_λ inside a, §16) and the vertex
  `EntranceVertex(LevelFunction, w, w'(a), ε, B)` with w the projection of the DWBA source on the channel; the
  zero-range prior weight M(E) is `ZeroRangeOrbitalIntegral(RadialOverlap, in, out, β, 1, θ)` (§13);
- checks and the full DWBA where wanted: `FiniteRangeCrossSection(TransferReaction, in, out, TransferGrid, θ)` with
  the residual a `ScatteringState` (d²σ/dΩ dE, §15) and `interaction = PriorFormWithRemnant(...)` or the post form.

`ThmDistortion` becomes a caller of `ZeroRangeOrbitalIntegral` (no physics of its own), `ThmDwVertex` of
`EntranceVertex` with the distorted w; their tables stay owned by the stage. Not done here.

## 13. The engine's THM weight with this library

The engine's zero-range prior-form weight M(E) of `src/ThmDistortion.cpp` is Eq. (6.7) with c_i = β, c_f = 1, l = 0
and the Whittaker tail of unit ANC (`AnalyticTail`), M = √(4π) O_0. For the case `optical` of
`tests/reference/thm_distortion_test.cpp` (¹²C(¹⁴N,d) at 30 MeV, x = ¹²C, E = 1.0 MeV, QF direction, Woods–Saxon
potentials in both channels; the masses, energies and potentials of that test, our constants):

| | |M|² (fm³ in the engine's units) |
|---|---|
| engine (`thm_distortion_test -v`, build of this branch) | 2.705443758408 × 10⁻⁴ |
| this library, h = 0.01 fm, l ≤ 40 | 2.705443507050 × 10⁻⁴ (−9 × 10⁻⁸) |
| this library, h = 0.005 fm | 2.705444832817 × 10⁻⁴ (+4.0 × 10⁻⁷) |
| the test's mpmath/scipy reference (older ħc, u) | 2.705453160181 × 10⁻⁴ |

The library covers the same integral; the replacement needs no new physics, only the THM kinematics (k_sF, k_aA and
their angle) from the engine, which stay with the THM observable of Phase 5. Scratch driver:
`work/dwba/drivers/engine_m.cpp` (not in the repository).

## 14. Scattering states and bins of x + A

**Scattering state.** In a real potential V (Coulomb, Woods–Saxon, spin-orbit times ⟨2 l·s⟩ for the (l, j) asked) the
regular solution of (3.1) at −B = E > 0 is, beyond the nuclear range, a combination of F_l and G_l. Of the distorted
wave (4.2), u = (i/2)(H⁻ − S H⁺) with S = e^{2iδ} (|S| = 1 for a real V),

```
(i/2)(H⁻ − e^{2iδ} H⁺) = e^{iδ} [cos δ F_l + sin δ G_l],   so   u_E ≡ e^{−iδ} u = cos δ F_l + sin δ G_l,          (14.1)
```

real, with asymptotic amplitude 1: u_E → sin(kr − lπ/2 − η ln 2kr + σ_l + δ). `SolveScatteringState` returns u_E,
δ in (−π/2, π/2] and V on the grid (it calls the distorted-wave solver of §4 and multiplies by e^{−iδ}).

**Normalisations.** For two such states of one Hamiltonian, (k² − k'²)∫₀ᴿ u_k u_k' dr = [u_k u'_k' − u'_k u_k']_R (the
radial equations multiplied crosswise and integrated). With the asymptotic forms the bracket is
½[(k + k') sin(a − b) + (k' − k) sin(a + b)] (a = kR + φ, b = k'R + φ'), and
sin((k − k')R + ...)/(k − k') → π δ(k − k') as R → ∞, so ∫₀^∞ u_k u_k' dr = (π/2) δ(k − k'). Hence

```
δ(k − k'):  u_k = √(2/π) u_E (fm^1/2);      δ(E − E'):  u^(E) = √(2μ/(πħ²k)) u_E (MeV^−1/2 fm^−1/2),           (14.2)
```

the second from δ(k − k') = (dE/dk) δ(E − E') = (ħ²k/μ) δ(E − E'). `NormalisationFactor` gives the factors.

**Bins** (FRESCO's `&OVERLAP` with BE < 0; its namelist manual, "CONTINUUM BINS", and the bin listing of its output):

```
u_bin(r) = √(2/(πN)) ∫_{k_min}^{k_max} dk w(k) u_E(k)(r),   N = ∫ dk w(k)²,                                     (14.3)
```

w = 1 (ISC = 2: FRESCO weights its wave e^{iδ}u_E by e^{−iδ}), w = sin δ (ISC = 4: by T* = e^{−iδ} sin δ), times k for
ISC = 12, 14; E_min = max(E − |ER|/2, 10⁻⁵ MeV), E_max = E_min + |ER|. By (14.2), ∫|u_bin|² = (2/(πN))(π/2)∫w² = 1,
and two bins of one potential over disjoint intervals, or a bin and a bound state of it, are orthogonal. δ is
continued across the bin from its principal value at k_min (a jump of more than π/2 between nodes adds π, as
FRESCO does), so that u_E(k) is continuous in k. The k integral: FRESCO's rule (NK equally spaced nodes including
both ends, each weighted (k_max − k_min)/(NK − 1), in the integral and in N; it converges as 1/NK) or Simpson's rule.
FRESCO builds the bin on [0, RSP] (RSP = rmatch unless set) and leaves the value at its last node unnormalised
(ours is normalised there; the listing of §18 compares the other nodes).

**The 1/R tails.** Beyond the potential (neutral, l = 0, w = 1) integrating (14.3) by parts gives
u_bin ≈ √(2/(πΔk)) [cos(k_min r + δ₁) − cos(k_max r + δ₂)]/r + O(r⁻²), so ∫_R^∞ u_bin² dr = 2/(πΔk R) + O(R⁻²) and,
for two adjacent bins sharing k₂, ∫_R^∞ u_b u_b' dr = −1/(π√(ΔkΔk') R) + O(R⁻²). The test of orthonormality uses
these (§18). The same tails make the transfer to a bin depend on how far the bin is kept (§15).

## 15. Transfer to the continuum, FRESCO's bins and the narrow-bin limit

**Kernel.** The kernel (6.4) reads only u_F and V on a grid: the residual overlap is a row of the table
`OverlapKinds()` (bound state, bin, scattering state) and enters through `RadialOverlap`; the bound state and the bin go
through the same code. The exit channel s + F* has E_f = E_i − B_a − E (B_a the binding of x in a, E the x + A
energy; FRESCO: Q of the ground state and the state energy E_x = S_x + E).

**d²σ/dΩ dE.** Let the final x + A states be normalised to δ(E − E') per partial wave (l, j, J_F). The sum over final
states of the golden rule becomes ∫dE over each partial wave, and the projections of different partial waves are
orthogonal when the x + A direction is not observed, so they add incoherently:

```
d²σ/dΩ_s dE = (μ_i μ_f/(2πħ²)²) (k_f(E)/k_i) Σ|T(E)|² / ((2J_a + 1)(2J_A + 1)),   u_F = u^(E) of (14.2),           (15.1)
```

the formula (6.6) with u^(E) in place of a bound state: mb/(sr MeV). The phase e^{−iδ} of the incoming-wave final
state is a constant per partial wave and drops out of |T|². **Narrow-bin limit:** for a bin of width ΔE around E
(w = 1), u_bin = √(2/(πΔk))∫u_E dk; with u_E = u^(E)/√(2μ/(πħ²k)) and dk = (μ/(ħ²k))dE, u_bin ≈ √ΔE u^(E)(1 + O(ΔE²)),
so σ_bin/ΔE → d²σ/dΩ dE with an O(ΔE²) error (checked: ratio of errors 4.00 between ΔE = 50 and 25 keV, §18).

**Relation to THM.md §1.** With plane-wave normalised continuum states ψ_k (∫d³k/(2π)³|ψ_k⟩⟨ψ_k| = 1 on the x + A
continuum) and d³k_xA = k_xA²dk_xA dΩ_xA = (μ_xA k_xA/ħ²) dE dΩ_xA, the three-body cross section of A(a, s)[x + A] is

```
d³σ/(dΩ_s dE dΩ_xA) = (μ_i μ_f/(2πħ²)²)(k_f/k_i) (μ_xA k_xA/((2π)³ħ²)) Σ|T(k_f, k_xA)|² / ((2J_a + 1)(2J_A + 1)),   (15.2)
```

whose integral over Ω_xA is (15.1) (the partial-wave expansion of ψ_k gives the δ(E) states of (14.2)). The measured
variables of THM are the lab energy and angles of the two fragments c, C of x + A → c + C (for elastic x + A, c = x,
C = A) with s undetected. The phase space of three free particles, ∫δ⁴(P_f − P_i) d³p_c d³p_C d³p_s, integrated
over p_s (momentum) and over |p_C| (energy, with p_s = P − p_c − p_C), with d³p_c = m_c p_c dE_c dΩ_c, gives

```
d³σ/(dE_c dΩ_c dΩ_C) ∝ |T|² · m_c p_c p_C² / |p_C/m_C − (p_s·p̂_C)/m_s|   (the kinematic factor KF of THM.md §1),   (15.3)
```

and with the same |T|² the ratio of (15.3) to (15.2) is the Jacobian between (E_c, Ω_c, Ω_C) and (Ω_s, E, Ω_xA). In
the plane-wave impulse limit T factorises into the momentum distribution φ(p_s) of x in a and the x + A → c + C
half-off-shell amplitude, which is THM.md §1's KF |φ(p_s)|² dσ^HOES/dΩ; (15.1) is the d²σ/dΩ dE at the vertex the
library computes with distortions, and its interior part is §16.

**FRESCO's bins.** FRESCO builds the bin on [0, RSP]. The transfer to a bin then depends strongly on RSP when the bin
reaches far (15N(d,n), narrow resonance, ER = 1 MeV: 88 mb/sr at 0° with RSP = rmatch = 40 fm, 58 mb/sr at 60 fm, in
both codes): the distorted waves are not absorbed outside the nucleus and the bin's tail (§14) is momentum-matched over
tens of fm. The benchmark keeps RSP = 40 fm fixed (the `&FRESCO` variable `rsp`) and converges rmatch separately, and
chooses the bins' widths so that the bin's envelope has a node near RSP (its value at the cut is 1–4 % of its
maximum): a bin cut where it is large makes the integrand discontinuous and both codes' steps then matter at 1 %.
The three bins (p + 15N, kp 5, no spin-orbit; 15N(d,n) at E_lab = 10 MeV):

| Bin | l, V | E (MeV) | Resonance | Weight |
|---|---|---|---|---|
| non-resonant | 1, 50 MeV | 1.055–2.945 | none below 8 MeV | ISC = 2 |
| narrow | 2, 53 MeV | 1.445–2.355 | E_r = 1.9015 MeV, Γ = 95 keV (2/(dδ/dE) at δ = π/2) < width | ISC = 4 |
| wide | 0, 47 MeV | 1.000–2.160 | E_r ≈ 1.58 MeV, Γ ≈ 3.3 MeV (2/(dδ/dE), 10 keV scan) > width | ISC = 4 |

**FRESCO's limits met.** (1) Zero range (KIND = 5, INH = 2) to a bin is not consistent with FRESCO's own finite range
(narrow bin: 554 against 78.5 mb/sr at 0°; non-resonant: an oscillating shape), while ours is within 4 % of our finite
range, as for the bound state; with INH = 0 FRESCO gives a smooth distribution again (no recoil, a different
quantity). The zero-range bin data are kept, not compared. (2) For ¹²C(⁶Li,d) to α + ¹²C bins FRESCO's ERWIN loses 14
or more digits in the d + ¹⁶O* channel at high L (exit energies 8–10 MeV, turning points beyond rmatch): its results
move by 5–40 % with cutl at RSP = 20–30 fm. With RSP = 12 fm it is stable to cutl to 0.5 % (`li6d_bin_wide_fr`, data
only) and agrees with ours to 0.1–1.7 %, where ours moves by up to 1.6 % with the step because the cut is a step of
the integrand. Not a benchmark; numbers in §18.

## 16. The interior part and the level function

Inside the channel radius a, a scattering state of a potential with an isolated resonance is, up to the energy
dependence, the level function of §8. Let u_λ be u_E at E_λ = E_r on [0, a], normalised there; its boundary value is
B = a u_λ'(a)/u_λ(a) (so the shift vanishes at E_r, THM.md §3 item 2), γ_λ = √(ħ²/(2μa)) u_λ(a). For any E,
multiplying the radial equations of u_E and u_λ crosswise and integrating over [0, a] (Green's identity, as for (8.2)):

```
A_λ(E) ≡ ∫₀ᵃ u_λ u_E dr = (ħ²/2μ) u_λ(a) [u_E'(a) − (B/a) u_E(a)] / (E − E_λ),                                   (16.1)
∫₀ᵃ w V u_E dr = Σ_μ A_μ(E) W_μ,   W_μ = ∫₀ᵃ w V u_μ dr   (the vertex (8.1) of the level μ),                           (16.2)
```

(16.2) from the completeness of the levels of B on [0, a]. Near E_r the term μ = λ dominates; at E_r it is the whole
sum (u_E ∝ u_λ inside a). If λ alone is kept, u_E(a) = Σ_μ u_μ(a)A_μ gives A_λ ≈ √(ħ²/(2μa)) u_E(a)/γ_λ, and with
(14.1) u_E(a) = |H⁺(ka)| sin(δ − φ) (φ the hard-sphere phase), i.e. the Breit–Wigner amplitude
√(Γ(E)/2)/|E_λ − E − iΓ/2| in the δ(E − E') normalisation up to √(1/π). The table: p + 15N, l = 2, V = 53 MeV
(E_r = 1.901482 MeV, Γ = 95.19 keV), w = r j₂(pr) with ε = ħ²p²/2μ = E + 2.2246 MeV (the d → p + n Trojan horse at
zero spectator energy), u_E in the δ(E − E') normalisation, h = 0.01 fm:

| a (fm) | E | ∫₀ᵃ wVu_E dr | A_λ W_λ | remainder (1 − A_λW_λ/∫) | exterior ∫_a^∞ w V_N u_E / interior | A_λ | one-level A |
|---|---|---|---|---|---|---|---|
| 5 | E_r − Γ | −8.8469 | −8.8581 | −1.26 × 10⁻³ | 9.29 × 10⁻² | 0.9613 | 0.9815 |
| 5 | E_r | −21.4201 | −21.4201 | 5 × 10⁻¹⁶ | 8.76 × 10⁻² | 2.2786 | 2.2786 |
| 5 | E_r + Γ | 10.2293 | 10.2162 | 1.28 × 10⁻³ | 8.24 × 10⁻² | −1.0658 | −1.0434 |
| 6 | E_r − Γ | −8.6247 | −8.6805 | −6.47 × 10⁻³ | 2.32 × 10⁻² | 1.0227 | 1.0793 |
| 6 | E_r | −20.9120 | −20.9120 | 2 × 10⁻¹⁵ | 2.06 × 10⁻² | 2.4145 | 2.4145 |
| 6 | E_r + Γ | 10.0007 | 9.9363 | 6.44 × 10⁻³ | 1.81 × 10⁻² | −1.1249 | −1.0623 |
| 7 | E_r − Γ | −7.9803 | −8.1297 | −1.87 × 10⁻² | 5.84 × 10⁻³ | 1.0602 | 1.1879 |
| 7 | E_r | −19.5028 | −19.5028 | 2 × 10⁻¹⁵ | 4.71 × 10⁻³ | 2.4904 | 2.4904 |
| 7 | E_r + Γ | 9.3979 | 9.2273 | 1.82 × 10⁻² | 3.65 × 10⁻³ | −1.1544 | −1.0139 |

(MeV fm^3/2 MeV^−1/2 for the vertices.) What it shows about the engine's vertex: the interior continuum amplitude is the
level's vertex times A_λ(E), exactly at E_r and to 0.13–1.9 % at E_r ± Γ (the other levels of B: they grow with a);
the exterior nuclear part, which the engine's surface/full form factor leaves to the channel, is 9 % at a = 5 fm and
0.5 % at 7 fm (V_N at a: 2.6 and 0.12 MeV). The one-level Breit–Wigner amplitude differs from the projection A_λ by
2–14 % at E_r ± Γ for the same reason. So the engine's form factor is the right limit when a lies beyond the nuclear
potential and one level dominates; the remainder is the background the R-matrix puts in other levels.

## 17. Prior form and the remnants

With H = T + V_sx + V_xA + V_sA, the exact final state ⟨χ_f I_F| satisfies (E − T_f − U_sF − H_F) = 0 on the left and the
initial state (E − T_i − U_aA − H_a)|I_a χ_i⟩ = 0, so T = ⟨χ_f I_F| H − T − U_sF − H_F + ... ⟩ gives the two forms

```
V_post = V_sx(r_sx) + [V_sA(r_sA) − U_sF(r_sF)],     V_prior = V_xA(r_xA) + [V_sA(r_sA) − U_aA(r_aA)],                (17.1)
```

with r_sA = r + ρ (vectors of §5); they give the same DWBA amplitude when both remnants are kept (the difference
⟨χ_f I_F|(H − E)... − (H − E)|I_a χ_i⟩ vanishes because each state solves its own equation; with complex optical
potentials this needs χ⁽⁻⁾* to solve the equation with U, not U*, which is how §4 defines it). In the kernel the
interaction is a sum of `InteractionTerm`s (`Interaction.h`), one row each in `InteractionKinds()`: the binding term
of the projectile (post) or the residual (prior) is tabulated as V u on its overlap's grid, the others (V_sA,
−U_sF, −U_aA, complex) are evaluated at every kernel node from (r, ρ) and multiply u_F(r) u_a(ρ). Nothing beyond a
potential difference on the kernel's grid was needed. The prior form and the remnants multiply u_a, not V_sx u_a: the
deuteron's tail needs ρ ≤ 40 fm (post = prior to 1.1 × 10⁻³ there, 4 × 10⁻⁵ at 0°; at ρ ≤ 20 fm they differ by 5 %).

FRESCO (KIND = 7, IP1 = 1 prior, IP2 = −1 complex remnant, IP3 = KPCORE) does not converge in its non-local range
for these runs: with the remnant the 0° cross section is 0.29, 0.37, 1.10 mb/sr for rnl = 20, 26, 32 fm (ours,
converged in ρ: 0.3833 post, 0.3833 prior); the prior form without remnant (ours, converged: 9.617 mb/sr at 0°) is 9.685, 9.590, 9.815 mb/sr
for rnl = 20, 26, 32 fm: within 0.3–2.1 % of ours at 0–90°, but FRESCO itself moves by 2 % with rnl (ERWIN accuracy
warnings in every run). Data in `dn_prior*`, `dn_post_remnant`; not compared
by the test.

## 18. Validation of part 2

Tests: `transfer_continuum`, `transfer_continuum_transfer`, `transfer_continuum_fresco` (and `transfer_fresco`
unchanged in its numbers). Each was run once with an injected error and failed (HANDOFF entry of 2026-11-12).

| Test | Compared with | Asked | Reached (worst, where) |
|---|---|---|---|
| Scattering state, no potential, l ≤ 6 | F_l(η, kr), δ = 0 (η = 0, 0.49) | 1e-10 | 1.2e-11 (h = 0.004 fm) |
| Neutral square well, l = 0 | δ and u from K cot KR = k cot(kR + δ) | 1e-6 | 4.5e-8 (δ), 6.0e-8 (u) at h = 2e-4 fm |
| δ(E − E') normalisation (14.2) | bin of ΔE = 0.1 keV / √ΔE | 1e-6 | 1.6e-10 |
| ⟨bound \| bin⟩, one well | 0 | 1e-6 | 4.7e-9 |
| ∫u_bin² on [0, R], adjacent and separated bins | 1 − 2/(πΔkR), 1/(πΔkR), 0 (§14) | 1e-5 | 2.1e-6, 4.7e-6, 1.5e-6 (R = 1000 fm) |
| FRESCO bin listing (3 bins) | energies, δ at the nodes, Norm, RMS, D0, u(r) | its digits | equal; δ to 1.7e-3 deg (3 decimals printed); u 1.6e-4 of max |
| FRESCO ¹⁵N(d,n) → non-resonant bin, finite range | dσ/dΩ, 0–90°, > 1 % of peak, NK = 1600 / 3200 | 1 % | 4.3e-3 / 1.1e-3 (25°, a minimum); FRESCO steps 3.2e-3, rmatch 2.9e-4, rnl 2.9e-4, NK 1600→3200 8e-3 |
| FRESCO ¹⁵N(d,n) → narrow-resonance bin | same | 1 % | 1.8e-3 / 8.7e-4 (85°); FRESCO steps 8.6e-4, rmatch 3.0e-4, rnl 2.7e-4 |
| FRESCO ¹⁵N(d,n) → wide-resonance bin | same | 1 % | 5.2e-3 / 2.6e-3 (85–90°); FRESCO steps 4.7e-3, rmatch 3e-5, rnl 1.3e-4 |
| FRESCO zero range to bins | — | 1 % | not reached: FRESCO's INH = 2 to a bin is inconsistent with its own finite range (×7, §15) |
| FRESCO ¹²C(⁶Li,d) → wide-resonance bin (RSP 12 fm) | its refined run | 1 % | 1.1e-3, 6e-3, 1.5e-2, 1.1e-2 at 0, 30, 60, 90°: not reached; FRESCO moves by 0.3–0.6 % with hcm and cutl, ours by up to 1.6 % with h (the cut at RSP) |
| Narrow-bin limit (15) | d²σ/dΩ dE from u^(E) | 1e-3 | 6.3e-4 at ΔE = 25 keV; error ratio 4.00 between 50 and 25 keV |
| Interior vertex (16.2) | A_λ W_λ, table of §16 | 1e-10 at E_r; table 1e-3 | 2e-15; reproduced |
| Post = prior with remnants (17.1) | each other | numerics | 1.1e-3 (35°), 8e-5 at 0° (ρ ≤ 40 fm) |
| FRESCO prior form / remnant | dσ/dΩ | 1 % | not reached: FRESCO does not converge in rnl (§17) |

FRESCO's convergence asked was 0.2 %: its step refinement (hcm 0.025 → 0.0125) moves the bin cross sections by
0.09 %, 0.32 % and 0.47 % at their worst angles (minima near 1 % of the peak); its k rule moves them by 1/NK (0.8 %
from NK 1600 to 3200 at the non-resonant minimum, the same in ours with the same rule: ours − FRESCO 1.1e-3 at NK
3200). Simpson's rule (also in the library) converges as 1/NK⁴; it is used for the narrow bins of §15.


## 19. Spin-orbit in the distorted waves of a transfer

The distorted waves of §4 with a spin-orbit term depend on (l, j). The projectile a has spin J_a (the composite of
the projectile overlap) and the ejectile s spin J_s (its core); the target's spin stays a spectator, as in §4. The
transfer amplitude then mixes the projections of J_a and J_s between the asymptotic states and the vertex.

**Entrance wave.** The incident wave e^{ik·r}|J_a m⟩, k = k_i ẑ, expanded in partial waves (§4), with Y_{l m_l}(r̂)|J_a m⟩
coupled to |(l J_a) j⟩, each component multiplied by its u_lj and decoupled again:

```
χ⁽⁺⁾_{m'm}(r) = (4π/(k r)) Σ_lj i^l e^{iσ_l} √((2l+1)/4π) u_lj(r) ⟨l 0 J_a m|j m⟩⟨l m−m' J_a m'|j m⟩ Y_{l,m−m'}(r̂),   (19.1)
```

the amplitude of spin state m' at r when m comes in (Y*_{l m_l}(ẑ) = √((2l+1)/4π) δ_{m_l 0}).

**Exit wave.** The bra ⟨χ⁽⁻⁾_{k_f μ}|r μ'⟩ = [χ⁽⁻⁾_{μ'μ}(k_f, r)]*. χ⁽⁻⁾ has the expansion of (19.1) with k̂_f, i^l e^{−iσ_l}
and the incoming-wave radial function, whose complex conjugate solves the equation with U (§4: "U, not U*"), so it is
u_lj:

```
χ̃_{μμ'}(r) = (4π/(k_f r)) Σ_lj Σ_{n n'} i^{−l} e^{iσ_l} u_lj(r) Y_{l n}(k̂_f) ⟨l n J_s μ|j N⟩⟨l n' J_s μ'|j N⟩ Y*_{l n'}(r̂),   (19.2)
```

N = n + μ = n' + μ'. Without spin-orbit the sums over j of (19.1), (19.2) are δ_{m m'}δ_{m_l 0} and δ_{μμ'}δ_{nn'} by the
orthogonality of the Clebsch–Gordan coefficients (Σ_j ⟨l a s b|j c⟩⟨l a' s b'|j c⟩ = δ_aa'δ_bb'), and §4's waves return.

**Angular reduction.** The vertex of (6.2) is (−1)^{m_F} Σ_LM ⟨l_F −m_F l_a m_a|L M⟩ G_LM with G_LM a rank-L tensor of the
pair (r̂_i, r̂_f). By rotational invariance (Wigner–Eckart for the coupled product Y_{l_i m_i} G_LM, projected on Y_{l_f n'}),

```
∫∫ dΩ_i dΩ_f Y*_{l_f n'}(r̂_f) Y_{l_i m_i}(r̂_i) G_LM = ⟨l_i m_i L M|l_f n'⟩ R^L_{l_f l_i},                                (19.3)
```

R independent of the projections. At m_i = 0 this is the quantity of (6.3)–(6.4); with
⟨l_f M l_i 0|L M⟩ = (−1)^{l_i} √((2L+1)/(2l_f+1)) ⟨l_i 0 L M|l_f M⟩ (symmetries of the coefficients; the reduction test
below checks it to 7 × 10⁻¹⁶),

```
R^L_{(l_f j_f)(l_i j_i)} = (−1)^{l_i} √((2L+1)/(2l_f+1)) I^L_{(l_f j_f)(l_i j_i)},                                         (19.4)
```

with I^L the kernel integral (6.4) of the radial waves u_{l_f j_f}, u_{l_i j_i}: computed once per pair of waves (the
kernel's angular part h^L_{l_f M} depends on l_f only and is shared by the j's).

**Amplitude.** With m_i = m − m', M = m_a − m_F, n' = m_i + M:

```
Q_{m_F m_a}(m', m; μ', μ; θ) = P (−1)^{m_F} Σ_{L, (l_f j_f), (l_i j_i)} ⟨l_F −m_F l_a m_a|L M⟩ i^{l_i−l_f} e^{i(σ_{l_i}+σ_{l_f})} √((2l_i+1)/4π) R^L
     × ⟨l_i 0 J_a m|j_i m⟩⟨l_i m_i J_a m'|j_i m⟩ ⟨l_i m_i L M|l_f n'⟩ Σ_n ⟨l_f n J_s μ|j_f N⟩⟨l_f n' J_s μ'|j_f N⟩ Y_{l_f n}(θ, 0),      (19.5)
T(μ, M_F; m, M_A) = A_a A_F Σ_{m' μ' σ} ⟨l_a m_a s_x σ|j_a μ_a⟩⟨j_a μ_a J_s μ'|J_a m'⟩⟨l_F m_F s_x σ|j_F μ_F⟩⟨j_F μ_F J_A M_A|J_F M_F⟩
     × Q_{m_F m_a}(m', m; μ', μ),   μ_a = m' − μ', μ_F = M_F − M_A,                                                         (19.6)
```

P = (4π)²/(k_i k_f (1 − αβ)³); dσ/dΩ is (6.6) with Σ|T|² over all four projections. Without spin-orbit the j sums of
(19.5) reduce it to

```
Q = δ_{m m'} δ_{μ μ'} O_{m_F m_a}(θ)   (Eq. (6.3)),   and (19.6) to (6.5).                                                  (19.7)
```

One function (`TransferAmplitudes.cpp`) does (19.5)–(19.6) for every case: a wave set without spin-orbit enters with all
j of (l, J) on its one radial wave, so the identity (19.7) is computed, not assumed. The θ-independent part of (19.5)
(everything but Y_{l_f n}) is formed once per reaction. Zero range is the same table with l_a = 0, L = l_F and
I = D0 g_{l_f l_i} R_{l_f l_i} (1 − αβ)³ (g, R of (6.7)), so finite and zero range share (19.5)–(19.6).

**Analysing power.** With T(μ, M_F; m, M_A) for every projection the projectile's polarisation observables need no
further integral. In the Madison frame (ŷ along k_i × k_f: our k_f at (θ, 0) gives +ŷ for θ > 0), with spin tensors
normalised to Tr τ_kq τ_kq† = 2J_a + 1, τ_11 = −√(3/(2s(s+1))) (S_x + i S_y), s = J_a, and
T_kq = Tr(T τ_kq T†)/Tr(T T†) (traces over the final projections, the matrix in the initial m):

```
iT11 = √(3/(2s(s+1))) Tr(T S_y T†)/Tr(T T†)   (Tr(T S_x T†) = 0 by the reflection symmetry below),                      (19.8)
```

which is (√3/2)·A_y for a deuteron and A_y/√2 for spin ½. **Reflection in the reaction plane** (parity times a rotation
by π about ŷ leaves k_i, k_f in the xz plane unchanged; |J M⟩ → (−1)^{J−M}|J −M⟩ times the intrinsic parities, whose
product is (−1)^{l_a+l_F}):

```
T(−μ, −M_F; −m, −M_A) = (−1)^{l_a + l_F} (−1)^{(J_s−μ)+(J_F−M_F)+(J_a−m)+(J_A−M_A)} T(μ, M_F; m, M_A),                      (19.9)
```

a property of the result the code never uses, so a test of its phases.

**Cost.** Per pair of waves the kernel accumulates one more product (3 × 2 per l pair for d and n); ¹⁵N(d,n) with
spin-orbit in both channels, distorted waves, finite and zero range at 37 angles: 1.0 s at one thread (scratch driver
`work/dwba3/drivers/dnso.cpp`).

## 20. The ³He Trojan horse: ¹⁸O(³He,d)¹⁹F

Proton transfer from ³He = d + p to ¹⁸O at E_lab = 25 MeV (E_cm = 21.43 MeV; a typical ³He energy of THM runs with the
proton as participant). Inputs (the benchmark's, not claims about the nuclei):

| | |
|---|---|
| ³He = d + p | l = 0, j = ½, no node, B = 5.4935 MeV; Woods–Saxon R = 1.5 fm, a = 0.5 fm (depth searched: 57.5764 MeV), Coulomb sphere R_C = 1.5 fm; D0 = −152.53 MeV fm^3/2 (ours −152.5335, FRESCO −152.5), ANC 1.631716 fm^-1/2 |
| ¹⁹F = ¹⁸O + p | g.s. ½⁺ 2s½ (one node, B = 7.9937 MeV) and 5/2⁺ 1d5/2 at 0.1972 MeV (no node, B = 7.7965 MeV); WS r0 = 1.25, a = 0.65, spin-orbit 6.0 MeV (r0 1.25, a 0.65), R_C = 1.25 A^1/3 fm |
| ³He + ¹⁸O | V = 155 MeV (1.2, 0.72), W = 10 MeV (1.4, 0.88), spin-orbit 2.5 MeV (1.2, 0.72) when on, R_C = 1.3 A^1/3 fm |
| d + ¹⁹F | V = 88 MeV (1.17, 0.75), W = 2 MeV (1.32, 0.65), W_D = 12 MeV (1.32, 0.65), spin-orbit 6.0 MeV (1.07, 0.66) when on |
| Q | 2.5002 MeV (g.s.), 2.3030 MeV (5/2⁺); masses 3, 18, 2, 19 u; spectroscopic amplitudes 1 |

Radii r0 A_t^{1/3} unless stated; spin-orbit strengths in the convention of (2.1) (FRESCO's TYPE=3). The ³He spin ½ is
the entrance waves' spin, the deuteron's 1 that of the exit: the case exercises (19.5) with J_a = ½, J_s = 1, the
reverse of ¹⁵N(d,n). Results (FRESCO data `he3d_<gs|d52>_<finite|zero>_range[_so]`, test `transfer_he3_fresco`), the
worst relative difference from 0 to 90° where σ > 1 % of its maximum; the iT11 of ³He in absolute units:

| Case | finite range σ; iT11 | zero range σ; iT11 | FRESCO ref − refined (FR σ) |
|---|---|---|---|
| g.s., no spin-orbit | 3.5e-4 (40°); 1.0e-5 | 9.2e-5; 1e-16 (both ≈ 0) | 1.8e-3 |
| g.s., spin-orbit | 3.4e-4 (40°); 7.5e-5 | 6.4e-5; 2.2e-5 | 1.5e-3 |
| 5/2⁺, no spin-orbit | 8.9e-4 (80°); 2.2e-4 | 5.1e-5; 1.2e-5 | 7.9e-4 |
| 5/2⁺, spin-orbit | 7.2e-4 (50°); 4.3e-4 | 2.7e-5; 2.0e-5 | 6.8e-4 |

Asked 1 %: reached everywhere. Numerical findings:

- **The vertex range.** V_pd u_³He keeps the p + d Coulomb tail e²/ρ × u (κ = 0.417 fm⁻¹): our result moves by 1.2e-3
  at 0° and 2.1e-3 at 20° between ρ_max = 10 and 16 fm and by < 2e-4 beyond (16, 20 fm); the test uses 16 fm. At
  h = 0.025, h_K = 0.05 fm, ρ_max 18 fm, 64 x nodes, r_max 28 fm, l ≤ 39 ours moves by ≤ 2.1e-4. 0.5–2.0 s per
  case at one thread.
- **FRESCO's partial waves.** FRESCO's ERWIN loses 14 digits at L ≥ 35 without a radial cut, and one such wave
  spoils the transfer: at hcm = 0.02 fm the g.s. cross section came out ×16 (2745 against 174 mb/sr at 0°), and with
  jtmax = 45 instead of 36 by 3–7 %. With cutl = −3.2 (R_min = 0.08 L fm) jtmax = 36, 40, 45, 55 agree to 3e-6; all
  finite-range runs use cutl and jtmax = 36. cutl = −1.6 is not enough for the 5/2⁺ (its cross sections ×20 with
  and without spin-orbit).
- **FRESCO's convergence**: its refined run (hcm 0.0125, rintp 0.025, hnl 0.0125, cutl −6.4) moves the finite-range
  cross sections by up to 1.8e-3 (g.s., at 10–20°), above the 0.2 % asked; its larger radii (rmatch 35, rnl 12, nnu
  96) by ≤ 2.4e-4; cutl −1.6 against −3.2 by 0.9–1.4e-3 (g.s.). Ours is converged to 2e-4 and lies 0.04–0.09 % from FRESCO's
  reference, 0.13–0.22 % from its refined run. The remaining 0.1–0.2 % is FRESCO's step and cut dependence.
- **The p + d ANC**: ours 1.631716 fm^-1/2 at h = 0.05, 0.025, 0.0125 fm (depth 57.57637–57.57639 MeV); FRESCO's
  u/W plateau 1.63168 at 8, 10, 12 fm and at both its steps: 2.2e-5 apart (part 1's p + ¹⁵N agreed to 4e-7). Not
  resolved; the test allows 3e-5 for this one ANC; it moves σ by < 5e-5.
- iT11 of ³He: the reflection test of §22 and FRESCO agree on the sign and normalisation for spin ½ (iT11 = A_y/√2).

## 22. Validation of part 3

Tests: `transfer_spin_orbit` (and `transfer_fresco`, `transfer_continuum*` unchanged in their printed lines: their
outputs before and after part 3 were compared line by line; the cross sections of part 1's four cases agree with
the previous code to 8 × 10⁻¹⁵). Each test was run once with an injected error and failed (HANDOFF entry of
2026-11-12, part 3).

| Test | Compared with | Asked | Reached (worst, where) |
|---|---|---|---|
| Reduction (19.7): waves flagged with a spin-orbit term of 10⁻³⁰⁰ MeV against waves without | each other, FR and ZR, σ and iT11 | 1e-12 | 6.7e-16 (σ), 1.9e-16 (iT11) |
| Reflection (19.9), spin-orbit in both channels, θ = 10, 35, 80° | T itself; Tr(T S_x T†) = 0 | 1e-10 | 2.0e-16 of max \|T\|; 4.1e-17 |
| FRESCO ¹⁵N(d,n)¹⁶O, spin-orbit in d + ¹⁵N and n + ¹⁶O, finite range | dσ/dΩ, 0–90°, > 1 % of peak | 1 % | 2.0e-4 (75°); FRESCO ref − refined 1.0e-4 |
| same | iT11 of the deuteron | 1e-3 abs. | 1.0e-4 (at 0°: ours at 0.001°, FRESCO at 0°; elsewhere ≤ 8e-5); FRESCO 4e-5 |
| FRESCO ¹⁵N(d,n)¹⁶O, spin-orbit, zero range | dσ/dΩ; iT11 | 1 %; 1e-3 | 6.2e-6 (0°); 7.8e-5 (0°) |

FRESCO's convergence for `dn_finite_range_so` (data in its `reference.txt`): reference hcm 0.025, rintp 0.05, hnl
0.025, rnl 10, rmatch 30; refined hcm 0.02 (1.0e-4); rmatch 35, rnl 12, nnu 96 (≤ 4e-5 at 0, 30, 90°). At hcm =
0.0125 without a radial cut FRESCO's ERWIN loses 14 digits at L ≥ 28 and its 0° cross section moves by 3 %
(0.4778 against 0.4626 mb/sr); with cutl = −3.2 (hcm 0.025) or −6.4 (hcm 0.0125) it gives 0.46314 and 0.46364 mb/sr
at 0°, 0.1–0.2 % above its uncut runs and ours, while ours moves by 1.2e-5 between (h, h_K) = (0.05, 0.1) and
(0.025, 0.05) fm and by 3.5e-6 with ρ_max 12 fm, 64 x nodes, r_max 28 fm. The cut runs are data, not compared.
