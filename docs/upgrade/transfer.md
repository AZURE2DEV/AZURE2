# physics/transfer: optical potentials, bound states, distorted waves and DWBA transfer (part 1)

The Trojan-horse vertices of `THM.md` §3 need distorted waves and a DWBA transfer amplitude that are checked against an
independent reaction code. This note describes the leaf library `physics/transfer` (namespace
`azure2::physics::transfer`, CMake target `azure2_physics_transfer`): an optical potential as a value, the bound state
of two clusters, distorted waves with spin 0, ½ and 1, the transfer A(a, s)F of a cluster x (a = s + x, F = A + x) to
a bound state of F in zero range and in finite range (post form), and the plane-wave limit of the THM entrance
vertex. Every part is compared with a closed form or with FRESCO 3.5 (§9). Nothing of the engine uses the library
yet; the transfer to the x + A continuum and the wiring are the next parts (§11, §12).

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
| `.../Transfer.h`, `src/transfer/Transfer.cpp` | Eqs. (5.1)–(6.7) |
| `.../ThmVertex.h`, `src/transfer/ThmVertex.cpp` | Eqs. (8.1)–(8.3) |
| `physics/tests/transfer_*_test.cpp`, `physics/tests/data/fresco/`, `physics/tests/data/make_fresco_benchmarks.py` | §9 |

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
- Not implemented in this part: spin-orbit in the distorted waves of a transfer (refused), the remnant term, prior
  form, transfer to unbound states (§11), global optical potentials, relativistic kinematics, target spin in the
  distortion (spin-spin), tensor forces, l > 40.

## 11. The continuum part (next)

THM needs the transfer to x + A at E > 0 inside the channel radius: in (6.2) the bound overlap φ_{l_F}(r) of F is
replaced by the x + A scattering function at E (or its interior part, the level function of §8), and the cross
section becomes the triple-differential one of `THM.md` §1. In the code: `Overlap` gains a continuum kind (u at E with
its normalisation), the kernel (6.4) is unchanged (it only needs u_F on a grid), and the vertex of §8 becomes the
projection of the kernel on the x + A channel inside a. FRESCO benchmarks there: continuum bins (`&OVERLAP` with BE < 0,
ISC = 2 or 4, ER, NK) for transfer to x + A, which are normalised averages of scattering states over an energy bin; our
amplitude with the same bin average must reproduce FRESCO's bin cross sections (first for a non-resonant bin, then a
resonant one), and the narrow-bin limit gives the transfer to a single energy. Also the prior form with remnant (to
compare post and prior, which agree when both remnants are kept).

## 12. Wiring (later branch)

The engine's `ThmDistortion` (zero-range prior M(E)) and `ThmDwVertex` (prior surface vertex) become callers of this
library: M(E) is (6.7) with c_i = β, c_f = 1 (§13); the optical potentials of `ThmOptical` become `OpticalPotential`
values filled by a global-set table (Koning–Delaroche, An–Cai, GDP08, Avrigeanu: one table, one row each, a test per
row); the vertex tables are built once per fit stage from `EntranceVertex`. That needs Phase 5's σ_HOES observable; it
is not done here.

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
