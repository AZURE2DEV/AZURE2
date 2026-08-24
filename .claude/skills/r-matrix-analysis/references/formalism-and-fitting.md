# Formalism and fitting — the equations behind the picture

Companion to `SKILL.md`. This file collects the working equations and points at the
exact source equations. Notation follows Lane & Thomas (LT) [Rev. Mod. Phys. 30,
257 (1958)] and is kept by Brune and deBoer. Symbols: `λ,µ` levels; `c = (α,s,l)`
channels; `E` energy; `a_c` channel radius; `B_c` boundary constant.

## 1. The R-matrix and the collision matrix

Reduced-width amplitude `γ_λc` = amplitude of level `λ`'s eigenfunction at the
surface of channel `c`. The R-matrix (LT; deBoer Eq. 17):

```
R_{c'c}(E) = Σ_λ  γ_{λc'} γ_{λc} / (E_λ − E)
```

R encodes the logarithmic derivative of the internal radial wavefunction at
`r = a_c`. The collision (scattering) matrix `U` comes from matching to the external
Coulomb functions; in level-matrix form (Brune 2005 Eqs. 5–7):

```
[A^{-1}]_{λµ} = (E_λ − E) δ_{λµ} − Σ_c γ_{λc} γ_{µc} (S_c + i P_c − B_c)
U_{c'c}       = Ω_{c'} Ω_c [ δ_{c'c} + 2 i (P_{c'} P_c)^{1/2} γ_{c'}^T A γ_c ]
```

`A` is the level matrix; `Ω_c` are hard-sphere phase factors from Coulomb functions.
Cross sections and angular distributions follow from `U`. `U` is unitary and
symmetric (flux conservation + time reversal).

## 2. Penetration, shift, phase — the energy-dependent surface functions

Built from regular/irregular Coulomb functions `F_l, G_l` evaluated at `ρ = k a_c`:

- **Penetration** `P_c(E)` — barrier transmission; controls widths. Rises steeply
  with `E`; near a threshold `P_c → 0`, throttling widths (source of the ghost
  anomaly).
- **Shift** `S_c(E)` — real part of the logarithmic derivative; shifts a resonance
  from `E_λ`.
- **Hard-sphere phase** `φ_c(E)` — smooth non-resonant (potential) scattering phase.

These are geometry (via `a_c`) plus Coulomb (`η`) — not fit parameters. AZURE2
computes them.

## 3. Widths: formal vs observed (the ~30% trap)

**Formal** partial and total widths (Brune 2005 Eqs. 9–11):

```
Γ_λc = 2 γ²_λc P_c(E)        Γ_λ = Σ_c Γ_λc
```

**Level shift** `Δ_λ = − Σ_c γ²_λc [S_c(E) − B_c]`. Choosing `B_c = S_c(E_λ)` makes
`Δ_λ(E_λ) = 0`, so `E_λ` is the resonance energy — but only for *one* level.

**Observed** (physical) partial width, from matching the single-level R-matrix to
Breit-Wigner (Brune 2005 Eqs. 15–16):

```
E_R = E_λ
Γ°_λc = Γ_λc(E_λ) / [ 1 + Σ_c γ²_λc (dS_c/dE)|_{E_λ} ]
      = 2 γ²_λc P_c(E_λ) / [ 1 + Σ_c γ²_λc (dS_c/dE)|_{E_λ} ]
```

The denominator `[1 + Σ_c γ²_λc dS_c/dE]` is why **formal ≠ observed**, often by
~30% for broad levels. When importing literature ("observed") widths into a
non-Brune fit, invert this; in Brune mode it is handled automatically.

## 4. Boundary-condition gauge freedom

`U` is **invariant** under a change `B_c → B_c'` if `E_λ, γ_λc` are transformed
consistently (Barker 1972; Brune 2005 §VI). Construct the real symmetric matrix and
diagonalize:

```
C = e − Σ_c γ_c γ_c^T (B_c' − B_c) = K^T D K   →   E_λ' = D_λ,   γ_c' = K γ_c
```

(`e = diag(E_λ)`.) Consequences: (1) fits are equivalent for any `B_c`; (2) you may
choose `B_c` freely for convenience. Analogy: an electromagnetic gauge
transformation leaving fields invariant. **Implication:** raw `E_λ, γ_λc` are not
physical by themselves — they carry a `B_c` convention.

## 5. Observed parameters directly — Brune's alternative parametrization

Solve the nonlinear eigenvalue problem (Brune 2005 §VII, Eqs. 22–25):

```
E a_i = Ẽ_i a_i ,   E = e − Σ_c γ_c γ_c^T [ S_c(E) − B_c ]
γ̃_ic = a_i^T γ_c   ("on-resonance" reduced-width amplitudes)
Γ°_ic = 2 P_c(Ẽ_i) γ̃²_ic / [ 1 + Σ_c γ̃²_ic (dS_c/dE)|_{Ẽ_i} ]
```

The eigenvalues `Ẽ_i` **are** the physical resonance energies; `γ̃_ic` map simply to
observed widths. One can compute `U` directly from `(Ẽ_i, γ̃_ic)` via an alternative
level matrix `Ã` (Brune 2005 Eqs. 34–35) in which **`B_c` never appears**:

```
[Ã^{-1}]_{ij} = (Ẽ_i − E) δ_{ij} − Σ_c γ̃_ic γ̃_jc (S_c + i P_c)  +  (level-shift term, Eq. 34)
U_{c'c}       = Ω_{c'} Ω_c [ δ_{c'c} + 2 i (P_{c'} P_c)^{1/2} γ̃_{c'}^T Ã γ̃_c ]
```

Full boundary-condition independence; channel radii still required. **This is the
parametrization AZURE2 uses** — enter and read *observed* energies/widths. To go
back to standard `E_λ, γ_λc`, solve the generalized eigenvalue problem
`N b_λ = E_λ M b_λ` (Brune 2005 Eqs. 27–33). Full derivation: Brune, Phys. Rev. C
66, 044611 (2002).

## 6. Single-level cross section (Breit-Wigner limit)

For `N_λ = 1` (Brune 2005 Eq. 8):

```
σ_{cc'}(E) = (π ω / k_c²) · Γ_λc Γ_λc' / [ (E_λ + Δ_λ − E)² + (Γ_λ/2)² ]
ω = (2J+1) / [ (2J_1+1)(2J_2+1) ]   (statistical / spin factor)
```

With `Δ_λ ≈ (E_λ − E) Σ_c γ²_λc dS_c/dE` near the pole. Multi-level: keep the full
`A`/`R` — the poles interfere (next section).

## 7. Interference

Two levels of the **same `Jπ`** add as amplitudes in `R`/`A`; the cross section
carries a cross term. Schematically the interference contribution scales as
(deBoer Eq. 85):

```
σ_interference ∝ 2 √(σ₁ σ₂) · (relative sign of γ amplitudes)
```

- **Same `Jπ`** → energy-dependent interference (valleys, tails). Governs the
  low-energy extrapolation. The relative **sign** of the `γ_λc` is a discrete fit
  choice; test both.
- **Different `Jπ`** → no interference in the total cross section, but interference
  in the **angular distribution** (odd Legendre terms) — used to separate multipoles.

Signs are generally not predictable from theory; capture data in the off-resonance
regions determine them.

## 8. Capture (photon) channels

Photon channels `p = (ε, L, λ_f)` (ε = 1 electric / 0 magnetic; `L` multipolarity;
`λ_f` final state) are treated by **perturbation theory** (deBoer Sec. IV; Eqs.
56–62):

- Transition matrix = **internal** part (photon reduced width `γ_λp`, a real
  parameter) + **external** part (depends only on nuclear parameters and the
  **final-state ANC**), via the Siegert long-wavelength form.
- External capture splits into **hard-sphere** (non-resonant, entrance channel only)
  and **channel** (resonant) capture.
- Photon channels **do not enter** the level matrix `A` or the total width in the
  Breit-Wigner denominator: particle widths set the width, γ widths set the height.
- Effective charge `ē_α^L` (Eq. 62) can nearly cancel (e.g. `12C+α`), suppressing
  external capture.

## 9. Angular distributions

Differential cross section as a Legendre series (deBoer Eqs. 33–36):

```
dσ/dΩ = Σ_L B_L P_L(cos θ)
```

with `B_L` built from products of transition-matrix elements and Racah/Z angular-
momentum coefficients. Same-`Jπ` amplitudes contribute even-`L`; different-`Jπ`
(e.g. E1×E2) contribute odd-`L` → forward-backward asymmetry. Example (ground-state
E1/E2 capture, deBoer Eqs. 91–92): terms `P_1…P_4` with attenuation factors `Q_n`
and an E1–E2 phase `φ`, `cos φ = cos[δ_{α1} − δ_{α2} − tan^{-1}(η/2)]` (Watson's
theorem). Measuring at 90° isolates a single multipole.

## 10. Fitting machinery (deBoer Sec. VI–VII)

Objective (deBoer Eq. 87), per-dataset normalization `n_i` and a systematic term:

```
χ² = Σ_i Σ_k [ (n_i y_ik^theory − y_ik^data) / (n_i δy_ik) ]²  +  Σ_i [ (n_i − 1)/δn_i ]²
```

- The systematic normalization term addresses **Peelle's Pertinent Puzzle** (fits
  biased low when common normalization uncertainty is mishandled).
- Optionally reduced-`χ²` weighting so a single large dataset does not dominate;
  constraint terms for well-known (e.g. subthreshold) parameters; robust
  **Sivia–Skilling** estimator to down-weight outliers (Eq. 95).
- Uncertainties by Monte Carlo, but trustworthy only if reduced-`χ²` ≈ 1.
- Yields → cross sections by an (approximate) convolution/deconvolution for
  target/beam energy spread.
- **Uniqueness:** particle parameters are constrained by unitarity; γ-width
  *magnitudes* are constrained but multipolarity assignments are tentative;
  background-pole values and some interference/cascade signs are **not** unique.
- Adopt the **most physical** fit, not merely the lowest `χ²`.

Bayesian / MCMC alternatives (BRICK on top of AZURE2): Odell, Brune & Phillips
et al. — see `key-papers.md`.

## 11. Quantity ↔ AZURE2 input cheat-sheet

| Concept | Symbol | AZURE2 `.azr` location |
| --- | --- | --- |
| Partition / pair | `α` | particle pairs (spins, parities, Q, charges) |
| Channel radius | `a_c` | particle-pair radius |
| Level (state or background pole) | `λ`, `E_λ`→`Ẽ_λ` | levels (Jπ, E_x) |
| Reduced-width amplitude | `γ_λc`→`γ̃_λc` | channel entry within a level (sign = interference) |
| Partial width (observed) | `Γ°_λc` | derived from `γ̃_λc` (Brune mode) |
| Photon reduced width | `γ_λp` | photon channel per multipole (E1/E2/M1) |
| ANC (bound/subthreshold) | `C_λc` | channel entry for the closed channel |
| Boundary condition | `B_c` | not exposed — AZURE2 uses Brune parametrization |
| Data (σ, dσ/dΩ, analyzing power) | — | segments (with normalizations) |

Mechanics of editing these live in the `azure2-eval` skill; the *meaning* of each
lives in `SKILL.md`.
