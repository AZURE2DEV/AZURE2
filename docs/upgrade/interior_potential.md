# The interior R-function of a hybrid potential

A design note for review; no engine code is changed by it. Asked for by the owner (R. J. deBoer, 2026-11-19) after a
test of the hybrid model on ¹²C + α. Code references are to `3a96c0a1`. The numbers of §2.5 were run with
`docs/upgrade/interior_potential_check.py` (scipy 1.18 DOP853, mpmath Coulomb functions; not the engine's code);
everything marked *unchecked* is believed, not verified.

## 1. Why

Today the potential of a hybrid pair acts only **outside** the channel radius a_c: `CoulFunc::computeHybrid`
(`src/CoulFunc.cpp:267`) replaces the Coulomb functions by the solutions of the radial equation with the potential,
integrated inward from r_m, and P, S, φ and the external-capture waves follow (`spin_orbit_potential.md` §1.1, §1.3).
Inside a_c the model is still a sum of poles.

The owner's test (¹²C + α, the l = 1 Woods–Saxon of Koonin, Tombrello and Fox, V₀ = 45.891 MeV, R = 3.55 fm,
a = 0.60 fm, at a_c = 5.35 fm; old engine `b60aed09`) found that this changes almost only a factor in P (1.14 to
1.21, constant to 5.6 % over 0.1 to 5 MeV) and a near-constant in S (+0.56 to +0.86). The reduced widths absorb the
first and Brune's S(E) − S(E_λ) the second: after a refit the E1 capture fit is the same to better than 1 % over the
data. So the present model cannot do what the hybrid papers use a potential for, which is to produce the broad
resonance and the smooth background from the potential **instead of background poles**:

- S. E. Koonin, T. A. Tombrello and G. Fox, Nucl. Phys. A 220 (1974) 221 (KTF);
- K. Langanke and S. E. Koonin, Nucl. Phys. A 410 (1983) 334; A 439 (1985) 384;
- J.-M. Sparenberg, Phys. Rev. C 69 (2004) 034601, Eqs. (1) to (4).

The owner's decisions for this note (2026-11-19): built on this branch's engine only; potentials per (pair, l);
a level that the potential already produces is **warned about**, not refused; capture from the potential's interior is
a **second step** (§6).

## 2. The physics

### 2.1 The R-function of the potential

For a channel c of a hybrid pair let u(r) be the **regular** solution of the radial equation (1) of
`spin_orbit_potential.md` (Coulomb sphere + V_N + centrifugal + spin–orbit term), u(0) = 0, and

  f_c(E) = a_c u′(a_c) / u(a_c)        (real, any E, above or below threshold).        (1)

The R-function of the potential with a boundary constant b_c is R_pot,c(E) = 1 / (f_c(E) − b_c); it is diagonal in the
channels. The model is (Sparenberg Eq. (1))

  R_cc′(E) = R_pot,c(E) δ_cc′ + Σ_λ γ_λc γ_λc′ / (E_λ − E).        (2)

### 2.2 What changes in the engine: three channel functions, nothing else

Write L_c = S_c + i P_c for the logarithmic derivative of the outgoing solution at a_c (today's hybrid S and P) and
φ_c for today's phase. With (2) and R_pot diagonal, the collision matrix keeps **exactly** today's form

  U_cc′ = Ω′_c Ω′_c′ [δ_cc′ + 2i P′_c^½ P′_c′^½ (γᵀ A γ)_cc′],  A⁻¹_λμ = (E_λ − E) δ_λμ − Σ_c γ_λc (S′_c − B_c + i P′_c) γ_μc,   (3)

with the primed functions

  L′_c ≡ S′_c + i P′_c = b_c + (L_c − b_c)(f_c − b_c) / (f_c − L_c),        (4)
  φ′_c = φ_c + arg[(f_c − L_c) / (f_c − b_c)]        (arg of a negative real number is π).        (5)

Derivation: 1 − R L⁰ = (1 − R_res L⁰′)(1 − R_pot L⁰) with L⁰ = L − b and L⁰′ = L⁰ / (1 − R_pot L⁰); the diagonal
factor D = 1 − R_pot L⁰ = (f − L)/(f − b) gives the phase of (5) and P′ = P / |D|². This is Sparenberg's Eq. (3) for
any number of channels. Consequences:

- **The level matrix and the collision matrix are not touched.** `src/AMatrixFunc.cpp` and the physics solvers read
  P, S and φ per channel; they get the primed ones. The change is where those three numbers are made
  (`CoulFunc`, `ChannelFunc`, `ChannelTable::Compute`) and in `physics/HybridPotential` (§3.1).
- **Without levels U_c = Ω′_c² = e^{−2iφ} (f − L*)/(f − L)** is the scattering of the potential alone, for any a_c and
  any b_c (§2.5, check 1).
- f → ∞ (no interior potential, a hard core) gives L′ = L and φ′ = φ: today's model.
- At a pole of R_pot (f = b) P′ = 0 and D⁻¹ changes sign: the factor P′^½ Ω′ of (3) passes through zero smoothly; the
  π of (5) is that sign. An implementation that carries the real factor (f − b)/|f − L| instead of |·| and a phase
  avoids the branch.
- Below threshold P = 0, L = S is real and (4) gives a real S′, with a pole where f = S: that is a **bound state of the
  potential**. The decaying solution and S(E < 0) of a hybrid pair are already the potential's
  (`physics/` `HybridBoundShift`, `DecayingSolution`), so S′ is continuous through threshold where S is. *Unchecked*:
  the ANC ↔ γ conversion with (4) (it needs dS′/dE and the level's own normalisation; §4, step 4).

### 2.3 The boundary constant b_c, and Brune's parameters

The split (2) depends on b_c: under a change b → b + β the pole part does **not** stay a sum of poles with constant
widths (in one channel γ² → γ² / (1 + β R_pot(E))², energy dependent). So b_c is part of the model, an input of the
potential record, not something a transformation removes. KTF use b = 0.

The level matrix's own boundary condition B_c is a different thing and keeps all its properties: (3) has the standard
form with the real function S′_c(E) in place of S_c(E) and P′_c in place of P_c, and the Barker and Brune
transformations are algebraic identities of that form for any real S′(E). So the alternative (Brune) parametrisation,
Park's, and `--no-brune` with B_c at the lowest level all work with S′, P′ and dS′/dE, and "observed" energies and
widths keep their definitions (Sparenberg's Eq. (4) is the one-level case). *Unchecked* beyond the one-level
identity of §2.5, check 3: no multi-level run was made.

### 2.4 The energy derivative

dS′/dE needs df/dE, which costs nothing extra:

  df/dE = −(2μ a_c / ħ²) ∫₀^{a_c} u² dr / u(a_c)²        (6)

(§2.5, check 4). The engine takes dS/dE by a Richardson difference today (`ChannelFunctions.h`); either route works,
and (6) is the test of the other.

### 2.5 Checks run for this note

`python3 docs/upgrade/interior_potential_check.py` (¹²C + α, l = 1, the KTF potential with R_C = 3.55 fm, r_m = 15 fm):

| # | check | result |
|---|---|---|
| 1 | no levels: the phase of U_pot of §2.2 against the phase shift of the potential matched to Coulomb functions at 15 fm, E = 0.5 to 5 MeV, a_c = 4.0, 5.35 and 6.5 fm | worst difference 3.3 × 10⁻¹¹ degrees; \|U_pot\| = 1 to 2 × 10⁻¹⁶ |
| 2 | the pole of R_pot at a_c = 5.35 fm, b = 0 | E = 1.9245 MeV, γ² = 0.6937 MeV, θ² = 0.951 in units of 3ħ²/(2μa²); KTF Fig. 1: 1.926 MeV, θ² = 0.94 |
| 3 | one level (E_λ = 3.2 MeV, γ² = 0.15 MeV) plus R_pot: U from (2) directly against (3) with (4), (5), E = 1 to 4 MeV | agree to 5 × 10⁻¹⁶ |
| 4 | Eq. (6) against a central difference of f at 2 MeV | −1.479213 MeV⁻¹ both |

The owner's independent scripts (his fit archive, directory `12C+a_onefile/11-19-26_hpot_test/partA`) give the
potential's l = 1 resonance at E_cm = 2.435 MeV, Γ = 479 keV with R_C = 3.5 fm; they are the reference for the engine
tests of §5.

## 3. Design

### 3.1 `physics/`

`physics/HybridPotential` gains the outward integration: `RegularLogDerivative(p, z1, z2, mu, l, two_l_dot_i, a, E)`
returning f and ∫₀^a u² dr / u(a)² (for (6)); Numerov on the same grid as `IntegrateInward`, started from u ∝ r^{l+1}
(outward is the stable direction for the regular solution). A free function, no state. `PrimedChannelFunctions(L, φ,
f, b)` returns S′, P′ and the amplitude factor of (5) (the engine's e^{−iφ} with the sign of f − b), with the limits
of §2.2 as unit tests.

### 3.2 The engine

- `CoulFunc`/`ChannelFunc`: for a channel whose potential has `interior` on, P, S, e^{−iφ}, dS/dE and S(E < 0) are the
  primed ones. Everything downstream (the channel table, the level-parameter transformations, the adjoint) calls
  these and needs no change. *Unchecked*: that no caller rebuilds P or φ from F and G itself. `ECIntegral` uses
  F and G at the quadrature radii **and** √P, e^{−iφ} at a_c (`spin_orbit_potential.md` §1.3); the outgoing wave
  O_c(r) for r > a_c is unchanged by R_pot, only its coefficients are, so the expectation is that it needs the primed
  √P and phase and the unprimed F, G. This is the first thing to verify (§4, step 3).
- The memo and table keys already compare the potential by value (`PotentialKey`); `interior` and b join the key.
- No added line in `src/AMatrixFunc.cpp` or `src/WorkingCopy.cpp` is expected (the size rule).

### 3.3 The record (format 3)

```
potential a12C_p
    pair            12C+a
    l               1                          # new, optional: the orbital angular momenta this record applies to
    coulomb_radius  3.55 fm
    central         woods_saxon depth 45.891 MeV radius 3.55 fm diffuseness 0.60 fm
    interior        boundary 0                 # new, optional: R_pot with this b_c; absent = today's model
```

- `l` (a list): a pair may have several records with disjoint l; a record without `l` applies to every l that no
  other record of the pair names (today's meaning, so present projects keep their bits). KTF and Langanke–Koonin
  potentials are l-dependent, and with R_pot the l-dependence matters.
- `interior boundary <b>`: b is required when `interior` is given; no default (it changes numbers).
- The old `<potential>` block has neither word; the converter writes neither.
- pyazr and the web page get both words from the schema rows.

### 3.4 Diagnostics

| code (proposed) | when |
|---|---|
| `W_POT_LEVEL_DOUBLE` | a level with a non-zero width in a channel with `interior` lies within one width of a pole or a resonance of that channel's potential: the potential already produces it (the owner: warn, do not refuse) |
| `E_POT_L_OVERLAP` | two records of a pair name the same l |
| `E_POT_INTERIOR_THM`, `E_POT_INTERIOR_RELATIVISTIC` | as the hybrid model's existing refusals |
| `W_POT_INTERIOR_CAPTURE` | capture data with an `interior` entrance channel, until §6 is built: the potential's interior does not contribute to the capture |

"Within one width" needs the potential's resonances; they come from f(E) = b (poles of R_pot) and from the phase of
U_pot, computed once per record at set-up over the project's energy range.

## 4. Steps

1. `physics/`: `RegularLogDerivative`, `PrimedChannelFunctions`, unit tests with the table of §2.5 as reference data
   (regenerated by the script, stored as CSV). No engine change; every existing output keeps its bits.
2. The record words `l` and `interior` (schema, reader, converter untouched, diagnostics of §3.4); `l` alone changes
   nothing for records without it.
3. The engine: primed functions for `interior` channels; the external-capture audit of §3.2.
4. Level parameters: Brune, Park and formal with S′; ANC ↔ γ below threshold; the tests of §5.
5. The user guide (`particle_pairs.rst`), the conventions page (Eq. (3) to (5)), `FOR_THE_OWNER.md` §3 row (numbers
   move only for projects that write `interior`).

## 5. Tests with independent answers

- No levels, KTF potential: the elastic phase shift equals the potential's for a_c = 4, 5.35, 6.5 fm (check 1), and the
  resonance is at 2.435 MeV with Γ = 479 keV (R_C = 3.5 fm).
- The pole of R_pot of check 2.
- One level: check 3 through `AZURE2 calculate`, in the three parametrisations.
- A potential plus a narrow level against the same physics entered as poles only, where a pole expansion of R_pot over
  the fitted range is accurate (a consistency test of the split, not of the physics).
- The owner's E1 regression fits (`partB_*` of the same directory) as a physics suite once step 3 is in.
- `interior` absent: every existing reference unchanged, bit for bit.

## 6. Second step: capture from the potential's interior

KTF add a capture R-function from the potential: the direct E1 matrix element of the potential's scattering wave over
the **whole** radial range, scaled by an effective dipole strength. The counterpart here is the external-capture
integral run from 0 instead of a_c with the regular solution inside, and an effective-charge factor as a parameter.
Not designed in this note.

## 7. Questions

1. **b_c.** An explicit input with no default (as proposed), or a rule (b = −l; b = S′ at a chosen energy)?
2. **Attractive only?** V₀ > 0 is required today. Sparenberg's supersymmetric partner potentials are repulsive at
   short range and cannot be entered. Allow a signed depth, or a tabulated V(r), in the same change?
3. **Two spinning particles.** `spin_orbit_potential.md` refuses a spin–orbit term there (`E_SO_TWO_SPINS`) because the
   term couples the channel spins of one (l, J). R_pot would then not be diagonal in the channels and (4), (5) become
   matrix equations. Proposed: `interior` follows the same refusal.
4. **The width in `W_POT_LEVEL_DOUBLE`.** The level's formal width, its observed width, or the potential resonance's?
5. **Fit parameters.** The potential's depth, radius and diffuseness are not fit parameters today
   (`spin_orbit_potential.md` §1.4 item 3). With R_pot standing in for background poles the depth is what one wants
   to fit. In this change, or after?
