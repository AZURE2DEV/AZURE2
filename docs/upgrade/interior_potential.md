# The interior R-function of a hybrid potential

**Status (2026-11-20): built on `feature/upgrade-6-interior-potential`** — steps 1 to 4 of §4 and the bound levels;
§8 lists what was built and where it differs from the design below. §6 (capture from the potential's interior):
steps 1 to 3 of §6.5 are built (§9); unbound final states, the rate mode's check and the fit parameters are not. The text of §1 to §7 is the design note as
reviewed ("no engine code is changed by it" was true of the note).

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

**Status (2026-11-20): design, asked for by the owner on that day; steps 1 to 3 of §6.5 were built the same day
(§9, which also lists the differences).** Code references are to `78587e70`. The formulas with levels (§6.2) were
derived here; the checks of §6.6 have since confirmed them (§9).

### 6.1 Why

With `interior_boundary` the broad states and the background come from the potential, but the capture still has only
the levels' internal photon widths and the integral **outside** a_c. The owner's ¹²C(α,γ) E1 test (Dyer and Barnes,
`12C+a_onefile/11-19-26_hpot_test/rpot`, 2026-11-20): χ² 24.9 without the potential; 2073 to 3805 with R_pot standing
in for the 9.585 MeV level, missed at every point (4 to 12 % of the points within 1σ). The potential makes the 1⁻
resonance in the elastic channel and nothing gives it a γ width. KTF and Langanke–Koonin take that width from the E1
matrix element of the potential's wave over the whole radial range, with an effective dipole strength.

### 6.2 The physics

R of Eq. (2) is a sum over the explicit levels and over the (infinitely many) eigenstates of the potential with the
boundary constant b_c. The interior wave function is the same sum, and the part carried by the potential's states is
closed: with ψ_c the channel's radial wave at the surface and u_c the regular solution of §2.1,

  χ_c(r) = [u_c(r) / u_c(a_c)] · R_pot,c · (a_c ψ′_c − b_c ψ_c),        r < a_c,        (7)

because Σ_n φ_n(r) φ_n(a)/(E_n − E) is the Green function of the potential with that boundary condition. The surface
value splits the same way: ψ_c(a) = χ_c(a) + (the levels' part). Outside, ψ_c = I_c δ_cc₀ − U_cc₀ O_c, and U is the
no-level term plus the level term of Eq. (3). The two terms of today's external capture
(`ECIntegral::CalculateDirect`, `src/ECIntegral.cpp:332`) each get an interior piece:

- **The non-resonant term.** Without levels ψ is the potential's own scattering wave, and (7) is its continuation:
  the standing wave v(r) = cos φ′ F(r) − sin φ′ G(r) of today's term (the primed phase is already used,
  `EntrancePhase`) continues inside as v(a) u_c(r)/u_c(a). The term becomes the integral from 0 to ∞ of the
  potential's wave: direct capture of the potential model.
- **The channel term of a level.** The level's outgoing wave O_c(r)/O_c(a) outside has the interior part

    [(L_c − b_c) / (f_c − b_c)] · u_c(r)/u_c(a_c),        (8)

  complex above threshold; below threshold L = S and it is real (the entrance factor of the closed-channel branch).
  The wave is **not** continuous at a_c: the remainder, (f − L)/(f − b), is the level's own interior, which is the
  level's internal photon width, a parameter as today.
- **The final state.** A bound final state is C W(r) outside. Its potential part inside is, by the same argument,

    w_f(r) = C W(a) · [(S_f − b_f) / (f_f − b_f)] · u_f(r)/u_f(a),        (9)

  the term whose square is already in the final state's norm (§8, the second term of dS′/dE). At a bound state of
  the potential (f_f = S_f) the factor is 1 and the final state is the potential's wave function, KTF's case.
  A final channel whose l has no `interior_boundary` has no interior wave in this model: the interior term is zero.

The interior radial integral of a pathway is then

  ∫₀^{a} w_f(r) · 𝒪_L(r) · [u_i(r)/u_i(a)] dr        (10)

with the same operator 𝒪_L as outside (r^L, or the charge function and the convection term without the
long-wavelength approximation), added to the FW/GW combination with the coefficient v(a) (non-resonant term) or
(8) × O_i(a) (channel term). Everything after the radial integral (angular factors, the ANC conversion, the level
matrix, Brune's transformation of the photon widths) is unchanged.

**What the model leaves out**, and says so: the matrix element between the *level* part of the final state and the
*potential* part of the entrance wave. (The other cross term, level part of the entrance wave with the whole final
state, is the level's internal photon width.) It vanishes when the final state is a state of the potential.

### 6.3 The dipole strength

For N = Z pairs the E1 operator of the relative motion nearly vanishes (¹²C + α: Z₁m₂/M − Z₂m₁/M ≈ 10⁻³), so the
kinematic operator inside gives the potential's resonance almost no γ width; the real strength is isospin mixing in
the interior. KTF fit an effective dipole strength. Proposed: one number per potential record and multipolarity,

```
    interior_charge  E1 0.012                  # new, optional; units of e
```

an **additional** effective charge q acting inside a_c only: the interior integrand gains q r^L (long-wavelength form;
k_γ r < 0.3 inside the radius, the neglected terms are below 1 %). Absent or 0: the pair's own operator, continuous
through a_c. The record is the entrance channel's (the one of l_i). It enters the amplitude linearly.

### 6.4 Design

- `physics/HybridPotential`: `RegularSolution` returns u on its Numerov grid (today: f and the norm only);
  `InteriorOverlap(u_i, u_f, weight)` does (10) on that grid (Simpson). u_f is made once per final state, u_i once per
  (channel, energy): the outward integration that f already needs.
- `ECIntegral::CalculateDirect`: the interior term in the open and the closed entrance branch, for a bound final
  state. The effective-charge function is evaluated at r < a (it is regular at 0). `ECTable` interpolates the sum as
  it does today's; the rate mode's fixed grid (`ECRateGridKey`) keeps the outside part and adds the interior number.
- `PotentialKey` gains q; the potential's interior and b are already in it, so every cache separates the cases.
- `--ignore-externals` switches the interior term off with the external one (both are computed capture).
- No change in `ExternalCaptureFinalStates.cpp` (its factors are those of (9)), none in `src/AMatrixFunc.cpp` or
  `src/WorkingCopy.cpp` (the size rule).

| code (proposed) | when |
|---|---|
| `W_POT_INTERIOR_FINAL` | capture from a channel with `interior_boundary` into a final channel without one: no interior capture for this final state |
| `E_POT_INTERIOR_UNBOUND_FINAL` | an unbound final state (`linewindow`) with `interior_boundary` in either channel, until step 4 |
| `E_POT_INTERIOR_CHARGE` | `interior_charge` without `interior_boundary`, or a multipolarity other than E1, E2 |

### 6.5 Steps

1. `physics/`: u on the grid, `InteriorOverlap`, unit tests against `physics/tests/make_interior_capture_reference.py`
   (scipy; not the engine's code).
2. The engine: the interior term for bound final states, both entrance branches; checks 1 to 4 of §6.6.
3. `interior_charge` (schema, reader, pyazr view, the diagnostics); check 5; the ¹²C(α,γ) trial of §6.6.
4. Unbound final states (the outgoing final wave with the factor (L_f − b_f)/(f_f − b_f)), the rate mode.
5. The theory page (`hybrid_potential.rst`), the user guide, `FOR_THE_OWNER.md`.

q as a **fit parameter** belongs with the potential's depth, radius and diffuseness (§7, question 5): step 3 reads a
fixed number.

### 6.6 Tests with independent answers

1. No levels, the KTF potential, a final state at a bound state of the l = 0 potential: the capture amplitude equals
   the potential-model integral ∫₀^∞ of an independent script, and is the same for a_c = 4.0, 5.35 and 6.5 fm.
2. The same with a final state off the potential's bound state: the factor of (9) against the script.
3. M1 with l_i = l_f in one potential, no levels: the wave functions are orthogonal, the amplitude over 0 to ∞ is
   zero (today's external part alone is not).
4. One level on top (the case of §2.5, check 3): the capture amplitude built from (7) directly, with ψ(a) and ψ′(a)
   from U, against the engine's two terms with (8). This is the check of (8).
5. `interior_charge`: the amplitude is linear in q; q = 0 reproduces step 2 bit for bit.
6. `interior_boundary` absent: every existing reference unchanged, bit for bit.
7. The trial: the owner's ¹²C(α,γ) E1 runs (R_pot + 7.12 MeV level + background pole) with q scanned; reported per
   point. This says whether the model can describe the data; it is not a test of the code.

### 6.7 Questions

1. **`interior_charge`**: as proposed (additive, inside a_c only, per record and multipolarity), or a factor on the
   whole potential capture, inside and outside (KTF's "effective dipole strength")? The first keeps the outside
   operator the one the masses fix.
2. **The cross term left out** (§6.2): accept, or a second parameter per final state? Proposed: accept; it is zero in
   KTF's case.
3. **The final state without a potential of its own l.** Zero interior capture and a warning (proposed), or refuse?
4. **Order.** Steps 1 to 3 make the ¹²C(α,γ) trial possible; q as a fit parameter then decides whether fits can use
   it. Do the fit parameters (potential and q) come before step 4?

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

## 8. As built (2026-11-20)

The owner's answers to §7, as relayed 2026-11-19: b is an explicit input with no default; a signed depth is
allowed; a pair with two spinning particles is refused; the width of the warning is the larger of the level's
observed width and the potential state's; the potential's parameters as fit parameters come later.

| design | built |
|---|---|
| `RegularLogDerivative`, `PrimedChannelFunctions` | `RegularSolution` (f and the interior norm), `InteriorChannelFunctions` (`physics/HybridPotential`); test `physics_interior_potential` against `physics/tests/make_interior_potential_reference.py` |
| record word `interior boundary <b>` | **`interior_boundary <b>`** (one word, a number); `l <l> ...` as designed |
| `E_POT_L_OVERLAP` | `E_POTENTIAL_L` (an l named twice, not an integer >= 0, or an l list with a spin-orbit term) |
| `E_POT_INTERIOR_THM`, `E_POT_INTERIOR_RELATIVISTIC` | the hybrid model's existing refusals apply unchanged; plus `E_POT_INTERIOR_TWO_SPINS` |
| `W_POT_LEVEL_DOUBLE` as a record diagnostic | an engine set-up warning (`WarnLevelDoublesPotential`, `src/CompoundSystemChecks.cpp`): it needs the channel functions. States are the zeros of f − S (bound states and resonances; above threshold only zeros narrower than half their energy); a pole of R_pot (f = b) is not a state and is not looked for |
| `W_POT_INTERIOR_CAPTURE` | not built: capture runs with the primed entrance functions and the exterior integral; §6 is still open |
| depth > 0 | signed (`E_OUT_OF_RANGE` for 0) |

**Bound levels** (not worked out in §2). Below threshold S′ = b + (S − b)(f − b)/(f − S). A level's amplitude at the
surface is γ/D, 1/D = (f − b)/(f − S), and

    C² = (2μa/ħ²W(a)²) γ² D⁻² / (1 + Σ γ² dS′/dE),
    (ħ²/2μa) dS′/dE = D⁻² [ ∫ₐ^∞ W²/W(a)² dr + ((S − b)/(f − b))² ∫₀ᵃ u²/u(a)² dr ].

`ParticleChannelValues` carries D⁻² in the ANC conversion factor; an external-capture final state
(`ECFinalStateFactors`) carries |1/D| in its conversion and the bracket times D⁻² as its norm integral (an unbound
final state: √(P′/P)). At a bound state of the potential (f = S) the limit of a large γ is the potential's
single-particle ANC. The sign of an ANC stays the sign of γ; the surface amplitude has the sign of γ/D, which matters
only for a level with two closed particle channels (not carried: a known limit).

**Checks.** `upgrade_interior_bound` (tests/upgrade/interior_bound, ¹²C + α, the potential of §2.5, a = 5.35 fm):
the factor to 2 × 10⁻¹⁶; the norm identity to 1.3 × 10⁻⁵ over l = 0 to 3 and bindings 0.05 to 30 MeV; the
single-particle ANC² at the potential's bound states (l = 0: 2.4276 MeV, l = 1: 9.6065 MeV, l = 2: 2.6152 MeV below
threshold) to 2.3 × 10⁻⁵; the phase without levels against the regular solution integrated to 14 fm, 2.3 × 10⁻⁸ rad
at 80 points; the potential's l = 1 resonance at 2.3988 MeV, 412 keV (f = S; the owner's scripts with R_C = 3.5 fm:
2.435 MeV, 479 keV by the phase), and the warning for the level placed on it.

**Still open.** Capture from the potential's interior (§6); the potential's parameters as fit parameters; the
physics suite of the owner's E1 fits (§5, last but one item); the signed surface amplitude above.

## 9. Capture from the interior, as built (2026-11-20)

The owner's "build it" (2026-11-20) without answers to §6.7: built with the proposals, each a small change to undo —
an additive interior charge inside a_c; the cross term of §6.2 left out; no interior capture and a warning for a
final channel without `interior_boundary` of its own l; steps 1 to 3 before step 4 and the fit parameters. They are
listed in `FOR_THE_OWNER.md` §8.

| design (§6) | built |
|---|---|
| `RegularSolution` returns u; `InteriorOverlap` | `RegularWaveOnGrid` (u/u(a) on the Numerov grid, f, the norm; `RegularSolution` is its two numbers, bit for bit) and `InteriorOverlap` (`physics/HybridPotential`) |
| the interior term in `ECIntegral::CalculateDirect`, both entrance branches, bound final states | `ECIntegral::Interior` and the three places of `CalculateDirect` (non-resonant, channel term above and below threshold); the final state's wave is kept per (l_f, J_f, binding) |
| `interior_charge E1 <q>` | **`interior_charge_e1 <q>`, `interior_charge_e2 <q>`** (two words, a number each); in `PotentialKey` after b, so every cache and the amplitude file's signature separate the cases |
| `E_POT_INTERIOR_CHARGE` | built (the word without `interior_boundary`); other multipolarities have no word |
| `W_POT_INTERIOR_FINAL`, `E_POT_INTERIOR_UNBOUND_FINAL` | one set-up warning for both cases (`WarnInteriorCaptureFinal`, `src/CompoundSystemChecks.cpp`); an unbound final state is **warned about, not refused**: it keeps the capture it had (√(P′/P) and the exterior integral) |
| `ECTable` interpolates the sum | **not interpolated**: a term whose entrance channel has `interior_boundary` is integrated at every energy (`ECTable::Integral`). The table divides by √P′ or e^{−iφ′}, which vanish or change sign at a pole of R_pot (f = b), and the terms carry the potential's resonances (1/(f − L)): measured 1 % off on the 412 keV wide 1⁻ resonance and of the order of the value near f = b. **This was also wrong at `78587e70`** for the exterior term of such a channel (its primed √P and phase), with the tables on, which is the default |
| the rate mode | untested (step 4) |

**Checks** (§6.6).

- 1, 2, 4 and the closed branch: `upgrade_interior_capture` against `tests/upgrade/interior_capture/reference.csv`
  (`make_reference.py`: scipy DOP853, mpmath). Each engine term with the interior over the same term with the final
  channel's potential without its boundary constant, 162 ratios at a = 4.0, 5.35, 6.5 fm, b = 0 and −1.5, final
  state on and off the potential's bound state: worst 2.7 × 10⁻⁵ at the default step, 9 × 10⁻⁷ at 0.01 fm. The
  script asserts the identity behind (8) (one level on R_pot, the potential's part of ψ(a) from U) to 10⁻¹⁰. The
  non-resonant term for the potential's bound state is the same at the three radii to 1.2 × 10⁻⁵ (exact form:
  1.7 × 10⁻⁵).
- 3: M1 between two l = 1 states of the potential, no level: the whole amplitude is 5.7 × 10⁻⁷ of its exterior part.
- 5: linear in q exactly; q equal to the pair's own E1 charge number doubles the interior part to 1.6 × 10⁻¹³.
- The whole engine (`interior_capture`, `tests/interior_capture`): the cross section of the potential model's
  direct capture at 51 energies across the 1⁻ resonance is the same at a = 4.0, 5.35 and 6.5 fm to 2.3 × 10⁻⁴ (the
  direct integrals' tolerance is 10⁻⁴); without the final state's interior the radii differ by a factor 7000.
  `interior_charge_e1 0.01` raises it by up to a factor 140: the pair's own E1 charge number is −9.7 × 10⁻⁴.
- 6 (`interior_boundary` absent: every reference bit for bit) and 7 (the ¹²C(α,γ) trial): the gate and the owner's
  runs, recorded where they are made.

**For the ¹²C(α,γ) trial.** The interior term needs a final state with an ANC in the α + ¹²C channel (the ground
state as a 0⁺ level, as external capture always did) and a potential record with `interior_boundary` that covers
l = 0 as well as l = 1; a project without that final state has no computed capture at all.

**Still open.** Unbound final states and the rate mode (step 4); q and the potential's parameters as fit
parameters; the cross term of §6.2; the theory page's worked example.

