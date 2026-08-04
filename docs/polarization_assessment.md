# Is the polarization implementation correct?

An honest assessment of the vector analyzing power in AZURE2: what was built,
what has actually been demonstrated, and what has not. Written to be useful when
deciding how far to trust a number, so the weak points are given the same space
as the strong ones.

Short version: **the machinery that builds the amplitudes is verified to
9–11 digits against a production code path, and the observable extracted from it
is verified against measurement for a spin-0 target. For a target that carries
spin, the extraction has no experimental validation at all — it rests on a
reduction argument.**

---

## 1. What was built, and why this way

### The amplitudes

Everything follows from the channel-spin amplitude matrix `M_{s'ν' sν}(θ)`, taken
from Seyler, *Nucl. Phys.* **A124** (1969) 253, Eq. (4). Seyler writes `M`
directly in terms of the Lane–Thomas collision matrix `U` — the object an
R-matrix code already builds — rather than in terms of phase shifts.

That choice is the reason the implementation is small. The decisive discovery was
that AZURE2's existing `tmatrix` (`src/AMatrixFunc.cpp:408`) *already is* Seyler's
bracket with its Coulomb phases applied:

```
tmatrix = exp{i(ω_l + ω_l')} · (δ_ss' δ_ll' − U^J_{s'l'sl})
```

So no phase is reconstructed and no convention is re-derived. More importantly,
the polarization code sees **the same `U`** the cross-section code sees, with
every boundary-condition and Brune-parametrization subtlety already in it. Had
`U` been rebuilt independently, the two paths could disagree without either being
obviously wrong.

The one genuinely new ingredient was angular: Seyler's expression needs
`Y_l'^{ν−ν'}(θ,0)`, and AZURE2 had only `P_L(cos θ)`. This is structural, not an
oversight — spin flip means `ν ≠ ν'`, hence `μ ≠ 0`, hence associated Legendre
functions. A code with only Legendre polynomials **cannot** produce a non-zero
analyzing power, because the interfering amplitudes are not representable in it.

### The observable

`A_y` is a density-matrix trace, `ρ_out = M ρ_in M†`, contracted with `σ_y`:

```
A_y = Tr(M σ_y M†) / Tr(M M†)
```

with `ŷ = k̂_in × k̂_out` (Madison convention). Written out:

```
A_y = 2 Σ Im[ M(+½) M*(−½) ] / Σ ( |M(+½)|² + |M(−½)|² )
```

**The subtlety that caused the 15N bug.** `σ_y` acts on the **projectile** spin
alone; the target spin is traced over. But `M` is stored in the *channel-spin*
basis, where projectile and target spin are already coupled. When the target is
spin-0 the two coincide — the channel spin *is* the projectile's — and the
formula can be read straight off the channel-spin amplitudes. When the target
carries spin they do not, and the entrance index must be decomposed first:

```
M_{out; m1 m2} = Σ_s ⟨j1 m1 j2 m2 | s, m1+m2⟩ · M_{out; s, m1+m2}
```

summed over the target projection `m2` alongside the exit configurations. For
p + ¹⁵N both particles are spin-½, so the channel spins are **0 and 1, never ½**.
The original code looked for channel spin ½, found none, and returned exactly
zero for a large, well-defined observable.

### Two consequences that are easy to get wrong

**Targets.** `A_y` is a ratio, so averaging it over a target thickness with equal
weight is meaningless. The measured quantity is the ratio of integrated yields,
`⟨A_y⟩ = ∫A_y σ dE / ∫σ dE`. This is implemented by running the *existing* yield
integrator twice — once on σ, once on σ·A_y — so quadrature, straggling and
energy loss are identical to the cross-section path by construction.

**Derivatives.** `M` is *linear* in `T`, so the coefficient the forward pass
multiplies `T` by is already the derivative. The adjoint walks the same loop over
the same `(k,m)` indices and contracts against the cotangents, which is why it is
short.

---

## 2. What has actually been checked

Graded by how much each one would catch.

### Strong

**(a) The amplitude matrix, via angle independence.** `M` must reproduce what
`GenMatrixFunc::CalculateCrossSection` produces by the Blatt–Biedenharn route.
The two differ by kinematic factors that are tedious to match, so the test used
is sharper and needs no bookkeeping: *at fixed energy, the ratio of the two must
be constant in angle*. Any error in coupling order, in the choice of `l` vs `l'`
inside the spherical harmonic, in the phase convention, or in the pathway
enumeration produces an angle-dependent ratio.

| system | target spin | angles | spread of the ratio |
|---|---|---|---|
| ¹²C + p | 0 | — | 1e-9 … 1e-11 |
| ¹⁵N + p | ½ | 33 | 1.3e-9 |

This is the single most informative check in the set. It tests the whole
amplitude construction at once against code that has been in production for
years — and it now covers a target with spin.

**(b) The observable and its sign, against measurement.** No internal check can
fix the overall sign: every bound and every zero is invariant under
`A_y → −A_y`. The sign depends on a chain of conventions (direction of `n̂`,
Condon–Shortley phase, Clebsch–Gordan ordering), and getting one backwards flips
it while disturbing nothing else. Baumann *et al.* (1992) mark four points where
`|A_y|` reaches unity; the ¹³N model reproduces all four in position **and in
sign pattern** (three positive, one negative). This is what pins the convention.

**(c) The analytic derivative.** Central differences on all 71 analyzing-power
rows × 14 parameters: twelve agree to 1e-9 or better. The two level energies sat
at 4e-6, which was chased rather than accepted — the residue falls as `h²`
(1.5e-5 → 1.5e-7 → 1.7e-9 for h = 1e-4, 1e-5, 1e-6) before hitting the roundoff
floor, so it is finite-difference truncation, not the adjoint. Capture widths
give identically zero, correctly. Re-verified on p + ¹⁵N: 4e-7 … 3e-11.

**(d) The generalization did not disturb the old path.** After rewriting for
arbitrary target spin, `tests/13N` returns a **bit-identical** chi-squared, so
the spin-0 case provably collapses to what was already validated in (b).

### Weak, or absent

**(e) There is no experimental validation of `A_y` for a target with spin.**
This is the significant gap. Correctness for p + ¹⁵N rests on three things, none
of which is a measurement:

- it reduces *exactly* to the spin-0 formula when `j2 = 0`, and that formula is
  validated by (b);
- the amplitudes it consumes are validated by (a), including for this system;
- `|A_y| ≤ 1` holds on real data.

That is a reduction argument plus a necessary condition. It is not the same as
having reproduced a measured curve.

**(f) `|A_y| ≤ 1` is necessary, not sufficient.** It would not catch a wrong but
bounded answer.

**(g) The target-thickness weighting is derived, not measured.** The formula is
not in doubt, but no comparison against a thick-target analyzing-power
measurement has been made.

**(h) No independent code cross-check.** Nothing has been compared against
another R-matrix or optical-model code.

**(i) The ¹⁵N comparison does not yet agree with data**, and this is unresolved.
On James's file the calculation gives ±0.15 where the Darden data is ~±0.30. The
most likely explanation is simply that his model has never been fitted to these
data. Weak evidence for that reading: the ratio of data to calculation is *not
constant* (5.9, 4.2, 3.4, 3.3, 3.1 over the first five points), so it does not
look like a missing scale factor, which is what a formula error would most
plausibly produce. But this is an argument from shape, not a resolution — the
honest statement is that the ¹⁵N case is **not yet confirmed against
experiment**.

---

## 3. Where it could still be wrong

**Inverse kinematics.** `A_y` is computed with respect to `GetJ(1)`, which AZURE2
defines as the **light** particle of the pair (`CNuc.cpp:419` prints "Light
Particle J"; the GUI labels the fields "Light Particle" / "Heavy Particle"). For
normal kinematics the light particle is the beam and this is right. For inverse
kinematics — heavy beam on light target — it is not, and the number returned
would be the analyzing power of the *target* rather than the beam. Nothing warns
about this.

**The Clebsch–Gordan ordering — resolved, and worth recording how.** The concern
was that if the ordering in the decomposition were transposed, the reduction to
`j2 = 0` would still hold (only one term survives), while a spin-carrying target
would pick up a wrong relative phase between the `s = 0` and `s = 1`
contributions. Check (a) cannot catch this: it tests `Σ|M|²`, which is blind to
how the entrance index is recombined afterwards. Three steps settled it.

*It matters, but modestly.* Computing both orderings side by side on p + ¹⁵N,
the difference is 1e-3 to 1e-2 in absolute `A_y` — around 5–7% where `A_y` is
largest, and 100% only where it is near zero anyway. On a spin-0 target the two
are identical to the last bit, as the reduction argument predicts. So this is a
real effect, but **it cannot explain a factor of 2–6**, and the ¹⁵N data
disagreement (i) is therefore a separate matter.

*Nothing inside AZURE2 fixes the order.* `CNuc.cpp:258` only enumerates the
allowed `s` values; every other `ClebGord` call in the engine involves orbital
momenta. The unpolarized cross section adds channel spins incoherently, so no
existing calculation is sensitive to it. The convention is one this code
introduces.

*Lane and Thomas fix it.* AZURE2's formalism is theirs, and they are explicit
(RMP **30** (1958) 257, sec. III.2a): *"this channel spin s is formed by coupling
I₁ and I₂ together: s = I₁ + I₂"*, with coefficients *"(I₁I₂i₁i₂|sν) … elements
of the matrix of the orthogonal transformation from the (I₁i₁, I₂i₂) scheme to
the (I₁I₂, sν) scheme … as discussed by Condon and Shortley."* Particle 1 first,
Condon–Shortley phases. `AngCoeff::ClebGord` was verified against the singlet and
triplet decompositions to be exactly `⟨j₁m₁j₂m₂|JM⟩` in that convention, and the
implementation calls it with particle 1 — which AZURE2 defines as the light
particle — first. **The implementation matches its own stated formalism.**

One residual caveat, which is a modelling matter rather than a code one: the
signs of the user's reduced width amplitudes for `(l, s)` channels must have been
taken in the same convention. Nothing can check that from inside the code.

**Identical particles are not handled, and this one is not gated.** For an
identical pair the amplitude must be symmetrised. The Coulomb side already is --
`EPoint::GetCoulombAmplitude` returns the Mott amplitude
`f_C(θ) + ε f_C(π−θ)` (`EPoint.cpp:1267`) and the polarization code consumes it
unchanged. The nuclear side is not: `PolarizationFunc.cpp` contains no reference
to `IsIdentical`, while the cross-section path multiplies the resonant term by 4
and the interference term by 2 (`GenMatrixFunc.cpp:285`) -- equivalent to
doubling the nuclear amplitude. So for an identical pair the nuclear amplitude
sits a factor of two low against the Coulomb one and `A_y` comes out wrong
rather than zero, which is worse, because it looks plausible. For α+α it is
harmless (`I₁ = 0`, so `A_y = 0` regardless); for polarized p+p elastic
scattering it is not. The fix appears to be one line, but is deliberately not
applied while untested — and the test that would settle it already exists, since
check (a) should currently *fail* for an identical pair and pass once the factor
is right.

**Not implemented, and gated rather than approximated:** tensor observables
(need spin-1 projectile and rank-2 operators); capture channels (need Seyler &
Weller, *PRC* **20** (1979) 453); and the analytic derivative of a
target-integrated `A_y`, which is a ratio of two integrals and falls back to
numerical differentiation for the whole fit rather than returning something
approximate.

**Unrelated but adjacent:** beta-delayed particle emission (`PType 20`) has no
analytic adjoint at all (`AMatrixFunc.cpp:525`), and because the fallback is
all-or-nothing, *one* such point disables the analytic Jacobian for an entire
fit. This is why the new fitting routine fails on James's ¹²C(α,γ) file — the
energy shifts are not the cause.

---

## 4. What would settle the open question

One measured analyzing power on a target that carries spin, reproduced in shape
and sign. Candidates:

- **Darden's ¹⁵N(p,p)**, already in James's file, once the model is fitted to it.
  This is the cheapest route and would directly close (e) and the sign concern.
- **Brune *et al.* (1998)**, ⁹Be(p⃗,d)⁸Be and ⁹Be(p⃗,α)⁶Li — ⁹Be is spin 3/2, so
  it tests a larger channel-spin decomposition, and the data were already
  R-matrix analyzed.

Until one of those is done, the honest position is: **the spin-0 target case is
verified end to end; the spin-carrying target case is verified up to the
amplitudes, its coupling convention is now pinned to the formalism the code is
built on, and what remains unconfirmed is the result itself against a
measurement.**
