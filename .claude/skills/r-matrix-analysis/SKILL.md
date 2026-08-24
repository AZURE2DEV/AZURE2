---
name: r-matrix-analysis
description: The "global picture" of phenomenological R-matrix theory for guiding resonance selection, understanding interference, and reasoning about R-matrix fits (e.g. in AZURE2). Explains — intuition first, then formalism — what channels/levels/boundary conditions are, how physical resonances map to R-matrix levels (energies, partial widths, Jπ, ANCs, subthreshold/bound states), how to decide which levels to include, how levels of the same Jπ interfere (constructive/destructive, ghost anomaly), how capture (γ) vs particle channels are treated, background poles, boundary-condition energy, and Brune's parametrization. Use whenever the task needs the theory/methodology behind an R-matrix evaluation — choosing levels, diagnosing an interference pattern, understanding what a fit parameter means, or connecting known level schemes to AZURE2 inputs. The companion `azure2-eval` skill runs the fits; this skill explains what the fit means and why.
---

# The global picture of R-matrix analysis

This skill is the **theory and methodology counterpart** to `azure2-eval` (which
runs the fits) and `nds-explorer` (which fetches the data and level schemes). It
does not run code. Its job is to let you reason about an R-matrix evaluation: to
choose which resonances belong in a fit, to read an interference pattern, to know
what each fit parameter physically means, and to translate a known level scheme
into sensible R-matrix inputs. Read this before deciding *what* to fit; use
`azure2-eval` to actually fit it.

The presentation is intuition-first. Formulas appear only after the picture is in
place, and the detailed formalism, equation-by-equation, lives in
`references/formalism-and-fitting.md`. Annotated pointers to the primary
literature (Lane & Thomas, Brune, deBoer, Descouvemont & Baye, Azuma) are in
`references/key-papers.md`. When you need a number, an equation, or a citable
sentence, go to those files rather than reconstructing it from memory.

## 1. What R-matrix theory is, in one picture

Imagine a nuclear reaction — say `12C + α → 16O + γ`, or `p + 14N` scattering. A
projectile and target approach, spend some time as a complicated many-body
compound system, and then fly apart into some exit channel. R-matrix theory is a
bookkeeping scheme for exactly this, and it rests on a single geometric idea:
**divide space into two regions at a sphere of radius `a_c` (the channel radius),
one per channel.**

- **Internal region** (`r < a_c`): the many-body nuclear soup. Too complicated to
  solve from first principles, so we do *not* try. Instead we parametrize whatever
  happens in here with a handful of numbers — a set of **levels** with energies
  `E_λ` and **reduced-width amplitudes** `γ_λc` that say how strongly each level
  couples to each channel.
- **External region** (`r > a_c`): the projectile and target are two separate
  clusters feeling only the Coulomb force (plus centrifugal). Here the physics is
  *known exactly* — the wavefunctions are Coulomb functions.

The whole method is: solve the known external physics, describe the unknown
internal physics with a few parameters, and **match the two at `r = a_c`**. The
matching condition is the logarithmic derivative of the radial wavefunction at the
surface, and the R-matrix is precisely the function that encodes that logarithmic
derivative versus energy:

```
R_{c'c}(E) = Σ_λ  γ_{λc'} γ_{λc} / (E_λ − E)
```

Each level `λ` is a pole. Near `E ≈ E_λ` one term blows up and dominates — that is
a resonance. Away from any pole, the smooth sum of tails is the "background." From
`R` you build the collision (scattering) matrix `U`, and from `U` come all cross
sections and angular distributions. **You never need the internal wavefunction
itself** — only its behavior at the surface. That is the entire trick.

Two flavors exist. The **calculable** R-matrix computes `E_λ` and `γ_λc` from a
model Hamiltonian. **Phenomenological** R-matrix — the kind used in AZURE2 and the
kind this skill is about — treats `E_λ` and `γ_λc` as *fit parameters* adjusted to
reproduce measured data. Its central use is **extrapolation**: measure a cross
section where you can, fit it, and evaluate it where you cannot (typically at
astrophysical energies far below any data). The R-matrix is the physically
constrained interpolating/extrapolating function that makes that honest.

## 2. The vocabulary: channels, levels, boundary conditions

**Channel** `c = (α, s, l)`: a way in or out. `α` is the pair of nuclei (the
partition, e.g. `12C+α` or `16O+γ`), `s` is the channel spin (vector sum of the two
nuclear spins), `l` is the orbital angular momentum. `p + 14N` with the two spins
coupled to `s=1` in a `d`-wave is a different channel from the same pair in an
`s`-wave. Parity and total `J` must be conserved across channels that connect.

**Level** `λ`: an internal-region basis state, with energy `E_λ` and one reduced-
width amplitude `γ_λc` per channel. In phenomenological work a level *usually*
corresponds to a physical excited state of the compound nucleus — but not always
(see background poles, §7). A level of a given `Jπ` only couples to channels with
that same `Jπ`.

**Reduced width `γ_λc` / partial width `Γ_λc`**: `γ_λc` is the amplitude of the
level's wavefunction at the channel surface — how strongly level `λ` opens into
channel `c`. Its *sign matters* (that is what produces interference, §6). The
observable **partial width** is `Γ_λc = 2 P_c(E) γ²_λc`, where `P_c` is the
penetration factor — the probability of getting through the Coulomb + centrifugal
barrier. Because `P_c` is strongly energy-dependent, partial widths are energy-
dependent, and a low-energy resonance is narrow not because it couples weakly but
because the barrier is nearly shut.

**Penetration `P_c`, shift `S_c`, hard-sphere phase `φ_c`**: three real functions
of energy built from Coulomb functions at `a_c`. `P_c` controls widths, `S_c`
controls the energy shift of a resonance, `φ_c` is the smooth non-resonant
(potential/hard-sphere) scattering phase. You rarely compute these by hand — AZURE2
does — but knowing that `P_c` rises steeply with energy explains most qualitative
behavior of excitation functions.

**Boundary condition `B_c`**: a real constant you must fix per channel; it defines
the logarithmic derivative the internal basis states satisfy at the surface. The
awkward truth: `E_λ` and `γ_λc` **depend on the arbitrary choice of `B_c`**, so the
raw fit parameters are not directly the physical resonance energy and width. This
is the single biggest conceptual trap in R-matrix work, and §8 (Brune) is how the
field escapes it. A common convenient choice is `B_c = S_c(E_λ)`, which makes the
level shift vanish so `E_λ` coincides with the resonance energy — but it can only be
set to satisfy *one* level at a time.

## 3. How physical resonances map to R-matrix levels

This is the daily translation task. You have a level scheme (from ENSDF/NDS via
`nds-explorer`); you need R-matrix inputs.

- **A resonance above threshold** → a level with `E_λ` at the resonance energy and
  `γ_λc` set so that `Γ_λc = 2 P_c γ²_λc` reproduces the tabulated partial widths.
  Its `Jπ` fixes which channels it couples to.
- **Energies**: use the *observed* (physical) resonance energy, and be aware it may
  differ from the peak of the excitation function because `P_c` and the phase-space
  factor `1/k²` also vary with energy.
- **Partial widths**: literature widths are "observed" widths. The R-matrix "formal"
  width `Γ_λc = 2P_cγ²_λc` differs from the observed width — often by ~30% for broad
  levels — through the factor `[1 + Σ_c γ²_λc dS_c/dE]`. Do not paste an observed
  width in where a formal one is expected without accounting for this (Brune's
  parametrization, §8, removes the hazard by letting you use observed values
  directly).
- **Spins/parities**: `Jπ` decides everything about which channels and which
  electromagnetic multipoles participate. Getting `Jπ` wrong is the most common way
  to build a physically impossible fit.
- **Bound states and subthreshold levels**: a state *below* a particle threshold
  cannot decay into that channel, so its partial width there is zero. It is still a
  genuine **pole of the scattering matrix** and it still shapes the cross section
  above threshold through its tail. You parametrize it not by a width but by an
  **ANC (asymptotic normalization coefficient)** — the amplitude of the bound-state
  wavefunction in the external region. Subthreshold states frequently *dominate* the
  low-energy cross section (in `12C(α,γ)` the −45 keV `1⁻` and −245 keV `2⁺`
  subthreshold states are the largest contributors at the 300 keV Gamow energy), so
  their ANCs are often the most important — and hardest to pin down — parameters in
  an astrophysical fit. ANCs come from transfer reactions or `β`-delayed emission,
  not from the resonance region itself.

## 4. Selecting which levels to include

Truncation is unavoidable: the exact sum over levels is infinite, and you must
choose a finite set. The guiding principle is **physical necessity first, then
completeness, then background.**

**Physically required levels** — always include:

- Every known resonance of the right `Jπ` that sits *inside or near* the energy
  range of your data. These produce the peaks you are fitting.
- Every subthreshold/bound state of the right `Jπ` whose tail reaches into your
  energy range — especially for low-energy extrapolation. Omitting a subthreshold
  state that dominates the Gamow window makes the extrapolation meaningless.
- Any state that is known to interfere strongly with the above (same `Jπ`,
  comparable or broad width).

**Optional / judgement-call levels**:

- Narrow, weakly-coupled resonances far from your data range: often safely fixed at
  literature values or omitted, but include them if a nearby measurement constrains
  them.
- High-`l` channels and energetically strongly-closed channels: typically dropped.

**Practical rules of thumb**:

- Start from the known level scheme up to a chosen excitation energy; model
  everything above with background poles (§7) rather than dozens of real levels.
- Only levels of a given `Jπ` contribute to an observable that selects that `Jπ`
  (e.g. for E1 ground-state capture in `16O`, only `1⁻` levels contribute; for E2,
  only `2⁺`). Use this to keep each partial cross section's level list minimal.
- Fix well-known parameters at literature values; free only what the data can
  actually constrain. Freeing everything invites overfitting.
- Before adding a level to improve `χ²`, ask whether the improvement is physical or
  just more free parameters. Use the significance test in `azure2-eval` (an F-test /
  Δχ² criterion) to justify each added level.

## 5. Capture (γ) channels vs particle channels

Particle channels (elastic, transfer, breakup) and photon channels behave
differently, and mixing them up causes errors.

- **Particle channels** are unitary and enter the collision matrix `U` fully. Their
  widths are large and they set a resonance's **total width** (its horizontal extent
  in energy).
- **Photon (capture) channels** are weak and are treated by **perturbation theory**,
  *outside* the unitary particle framework. Consequences you must remember:
  - γ channels **do not contribute to the total width** in the Breit-Wigner
    denominator. The particle width(s) set the resonance's width; the γ width sets
    only its *height*.
  - Capture has an **internal** part (the photon reduced width `γ_λp`, a real fit
    parameter per level per multipole) and an **external / channel** part that
    depends only on nuclear parameters and the **final-state ANC**. At low energy,
    external (including "hard-sphere") capture and cascade transitions are often
    dominated by the external term, so they are relatively insensitive to internal
    parameters.
  - Multipolarity matters: E1, E2, M1 each select different `Jπ` combinations and
    each has its own set of contributing levels and its own interference pattern.
    Separating E1 from E2 is usually done through **angular distributions**, not the
    total cross section (§9).
  - Because γ widths ≪ particle widths, capture data mainly constrain the *products*
    of γ amplitudes with particle amplitudes and the interferences among them —
    which is exactly why capture cross sections are so sensitive to relative signs.

## 6. Interference — the heart of the analysis

Interference is where R-matrix intuition pays off, and where most fitting mistakes
hide.

**Same-`Jπ` (energy-dependent) interference.** Two levels of the *same* `Jπ`
contribute *amplitudes* to the same channel, and amplitudes add before squaring.
The cross section carries a cross term whose size goes like `2√(σ₁σ₂)` — so **even a
weak level can strongly modify the cross section far from its own peak** through
interference with a strong one. The **relative sign of the reduced-width
amplitudes** decides whether they add (constructive) or cancel (destructive) in the
valley between and in the tails. This is not a small effect: in `12C(α,γ)` the
constructive vs destructive E1 choice moves the extrapolated S-factor by a large
factor and the two solutions differ by hundreds in `χ²`. **Theory usually cannot
predict the sign** — you determine it by which sign fits the data, and capture data
between resonances are the decisive constraint. Practically: always test both
relative signs of interfering same-`Jπ` levels; the off-resonance / valley regions
and the low-energy tail are where the sign shows itself.

**Different-`Jπ` (angle-dependent) interference.** Amplitudes with different `Jπ`
cannot interfere in the *total* cross section, but they *do* interfere in the
**angular distribution**, producing odd Legendre terms (e.g. E1–E2 interference
gives a forward-backward asymmetry). This is how multipoles are separated (§9).

**The ghost anomaly.** A broad resonance sitting just above a threshold has a width
`Γ_c(E) = 2P_c(E)γ²_c` that is throttled to near-zero at threshold by the
penetrability and then rises steeply. The result is a strongly skewed, asymmetric
line shape: the resonance is pushed up and a secondary "ghost" bump can appear
displaced from the nominal energy. The textbook case is the `8Be` ground state seen
in `α+α`; it is a real, calculable consequence of the energy dependence of `P_c` and
`S_c`, not an artifact. Recognize it so you do not "fix" it by inventing an extra
level. (Barker & Treacy 1962; discussed in Lane & Thomas — see
`references/key-papers.md`.)

## 7. Background poles

The levels you do not include explicitly — everything high above your data range —
do not vanish; their low-energy tails add up to a smooth background. You represent
that sum with one or more **background poles** per `Jπ`, placed at an energy well
*above* your highest data point.

- Their energies are **arbitrary** and their individual parameters are **not
  unique** — what is physically meaningful is the roughly energy-independent
  background contribution they produce, which can be reproduced many equivalent
  ways. Do not over-interpret a background pole's `E_λ` or width as a real state.
- More background poles → lower `χ²` almost automatically. Always ask whether an
  added pole is physically warranted or just soaking up free parameters. If explicit
  known levels already describe the capture data, you may need *no* capture-partition
  background pole even while the scattering partition still needs several.
- Background poles and the **channel radius are correlated**: a smaller radius
  encloses less of the interaction and leans harder on background poles; a larger
  radius reduces that need but increases the density of background poles and can make
  them behave badly. Typical phenomenological radii (~4–7 fm) enclose most but not
  all of the nuclear interaction. The channel radius is *not* a physical nuclear
  radius even though the fitted value often lands close to one.

## 8. Boundary conditions and Brune's parametrization

The `B_c` problem from §2: standard R-matrix parameters `E_λ, γ_λc` depend on the
arbitrary boundary constants, so they are not directly the physical resonance
energies and widths, and you cannot in general set `B_c = S_c(E_λ)` for every level
at once. Historically you converted formal ↔ observed parameters level by level,
iteratively — error-prone and confusing.

**Brune's alternative parametrization (Brune 2002)** rewrites the theory so the
scattering matrix is expressed directly in terms of the **observed** resonance
energies `Ẽ_λ` and **on-resonance** reduced-width amplitudes `γ̃_λc` of *all*
levels. In this formulation:

- The boundary constants `B_c` **never appear** — the results are automatically
  boundary-condition independent.
- Every parameter has a **direct physical meaning**: `Ẽ_λ` is the resonance energy,
  `γ̃_λc` maps simply onto the observed partial width.
- You can **drop literature values straight in** as starting parameters and read
  fitted values straight out, without formal↔observed gymnastics.

It is mathematically equivalent to the original theory (same cross sections) — only
the parametrization changes. The channel radii still must be chosen. **Prefer
Brune's parametrization** for phenomenological fits that draw on tabulated level
data; it is what modern analyses (and AZURE2) use. See
`references/key-papers.md` (Brune) and `references/formalism-and-fitting.md` for the
transformation itself.

## 9. Angular distributions

R-matrix predicts full differential cross sections, not just angle-integrated ones.
The differential cross section expands in Legendre polynomials `P_L(cos θ)` with
coefficients built from the same transition-matrix elements. Two uses dominate:

- **Multipole separation.** Different-`Jπ` amplitudes (e.g. E1 and E2) interfere
  only in the angular distribution, generating specific Legendre terms. Fitting the
  angular distribution — often via the E1–E2 phase difference — is how you
  disentangle multipoles that the total cross section cannot separate. Measuring at
  90° can isolate a single multipole because interference and some terms vanish
  there.
- **Spin/parity assignment.** The shape of the angular distribution constrains the
  `Jπ` of a resonance, feeding back into level selection (§4).

Analyzing powers and polarization observables extend this further and are handled by
AZURE2's segment machinery (see `azure2-eval`).

## 10. Common pitfalls

- **Using observed widths as formal widths** (or vice versa) without the conversion.
  Use Brune's parametrization to sidestep it.
- **Wrong `Jπ`**: coupling a level to channels or multipoles it cannot reach.
- **Forgetting a subthreshold state** that dominates the low-energy extrapolation.
- **Not testing both interference signs** for same-`Jπ` levels — the fit can look
  fine on-resonance and be wildly wrong in the extrapolation.
- **Over-parametrizing with background poles** and reading physical meaning into
  them; ignoring the radius–background correlation.
- **Chasing the lowest `χ²`** rather than the most physical fit. The best fit is not
  always the minimum-`χ²` fit; data uncertainties are rarely fully correct, and
  normalization/systematic treatment (Peelle's pertinent puzzle) matters.
- **Treating the channel radius as a physical radius**, or changing it without
  re-examining the background poles.
- **Overlooking that γ widths set peak height, not total width.**

## 11. Connecting to AZURE2 inputs

The theory maps onto the `.azr` project as follows (details and mechanics in
`azure2-eval`):

- **Particle pairs** ↔ the partitions `α` (entrance/exit pairs, with their spins,
  parities, `Q`-values, channel radii). The channel radius (§2, §7) is set here.
- **Levels** ↔ the R-matrix levels `λ`: one entry per compound-nucleus state, with
  `Jπ` and excitation energy. Subthreshold states are levels below the relevant
  threshold. Background poles are extra high-lying levels of a given `Jπ`.
- **Channels within a level** ↔ the `γ_λc` (or, in Brune mode, `γ̃_λc`): one
  reduced-width amplitude per open channel, plus photon channels per multipole. Their
  **signs** are the interference knobs (§6). ANCs enter for bound/subthreshold
  channels.
- **Segments** ↔ the data being fit (integrated or differential cross sections,
  analyzing powers), each with its own normalization.
- **Boundary condition** ↔ AZURE2 uses Brune's parametrization, so you enter and read
  **observed** energies and widths directly.

Workflow: use `nds-explorer` to pull the level scheme and data; use *this* skill to
decide the level list, the interference signs to test, the subthreshold states and
background poles needed, and what each parameter means; then use `azure2-eval` to
build the `.azr`, fit, run significance tests on candidate levels, decompose the
cross section into level/interference/external contributions, and extrapolate.

## Reference files

- `references/key-papers.md` — annotated digest of Lane & Thomas (1958), Brune
  (2002 + 2005 notes), deBoer et al. (2017), Descouvemont & Baye (2010), Azuma et al.
  (2010), and others: what each contributes, and when to reach for it.
- `references/formalism-and-fitting.md` — the equations behind the picture: R-matrix
  and collision matrix, penetration/shift/phase functions, formal vs observed
  parameters, the Brune transformation, capture formalism, and the `χ²`/fitting
  machinery — with pointers to the exact source equations.
