# Key papers — annotated digest

Grounding literature for the `r-matrix-analysis` skill. Each entry says **what the
paper gives you** and **when to reach for it**. Items marked *(in library)* are in
the user's Zotero library with the item key shown, so they can be pulled with the
Zotero tools.

---

## Lane & Thomas (1958) — the foundation
**A. M. Lane and R. G. Thomas, "R-Matrix Theory of Nuclear Reactions,"
Rev. Mod. Phys. 30, 257 (1958).** *(in library, Zotero key `RPC23CVA`)*

The canonical, exhaustive formulation of R-matrix theory. Everything downstream —
channels, the internal/external division, the R-matrix and collision-matrix
equations, penetration/shift factors, boundary conditions, formal vs observed
widths, the single-level (Breit-Wigner) reduction — is defined here, and modern
papers keep its notation (`E_λ`, `γ_λc`, `B_c`, `P_c`, `S_c`, `A`, `U`).

- **What to take from it:** the definitions and the internal/external picture; the
  single-level formula and the relation between formal and observed widths
  (Sec. XII.3); the origin of the boundary-condition dependence. Its Sec. on the
  "ghost" of a level near threshold underlies the ghost-anomaly discussion.
- **When to reach for it:** you need the *authoritative* definition of a quantity or
  equation, or you need to check exactly how a width/shift/penetrability is defined.
  It is a reference work, not a tutorial — huge and dense; do not read it front to
  back, look up the specific section.

## Descouvemont & Baye (2010) — the pedagogical review
**P. Descouvemont and D. Baye, "The R-matrix theory," Rep. Prog. Phys. 73,
036301 (2010).** *(in library, Zotero key `5EHXWJ2N`)*

The best modern *pedagogical* entry point. Presents both facets in one framework:
the **calculable** R-matrix (a numerical tool to solve the Schrödinger equation)
and the **phenomenological** R-matrix (fitting data). Clean derivations of the
penetration and shift functions, boundary conditions, and the collision matrix.

- **What to take from it:** intuition and clean derivations; the clearest statement
  of the calculable-vs-phenomenological distinction; a readable account of `P_c`,
  `S_c`, `φ_c` and channel radii.
- **When to reach for it:** you want to *understand* a piece of the formalism (not
  just cite it), or explain it accessibly to a non-expert.

## Brune (2002) — the alternative parametrization
**C. R. Brune, "Alternative parametrization of R-matrix theory,"
Phys. Rev. C 66, 044611 (2002).** *(primary source; cited by the 2005 notes below,
which are in the library)*

Recasts R-matrix theory so the collision matrix is expressed **directly** in terms
of the observed resonance energies `Ẽ_λ` and on-resonance reduced-width amplitudes
`γ̃_λc`, with the **boundary constants `B_c` eliminated entirely**. Mathematically
equivalent to Lane & Thomas (same observables), but every parameter has a direct
physical meaning and literature values can be used as-is.

- **What to take from it:** *why the raw R-matrix parameters are not physical* (they
  depend on arbitrary `B_c`) and how the alternative parametrization fixes it; this
  is what AZURE2 uses and what makes tabulated level parameters usable directly.
- **When to reach for it:** any time you are entering or interpreting fit parameters,
  or explaining why AZURE2's energies/widths are the "observed" ones.

## Brune (2005) — the accessible lecture notes on the same idea
**C. R. Brune, "Formal and Physical R-matrix parameters," notes from the 2004 JINA
R-Matrix School, arXiv:nucl-th/0502087 (2005).** *(in library, Zotero key
`WZ6NKJCU`)*

A short, teaching-oriented walk-through of the relationship between **formal**
R-matrix parameters (`E_λ`, `γ_λc`, tied to `B_c`) and **physical/observed**
resonance parameters. Contains the essential results in digestible form:

- How to define a resonance energy (peak, phase = π/2, S-matrix pole, Breit-Wigner)
  and why prescriptions differ most for **broad** resonances.
- The single-level R-matrix → Breit-Wigner reduction, and the identifications
  `E_R = E_λ` (with `B_c = S_c(E_λ)`) and the observed partial width
  `Γ°_λc = 2 P_c γ²_λc / [1 + Σ_c γ²_λc dS_c/dE]`. Note the caution that **formal
  vs observed widths often differ by ~30%** — large next to experimental errors.
- The **invariance of `U` under changes of `B_c`** (with a compensating `E_λ, γ_λc`
  transformation) — the "gauge freedom" of R-matrix theory.
- The eigenvalue route to observed parameters, converting observed↔formal, and
  computing observables directly from physical parameters (the alternative
  parametrization) without ever choosing `B_c`.
- A concrete `16O` `1⁻` (`12C+α` / `16O+γ`) worked example.

- **When to reach for it:** you need the actual equations for formal↔observed
  conversion, or a citable accessible explanation of the `B_c` gauge freedom. Start
  here before the denser 2002 paper.

## deBoer et al. (2017) — the modern practical/global analysis
**R. J. deBoer et al., "The 12C(α,γ)16O reaction and its implications for stellar
helium burning," Rev. Mod. Phys. 89, 035007 (2017).** *(in library, Zotero key
`6CPV5J8A`)*

The best worked demonstration of a **large, global, phenomenological R-matrix
evaluation** and a goldmine of practical methodology. Even if you never touch
`12C(α,γ)`, its R-matrix section (Sec. IV) and fitting discussion are the template.

Practical points it establishes (with representative findings):

- **Role & extrapolation:** phenomenological R-matrix is the tool of choice for the
  resolved-resonance region; the whole point is extrapolating a fit to energies well
  below the data (the 300 keV Gamow window here).
- **Channel radius:** must enclose *most but not all* of the nuclear interaction; too
  large increases background-pole density and makes them misbehave; it is **not** a
  physical nuclear radius even when the fitted value lands near one. (Best-fit values
  ~5 fm here.)
- **Level truncation & selection:** include known levels up to some `E_x`, model the
  rest with background poles; drop strongly-closed channels and high-`l` channels;
  for a given multipole only the matching `Jπ` levels contribute (E1 ← `1⁻`, E2 ←
  `2⁺`).
- **Background poles:** placed *above* the highest data; individually **non-unique**
  (only the net background they produce is meaningful); more poles lower `χ²`
  trivially, so justify them; radius↔background-pole correlation.
- **Interference:** implemented via the **relative signs of reduced-width
  amplitudes**; interference `∝ 2√(σ₁σ₂)`, so weak levels matter off-resonance;
  same-`Jπ` → energy-dependent effects (critical for extrapolation), different-`Jπ`
  → angle-dependent effects (separate multipoles). Signs are not predictable by
  theory — capture data decide them; here the destructive E1 solution is ruled out
  by Δχ² ≈ 324.
- **Capture treated perturbatively:** photon channels do **not** enter the total
  width; internal (photon reduced width) + external/channel capture (depends on
  final-state **ANC**); external capture can be suppressed by an effective-charge
  cancellation; cascades dominated by external capture at low energy.
- **Subthreshold states & ANCs:** the −45 keV `1⁻` and −245 keV `2⁺` subthreshold
  states dominate the low-energy cross section; parametrized by ANCs (from transfer
  and `16N(βα)`), which are the dominant extrapolation uncertainty.
- **Brune parametrization used exclusively**, "in order to more conveniently utilize
  level parameters from the literature."
- **Fitting realities:** the adopted "best fit" was *not* the lowest-`χ²` fit but the
  most physically reasonable; uncertainties in the data are rarely fully correct;
  normalization factors and a systematic term address Peelle's pertinent puzzle;
  robust estimators (Sivia–Skilling) down-weight outliers; Monte Carlo for
  uncertainties needs reduced-`χ²` ≈ 1 to be trustworthy.
- **Angular distributions:** Legendre expansion; E1–E2 phase difference separates
  multipoles; 90° isolates a single multipole.

- **When to reach for it:** designing or critiquing a real fit — level lists,
  background poles, interference-sign tests, capture treatment, subthreshold/ANC
  handling, error estimation. This is the "how it's actually done" reference.

## Azuma et al. (2010) — the AZURE code paper
**R. E. Azuma et al., "AZURE: An R-matrix code for nuclear astrophysics,"
Phys. Rev. C 81, 045805 (2010).** *(in library, Zotero key `BHW47FJJ`)*

Describes the multi-channel, multi-level R-matrix code (predecessor of AZURE2).
Documents how the theory becomes software: partitions/particle pairs, levels and
channels, data segments, the level matrix, and capture treatment.

- **When to reach for it:** connecting formalism to the `.azr` input structure, or
  citing the code. For actually driving AZURE2, use the `azure2-eval` skill.

---

## Supporting / related items in the library

- **I. J. Thompson et al., "Verification of R-matrix calculations ... 7Be system,"
  Phys. Rev. C 100, 015805 (2019)** *(Zotero key `RSHCS7N3`)* — code
  cross-verification (AMUR, AZURE2, EDA, SAMMY, …); useful for trusting/benchmarking
  results and for standard definitions.
- **D. Odell, C. R. Brune, D. R. Phillips et al., "Performing Bayesian Analyses With
  AZURE2 Using BRICK ... 7Be,"** *(Zotero key `VF28TJ8S`)* and **"How Bayesian
  methods can improve R-matrix analyses ... d+t,"** *(Zotero key `8PFWKX3B`)* —
  modern Bayesian/uncertainty-quantification workflows layered on AZURE2 (the BRICK
  toolkit). Reach for these when the task is uncertainty estimation or MCMC, not just
  a point fit.
- **Leeb, Dimitriou & Thompson, IAEA R-matrix code reports** *(Zotero keys
  `H8R39NBY`, `UZCCXHAS`, `6NEQH6NW`)* — the IAEA code-comparison exercises; context
  for evaluation standards.
- **Acharya et al., "Solar Fusion III" (2025)** *(Zotero key `XNNCPMAM`)* — recent
  evaluated reaction data and methodology for hydrogen-burning reactions; good for
  current recommended values.

## Not in the library (cite from primary source)

- **F. C. Barker and P. B. Treacy, Nucl. Phys. 38, 33 (1962)** — origin of the
  **ghost anomaly** treatment (broad near-threshold levels; classic `8Be` case).
- **E. P. Wigner and L. Eisenbud, Phys. Rev. 72, 29 (1947)** — the original R-matrix
  paper that Lane & Thomas systematized.
- **F. C. Barker, Aust. J. Phys. 25, 341 (1972)** — invariance of `U` under `B_c`
  changes for finite level number (the result Brune builds on).
