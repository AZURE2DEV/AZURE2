# The R-matrix picture for building a model

This reference states, in the order needed to build a project, what the R-matrix parameters mean physically and
how a known level scheme becomes a list of levels and channels. Interference, background poles, boundary conditions
and parametrization bases (formal, Brune, Park) are in azure2-fit; the equations are in azure2-fit's formalism
reference; the primary literature is in references/key-papers.md.

## 1. The picture

Space is divided at the channel radius a_c (fm), one per channel.

- **Internal region** (r < a_c): the many-body compound system, not solved. It is parametrized by levels λ with
  energies E_λ (MeV) and reduced-width amplitudes γ_λc (MeV^½) that say how strongly level λ couples to channel c.
- **External region** (r > a_c): two clusters interacting through the Coulomb and centrifugal potentials only;
  the radial functions are Coulomb functions, known exactly.

Matching the two at r = a_c through the logarithmic derivative gives the R-matrix

```
R_c'c(E) = Σ_λ γ_λc' γ_λc / (E_λ − E)                                         (1)
```

Each level is a pole; near E ≈ E_λ one term dominates (a resonance); away from all poles the sum of tails is the
background. The collision matrix U follows from R, and every cross section and angular distribution from U. Only
the surface behaviour of the internal wave function enters.

The *calculable* R-matrix derives E_λ and γ_λc from a Hamiltonian. The *phenomenological* R-matrix (AZURE2) fits
them to data. Its main use is extrapolation: fit where the cross section is measured, evaluate it where it is not
(typically astrophysical energies below the data).

## 2. Vocabulary

- **Channel** c = (α, s, l): partition α (the pair, e.g. ¹²C+α or ¹⁶O+γ), channel spin s (the vector sum of the
  two spins), orbital angular momentum l. Total J and parity are conserved; a level of a given J^π couples only to
  channels of that J^π.
- **Level** λ: an internal basis state with E_λ and one γ_λc per channel. Usually a physical state of the compound
  nucleus, but not always (background poles).
- **Reduced and partial width.** γ_λc is the amplitude of the level at the channel surface; its sign matters
  (interference). The partial width is Γ_λc = 2 P_c(E) γ²_λc, with P_c the penetration factor. P_c rises steeply
  with energy, so a low-energy resonance is narrow because the barrier is nearly closed, not because it couples
  weakly.
- **P_c, S_c, φ_c**: penetration factor, shift function and hard-sphere phase, real functions of E built from the
  Coulomb functions at a_c. P_c controls widths, S_c the resonance shift, φ_c the non-resonant phase.
- **Boundary condition B_c**: a real constant per channel defining the logarithmic derivative of the basis states at
  the surface. E_λ and γ_λc depend on the choice of B_c, so the raw parameters are not the physical resonance energy
  and width. B_c = S_c(E_λ) removes the shift of one level only. AZURE2 parametrizes by observed (Brune) or by
  Park parameters, so a project carries physical energies and widths (azure2-fit, azure2-exchange).

## 3. From a level scheme to levels

The level scheme comes from ENSDF (references/nds-exfor.md §5).

- **A resonance above threshold** → a level with E_λ at the resonance energy and widths reproducing the tabulated
  partial widths. Its J^π fixes the channels.
- **Energies.** Use the observed resonance energy; the peak of an excitation function can differ from it because
  P_c and 1/k² vary with energy.
- **Widths.** Literature widths are observed widths. The formal width 2P_cγ² differs from the observed one by the
  factor [1 + Σ_c γ²_λc dS_c/dE], often about 30 % for broad levels. AZURE2's `<levels>` block takes observed
  values (Γ in eV), so they can be entered directly.
- **J^π** decides which channels and multipoles take part. A wrong J^π is the most common way to build a physically
  impossible fit.
- **Bound and sub-threshold states.** A state below a threshold has no width in that channel but is still a pole of
  the collision matrix and shapes the cross section above threshold through its tail. It is parametrized by an ANC
  (fm^-1/2), the amplitude of its external wave function. Sub-threshold states often dominate low-energy cross
  sections: in ¹²C(α,γ) the 1⁻ state 45 keV and the 2⁺ state 245 keV below the α threshold are the largest
  contributions at the Gamow energy E_c.m. ≈ 300 keV. ANCs come from transfer reactions or β-delayed emission, not
  from the resonance region.

## 4. Which levels to include

The exact sum over levels is infinite; the truncation follows physical necessity first, then completeness, then
background.

Always include:
- every known resonance of a contributing J^π inside or near the energy range of the data;
- every bound or sub-threshold state of a contributing J^π whose tail reaches the range, above all for low-energy
  extrapolation;
- states known to interfere strongly with these (same J^π, comparable or large width).

Judgement calls:
- narrow, weakly coupled resonances far from the data: fix at literature values or omit, unless a nearby
  measurement constrains them;
- high-l and strongly closed channels: usually dropped.

Rules of thumb:
- Start from the known scheme up to a chosen E_x; represent everything above by background poles (azure2-fit).
- An observable that selects a J^π sees only those levels (E1 ground-state capture in ¹⁶O: 1⁻ only; E2: 2⁺ only).
  Keep each partial cross section's level list minimal.
- Fix well-known parameters at literature values; free only what the data constrain.
- Justify each added level by a significance test, not by Δχ² alone (azure2-check, azure2-resonance-search).

## 5. Photon channels versus particle channels

- Particle channels enter the unitary collision matrix fully and set the total width of a resonance.
- Photon channels are treated in first-order perturbation theory, outside the unitary particle framework:
  - Γ_γ does not enter the total width in the denominator; it sets the peak height, not the width.
  - Capture has an internal part (a photon reduced width per level and multipole, a fit parameter) and an external
    (channel) part depending on the nuclear parameters and the final-state ANC. At low energy, external capture and
    cascade transitions are often dominated by the external term.
  - E1, E2, M1 select different J^π combinations, each with its own levels and interference pattern. E1 and E2 are
    usually separated through angular distributions (§6).
  - Because Γ_γ ≪ Γ_particle, capture data constrain products of γ and particle amplitudes and their interference,
    which makes capture sensitive to relative signs.

## 6. Angular distributions

The differential cross section expands in Legendre polynomials P_L(cos θ_c.m.) with coefficients built from the same
transition-matrix elements.
- **Multipole separation**: amplitudes of different J^π (E1 and E2) interfere only in the angular distribution.
  Fitting it, often through the E1–E2 phase difference, separates multipoles the total cross section cannot. At
  θ_c.m. = 90° some interference terms vanish and a single multipole can be isolated.
- **J^π assignment**: the shape constrains the J^π of a resonance and feeds back into level selection (§4).
- Analysing powers and polarization observables are segment observables (azure2-calculate).

## 7. Mapping onto the `.azr`

- **Particle pairs** ↔ partitions α: spins, parities, masses, separation energies, channel radii a_c.
- **Levels** ↔ R-matrix levels λ: one per compound-nucleus state, with J^π and E_x. Sub-threshold states are levels
  below their threshold; background poles are extra high-lying levels of a J^π.
- **Channels of a level** ↔ γ_λc (in the file, the observed Γ or ANC): one per open channel plus one per photon
  multipole. Their signs are the interference choices; ANCs for bound or sub-threshold channels.
- **Segments** ↔ the data (integrated or differential cross sections, analysing powers), each with its own
  normalization.

Order of work: data and level scheme (azure2-setup) → level list, signs to test, sub-threshold states and background
poles (this reference and azure2-fit) → fit (azure2-fit) → significance of candidate levels (azure2-check) →
decomposition and extrapolation (azure2-calculate).
