# Background poles

The levels not included explicitly (everything high above the data range) do not vanish: their low-energy tails add up
to a smooth background in each J^π. A fit represents that sum by one or more background poles per J^π, levels placed
at an excitation energy E_x well above the highest data point. This page states what their parameters mean, how to
restrict a pole to the channels that need it, and what to fix. A pole that runs to saturation is a trust question:
`azure2-check`.

## What a background pole means

- Its energy is arbitrary and its individual parameters are not unique. What the data determine is the roughly
  energy-independent contribution it produces, which many parameter sets reproduce. Do not interpret a background
  pole's E_λ or width as a state.
- More background poles lower χ² almost automatically. Ask whether an added pole is warranted or only absorbs free
  parameters. If explicit known levels already describe the capture data, no capture-channel background pole may be
  needed even while the scattering channels need several.
- Background poles and the channel radius a_c (fm) are correlated: a smaller radius encloses less of the interaction
  and leans harder on background poles; a larger one reduces that need but raises the density of background poles
  and can make them behave badly. Typical phenomenological radii, 4–7 fm, enclose most but not all of the nuclear
  interaction. The channel radius is not a physical nuclear radius, even when the fitted value lands near one. A
  radius change calls for re-examining the background poles (radius sensitivity: `azure2-check`).

## Restricting a pole to the channels that need it

A background pole written with every channel of its J^π couples all of them: strength added for one reaction is added
for every reaction that shares those channels. When the extra strength is wanted in one data set only, give the pole
widths only in that data set's channels and leave every other channel at exactly zero, which the engine then holds
fixed (seed the ones you do want with a nonzero value).

Worked case (¹⁶O): ¹⁵N(p,γ₀) could not be described together with ¹²C(α,γ₀), because the shared photon widths of the
E_x = 12–13.3 MeV levels needed amplitudes about 2.4× larger for (p,γ₀) than (α,γ₀) allows. With a free, unpenalized
segment normalization the fit hid this (the normalization slid to 0.17, and the data set was effectively dropped while
still adding leverage); with the normalization anchored and only ordinary E_x = 40 MeV background poles free, the fit
reached for a 3⁻ pole at 10⁴ W.u. The solution of deBoer et al., Rev. Mod. Phys. 89 (2017) 035007, is a 1⁻ pole at
E_x = 17.09 MeV carrying only a p+¹⁵N width (l = 0, s = 1, fixed at Γ = 500 keV) and a free γ₀ E1 width (about 500 eV,
0.23 W.u.), with every α and cascade channel at zero. Having no α channel, it cannot change ¹²C(α,γ₀).

- **Fix one of the two widths.** A pole feeding one reaction through one entrance and one exit channel determines only
  the product γ_in γ_out, so fitting both is degenerate. Hold the particle width at a plausible value and fit the
  photon width.
- **Adding a pole renames the parameters.** Level parameter names are `width_<levelIndex>_<channelIndex>` in the
  engine's own ordering, so inserting a level shifts the indices of the others. A saved vector (`param.sav`, an
  `.npz`, an MCMC chain) cannot be carried into the edited model by name; map by physical identity
  (`references/parameters-and-sav.md`).
- **Seed in Wigner units.** A freed pole width that is exactly 0 stays fixed; seeding six E_x = 40 MeV capture
  background poles at 0.01 eV was what raised one model's free-parameter count from 366 to 372 and let the fit use
  them.

- **Seeding a dead channel takes two steps.** A physical seed of 20–200 keV is a reduced width of thousands of MeV in
  a high-l channel (P ~ 10⁻⁶); it moved one model by +137,000 to +167,000 in χ². Write a tiny physical value so the
  channel exists as a parameter (Γ = 10⁻¹² eV for an open channel, ANC = 10⁻⁶ fm^-1/2 for a closed one; a value of
  "1" in a closed channel's field is an ANC and no small seed), then set the amplitude to frac · γ_W in the engine
  (θ² = frac², 0.04 is enough) and write the model back. A seed should move the objective by less than the first
  fit step, not by 50 %.
- **A seed too large for the penetrability breaks the Brune transform silently.** At P ~ 10⁻¹² even 1 meV exceeds
  γ_W; the log shows `WARNING: Denominator less than zero while transforming J^pi=..., E=...`, the level's energy
  comes out wrong, and writing the model back rescales its unseeded widths (×2.08 in one case) with a perfect round
  trip. After seeding, grep the log for `Denominator less than zero` (it must be empty) and compare every unseeded
  width of the written file with the input (relative change ≲ 10⁻³).
- **An inert pole looks free.** A level with `levelFix = 1` (field 4) and Γ = 0 in every channel creates no parameter,
  although its `channelFix` fields (field 11) read 0. A level whose widths are all exactly zero is invisible to a
  gradient fit even when free: σ is quadratic in the amplitude, so its derivative vanishes at γ = 0. Count the free
  parameters per level and list the levels with none before concluding that a region needs another level; to
  activate such a pole, free it and seed physical widths in the channels the data observe.

## Caps

A Wigner-unit cap on background-pole amplitudes is a modelling decision (how large a pole is physically reasonable),
set per model from where the baseline amplitudes already sit, not a physics default. When a fit stops at the cap,
read how it stopped before raising it (`references/interference.md`, "Guard or runaway").
