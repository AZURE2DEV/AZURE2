---
name: azure2-targets
description: Model the apparatus between the AZURE2 cross section and the measured yield - target integration (thin, thick, compound, implanted, gas), stopping power and straggling, Gaussian beam spread and TOF resolution, SAMMY UDR, beam-profile kernels, free-gas Doppler, transmission and self-shielded yield, electron screening, fitted density or thickness, integration windows, Q attenuation coefficients. Use for any <targetInt> line or experimental effect. Building the project: azure2-setup; fitting: azure2-fit; stopping-profile inference from excitation curves: target-analysis.
---

# Targets, resolution and other experimental effects in AZURE2

An experimental effect is a `<targetInt>` line that folds the calculated c.m. cross section σ(E_c.m.) of the
segments it names with a kernel: the energy loss in the target, the beam energy spread, the detector or
time-of-flight resolution, the thermal motion of the target atoms, the sample attenuation, or the screening of the
nuclear charge. This skill says how to choose and write those lines and how to check them. What each kernel
computes, with its equations, is in `docs/source/user_guide/experimental_effects.rst`; the keyword syntax is in
`docs/source/reference/data_formats.rst`; the implementation (sub-point grids, adaptive lattices) is in
`docs/source/developer/edata_and_target_effects.rst`.

## When to use, and when a sibling skill is the owner

- Use this skill to decide which effect a data set needs, to write or audit its `<targetInt>` line, and to fit an
  apparatus parameter (target density, sample thickness, Doppler temperature, screening potential U_e).
- Project construction, masses, frames of the data file: `azure2-setup`.
- Running calculations and extrapolations: `azure2-calculate`. Freeing, seeding and staging parameters, energy
  shifts and normalizations: `azure2-fit`.
- Inferring a stopping profile or target composition from the excitation curves themselves: the user-level skill
  `target-analysis`.

## 1. Choosing the kernel

| Measurement | Effect on the `<targetInt>` line | Docs section |
|---|---|---|
| Charged beam on a solid target, ΔE comparable with a resonance width Γ | target integration, density N > 0 (thin) | "Target Integration" |
| Yield per incident particle from a target thicker than the range | target integration, density 0 (thick) | "Target Integration" |
| Beam energy spread or detector resolution, Gaussian | Gaussian convolution, σ_b in MeV, lab | "Gaussian Beam Spread" |
| Neutron TOF with an energy-dependent width | Gaussian with an expression of E_lab and `sigmalab` | "Gaussian Beam Spread" |
| Neutron TOF with a tabulated response (n_TOF) | `udr "<file>" L burstFWHM channel centred` | "User-defined numerical resolution function" |
| Broad, skewed γ beam with event-by-event energy (HIγS TPC) | `beamprofile` block | "Beam-profile kernel" |
| Narrow resonances (Γ of order 100 eV) in a room-temperature sample | `doppler <T>` | "Free-Gas Doppler Broadening" |
| Neutron transmission T = exp(−nσ) | segment keyword `transmission <n>` | `docs/source/user_guide/segments.rst` |
| Capture or reaction yield of a thick sample | segment keyword `yield <n>` | `docs/source/user_guide/segments.rst` |
| Electron screening of a low-energy charged-particle reaction | `screening U_e vary σ_U` | "Electron Screening" |
| Angular attenuation of a finite detector | Q coefficients | "Attenuation Coefficients" |

A thin-target point is the depth average Y(E_0) = (1/N) ∫_{E_0−Δ}^{E_0} σ(E′)/ε(E′) dE′, with N the areal density
(atoms cm⁻²), ε the stopping cross section and Δ the energy loss defined by ∫ dE′/ε = N (exact range integral, not
ε(E_0)N). A thick-target point (density 0) is ∫_0^{E_0} σ/ε dE without division by N, integrated down to 1 keV.

**Which kernel the paper describes.**
- **A foil thickness (μg cm⁻²) or "the energy loss in the target was X keV" means target integration**, with an areal
  density and the stopping-power formula, not a Gaussian. Add straggling with the straggling flag on the same line.
- **A detector resolution is not a beam spread.** "7.5 keV FWHM for protons scattered to 144°" is the resolution of
  the detected energy (used to separate peaks); it says nothing about the reaction energy. The Gaussian is for the
  accelerator's beam-energy spread, usually quoted separately (≤ 1 keV on an electrostatic machine).
- **The two kernels are not interchangeable.** Target integration is a one-sided window below E_0, so its first-order
  effect on a feature narrower than Δ is a shift of about Δ/2 to higher nominal energy; broadening is second order.
  A point-centred Gaussian shifts nothing. On a ¹²C(p,p) excitation function through the 1/2⁺ dip (7.5 μg cm⁻² C,
  ΔE_lab = 2.9 keV at 0.45 MeV) target integration moved the model by ±16 % on the two flanks and the dip by
  +1.29 keV (ΔE_c.m./2 = 1.32 keV); a Gaussian with σ_b = Δ/2 moved the points by ≤ 1.5 %. This is the expected
  result, not a defect.
- **"The loss was held constant" with foils of different thickness means N varied with beam energy.** One density
  with an energy-dependent ε does not reproduce it (7.5 μg cm⁻² gives 2.8 keV at 0.46 MeV and 1.1 keV at 2 MeV).
  Match the quoted loss at the energy that matters (2.0 × 10¹⁷ atoms cm⁻² gives 1.5 keV at the 1/2⁺ resonance) or use
  per-window `<targetInt>` lines with lab-energy ranges.
- **Before calling a target-integration discrepancy physics**, run the independent check in
  `references/target-integration-check.md`.

## 2. Frames and units on a `<targetInt>` line

- **Energies on the line are lab, entrance channel, MeV.** σ_b, the `sigmalab` expression (a function of E_lab in
  MeV), the beam profile, integration windows: all lab. AZURE2 converts to c.m. for every point with that point's
  own entrance pair, so one line may serve segments of different pairs (natural targets).
- **σ_b is a standard deviation, not a FWHM** (FWHM = 2.355 σ_b). Write `0.001` for 1 keV.
- **The straggling coefficient c_s is a c.m. quantity** (keV^½): σ_strag = c_s √(ΔE_c.m. [keV]) keV. A coefficient
  quoted for the beam in the lab enters as c_lab √(m₂/(m₁+m₂)) (0.85 for α + ¹¹B, 0.96 for p + ¹²C).
  Straggling switched on with c_s = 0 computes Bohr's width from Z₁, Z₂ and the line's ε. Bohr's value is
  0.3–0.5 keV^½ for light ions in light targets (p + ¹⁴N at 270 keV: 0.31), is the high-velocity limit (too large
  below about 1 MeV/u) and counts only the target nucleus's electrons. For a compound target with ε per molecule or
  per active atom enter c_s = √(f · 4π Z₁² e⁴ Z_entity / ε) yourself, with ε in MeV cm², e² = 1.44 × 10⁻¹³ MeV cm,
  f = m₂/(m₁+m₂). A coefficient of 0.04 keV^½ is 5–10 times below Bohr, i.e. no straggling.
- **Expressions without `sigmalab` are the legacy convention**: evaluated at E_c.m. and taken as the c.m. width.
  A lab TOF derivation written that way is off by the mass ratio on both axes. Convert to the lab form,
  σ_lab(x) = σ_cm(F x)/F with F = m_t/(m_p+m_t), or write new lines with `sigmalab`.
- **Equation coefficients are bound in ascending order**: the first number is a0, the coefficient of the lowest
  power. Written in descending order the line applies a different resolution without any message (a 20 % wider σ
  changed a fitted Γ_n from 0.61 to 0.44 keV in one case). Evaluate σ_b at one known energy before fitting.

## 3. Compound and implanted targets

The yield uses ε and N only through Δ = εN and 1/(Nε), so ε and N must count the same entity:

| Normalization of N | ε for A_xB_y with A active |
|---|---|
| active atoms cm⁻² ("active density") | ε_A + (y/x) ε_B |
| molecules cm⁻² | x ε_A + y ε_B |
| all atoms cm⁻² | (x ε_A + y ε_B)/(x+y) |

Pure-element ε_A with the active density misses (y/x) ε_B; the average-atom ε with the active density is too small
by (x+y)/x. When a paper
quotes the energy loss in keV, choose ε(E_0) N to reproduce it and skip the stoichiometry; for a 10–20 % target
N ≈ Δ/ε(E_0) holds only to first order (a few per cent).

- builds before 2026-09-28 (expires: 2027-02-01): the GUI's compound ERYA fetch returned the average-atom ε. An
  equation written by such a GUI is consistent only with all atoms cm⁻²; scale it by (x+y)/x for the active density.
- A thin target whose ε evaluates to ≤ 0 at some point is refused with an error naming the segment and energy:
  fix the stopping-power polynomial's range.

## 4. Narrow resonances and sub-point grids

- **A resonance narrower than the kernel needs sub-points on it.** The adaptive grid places a fine lattice around
  each level (`resonanceWidthMultiplier`, default 20 widths). A 200 eV resonance under a 100 eV Gaussian needed
  1000 sub-points for 10⁻⁴ accuracy (4 × 10⁻⁴ with 400). Check convergence by raising the sub-point count until the
  folded curve stops moving.
- **Grid cost against accuracy.** The last two optional numbers of a `<targetInt>` line are resonanceWidthMultiplier
  and pointsPerWidth (defaults 20 and 50). For one narrow level under target integration, measured at its fitted width:
  20/50 113 s (reference), 5/50 37 s (+0.005 % in χ²), 5/20 17 s (+0.03 %), 2/10 6 s (−0.65 %, 6 units of χ², more
  than the Δχ² = 1 interval). Screen at 5/20 and fit and report at 20/50; never use a multiplier of 2. Re-measure per
  project. A width floor on a template buys at most about 2.6× in cost; the number of patterns dominates.
- **Gaussian kernels reach ±5 σ_b** and are normalised by their own integral; a constant σ survives the fold
  exactly. builds before 2026-10-26 (expires: 2027-02-01): kernels stopped at ±3 σ_b, so points 2–3 σ_b from a
  resonance 10³ times above its surroundings came out 2–4 % low; refit such segments on a current build.
- **Channel thresholds.** Close to a charged-particle threshold the bare σ used to overflow (values up to 10³⁷ b
  within about 20 eV). A sub-point grid that lands there gives isolated one-point spikes in a folded curve, at
  energies that move with the sub-point count. builds before 2026-10-11 (expires: 2027-02-01) have this; on a
  current build set `AZURE_DUMP_NONFINITE=1` to locate any remaining non-finite point
  (`docs/source/reference/command_line.rst`).
- builds before 2026-09-11 (expires: 2027-02-01) mis-estimated level widths when building the narrow-resonance
  lattice, which can invert which Γ the data prefer; redo any fit of a narrow level through a `<targetInt>` kernel
  made with them.
- **A long in-process job keeps the extension it loaded.** Rebuilding `pyazr/_azure2*.so` during a campaign makes
  later evaluations incomparable with earlier ledger values; compare the `.so` time stamp with the job start
  before comparing objectives across a rebuild.

## 5. Integration windows (`"lo-hi,..."`, `"auto"`)

Windows restrict the full depth integral of a thin-target line to lab-energy ranges. Outside them a point is
σ(Ē), with Ē = ∫E dE/ε / ∫dE/ε over [E_0−Δ, E_0]; the error is |σ″| Δ²/24 (12.6 % on a 25 keV target where σ
changes by a factor e per 16 keV). Therefore:

- **Cover strongly curved regions, not only narrow resonances**: sub-barrier tails and every particle threshold
  of the compound nucleus (omitting the ¹⁶O(6.05 MeV) n₁ threshold cost 15 % in one case).
- **`"auto"`, `"auto:tol"`, `"auto:tol:margin"`** (defaults 0.01 and 3) place windows from the levels of the
  segment's entrance pair: a level with Γ < Δ/√(3 tol) gets [E_R − w, E_R + Δ + w], w = margin · max(Γ, Δ); Δ
  includes twice the beam spread. The run log prints, per segment, the windows and "I of M points integrated".
  Tighter tolerances integrate more points (one 13C+α model: χ² full 400,616; auto 400,881 with 58 % of points;
  auto:0.003 400,656 with 81 %; auto:0.001:5 400,621 with 99 %).
- **Check a windowed fit against the full integration at its end point** before trusting it: a model can react to
  the σ(Ē) points more than the data do and end worse.
- Windows are refused on thick targets and together with `doppler`.
- builds before 2026-11-06 (expires: 2027-02-01): points outside windows were σ(E_beam), and `"auto"` is not
  understood; rerun any windowed fit made with them.

## 6. Neutron TOF: UDR, transmission, yield, Doppler

- **UDR files.** The n_TOF response files already contain the proton burst: use `burstFWHM 0`; the HI (7 ns) and LI
  (14 ns) pulse types have separate files; check a new file by the width of its highest-energy blocks. There is no
  general EAR2 file: it is made per experiment, so obtain the one the authors' SAMMY fit used. If the data's TOF→E
  conversion fitted the flight path rather than re-centring, use `centred 0` with that L; use `centred 1` only when
  L already includes the mean moderation length. The kernel is evaluated at the point's nominal energy; one channel
  width per line (split a segment where the binning changes, e.g. 10 → 100 bins per decade at 300 keV). A positive
  delay maps to a higher true energy. Details: `references/udr-and-beam-profile.md`.
- **`transmission <n>` and `yield <n>`** (n in atoms/barn) model T = exp(−nσ) and Y₀ = (σ/σ_t)(1 − exp(−nσ_t)),
  applied at every sub-point before the kernel (transforming the broadened σ instead is 5–6 % off at nσ ≈ 1).
  Put transmission on a `total 1` line and fit T, not a σ_tot inverted from it. σ_t in `yield` is the
  optical-theorem total of the entrance pair; capture is not in it. Refused: both keywords on one line,
  non-angle-integrated data, RATIO, charged or non-particle entrance, terms on different sub-point grids; for
  `yield` also natural segments, the elastic/total line, other entrances or elastic exits, and `--use-rmc`.
  Multiple-scattering yields Y₁/Y₂ are not built. pyazr: `Segment.transmission`, `Segment.yield_thickness`.
- builds before 2026-11-04 (expires: 2027-02-01) ignore `transmission`/`yield` without a message and fit plain σ
  to T or Y data; check the build before trusting such a file.
- **`doppler <T>`** (T in K) broadens σ first, then applies `transmission`/`yield`, then the resolution (SAMMY's
  order). Full width Δ_D = √(4 m E kT/M): 63 eV for n + ¹²C at 300 K and 0.46 MeV. It matters for resonances a few
  hundred eV wide or narrower. A line may carry `doppler` alone:
  `1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 293.6`. Refused: `doppler 0`, a missing T, photon or β entrance.
  pyazr: `add_target_effect(..., doppler=T)`. builds before 2026-11-05 (expires: 2027-02-01) refuse the whole
  `<targetInt>` block.

## 7. Fitted apparatus parameters

Each is a parameter with a Gaussian prior about the value written on the line; the fitted value lives only in
`param.sav` (the line keeps its prior centre, `save_fit` does not move it). Parameter index n = place among the
active `<targetInt>` lines.

| Keyword | Parameter | Range and grids | pyazr |
|---|---|---|---|
| `doppler T vary σ [Tmax]` | `doppler_<n>` (K) | grids laid out for Tmax (default max(T+3σ, 1.21T); 2T without prior); MIGRAD keeps T in [T/4, Tmax] | `AzrModel.set_doppler(effect, T=, vary=, sigma=, tmax=)` |
| `transmission n vary σ`, `yield n vary σ`, `thicknessgroup <label>` | `segment_<key>_thickness` (atoms/barn) | group members follow the first segment | `AzrModel.set_sample_thickness(key, n=, vary=, sigma=, group=)` |
| `densityfit vary σ [max]`, `densitygroup <label>` | `density_<n>` (atoms cm⁻²; σ absolute or "2.3%") | grids laid out at max (default N₀+3σ, ≥1.1N₀; 1.5N₀ without prior) | `AzrModel.set_target_density(effect, density=, vary=, sigma=, max=, group=)` |
| `screening U_e vary σ_U` | `screening_<n>` (eV) | factor exp(πηU_e/E_c.m.) on σ at every point and sub-point | `add_target_effect(screening=(Ue, vary, sigma))`, `AzrModel.set_screening` |

- **A fit that ends at Tmax or at the density maximum wants a larger maximum**; beyond max the engine computes at
  max and warns once.
- **The fitted density changes only Δ and the straggling with depth**; the yield is still the depth average, not
  multiplied by N. Refused: no target integration, thick target, `beamprofile`/`udr`, a group without a
  `densityfit` leader.
- **Screening acts on cross sections only** (phase shifts, A_y and ratios are untouched; neutral pairs f = 1). The
  adiabatic limit is B(Z₁+Z₂) − B(Z₁) − B(Z₂) from NIST total binding energies (α + ¹³C 934.7 eV,
  p + ¹²C 442.4 eV). Check the energy frame of screening data: EXFOR entry O1282 tags c.m. energies as lab (EN),
  which, read as lab, puts the metal-target U_e 26–27 % low.
- rmfit runs these kinds (doppler, thickness, density, screening) under the single-session evaluator only
  (`rmfit-campaign`).
- builds before 2026-11-05 (expires: 2027-02-01): a build that knows `doppler T` or `transmission n` but not the fit
  form reads the first number and ignores the rest without a message.
- builds before 2026-10-29 (expires: 2027-02-01): the first evaluation of a session at a shift of exactly 0 used the
  `.azr`'s nonzero shift on segments with target-effect sub-points.

## 8. Every effect also folds the extrapolation segments

A `<targetInt>` line applies to `<segmentsTest>` segments with the same keys. For a bare σ or S factor, deactivate
the line (first field 0) in a copy of the `.azr` before mode 3. To see the folded shape on a fine grid, give a data
segment a dummy fine-grid data file and run mode 1. (Extrapolation procedure: `azure2-calculate`.)

## 9. Energy shifts and sub-points

Sub-points move rigidly with their data point: a constant, `sqrtshift` or `shiftfn` shift s(E_0) is a calibration
of the measured beam energy, so the integration window keeps its width and the resonance-anchored grid stays on the
resonance. Sub-points are built at the parent's c.m. angle.

- builds before 2026-10-26 (expires: 2027-02-01): (a) sub-points of lab-frame differential segments with any
  sub-point effect moved to the wrong angle (by the full lab − c.m. difference) whenever the segment's energy shift
  was active, and the shift derivatives were wrong with it; (b) under `sqrtshift`/`shiftfn` each sub-point moved by
  s(E_sub), stretching the window by 1 + ds/dE. Refit exposed segments.

## 10. Attenuation (Q) coefficients

For a particle exit, Q_K multiplies the lab Legendre coefficients and the result is returned to the point's c.m.
angle; photon, UPOS, A_y and charged-particle elastic data keep Q as entered. When g = v_c.m./v₃* is large, give
coefficients beyond the highest c.m. order. builds before 2026-10-25 (expires: 2027-02-01) applied Q per c.m. order
for particle exits, a few per cent off, and changed results even with every Q = 1.

## 11. Pitfalls by symptom

- **A shared line under-smoothed every segment it names**: builds before 2026-10-01 (expires: 2027-02-01) converted
  a shared σ_b lab → c.m. once per segment named, so a line naming N segments (or a convolved SUM/RATIO) applied
  σ_b (m_t/(m_p+m_t))^N. With such a build, write one line per segment.
- **A thick target with beam spread or straggling gives zero at its first point**:
  builds before 2026-10-27 (expires: 2027-02-01). On those builds a thin target whose ε went ≤ 0 also silently
  became thick.
- **Yields of thick (> 5 %) targets at low energy differ by 2–8 % between builds**:
  builds before 2026-10-27 (expires: 2027-02-01) used Δ = ε(E_0)N instead of the range integral; a flat σ came out 2 % low for a 5 % target,
  8 % for 20 %, more on a resonance tail.
- **A slow evaluation dominated by Whittaker functions** (`gstack` shows `RecalcEDependentValues → theWhitFunc`):
  builds before 2026-09-13 (expires: 2027-02-01) recomputed sub-point Coulomb functions every evaluation when
  external capture was on. To find such costs on any build, switch data groups off in copies of the model and time
  `calculate_chi2_rwa`, then sample stacks.
- **The Qt GUI wrote beam-profile numbers with 12 significant digits, pyazr writes 17**: harmless, but visible in a diff.

## Pointers

- `docs/source/user_guide/experimental_effects.rst`: every kernel, its equation and its dialog.
- `docs/source/user_guide/segments.rst`: `transmission`, `yield`, `total 1`.
- `docs/source/reference/data_formats.rst`: the `<targetInt>` line token by token.
- `docs/source/developer/edata_and_target_effects.rst`: sub-point grids, windows, adaptive lattice.
- `references/udr-and-beam-profile.md`: the beam-profile and UDR kernels in detail.
- Tests that show each effect against an independent result: `tests/convolution_window`, `tests/thick_target`,
  `tests/straggling_bohr`, `tests/target_thickness_exact`, `tests/shared_effect_pairs`, `tests/transmission_yield`,
  `tests/doppler_fgm`, `tests/electron_screening`, `tests/screening_6li`, `tests/target_density`,
  `tests/target_window_fallback`, `tests/target_auto_windows`, `tests/beam_profile_kernel`,
  `tests/q_attenuation`, `tests/shift_subpoint_angle`, `tests/shift_sqrt_subpoints`, `tests/target_effect_shared`.
