---
name: azure2-targets
description: Model the apparatus between the AZURE2 cross section and the measured yield - target integration (thin, thick, compound, implanted, gas), SRIM stopping (the stopping library) and straggling, Gaussian beam spread and TOF resolution, SAMMY UDR, beam-profile kernels, free-gas Doppler, transmission and self-shielded yield, electron screening, fitted density or thickness, Q attenuation coefficients. Use for any target, beam or detector record or experimental effect. Building the project: azure2-setup; fitting: azure2-fit; stopping-profile inference from excitation curves: target-analysis.
---

# Targets, resolution and other experimental effects in AZURE2

An experimental effect is a record of `<targets>` -- a `target` (energy loss, straggling, Doppler, screening, sample
thickness), a `beam` (energy spread, profile) or a `detector` (resolution, UDR, attenuation) -- that a data set or an
extrapolation names in its `target`, `beam` or `detector` field; it folds the calculated c.m. cross section σ(E_c.m.) of
that data set with a kernel: the energy loss in the target, the beam energy spread, the detector or time-of-flight
resolution, the thermal motion of the target atoms, the sample attenuation, or the screening of the nuclear charge. This
skill says how to choose and write those records and how to check them. What each kernel computes, with its equations,
is in `docs/source/user_guide/experimental_effects.rst`; the fields are in `docs/source/reference/format3.rst`
(generated from the table; `AZURE2 --describe` shows what a data set got); the implementation (sub-point grids, adaptive
lattices) is in `docs/source/developer/edata_and_target_effects.rst`.

## When to use, and when a sibling skill is the owner

- Use this skill to decide which effect a data set needs, to write or audit its target, beam or detector, and to fit an
  apparatus parameter (target density, sample thickness, Doppler temperature, screening potential U_e).
- Project construction, masses, frames of the data file: `azure2-setup`.
- Running calculations and extrapolations: `azure2-calculate`. Freeing, seeding and staging parameters, energy
  shifts and normalizations: `azure2-fit`.
- Inferring a stopping profile or target composition from the excitation curves themselves: the user-level skill
  `target-analysis`.

## Scripts

| Script | Does |
|---|---|
| `scripts/layer.py` | SRIM stopping of a material, the loss of a layer (range integral and first order), straggling, the record |
| `scripts/target_check.py` | the engine's thin-target fold of a data set against an integral outside it |
| `scripts/selftest.py` | the loss against a quadrature of 1/ε, a round trip, the record and the fold checked in T1 |

```bash
python3 .claude/skills/azure2-targets/scripts/layer.py --material TiN --projectile p --energy 0.4 \
    --areal-density 2e18 --active 15N        # 27.93 keV by the range integral, 2 % above eps(E0) N
python3 .claude/skills/azure2-targets/scripts/layer.py --material Ta2O5 --projectile p --energy 0.3 --loss-kev 10
```

## 1. Choosing the kernel

| Measurement | Record and field | Docs section |
|---|---|---|
| Charged beam on a solid target, ΔE comparable with a resonance width Γ | `target` of `kind layer` with `areal_density` N and `stopping` | "Target Integration" |
| Yield per incident particle from a target thicker than the range | `target` of `kind infinite` with `stopping` | "Target Integration" |
| Yield per particle or per µC from a layer of finite thickness (a resonance scan) | `quantity yield` on a `target` of `kind layer` | `docs/source/reference/data_formats.rst`, "Yields" |
| Beam energy spread or detector resolution, Gaussian | `beam` field `spread σ_b MeV lab` (or a `detector`'s `energy_resolution`) | "Gaussian Beam Spread" |
| Neutron TOF with an energy-dependent width | `beam` field `spread formula "<expression in x>" a0=... MeV lab` | "Gaussian Beam Spread" |
| Neutron TOF with a tabulated response (n_TOF) | `detector` field `udr file F flight_path L m burst B ns channel C ns centred yes` | "User-defined numerical resolution function" |
| Broad, skewed γ beam with event-by-event energy (HIγS TPC) | `beam` field `profile skewed_gaussian ...` | "Beam-profile kernel" |
| Narrow resonances (Γ of order 100 eV) in a room-temperature sample | `target` field `doppler free_gas T K fixed` | "Free-Gas Doppler Broadening" |
| Neutron transmission T = exp(−nσ) | `quantity transmission` on a `target` of `kind sample` (`areal_density n atoms/b`) | `docs/source/user_guide/segments.rst` |
| Capture or reaction yield of a thick sample | `quantity yield` on a `target` of `kind sample` | `docs/source/user_guide/segments.rst` |
| Electron screening of a low-energy charged-particle reaction | `target` field `screening U_e eV fit prior U_e σ_U` | "Electron Screening" |
| Angular attenuation of a finite detector | `detector` field `attenuation Q1 Q2 ...` | "Attenuation Coefficients" |

A thin-target point is the depth average Y(E_0) = (1/N) ∫_{E_0−Δ}^{E_0} σ(E′)/ε(E′) dE′, with N the areal density
(atoms cm⁻²), ε the stopping cross section and Δ the energy loss defined by ∫ dE′/ε = N (exact range integral, not
ε(E_0)N). A thick-target point (density 0) is ∫_0^{E_0} σ/ε dE without division by N, integrated down to 1 keV.
A layer yield (format 3, `quantity yield` on `kind layer`) is 10⁻²⁴ ∫_{E_out}^{E_0} σ/ε dE, not divided by N: it is the
data a resonance scan gives; N (fitted or fixed) sets the width of a narrow resonance's plateau, not its height.

**Which kernel the paper describes.**
- **A foil thickness (μg cm⁻²) or "the energy loss in the target was X keV" means target integration**, with an areal
  density and the stopping-power formula, not a Gaussian. Add straggling with the target's `straggling` field.
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
  one target record per energy window, each named by the data sets of its window.
- **Before calling a target-integration discrepancy physics**, run `scripts/target_check.py model.azr --dataset NAME`:
  it recomputes the depth average outside the engine (bare σ on a fine grid, SRIM's ε, the range integral) and
  compares point by point (agreement 1e-4 through a 33 keV resonance under an 8 keV loss in the self-test).

## 2. Stopping powers: SRIM through AZURE2's stopping library

A layer or infinite target needs a `stopping` field (`docs/source/reference/format3.rst`, "target"):

- `stopping evaluation` with `composition Ti 1 N 1` (and `phase gas` for a gaseous compound, `compound_correction`
  for SRIM's factor C of Bragg's rule): SRIM-2013's stopping, evaluated by the engine per atom of the active
  nucleus's element (ε of the material over that element's atom fraction). This is the default choice; every
  result that uses it cites SRIM (`pyazr.stopping.CITATION`, printed by `layer.py`).
- `stopping table <file> MeV eV/1e15atoms/cm2` (or `MeV*cm2/g`): a measured or published table.
- `stopping formula "<expression in x>" a0=... <unit>`: a fitted polynomial from a paper; check its range covers
  every energy the beam has in the layer (a stopping that is not positive there is refused, `E_LAYER_STOPPING`).
- `... scale 1 fit prior 1 4%` on any of them: a stopping-scale nuisance (parameter kind `apparatus`), the usual way
  to carry SRIM's uncertainty (`pyazr.stopping.stopping(...)["rel_uncertainty"]` gives SRIM's quoted accuracy for
  the ion). Two targets of one material and projectile share the scale and must write it alike.
- `straggling none | bohr | coefficient <c> keV^1/2 cm`; `epoch <name> areal_density ...` for a target whose
  thickness changed between runs, a data set naming `target T@<epoch>`.
- A data set that names a target, beam or detector must state its `quadrature points N` (§5).
- Python without a session: `pyazr.stopping` (`Material`, `stopping`, `energy_loss`, `energy_after`,
  `areal_density_between`, `range_of`, `straggling`, `provenance`); its units are per atom of the material and
  atoms/cm² of the material. Inferring a profile from the excitation curves themselves: the user-level
  target-analysis skill.

## 3. Frames and units on target, beam and detector records

- **Every number carries its unit and frame.** A spread or resolution says `lab` or `cm`; σ_b, a `formula` (a
  function of x = E_lab in MeV with `lab`) and the beam profile are entrance-channel lab energies.
  AZURE2 converts to c.m. for every point with that point's own entrance pair, so one record may serve data sets of
  different pairs (natural targets).
- **σ_b is a standard deviation, not a FWHM** (FWHM = 2.355 σ_b). Write `spread 1 keV lab` for 1 keV.
- **The straggling coefficient c_s is a c.m. quantity** (keV^½): σ_strag = c_s √(ΔE_c.m. [keV]) keV. A coefficient
  quoted for the beam in the lab enters as c_lab √(m₂/(m₁+m₂)) (0.85 for α + ¹¹B, 0.96 for p + ¹²C).
  `straggling bohr` computes Bohr's width from Z₁, Z₂ and the target's ε. Bohr's value is
  0.3–0.5 keV^½ for light ions in light targets (p + ¹⁴N at 270 keV: 0.31), is the high-velocity limit (too large
  below about 1 MeV/u) and counts only the target nucleus's electrons. For a compound target with ε per molecule or
  per active atom enter c_s = √(f · 4π Z₁² e⁴ Z_entity / ε) yourself, with ε in MeV cm², e² = 1.44 × 10⁻¹³ MeV cm,
  f = m₂/(m₁+m₂). A coefficient of 0.04 keV^½ is 5–10 times below Bohr, i.e. no straggling.
- **A `formula` with `cm` is evaluated at E_c.m. and taken as the c.m. width** (the converter writes an old
  line's c.m. width so). A lab TOF derivation written that way is off by the mass ratio on both axes:
  use `lab`, or σ_lab(x) = σ_cm(F x)/F with F = m_t/(m_p+m_t).
- **Formula coefficients are named** (`a0=... a1=...`): check that a0 is the coefficient of the lowest power.
  Swapped, the record applies a different resolution without any message (a 20 % wider σ
  changed a fitted Γ_n from 0.61 to 0.44 keV in one case). Evaluate σ_b at one known energy before fitting.

## 4. Compound and implanted targets

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

- A layer whose ε is not positive at an energy of its sub-point grid is refused (`E_LAYER_STOPPING`): fix the
  formula's range, or use `stopping evaluation`.

## 5. Narrow resonances and sub-point grids

- **A resonance narrower than the kernel needs sub-points on it.** The adaptive grid places a fine lattice around each
  level (the data set's `quadrature ... width_multiplier`, default 20 widths). A 200 eV resonance under a 100 eV
  Gaussian needed 1000 sub-points for 10⁻⁴ accuracy (4 × 10⁻⁴ with 400). Check convergence by raising the sub-point
  count until the folded curve stops moving.
- **Grid cost against accuracy.** The data set's field `quadrature points N width_multiplier W per_width P` sets
  the grid (W and P default 20 and 50). For one narrow level under target integration, measured at its fitted width:
  20/50 113 s (reference), 5/50 37 s (+0.005 % in χ²), 5/20 17 s (+0.03 %), 2/10 6 s (−0.65 %, 6 units of χ², more
  than the Δχ² = 1 interval). Screen at 5/20 and fit and report at 20/50; never use a multiplier of 2. Re-measure per
  project. A width floor on a template buys at most about 2.6× in cost; the number of patterns dominates.
- **Gaussian kernels reach ±5 σ_b** and are normalised by their own integral; a constant σ survives the fold
  exactly.
- **Channel thresholds.** Close to a charged-particle threshold the bare σ used to overflow (values up to 10³⁷ b
  within about 20 eV). A sub-point grid that lands there gives isolated one-point spikes in a folded curve, at
  energies that move with the sub-point count. Set `AZURE_DUMP_NONFINITE=1` to locate any remaining non-finite point
  (`docs/source/reference/command_line.rst`).
- **A long in-process job keeps the extension it loaded.** Rebuilding `pyazr/_azure2*.so` during a campaign makes
  later evaluations incomparable with earlier ledger values; compare the `.so` time stamp with the job start
  before comparing objectives across a rebuild.

## 6. Every point of a folded data set is folded

There are no integration windows, edge blends, automatic tolerance or σ(Ē) stand-in any more: every point of a data
set that names a target, beam or detector is folded with the whole kernel on the shared mesh of its compound nucleus.
A file that writes `fold_windows` is refused (`E_FIELD_DELETED`: remove the field); the converter of old projects
notes what it dropped. A windowed fit made with an older
build differs from today's model by at most what its windows left out (outside the windows of a beam spread, the whole
spread): refit it. The cost of a point is set by the mesh's settings (`quadrature`), not by which points are folded.

## 7. Neutron TOF: UDR, transmission, yield, Doppler

- **UDR files.** The n_TOF response files already contain the proton burst: use `burst 0 ns`; the HI (7 ns) and LI
  (14 ns) pulse types have separate files; check a new file by the width of its highest-energy blocks. There is no
  general EAR2 file: it is made per experiment, so obtain the one the authors' SAMMY fit used. If the data's TOF→E
  conversion fitted the flight path rather than re-centring, use `centred no` with that L; use `centred yes` only when
  L already includes the mean moderation length. The kernel is evaluated at the point's nominal energy; one channel
  width per detector (split a data set where the binning changes, e.g. 10 → 100 bins per decade at 300 keV). A positive
  delay maps to a higher true energy. Details: `references/udr-and-beam-profile.md`.
- **`quantity transmission` and `quantity yield` on a `kind sample` target** (`areal_density n atoms/b`) model
  T = exp(−nσ) and Y₀ = (σ/σ_t)(1 − exp(−nσ_t)), applied at every sub-point before the kernel (transforming the
  broadened σ instead is 5–6 % off at nσ ≈ 1). Fit T, not a σ_tot inverted from it. σ_t in a yield is the
  optical-theorem total of the entrance pair; capture is not in it. Refused: a differential quantity, a ratio, a
  charged or non-particle entrance, terms on different sub-point grids, a cross section on a sample target
  (`E_TARGET_KIND`); for a yield also natural samples, other entrances or elastic exits, and `--use-rmc`.
  Multiple-scattering yields Y₁/Y₂ are not built.
- **`doppler free_gas T K fixed`** on the target broadens σ first, then applies the transmission or yield, then the
  resolution (SAMMY's order). Full width Δ_D = √(4 m E kT/M): 63 eV for n + ¹²C at 300 K and 0.46 MeV. It matters
  for resonances a few hundred eV wide or narrower. A target may carry `doppler` alone. Refused: T ≤ 0, a photon or
  β entrance.

## 8. Fitted apparatus parameters

Each is a parameter of the target record, named after it, with a Gaussian prior about the value the field writes;
the fitted value lives in `param.sav` (`pyazr.Session.save_fit` writes it back into the field). One target named by
several data sets is one parameter for all of them (what the old `densitygroup` / `thicknessgroup` did). Edit the
fields with `pyazr.Project.set`.

| Field of the target | Parameter | Range and grids |
|---|---|---|
| `doppler free_gas T K fit [prior T σ] [grid_max Tmax]` | `<target>.doppler` (K) | grids laid out for Tmax (default max(T+3σ, 1.21T); 2T without prior); MIGRAD keeps T in [T/4, Tmax] |
| `areal_density n atoms/b fit [prior n σ]` on a `kind sample` target | `<target>.areal_density` (atoms/barn) | shared by the data sets on the target |
| `areal_density N atoms/cm2 fit [prior N σ] [grid_max Nmax]` on a `kind layer` target | `<target>.areal_density` (atoms cm⁻²; σ absolute or "2.3%") | grids laid out at grid_max (default N₀+3σ, ≥1.1N₀; 1.5N₀ without prior) |
| `screening U_e eV fit prior U_e σ_U` | `<target>.screening` (eV) | factor exp(πηU_e/E_c.m.) on σ at every point and sub-point |

- **A fit that ends at Tmax or at the density's grid_max wants a larger grid_max**; beyond it the engine computes at
  grid_max and warns once, naming the parameter.
- **The fitted density changes only Δ and the straggling with depth**; the yield is still the depth average, not
  multiplied by N. Refused: an infinite target, a beam profile or a UDR on the same data set.
- **Screening acts on cross sections only** (phase shifts, A_y and ratios are untouched; neutral pairs f = 1). The
  adiabatic limit is B(Z₁+Z₂) − B(Z₁) − B(Z₂) from NIST total binding energies (α + ¹³C 934.7 eV,
  p + ¹²C 442.4 eV). Check the energy frame of screening data: EXFOR entry O1282 tags c.m. energies as lab (EN),
  which, read as lab, puts the metal-target U_e 26–27 % low.
- rmfit runs these kinds (doppler, thickness, density, screening) under the single-session evaluator only
  (`rmfit-campaign`).

## 9. Extrapolations are folded only by the records they name

An extrapolation record is folded by the `target`, `beam` or `detector` it names, and by nothing else (the converter of
an old project writes an old line's folding as the extrapolation's own field, with a `# CONVERTED:` note). For a bare σ
or S factor, remove the field (`pyazr.Project.set` with `None`). To see the folded shape on a fine grid, name the target
in an extrapolation. (Extrapolation procedure: `azure2-calculate`.)

## 10. Energy shifts and sub-points

Sub-points move rigidly with their data point: a constant, `sqrt` or function shift s(E_0) (the data set's `shift` or
its `energy_scale`) is a calibration of the measured beam energy, so the integration window keeps its width and the
resonance-anchored grid stays on the resonance. Sub-points are built at the parent's c.m. angle.


## 11. Attenuation (Q) coefficients

For a particle exit, Q_K multiplies the lab Legendre coefficients and the result is returned to the point's c.m.
angle; photon, UPOS, A_y and charged-particle elastic data keep Q as entered. When g = v_c.m./v₃* is large, give
coefficients beyond the highest c.m. order.

## Pointers

- `docs/source/user_guide/experimental_effects.rst`: every kernel, its equation and its dialog.
- `docs/source/reference/format3.rst`: every field of target, beam, detector and data set, with its grammar.
- `docs/source/user_guide/older_projects.rst`: converting an old project's target lines with `AZURE2 --convert`.
- `docs/source/developer/edata_and_target_effects.rst`: sub-point grids, windows, adaptive lattice.
- `references/udr-and-beam-profile.md`: the beam-profile and UDR kernels in detail.
- Tests that show each effect against an independent result: `tests/convolution_window`, `tests/thick_target`,
  `tests/straggling_bohr`, `tests/target_thickness_exact`, `tests/shared_effect_pairs`, `tests/transmission_yield`,
  `tests/doppler_fgm`, `tests/electron_screening`, `tests/screening_6li`, `tests/target_density`,
  `tests/beam_profile_kernel`,
  `tests/q_attenuation`, `tests/shift_subpoint_angle`, `tests/shift_sqrt_subpoints`, `tests/target_effect_shared`.
