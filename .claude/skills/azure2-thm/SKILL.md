---
name: azure2-thm
description: R-matrix analysis of Trojan Horse Method (THM) data with AZURE2, with the emphasis on assessing and quantifying the model (theory) dependence of THM results. Covers establishing what published THM points actually are (HOES vs penetrability-corrected, carrier, KF convention, resolution), building the HOES model (<thm> options, binding-energy field, identical particles), joint fits with direct data, the model-dependence protocol and its "theory dependence index", reproducing a published THM analysis, and the pitfalls met in 7Li, 6Li, 15N, 17O, 18O and 12C+12C. Use whenever the task involves THM/HOES data, a <thm> block, modified R-matrix, THM strengths or S factors, or comparing THM results with direct data or literature. Generic AZURE2 usage (running, fitting, file anatomy, the THM segment mechanics) is in azure2-eval.
---

# THM analysis with AZURE2, and how model-dependent the answer is

This skill is the procedure, not the mechanics. How a THM segment is declared
(observable code +10, field 32 binding energy, field 33 `gammaIsRWA`, the `<thm>`
keys, `AzrModel.set_thm_option`) is in **`azure2-eval`**, section "Trojan Horse
(THM) segments"; the equations behind each option are in
`docs/source/theory/thm_implementation.rst`. Level choice and Brune/formal
questions: **`r-matrix-analysis`**. Fetching data: **`nds-explorer`**. New
compound nucleus: **`azr-project-builder`**.
GUI: all of it (`<thm>` options, experiments, fields 32/33) is edited in *Configure > THM Workspace...*, not in the classic tabs.

Calibration comes from seven cases, all in `examples/` (`.azr` holds the best
fit, `data/` headers name the sources): `li7_tumino2006`, `li6_pizzone2011`,
`n15_lacognata2007`, `o17_guardo2017_fit`, `o18_lacognata2008` (narrow 20/90/144 keV),
`o18_lacognata2010` (1/2⁺ doublet), `c12c12_tumino2018`. Regression pins:
`tests/18O_p_a_thm`, `tests/7Li_p_a`, `tests/6Li_d`, `tests/17O`, `tests/thm_options`.

The one-sentence lesson: **THM data fix shapes, energies and ratios of γ_x²|M_l|²; the
absolute widths and the conversion to on-shell strengths depend on the reaction
model, by ≤ 20 % when overlapping direct data anchor the widths and by factors
1.5–200 when strengths come from HOES peak areas.** Quote that dependence with every
number.

## 1. Establish what the published points are — before any modelling

Published "THM data" are rarely the raw HOES excitation function. Settle, per
dataset, and write it down:

| question | what we found | consequence |
|---|---|---|
| Raw HOES, or penetrability-corrected / normalised S(E) or σ? | 7Li (Tumino 2006, EXFOR O1653002) and 6Li (Pizzone 2011, D0649002): penetrability divided out, normalised to direct data → on-shell-equivalent. 15N Table 3 (C1788004): deconvolved S(E), already normalised (fitted scale 1.02). 18O ApJ 723: PWIA S factor σ_HOES·P₀e^{2πη} + linear background. 18O PRL 2008, 17O Guardo 2017, 12C+12C Nature 2018: HOES, arbitrary units | on-shell-equivalent points are fitted with the **on-shell** observable and a free scale, not as a THM segment |
| Carrier nucleus and B (field 32) | 7Li and 6Li both used ³He (B = 5.49 MeV), while `tests/7Li_p_a` / `tests/6Li_d` carry 2.2246 / 1.4735 MeV. d carrier 2.224566; ¹⁴N → ¹²C+d 10.2723 | wrong B shifts p = √(2μ(E+B))/ħ and every vertex node |
| KF convention (what was divided out) | full three-body density (Typel & Baur eq. 16, Tumino 2021 eq. 15) → `kf3body`; λ₃/λ₂ with on-shell p_Ax (Pizzone 2011, Tumino 2021 eq. 22) → `lambda32`; only |φ|² → `triple`; La Cognata's working formula or unknown → `lacognata` | only the energy dependence matters (norm is free); `lambda32` leaves 1/k_i, a factor 7 over 0.02–1 MeV |
| Resolution σ (and σ vs FWHM, frame) | 18O: 40 keV FWHM = 17 keV σ (ApJ 708 says σ). 17O: σ 20/30 keV (Gulino) vs 35 keV FWHM (EXFOR); fitted 22.2 ± 1.2. 7Li: "80–120 keV", convention not stated. 12C+12C: σ 30 keV (scan minimum 30–35) | fit σ as a check; `<targetInt>` σ is **lab**: ×(m₁+m₂)/m₂ (0.017 → 0.017948 MeV for p+¹⁸O; 0.030 → 0.060 for ¹²C+¹²C) |
| Angular window | 7Li: θ_cm 50–70°; 18O: 4π via direct angular distributions | incoherent entrance-l sum is exact only for 4π (below) |
| Energy offset, background subtracted? | 18O narrow: per-level offsets −1…+9 keV; ApJ 723: linear background (digitise the dashed line); 17O energies read as neutron-lab (χ² 13.0 vs 17.0 as c.m.) | fit an energy shift; test lab vs c.m. reading of the energy column |
| Channel radius and boundary | usually **not stated** (table below) | radius becomes a model variation, not an input |

What the papers state (open copies only):

| paper | radius | boundary / vertex | KF | resolution |
|---|---|---|---|---|
| La Cognata PRL 2008, ApJ 708 (18O narrow) | symbol only ("5 % from the radius", double ratio) | Breit–Wigner with M_i(E) | three-body (Dolinsky), not written | 40 keV FWHM = σ 17 |
| La Cognata ApJ 723 (18O doublet) | 5.1 fm p, 5.7 fm α (from Yagi/Mak) | B = S(E₁) labels the representation; **no M_l**, a fitted complex ratio L21 = (γ₁/γ₂)m21e^{iφ21} | not written | 17 keV |
| La Cognata 2006/07 (15N) | symbol only | PWIA + penetrability | "KF", not written | not quoted |
| Gulino 2013 (17O) | 4.1 fm for the penetrability only | outgoing-wave (on-shell) Wronskian | MC KF|φ|² | σ 20/30 keV |
| Tumino 2006 (7Li) | none | PWIA + penetrability | eq. 4 (p_c³ misprint) | 80–120 keV |
| Pizzone 2011 (6Li) | none | PWIA | λ₃/λ₂, on-shell p_ax | none |
| Tumino 2018 + Reply (12C+12C) | none (¹²C+d 3 fm only) | per level (Tumino 2021 eq. 51) | not stated | 30 keV σ |

No primary analysis includes the external Coulomb term C_l. No paper tabulates its
fitted curve.

**Test provenance by fitting.** Fit the same model two ways — the points as a HOES
THM segment, and as an on-shell observable with a free scale — and compare THM χ²/N,
the direct norms, and the residual pattern near threshold. 7Li: on-shell χ²/N
**3.5** against HOES **12.8** (12.5 with B = 5.49 MeV); direct norms return from
1.25–1.5 to 0.94–0.98; HOES residuals of −10 to −14σ below 0.4 MeV (HOES rises
toward threshold where the data fall: penetrability already divided out). 6Li: the
free on-shell scale comes out **1.11**, i.e. the points are absolute. 15N: HOES ×
P₀e^{2πη} gives 2.85 against on-shell 1.25. A χ² that does not discriminate (6Li
0.98 vs 1.32) leaves the question to the paper text.

## 2. Building the model

- **Level scheme** from direct R-matrix analyses and ENSDF (`nds-explorer`), not from
  the THM paper alone. Check the windows: the 18O 3/2⁻ 597.6 keV (Γ = 2 keV) was cut
  out by an old excluded window although it dominates the direct data locally.
- **Parameters.** AZURE2 inputs are Brune (observed) parameters and the CLI runs
  Brune by default (no switch to turn it off; `pyazr.azure2(..., use_brune=False)`
  does). Convert a published formal set (ApJ 723 Table 3) by solving the Brune
  eigenproblem and enter the amplitudes with `gammaIsRWA` = 1; check the on-shell σ
  against an independent formal calculation (agreed to 1e-5). Brune fails when
  γ²dS/dE is large (12C+12C θ² = 26 at 6.41 fm is not representable): use a larger
  radius or smaller seeds.
- **THM segment.** Observable 10 (angle-integrated) or 11; B in field 32 of every
  entrance-pair channel line; one `<targetInt>` Gaussian per THM segment (60 uniform
  sub-points for 12C+12C; check convergence, section 6); THM norm free → AZURE2
  profiles it analytically (n* = S_mm/S_md, no penalty). Several channels of one
  experiment share **one** scale (12C+12C: free per-channel scales came out 0.96–1.07). Declare
  such segments one THM experiment, `experiment[<name>] segments=1,2,5-7` in `<thm>`:
  one shared profiled norm (n* from the summed S sums), optional
  `background=const|linear|quadratic` added to the folded model (e.g. the linear
  background of ApJ 723), both closed-form linear least squares; coefficients and
  covariance in `output/thm_experiments.out` and `session.thm_background(name)`.
  `AzrModel.set_thm_experiment(...)`; `theta` is reserved (refused).
  Details: thm_implementation.rst, "THM experiments".
- **Spectator-momentum window**: `ps=hulthen:0-40` (or `hulthen:a,b:..`, `gauss:FWHM:..`,
  `table:file`; `psNodes=`, needs kinematics, excludes `spectatorEnergy`) averages the HOES
  cross section (incoherently) over p_s with weight |φ|²p², T_s = p²/2μ_sx added to E + B at
  each node; fills vertex nodes, changes 12C+12C by 0.3–1.7× (χ² 61 → 1578 unrefitted), 7Li
  ~3 % in shape below 1 MeV. `session.thm_vertex(name, E)`; thm_implementation.rst,
  "Spectator-momentum window".
- **Coulomb line shape** (charged spectator): `lineshape=on` on the experiment line
  (needs `beam/target/spectator/Ebeam`, Brune) multiplies each level's exit
  amplitude, inside the coherent sum, by N_C = e^{πζ/2}(E_λ−E−iΓ_λ/2)^{−iζ},
  |N_C|² = exp[2ζ arctan(2(E_λ−E)/Γ_λ)] (Mukhamedzhanov 2020 eqs. 56–62), with
  ζ = η_sB − η_0 = Z_s α(Z_B μ_sB − Z_F μ_sF)/k_sF, E_sF = E_aA − B − E (paper's case 2,
  η_sb dropped; ζ < 0 → peaks move up). Off by default and 1 for a neutron
  spectator; published fits assume N_C = 1, so turning it on changes the shape a lot
  for 12C(14N,d) (ζ ≈ −0.13…−0.49). `session.thm_lineshape(name, E)` gives ζ, E_sF,
  η_sb (validity: ≪ 1) and |N_C|² per level. Details: thm_implementation.rst,
  "Coulomb line shape".
- **Distortion factor R(E)**: `distortion=coulomb|optical|table:file` on the experiment line (kinematics needed; `opticalAA/SF=plane|coulomb|V,R,a,W,RW,aW,WD,RD,aD,RC`, `spectatorAngle=qf|<lab deg>|cm:<deg>`, `distortionRatio=dwpw|dw`, `boundState=whittaker|yukawa[:rmin]`, `distortionRef=`) multiplies the model by the zero-range DWBA |M/M_PW|² ratio (= dividing PWA S* by R); reproduces Mukhamedzhanov's 2019 12C+12C curve with defaults, the 2026 one with `dw` and E_sF +50 keV; 12C+12C χ² 61 → 3247 unrefitted, 18O(d,n) R = 0.93–1.08. `session.thm_distortion(name, E)`; thm_implementation.rst, "Distortion factor R(E)".
- **Recommended `<thm>` defaults** and why:
  - `entranceL=incoherent`: exact for a 4π-integrated, spin-summed observable. In a
    restricted window the l ≠ l′ cross term survives (7Li 50–70°: +77 % for s = 1,
    −9 % for s = 2); regression pin 2138 vs 3197 coherent.
  - `vertex=constant`: B = S(E₁) of the lowest level of each Jπ, factors out of the
    level sum, representation-invariant. `perlevel` mixes representations once two
    levels of one Jπ interfere (18O band 23 % vs 5 % peak rms; 7Li χ²/N 19.0 vs 12.8;
    15N HOES 10.3 vs 2.85) — use only to reproduce pre-Sep-2026 numbers. `onshell`
    (S + iP) is the physically cleanest and a mandatory variation.
  - `kinematics=` per the KF convention of section 1.
  - `coulombIntegral=0` for comparison with published analyses (none uses it); `1`
    is physically part of the amplitude (|C_l| = 20–45 % of the surface term) and is a
    variation; costs ×3–4 run time; nothing for neutrons.
- **Identical particles.**
  - Two identical 0⁺ bosons form only even-J natural-parity levels. AZURE2 applies
    (1+δ) = 2 to every cross section out of an identical *entrance* pair (df8c2c8;
    older builds were ×0.5 low). An identical *exit* pair (α+α) gets no factor: data
    that count both particles (4π α yield: Mani, Jeronymo at 2.03×, 2.1× the model)
    need ×½.
  - THM-only levels (odd J or unnatural parity, populated by transfer but unable to
    fuse: 12C+12C 0.877 MeV 1⁻, which carries the largest published S* peak). AZURE2
    keeps l = J with a Bose-symmetry warning and cannot average l = J±1. Current
    workaround: a second, identical entrance pair for the on-shell segments; odd-J
    levels couple to the THM pair only, allowed levels to both with the same width.
    Cost: the second channel enters the level matrix and Brune denominator twice
    (χ² 110.69 → 111.65 for levels near the Wigner limit). Exact alternative: two
    sessions, forbidden amplitudes zeroed for on-shell observables.

## 3. Joint fits with direct data

THM alone fixes shapes and ratios, not absolute widths: THM-only fits reach low χ²
with widths wrong by large factors on shell (7Li stage B: THM χ²/N 3.1, Rolfs χ²/N 186;
18O doublet THM-only: Γ_p → 0; 15N: widths unconstrained). Only joint fits give
anchored numbers.

- **Loop.** Own scipy `least_squares` (TRF) over pyazr `residuals()`; pyazr does not
  own fitting. THM scale profiled; each direct norm profiled or fitted with a Gaussian
  penalty of its quoted systematic, so the absolute scale comes **only** from direct
  data. THM Jacobian rows are finite differences (engine-side, parallel since
  dcad73e); if slow, forward-difference `residuals()` yourself.
- **Stages.** Direct only → THM only (to see what it wants) → joint with the few
  parameters the data constrain → free more. Fixed widths in `<levels>` are fixed as
  amplitudes, not physical widths (they drift through the Brune denominator): keep all
  free in the engine and mask in Python; carry the amplitude vector between stages.
- **Priors** (12C+12C, 200+ parameters): Gaussian in ln Γ_c, σ = 1 (factor 2.7)
  around the table; E σ = 30 keV; penalty ln(θ²)/0.5 for θ² > 1; box bounds ±40 keV.
  A width whose posterior error is not below half the prior width is **not measured**
  — mark it (only 23 of 136 exit widths passed).
- **Checks.** Per-dataset χ²/N (not only the total); direct norms near 1 (norms of
  1.25–1.5 signal a wrong THM observable); θ² against 3ħ²/(2μa²) (use `pyazr.widths`);
  parameters at bounds; super-Wigner high-l entrance amplitudes used as HOES shape
  knobs (7Li p+⁷Li l = 3 θ² = 6.5 — no entrance penetrability in HOES, so on-shell data
  cannot veto it); the ½ for identical exit particles; a χ² profile of the THM scale
  (12C+12C: factor 10 costs Δχ² = 28 from the direct data; the Wigner limit caps it
  from above).
- **Errors** from JᵀJ, × √(χ²/ν) when > 1; derived Γ, ωγ, θ² by linear propagation
  through `transform_rwa` (signed: keep the sign). Fit errors exclude the model
  dependence of section 4, which is often larger by 1–3 orders of magnitude.

## 4. Model-dependence protocol (the core)

For every quantity you will quote (ωγ, Γ_x, S(0), a ratio of strengths, Σωγ below
some E), vary one ingredient at a time from the adopted fit:

| variation | values | what to redo |
|---|---|---|
| `vertex` | constant / onshell / perlevel | refit; perlevel only as a check (representation-dependent) |
| `kinematics` | lacognata / triple / kf3body (+ lambda32 only if the data were divided by λ₃/λ₂; list it, exclude from the range otherwise) | refit; for multi-channel data the convention moves the α/p weight by k_f² (≈ 5 for 12C+12C) → refit the exit widths |
| `coulombIntegral` | 0 / 1 | refit |
| entrance channel radius | adopted ± 1 fm | **full refit** — it also changes the on-shell model (7Li joint χ² 1498 / 1329 / 1143) |
| `spectatorEnergy` | 0 / 0.5 MeV (30–50 MeV/c cuts give 0.5–1.4 MeV for a d spectator) | refit, or fixed-area estimate; the 18O doublet fit breaks down (ρ near a zero of M₀) |
| Coulomb-distortion weight `weight[k]=` R(E) | none / published R(E) curves / spectator penetrability ratio | heavy or sub-barrier spectators only (12C+12C); **full refit** — a restricted refit gave χ²_THM 147, a full one 56 (the weight is absorbed into widths) |
| resolution σ | quoted ± 10 %, or free | refit (fitted σ: 18O 17.0 ± 0.8, 17O 22.2 ± 1.2 keV) |
| data treatment | HOES vs on-shell reading; carrier B; background on/off; anchor dataset | full refit |

**Refit vs estimate.** Re-profiling the THM scale at fixed parameters only measures
goodness of fit (7Li lambda32: 70.7 re-profiled, 65.5 refit), never parameter shifts.
The fixed-peak-area estimate — γ_x² scales as [|M_l|²K]_old/[|M_l|²K]_new per level,
the anchor removes a common factor — reproduced refits to ≤ 5 % for isolated 12C+12C
levels below 1.75 MeV, and **fails near vertex nodes** (1.955 MeV 4⁺; 2.339 MeV 4⁺
moved ×490) and for interfering levels. Plot |M_l(pa)|² over the window first
(nodes for l = 0 at 3.15/1.32/0.283 MeV for a_p = 4.1/5.1/6.1 fm in
p+¹⁸O; 12C+12C l = 4 at 1.98 MeV for 7.5 fm, l = 0 at 0.88 MeV for 6.5 fm). Mark
estimates with `*`.

**Theory-dependence index** = max/min of the quantity over the values of one option.
Report the table (reaction × quantity × option) and the model range (envelope over
the accepted variants). For tensions use z_tot = ln(r)/σ, r = this work/literature,
with half the log-range added in quadrature to the statistical errors.

Calibration (index):

| case | quantity | dominant option | index |
|---|---|---|---|
| 7Li(p,α) joint | S(0); ωγ(2⁺ 20.1) | radius; vertex | 1.10; 1.15–1.21 (S(0) ≤ 4 % for every THM option) |
| 6Li(d,α) joint | S(0) | radius, spectator | ≤ 1.01 |
| 15N(p,α₀) joint | ωγ(312 keV); S(0) | vertex | 1.11; 1.14 (radius 1.04–1.08) |
| 18O doublet joint | Γ_p, ωγ(0.60, 0.80 MeV) | vertex | 1.16–1.19 (radius 1.04–1.12) |
| 17O(n,α) anchored ratio | ωγ(2⁺)/ωγ(3⁻); ωγ(5⁻) | radius | 1.54; 1.76 (spectator 1.35*) |
| 18O narrow, peak areas | ωγ(20), ωγ(90 keV) | radius a_p 4.1–6.1 fm | **×25** (spectator 3.4*, C_l 1.7, vertex/KF ≤ 1.03) |
| 12C+12C, peak areas | Σωγ(E < 1.5 MeV); ωγ(0.985 0⁺) | radius 6.5–8.5 fm | **×180–200** |
| | same | Coulomb weight R(E) | **×15–32** (spectator 30* for the sum, KF 1.6–2.6*, vertex ≤ 1.2) |

Pattern: when overlapping direct data fix the widths and scale, every THM option
moves results by ≤ 20 %; when strengths come from HOES peak areas relative to an
anchor, |M_l(pa)|² enters directly and the radius (and for heavy systems the Coulomb
weight and spectator energy) dominate. Quoted THM errors (10–40 %) contain none of
this. Note also that an anchored joint fit sees only R(E)/R(E_anchor): 12C+12C
strengths fell ×12–45, not the ×250 applied to the published S*.

**Figures that communicate it.** (i) Forest plot: one row per quantity, ratio this
work/literature on a log axis, statistical error bar, grey bar for the model range,
literature sources as marker shapes, a vertical line at 1. (ii) Index chart: bars (or
a colour grid, log₁₀) of max/min per option for each reaction/quantity, estimates
hatched or starred, "not applicable" shown as such. (iii) |M_l|² vs E for the adopted
and ±1 fm radii with nodes marked, next to the level energies.

## 5. Reproducing a published THM analysis

1. **Triage** (open sources only; do not bypass paywalls or captchas): full parameter
   set? resolution? data? fitted curve or only a band? formula? Most papers fail one;
   ApJ 723 was the only complete multilevel case.
2. **Collect**: formula (and its variants), formal vs observed tables, radii, B,
   carrier, KF, σ, windows, normalisation region, background.
3. **Digitise.** Prefer PDF vector path data (exact); otherwise tick-calibrated
   raster (check every major tick to ≤ 1 px; 18O: 0.7 keV, 9 MeV b). Record the
   digitisation error (≈ 1 %, 0.03 dex) and cross-check against any vector-exact
   replot by others. Cap-sized error bars: add ~5 % in quadrature.
4. **Compare** with exactly one free normalisation over a stated range: rms and max %
   (log data: rms dex), peak-region numbers, fraction of points inside the band.
   Compute the floor (e.g. free Lorentzian per level × channel: 12C+12C 0.8–1.9 %
   rms, 0.11 dex) so you know what "reproduced" can mean. Run an independent Python
   version of the paper's formula to separate engine from interpretation (18O:
   ≤ 0.04 %).
5. **Diagnose non-reproduction.** Internal consistency first (ΣΓ_c = Γ; formal vs
   observed tables: ApJ 723 Γ_p1 8.2 vs 11.1 keV). Then column permutations — all 24
   for four width columns; the 12C+12C Γp0/Γp1 swap beat the next by ×11 in χ², and
   per-level swap tests showed a transposition, not scattered typos. Formula variants
   (ApJ 723 literal eq. 14 with an extra E vs the OES-shaped reading; φ21 radians vs
   degrees). Missing background. Representation (Brune vs formal B; perlevel vs
   constant). Frames (E_lab = 2E_cm for 12C+12C files).
6. **Verdict wording**: REPRODUCED / PARTIALLY REPRODUCED / NOT REPRODUCED, then the
   exact condition ("only with the Γp0/Γp1 columns exchanged"; "with the constant
   vertex at B = S(E₁), 5 % peak rms, while the authors used a fitted L21") and a
   numbered list of what the authors would need to publish.

## 6. Pitfalls met

- Lab vs c.m.: AZURE2 data input is lab (E, angle, dσ/dΩ); Amsel 1967 dσ/dΩ(165°) is
  lab — enter it as such, convert once, compare in c.m. outputs.
- Folding eV–keV-wide levels at σ ≈ 17–30 keV: raise `<targetInt>` sub-points until
  the folded curve stops changing, or fold in Python on a tan-spaced grid around each
  level; a sub-point at E = 0 can return 1e18 spikes (clip, move the grid).
- Data windows that cut real narrow levels (18O 597 keV 3/2⁻: its tail changes the
  direct χ² far outside ±20 keV). The same level put into the HOES model predicts a
  THM peak 3× the broad structure that the published spectrum does not show — ask.
- Unfolded models against resolution-broadened data: Koehler & Graff 17O(n,α) TOF has
  ΔE ≈ ±16 keV (E/100 keV)^{3/2}; compare the folded, bin-averaged model.
- Odd-J / unnatural-parity levels in identical-boson systems: THM-only (section 2).
- Plane-wave vertex for heavy systems: pa ≈ 14 for 12C+12C, M_l has nodes in the
  window; plane-wave Coulomb neglect of a sub-barrier spectator (d+²⁴Mg η ≥ 1.55).
- pyazr units: `calculate_rwa` takes reduced-width amplitudes (Brune session);
  `calculate` takes physical widths (eV); `transform_rwa` maps; never convert a
  `gammaIsRWA` channel by hand. `Parameter.wigner_limit` is 3ħ²/(2μa²).
- Extrapolation grids: `add_extrapolation(..., frame="lab")` is the default; pass
  `frame="cm"` for c.m. energies (before 7a35b9a the c.m. doc was wrong).
- Regression-test projects are pins, not physics (`tests/7Li_p_a`, `tests/6Li_d`:
  wrong carrier, on-shell points fitted as HOES).

## 7. Plot conventions

PRC style: serif, inward ticks on all sides, Okabe–Ito colours with distinct markers,
residual panels under fits, PDF + 300 dpi PNG drawn at final size (full 6.5 in, half
3.25 in; labels 11 pt, ticks 10 pt — large enough to read in print). Every panel
carries the **reaction as a bold label inside the axes**
(`r'$\mathbf{{}^{18}O(p,\alpha){}^{15}N}$'`). Show data and our fits only — no
published curves stacked on top unless the figure is a reproduction. Scale data by
fitted norms (say so in the caption), give χ² per set in the caption, and state the
fold (σ, frame) of every model curve.
