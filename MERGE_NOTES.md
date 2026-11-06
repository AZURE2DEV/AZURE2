# thm -> dev: differences for classic projects

For the merge of `thm` into `dev`. "Classic" means a project without a `<thm>`
block and without a THM segment (isDiff >= 10). Everything listed here changes
something for such a project. The THM machinery itself (HOES cross section,
`<thm>` keys, experiments, R(E), the DW vertex, the ps and angle windows,
cbackground=, thm_experiments.out) is new and does nothing unless it is used.

`origin/dev` 3264896d is merged into `thm` (November 2026; see "Merge of
origin/dev 3264896d" below; the previous merges were 2d97b1f, 0919dcb,
8673c5c, 7e6dba8, 64ed575, 4a93bc1, and 6f3228e, 3628a87), so all of dev is
already in `thm`.
The dev changes merged by 18c49ea, f56c150, da72def, 497548a, e191991 and
3628a87 are not differences: the Wigner-limit bound, the elastic identical-pair factor
(7c34992), the adaptive-grid convergence work, target-integration anchoring,
the beam profile, the LM energy-shift Jacobian, CR line endings and the EC
cache guard.

## A. Results that move

| commit | change | reason | re-pinned (tests/*/expected) |
|---|---|---|---|
| b6cc41b | ShftFunc below threshold: exact S_l(0) at E = 0; an open channel with eta > 100 is closed (P = 0, S = 2 S(0) - S(-e)); near threshold (outside eta <= 1000 and (eta <= 50 or z >= 0.05)) S is interpolated linearly between S(-bindingMin) and S(0).  Since the merge of 64ed575 the exact value itself is dev's U ratio (`CoulombTable::ExactBound`, dev 441d14f) and the neutral closed form (dev af3075e); dev uses the same values where its table or the U ratio is finite, and the numerical radial derivative (0 or NaN close to threshold) elsewhere | NaN and underflow near threshold; ~1e-9 noise in S became 0.03-0.5 % in dS/dE (Brune transform) | 22c3d03: 13N, 13N_capture_ay, 7Li_p_ay, beam_profile_kernel, data_file_comments, hybrid_potential, polarization_product, target_effect_ranges, thick_target_resonance (all < 1e-3) |
| 12fe484 | ChannelFunc: one helper for S, P and dS/dE at a channel energy (e < 0, e = 0, eta > 100 closed, Coulomb), used by the CNuc transforms, the boundary conditions, AZUREGrad, the adaptive grid, EPoint and ECIntegral; GSL U trusted only for eta <= 50 or z >= 0.1 | a level at or near threshold aborted the CLI or gave NaN; dS/dE was 0 for a level exactly at threshold | none; a width that cannot be converted (P = 0) becomes gamma = 0 with a WARNING  Since the merge of 2d97b1f, ECIntegral near a Coulomb threshold is dev's (3e9b529: channel capture interpolated linearly between the threshold value and eta = 100, continuous at both ends) instead of thm's mirrored energy |
| ec416de | dS/dE central difference with starting step \|e\|/4 within [1 eV, 1 keV] (was 1 eV).  Since the merge of 64ed575 a stencil wholly above threshold uses dev's Richardson difference (PEShift_dE, 441d14f), as dev; the step rule remains for closed channels (dev: 1 eV) and stencils that cross threshold | round-off dominated: 1.6e-3 in chi2 between compilers | 7Li_p_ay, beam_profile_kernel, data_file_comments, target_effect_ranges; 15bb0a2: 13N, hybrid_potential (< 7e-5) |
| 02ae97b | adaptive grid: geometric tails beyond the resonance core; anchors quantised (width to 1.25^n, energy to width/4); only the finest of overlapping lattices kept; a pure Gaussian convolution whose `<targetInt>` line omits the multiplier uses a 2-width core (an explicit value, and target integration, keep 20) | narrow resonances folded 16-38 % high; the union of lattices tripled the sub-points | thick_target_resonance 19.0249 -> 30.8444, hybrid_potential (+1.6e-3); tests/convolution_unify (dev's check; its reference_old/ holds dev's grid): plain target integration 6.1e-5 from the reference, the constant-sigma ratio new/old 1.00259-1.00272 instead of 1.0027100 at every energy.  thm's grid is converged there (3e-8 and 2e-7 from a 4x finer grid), dev's is off by 6.3e-5 and 1.1e-4 near the 3.50/3.55 MeV levels; held to 1e-4 and to a 1.5e-4 spread (the median still to 1e-5)  dev fixed the same tail error differently (d216517, graded zones moving with E_R); the merge of 2d97b1f keeps thm's grid: dev's tests/reference/adaptive_grid_reference_test then fails 1 of 14 checks (below) |
| (02ae97b, 528b672) | the thm sub-point grid refresh (EPoint::RefreshSubPointGrid) is gone since the merge of 64ed575: dev re-places every grid at every evaluation (EData::ReplaceSubPointGrids, a163c06, `--no-grid-replacement`), from the observed values of the current parameters (TransformOut on a copy), which serves Brune, formal parameters and Park alike.  No difference from dev left; the quantized anchors of the row above make the re-placed grid piecewise constant | -- | -- |
| df8c2c8 | (1 + delta_12) = 2 for every cross section out of an identical entrance pair (reactions, capture, differential); none for an identical exit pair; the adjoint gradient follows; a rate run prints a note | dev applied it to elastic only: 12C(12C,a), d(d,p), d(d,g) were a factor 2 low | none (new tests/identical_entrance_reaction) |
| 76888ad | a `<targetInt>` line listing several segments gives each segment its own copy, so sigma lab -> c.m. is converted once | the shared object was converted once per listed segment | target_effect_ranges 3861.5 -> 3838.66. dev fixed the same bug with a once-only guard (f0444ab, pin 3838.73); both are kept, and the 7e-5 between the pins is the shift-function rows above. A line shared by segments of different entrance pairs still differs: thm converts each copy with its own pair's factor, dev uses the first segment's: tests/natural_neon ("5,6": n+22Ne and n+20Ne) segment 6 247.1 on thm, 246.304 on dev (dev with separate lines: 247.15), total 3346.1 (dev 3345.28); examples/natural_neon shares "10,11", "6,8" and "7,9" the same way (its documented table is dev's; not run here)  Since the merge of 2d97b1f dev converts every point with its own pair (170ab0c, EPoint::SetTargetKinematics) and so does this branch; the copies stay (screening/density parameters are per line), no conversion is done on them, and the natural_neon difference should be gone (not measured: memory) |
| 79358bb | EPoint::CalcEDependentValues clears its L_o, P and phase tables before refilling | a second call appended copies; lookups read the stale first copy | none. dev fixed the same doubling in RecalcEDependentValues (797c2ea, rows emptied in place); since the merge both functions use dev's ClearEDependentRows, so no output difference remains |
| 953486f | nuisance priors looked up by full Minuit index in chi2, gradient and LM setup | a fixed parameter before the prior's moved the prior to another parameter (also on dev) | none (new tests/nuisance_prior) |
| 178ae16 | AZURECalc::FixedMask from the minimizer: LM, GSL-LM, the band covariance and the FD gradient honour the fixed flags | --use-lm / --use-gsl-lm fitted fixed parameters | none (new tests/fixed_param_par); covariance.dat has one row per free parameter |

## B. Behaviour and format (no physics change)

- a0cd873: the CLI exits with -1 when a run fails (fill, initialize or transform abort); it always returned 0.  Since the merge of 2d97b1f it exits with 1, as dev (2b9ea26).
- d284224: a malformed `<potential>` block makes ReadConfigFile return -2; the binding refuses the file instead of running with defaults.
- 953486f: Calculate prints "Total Chi-Squared:" (the total with priors) to 12 significant digits, not 6.
- a2b8c81: "Denominator less than zero while transforming" names the channel and its saturation width.
- c465277: GslQuiet replaces the process-wide GSL handler swaps (WhitFunc, ShftFunc), which were not thread-safe. Values unchanged.
- e3f43df: segment norms and energy shifts are recognised by their name form, not a substring (no classic name contains the words).
- GUI:
  - 5832bee, ee092f9, 2887d99: level lines have 33 fields (32: THM binding energy, 33: gammaIsRWA; both 0 for classic). Dev's engine ignores trailing tokens.
  - d2d30a9, 2887d99, 3e7bd97: doubles are saved with round-trip precision (was 6 digits) and fixed-width columns keep a space; opening no longer rounds free level parameters. A GUI-saved classic project can therefore give a different chi2 than on dev (7Li_p_a: 2137.82733857 vs 2137.82673527 after open + save on dev).
  - f3210c7: `<parameterSettings>` level rows use AZURE2's numbering.  Since the merge of 7e6dba8 this is dev's implementation (be817bf, e36d485: the `# level numbering: engine` marker, unmarked blocks read by their value column); left as differences: a reduced-width-amplitude channel's row says `Width (MeV^(1/2))`, and a level without channel lines is not numbered (it writes no line).
  - 32701d7: pairs are identified by key; repeated segment lines read from a file are kept.
  - 57aba50, 065041b: the RWA -> physical conversion on a param.sav load runs on the live GUI state (a temporary snapshot) and leaves the open project alone.
  - 2b0755d, 0304f52: the plot validates cross section and S factor separately; the band follows the points kept.
  - this branch, October 2026: moving, adding or deleting a data segment line renumbers what refers to it by number: `<parameterSettings>` segment_N rows (settings and prior_centre), the Experimental Effects segment lists, and the `<thm>` experiments and weight[k]; a line left without segments is removed, with a status-bar note. Unedited projects save byte for byte.
- pyazr, visible on classic projects:
  - ac3611b: save_fit writes the fitted free norms and shifts by default (norms="fitted"; "nominal" is the old behaviour); write_output_files honours the `<config>` output directory; a Python list passed to calculate_chi2_rwa no longer crashes.
  - 8def2ab: save_fit writes fixed widths as the fit had them.
  - e7f7b16: save_fit writes prior_centre rows for the priors whose column it moves.
  - 783fef9: calculate_chi2_physical applies the norms of its vector; widths.teichmann_wigner() no longer multiplies by 1.5 twice.
  - 80cd478, d178634, 57aba50: observable codes 7/8 named; the final newline kept; the pair record has field 16 (binding_energy).
- Tests: tests/run_tests.sh and the check.sh scripts put a time limit on every run on every platform (tests/lib/guard.sh: timeout, gtimeout or a shell watchdog).

## C. New opt-in features (no effect unless used)

- e7f7b16: `<parameterSettings>` rows `segment_N_norm|segment_N_energy_shift prior_centre c`. Older engines and GUIs skip them.
- 5832bee, 0871c99: level-line field 33 gammaIsRWA (the channel value is an amplitude in MeV^1/2).
- pyazr: add_extrapolation(frame="cm") (7a35b9a; "lab" is the old default), tabulate() (7bbab4e), composite segments (f8bcd5e), transform_all_rwa(include_fixed=) (8def2ab), modelavg and scripts/thm_model_average.py (88fa280, 884571a).
- New classic test: tests/15N_p_a (0ee7ee7, pin 4601.32).
- GUI, October 2026: THM is opt-in, as the hybrid nuclear potential is
  (Configure > Runtime Options > Use Trojan Horse Method).  Off for a new or
  classic project: no THM Workspace entry, no THM tick in the segment
  dialogs, no THM Background sub-tab.  On by itself when an opened project
  has THM content; not stored in the file.  Switching it off hides the
  controls only.  The THM code of the main window and the Fitting tab is in
  AZURESetupThm.cpp and FittingTabThm.cpp.

## D. Speed and structure (values unchanged)

- this branch, October 2026: CoulFuncCache keeps the values of a key whose
  near-energy memo it gives up (4096 queries, under 5 % hits) in an
  exact-energy table (bits of E, capped at 32768, dropped if it fills without
  hits).  An exact hit is what recomputing gives, so every output is byte
  for byte the same; a second session in one process (pyazr, the GUI's THM
  diagnostics, a save_fit check) finds them (12C+12C THM example: 24.5 s ->
  2.6 s; the first session 49 -> 45.5 s, as such energies do recur).  Peak
  RSS +13 MB there.  The near-energy memo itself (energies within 1e-12 MeV
  taken as equal) is dev's and unchanged; it makes a warm-memo session differ
  from a fresh process in the last digits (6Li_d: 682.1465248015013 against
  682.1465248016414), on dev as here.
- this branch, October 2026: the classic files lost their THM blocks where
  they could move unchanged: EData.cpp -> EDataThm.cpp (second pass),
  Config.cpp -> ConfigThm.cpp and Config.h -> ThmOptions.h (Config::ThmOptions
  is a typedef of ::ThmOptions).  Config.cpp/Config.h now differ from dev by
  6 and 10 lines.

## Merge of origin/dev 6f3228e (3628a87, October 2026)

dev between d801ba5 and 6f3228e: the sqrt(E) energy-shift term (1c3e7e3),
Park's parametrization (c1ec777, 692e7f8, 397668d, `--use-park`,
`--no-brune`), isDiff 8 removed (330dc1e), the energy-dependent convolution
window from the c.m. energy (313c7d5), a shared `<targetInt>` sigma converted
once (f0444ab), the resident-memory fix (797c2ea), param.fit at the
evaluated point (67ac25d), EPoint::parentSegment_ dropped (588c584), the
`#parametrization` tag and name-keyed .sav readers (8cdc841), compound
stopping per active atom in the GUI (751bec5), skill notes.  All of it is in.

Decisions taken at the merge (they are differences from dev only where
marked):

- Parameter order: levels, norms, energy shifts, sqrt(E) coefficients (dev),
  then the THM cbkg block.  A classic project has dev's layout.
- `GetParameterInfo` type codes: 4 = sqrt(E) coefficient (dev), 5 = THM cbkg
  (was 4 on thm).  pyazr `_KINDS` and `bands._BAND_KINDS` follow;
  `bands._NFIELDS` is 16, the record length since thm added `input_is_rwa`
  (it was 15 on thm: a latent bug in `live_parameters`).
- `AZURELabel::IsEnergyShiftName` also matches `_energy_shift_sqrt`, so the
  LM penalties, the limits manager and `GetEnergyShiftIndices` treat the sqrt
  coefficient as dev's substring tests did.
- Park and thm's amplitude input (gammaIsRWA): the value in the file is
  Brune's amplitude in both modes (October 2026, with THM under Park, below);
  under Park it is read as gamma sqrt(J), J = denom / (2 numer) =
  1 / (1 + sum gamma_Brune^2 dS/dE), and written back
  (CNuc::GetTransformParams: save_fit, the GUI's parameter-file load) as
  gamma_Park / sqrt(J).  At the merge it had been taken as Park's amplitude,
  which made one file two different models.  Classic projects without the
  flag: as dev.
- Park's dS/dE (CalcShiftFunctions, ConvertAmplitudeBasis) is ChannelFunc's,
  as every other dS/dE on thm (row ec416de); d2S/dE2 is dev's.  **Differs
  from dev** in the last digits under `--use-park`; tests/park_formalism and
  pyazr_park_gradient pass.
- THM under `--use-park`: refused at the merge, allowed since October 2026
  (see "THM under Park" below).
- isDiff 18 (THM + P dsigma/dOmega) is refused like 8.
- Sub-grid refresh (02ae97b) skips a segment with a sqrt(E) shift (not a
  translation of the grid).
- The THM segment's unfloored shift (thm) uses dev's TotalEnergyShift, so a
  THM segment can carry the sqrt(E) term.
- EPoint: dev's in-place row clearing is used by thm's CalcEDependentValues
  too, so dev's memory fix is kept.
- The Park penalty is added in CalculateChi2RWA only (as dev), not in
  thm's shared EvaluateFilledChi2 (residuals, Chi2Physical).
- Dev's three new check.sh scripts take thm's time guard (tests/lib/guard.sh).
- c8d8db5: the GUI keeps a segment line's `sqrtshift` block on save.  Since the
  merge of 64ed575 this is dev's `extraTokens` (3209835: every keyword block,
  sqrtshift included), so no difference is left.

Verified after the merge (one heavy process at a time, OMP_NUM_THREADS=1):
classic test projects give the same files as 101e028 (thm before the merge)
byte for byte, except param.par's new `#parametrization` line and
`segment_N_energy_shift_sqrt` rows (dev's format) and the new
energy_shift_sqrt project; against dev they differ only by the rows of
section A (largest: thick_target_resonance, 02ae97b; hybrid_potential and 13N,
b6cc41b/ec416de/02ae97b; parameters.out may print a width at exactly 1 keV as
"1.000000 keV" where dev prints "1000.000000 eV").  Every THM example and THM
test project gives 101e028's files byte for byte, apart from the same two
param.par format rows.

dev's own suite (its tests/run_tests.sh and check.sh scripts) passes 16/16
with dev's binary and 15/16 with the merged one: thick_target_resonance moves
19.03 -> 30.84 (02ae97b, re-pinned on thm); the other pins agree within the
suite's tolerance.  All THM check.sh scripts of 101e028 give 101e028's files
byte for byte (3250 files), apart from param.par/.sav/.fit format rows, the
temporary directory names they print and dev's shifts.out format (sqrt_shift
column, complete rows only).  GUI: the merged GUI's code differs from 101e028
by exactly dev's GUI changes; an opened and saved classic project differs
from dev's save only by section B (fields 32/33, round-trip precision,
AZURE2's level numbering in `<parameterSettings>`, pairs by key) and every
classic and THM project saves byte for byte on a second open + save; the
`<thm>` blocks are kept verbatim.

Resolved from the old "to check" list:

- 79358bb: dev did recompute sub-points twice and fixed it the same way
  (797c2ea); no difference left.
- 02ae97b: the GUI's Add Experimental Effect dialog writes the multiplier
  (default 20), so GUI-made projects keep the 20-width core.

Still open:

- 12fe484: the narrower GSL U trust region could move S slightly for light
  pairs at large eta that no test covers.

## Merge of origin/dev 64ed575 (October 2026)

dev between 6f3228e and 64ed575 (137 commits): the universal Coulomb table
(stages 1-4: `azure2-coulomb-table`, `--coulomb-table[=FILE]`, on by default,
`--build-coulomb-table`, Configure > Build Coulomb Table; not distributed),
sub-point grids re-placed at every evaluation (EData::ReplaceSubPointGrids,
`--no-grid-replacement`), exact shift functions (U-ratio S for closed charged
channels, Richardson dS/dE for open ones), the closed form for closed neutral
channels, level widths matched to their J-group channel by (pair, l, s), the
Q-coefficient fix for UPOS points, the abundance prior in the LM/GSL cost and
MCMC, full precision for GUI level lines, the minuit_index fallback limited to
compatible names, the AZUREParams nan-error fix, several compound nuclei in
one project (samples, `<abundances>`, abundance parameters, `merge_projects`),
shift groups, normalization groups, energy-shift functions (`shiftfn`), the
SAMMY UDR resolution function, hard-sphere groups (`includeAllJPi`), one
Gaussian beam-spread quadrature (dad58dc: a normalised kernel, lab sigma(E)),
point mapping within a segment, PrecalculatePoints (one parallel pass), the
Jacobian G2 for `--use-gradient`, the target-integrated A_y adjoint, the EOF
prompt fix, ERYA table fix, the natural_neon suite and example, docs moved to
azure2dev.github.io, skills.  All of it is in, with these decisions:

- **Shift functions, one implementation.**  ShftFunc below threshold: S(0) at
  E = 0; dev's neutral closed form; dev's table bound sheet when a table is in
  use; else dev's U ratio (`CoulombTable::ExactBound`, U' = -a U(a+1, b+1, z)),
  with thm's linear continuation to S(0) where GSL's U is not trusted: eta >
  1000, or eta > 50 and z = 2 kappa a < 0.05 (b6cc41b had eta <= 300, z >= 0.1
  and the U(a-1, b)/U(a, b) form).  Against mpmath (p+p ... 12C+12C, l = 0-3,
  binding 1e-10-3 MeV) the largest error of the function is 8.7e-9 (was 2.2e-8
  with thm's form, and dev's U ratio alone is off by up to 1.3e-6 at p+p, eta
  111, z 1.2e-3, and NaN within ~1e-7 MeV of threshold, where dev falls back to
  the numerical derivative that underflows to 0).  WhittakerShift is that U
  ratio; `AZURE_BOUND_SHIFT_NUMERICAL=1` still selects dev's old route.  Open
  channels: dev's Richardson dS/dE (PEShift_dE) for every stencil wholly above
  threshold; ChannelFunc keeps thm's branches at and below threshold and for
  eta > 100.  ExactBound and the shift-function SVD of EData use GslQuiet
  (dev swapped the process-wide GSL handler, not thread-safe, and the SVD one
  was never restored).
- **Sub-point grids: dev's re-placement only.**  thm's
  EPoint::RefreshSubPointGrid, CurrentGridAnchors, the thm SubGridSpec and
  GridConfig::useFitParameters / formalParameters / parkAmplitudes are
  removed: ReplaceSubPointGrids re-places every grid each evaluation from the
  observed values of the current parameters (TransformOut on a copy), which
  covers Brune, formal parameters (thm_nonbrune_grid_test passes) and Park
  (e186e78's intent: the same grid in both parametrizations; thm_park passes),
  shifted segments (thm_energy_shift) and components.  The THM fold window
  (thmConvolutionRange, 5 sigma) is applied after dev's TargetEffect::
  GridExtent.  Note: `--no-grid-replacement` keeps the file-based grid, which
  for examples/f19_pag_thm is far off (chi2 626 instead of 77.64); thm's
  refresh used to replace it at the first evaluation.
- **Coulomb table and THM.**  The THM quantities that are the R-matrix's P and
  S (the vertex boundary S_c, ThmLevelWidth's P and dS/dE, R(E)'s
  penetrabilities through the point tables) use the table where it covers the
  channel, as the R-matrix does; the thread-local memos of
  ThmVertexBoundary::ShiftAtLevelEnergy and ThmLevelWidth carry the table in
  their keys.  The Coulomb waves of ThmCoulombTerm (N_C, coulombIntegral),
  the DW vertex and ThmDistortion (cwfcomp) bypass it, as does thm's
  CoulFuncCache exact-energy table (it holds F, G, which the table does not
  replace).  Table against exact on every THM example and test: models
  within 1e-8 relative, formal parameters of the transform-limited levels
  within 5e-6.
- **Parameter layout and type codes.**  Levels, abundances (dev), norms,
  shifts, sqrt(E), shift-function parameters (dev), THM cbkg (last).
  GetParameterInfo types: 5 shift_fn, 6 abundance (dev), **7 cbkg (was 5 on
  thm)**; kParamInfoFields 16 (field 15 input_is_rwa).  ParamKind order ...
  Abundance, ThmCoherent.  MCMC kinds: PARAM_ABUNDANCE 6 (dev: 5; internal),
  PARAM_CBKG 7, shift-function parameters as PARAM_SHIFT (dev's prior rule);
  the by-name classification knows abundance_* and segment_N_energy_shift_fnJ.
  pyazr `_KINDS`, `bands._BAND_KINDS = (0, 1, 7)`.
- **Pair record**: kPairInfoFields 20: dev's nucleus Z, A, abundance in
  fields 16-18 (dev's positions), thm's binding energy **moved to field 19**
  (was 16).  pyazr `Pair._NFIELDS = 20`.
- **Width matching by (pair, l, s)** (b70d1c8): ALevel::SetGammaFromLine and
  PadChannels carry thm's gammaIsRWA flag with the width (dev's code did not
  know it: a reordered level would have kept the flag on the wrong channel).
- AZURELabel::IsEnergyShiftName also matches `_energy_shift_fnJ`, so the LM
  penalty rows (dev's ForEachShiftPenalty priors), the limits manager and
  GetEnergyShiftIndices treat shift-function parameters as dev's substring
  test did.  The LM/GSL penalty rows: norms and shifts by name form (thm),
  shift priors from ForEachShiftPenalty and abundance rows (dev), nuisance
  priors by full index (thm, 953486f).
- minuit_index fallback (b10255d): `cbkg_*` settings, like `segment_*`, are
  exact Minuit names and never fall back by index; abundances are left out
  of the nuisance priors (they carry their own).
- Chi2On(thmOnly) / Chi2Value: dev's pooled copies and PrecalculatePoints
  (THM-only for the THM finite differences), thm's THM experiment profiling,
  the Park wall (thm dcb8064) and dev's abundance penalty (not in the
  THM-only part).  Gradient: dev's structure (Jacobian path, FdStep, shift
  groups), thm's FixedMask, the HOES adjoint rows, THM-only differences for
  THM norms; a shift parameter of a THM experiment segment is not local.
- CopyEDependentValuesFrom (dev, mapped points) copies the THM tables too.
- THM form factors (thm_jl_, rhodjl, coul) are stored for THM points only:
  for classic points they were three zeros per channel per sub-point, 0.6 GB
  with dev's hard-sphere groups on tests/convolution_unify (430k sub-points).
  The merged CLI still peaks 15 % above dev there (1.74 vs 1.52 GB; EPoint
  896 vs 744 bytes) -- the THM members of every point.
- EOF at a prompt: thm's readAnswer (the prompt named, exit 1, blank defaults
  taken) replaces dev's promptLine.  **Differs from dev**: at the main menu dev
  printed "End of input on stdin: exiting." and exited through the Exit
  option; here it is the error and exit 1 (tests/cli_eof).
- dev's banner left the log stream left-aligned (`std::left` before the
  version column): every later number of the log moved ("Iteration: 440
  Chi-Squared"); `std::right` restored after it.  **Differs from dev**.
- `--help` lists `--grid-replacement`; tests/cli_help reads `--flag[=FILE]`.
- pyazr: save_fit has thm's options (norms=, close_session=, cbkg, prior
  centres, chi2 check) and dev's verification (_snapshot_mismatch: amplitudes
  up to a level's sign) and abundances; the chi2 check is 1e-4 in the
  standard parametrization (dev's beta_sign_test: ~1e-5 from the iterative
  inverse transform), 1e-6 otherwise.  transform_all_rwa takes
  `include_fixed=` (thm) and `free_only=` (dev).  calculate_rwa in
  extrapolation mode accepts a data-mode vector (its free R-matrix and
  abundance block), as before dev's length check (extrap_layout_test).
  update_rwa_params_from_sav is dev's load_sav on the configured output dir.
  penalties(): nominal_norm (prior centres) with dev's normgroup rule.
- GUI: SegmentsData/TestData carry isTHM then dev's extraTokens (thm's
  sqrtShiftTokens is gone: dev keeps every keyword block); level lines keep
  thm's round-trip writer, in dev's per-nucleus order; PairsData column 15
  (THM binding) before dev's derived Compound Nucleus column (16); MCMC tab
  categories know abundance and shift_fn.
- Tests: tests/thm_coherent_background/toy's exit pair is n + 16O (one
  compound nucleus; it was n + 12C, which dev now refuses);
  park_chi2_fallback makes the adjoint bail with an auto-tolerance target
  effect (a convolved A_y no longer does, 2e726db); park_formalism's
  iteration regex takes either alignment; target_int_tab_test's long-number
  case is section 13.
- CI: thm's venv Python (pybind11, numpy, scipy, mpmath) for all three
  runners, dev's queued concurrency and `-LE slow`; dev's second Python
  pinning block (system python3, PYBIND11_FINDPYTHON) dropped as redundant.

Results that move (no table file; the table changes nothing beyond 1e-8):

- Classic against dev 64ed575 (dev's tests and examples, and the classic
  tests of thm, both binaries on the same inputs): identical apart from
  section A -- thick_target_resonance 19.0249 -> 30.8444 (02ae97b),
  hybrid_potential 2.46215e6 -> 2.46612e6 (1.6e-3), 13N segment 1 814.402 ->
  814.349 (6.5e-5), target_effect_ranges 3838.74 -> 3838.66 (2.1e-5),
  beam_profile_kernel 40.4474 -> 40.4393 (2.0e-4; thm's grid and dS/dE rows
  on top of dev's detailed-balance fix), parameters.out of 7Li_p_ay-type
  projects in the 7th digit.  All of these are the pins thm already had,
  except beam_profile_kernel (re-pinned 40.4393).  natural_neon (test and
  example) needs more than 2.4 GB with either binary and could not be run
  here.  Found by CI after the merge and measured with both binaries under a
  3.0 GB address-space cap: tests/natural_neon segment 6 247.1 (dev 246.304;
  row 76888ad, a "5,6" line over two entrance pairs; dev's binary needs more
  than 3.0 GB for the whole project and was run on segments 5 and 6 alone),
  and tests/convolution_unify's plain target integration and constant-sigma
  ratio (row 02ae97b: its reference holds dev's grid).
- THM against ab85222 (every THM example and test): dev's normalised Gaussian
  kernel moves every folded model by the kernel's missing tail, 2.70e-3
  (1 - erf(3/sqrt 2)) for direct data under a 3-sigma window and 5.7e-7 for
  THM segments (5 sigma); the shift-function and dS/dE changes move the rest
  by <= 1e-5 (formal parameters of transform-limited levels up to 1e-4, e.g.
  f19's 810.82 -> 810.90 MeV^1/2 at the same widths).  Data chi2:
  c12c12_tumino2018 111.651 -> 111.613 (-3.4e-4; total with priors 112.0758
  -> 112.0376), li7_tumino2006 439.104 -> 439.148 (+1.0e-4),
  li6_pizzone2011 156.991 -> 157.036 (+2.9e-4), f19_pag_thm 77.6398 ->
  77.6400, 7Li_p_a 2138.48 -> 2138.45 (-1.3e-5), 6Li_d 682.147 -> 682.146,
  18O_p_a_thm norm 2.9425e-5 -> 2.94251e-5; the others in the 7th digit.
  Re-pinned: tests/6Li_d, 7Li_p_a, 18O_p_a_thm, the examples' documented
  totals, pair_identity's 112.0376.  THM check.sh: all pass; the only pin
  moved is thm_options kinematics=lambda32, 5499.08 -> 5604.98 (1.9 %): the
  lowest point (E = 76 keV) folds the 1/k_i of lambda32 through E = 0, an
  integrable singularity that makes it depend on where the sub-points fall,
  and dev's re-placement keeps the anchor of level 5002 at -0.336 MeV (Gamma_alpha
  134 keV) that thm's refresh dropped; that point moves 18 %, every other by
  <= 1e-4.
- THM best fits (/home/almalinux/thm_bestfits/12C_12C, file alone): plain_6.5
  91.5788 -> 91.5550, plain_8.5 70.3198 -> 70.3222, NCxR_8.5 100.470 ->
  100.472 (the direct segments' kernel normalisation; single segments by up
  to 3 % of their small chi2).

Still open: the merged EPoint is 152 bytes larger than dev's; the three
convolution_unify models with target integration and a beam spread need more
than 3.0 GB of address space with either binary (checked here on 21 points).

## Merge of origin/dev 69826014 (November 2026)

dev between 3264896d and 69826014 (6 commits): thin-target integrations
restricted to energy windows evaluate the points outside the windows at the
depth-averaged energy E_bar (134f6ad1, tests/target_window_fallback), the EC
amplitude progress bars drawn one thread at a time (9fbb82c0,
tests/progress_bars), the light_excitation run-log compare, and the
cross-copy node cache with packed per-point tables (69826014).  Conflicts
only where thm extends the per-point tables:

- **Packed `EPoint::EDependentTables`**: dev's classic rows are now one flat
  array per kind with `rowStart`.  thm's THM vertex tables (`thm_jl`,
  `thm_rhodjl`, `thm_coul`) stay in row form inside the same struct: they
  are filled per entrance channel and read with bounds checks, and only THM
  points carry them.  `Clear()` and `ShrinkToFit()` cover them, so
  CalcEDependentValues no longer clears them separately (it rebuilds them
  right after ClearEDependentRows in any case).
- **NodeTableCache (dev)**: holds the tables and Legendre polynomials weakly;
  `thmPs` rides in the cached `EDependentShare` as a strong pointer (the
  spectator-window table is immutable, and an entry is swept once its tables
  expire).  The key (segment key, c.m. energy, c.m. angle) plus
  AdoptEDependentValues' exact check covers THM sub-points as well.
- tests/progress_bars (new from dev): its bar pattern now allows spaces
  before the segment number.  thm's stream prints the number right-aligned
  in the setw(12) field (`Segment #           1`), dev's left-aligned
  (`Segment #1           `); otherwise the same bars, whole at 8 threads.
  Failed the same way on thm 3ba23e6a + dev 02b0001e without the change.
- tests/shared_point_tables keeps thm's 109.250542525 (unchanged by this
  merge, 1 and 4 threads); dev's own pin is now 109.134 (dev's earlier
  109.136 was stale: dev 02b0001e gives 109.134).

Checked against thm 3ba23e6a merged with dev 02b0001e (everything but the
node cache and packing): the calculations of examples/ c12c12_tumino2018,
f19_pag_thm, li6_pizzone2011, li7_tumino2006, n15_lacognata2007,
o17_guardo2017_fit, o18_lacognata2008, o18_lacognata2010 and
tests/18O_p_a_thm are byte-identical; capped 4-thread fits (c12c12,
o17_guardo2017_fit, f19_pag_thm) give identical chi2.  Their param.fit and
fitStatus.out (call count) differ, but so do two runs of the reference: at
several threads thm's param.fit records a few parameters at a
finite-difference step of the last pooled evaluation (~1e-3 off the fitted
value).  Pre-existing on thm, not on dev (dev's 1- and 4-thread fit outputs
are identical, tests/shared_point_tables check c); open.

Tests of this merge: ctest -LE slow 91/92 (adaptive_grid_reference, the
known 1 of 14 below), run_tests.sh 128/128 (with progress_bars as above),
ctest -L slow 8/8.

## Merge of origin/dev 3264896d (November 2026)

dev between 2d97b1f0 and 3264896d (18 commits): free-gas Doppler broadening
and its fitted temperature (`doppler <T>`, `doppler T vary sigma [Tmax]`),
the fitted sample thickness of transmission / yield segments
(`transmission n vary sigma`, `thicknessgroup`), the prior parameter
covariance (`priorCovariance=<file>`, pyazr.priors), the light particle of a
pair in an excited state (`e1`), full-precision pi and CODATA 2018 constants,
the InputFiles cache fix, the GUI's tied systematics and macOS file dialogs,
the EData::Clone fix for segments with components (d2c10ac5), and the shared
per-point tables with the CoulFuncCache::TransientScope race fix (3264896d).
All of it is in, with these decisions:

- **`<levels>` trailing columns: told apart by their number.**  dev's
  optional 32nd column is `e1`; thm's columns 32 and 33 are the THM binding
  energy and the width flag, in every THM file.  NucLine (and pyazr's
  AzrChannel, the GUI's reader) now read the trailing fields by count:
  1 = `e1` (dev), 2 = binding energy and flag (thm), 3 = binding energy,
  flag and `e1`; more is refused.  Every existing dev and thm file reads as
  before.  Written: pyazr keeps a line's layout and adds `e1` (34th after the
  THM columns, 32nd without them) only when nonzero; the GUI always writes the
  two THM columns and then `e1` when nonzero.  A thm file whose line carried
  the binding energy alone (32 fields) would now read it as `e1`: none in
  tests/ or examples/ (all 33), and thm's writers always wrote both.  A
  trailing token that is not a number now fails the line (as dev; thm
  ignored it).  To decide when thm goes into dev: this layout, or the THM
  columns moved behind `e1`.
- **GetPairsInfo: 21 fields.**  19 = dev's light excitation, 20 = the THM
  binding energy (field 19 before); `kPairInfoFields` 21, pyazr `Pair`
  reads 21, dev's 20 (binding energy 0) and an older engine's 19.
- **Parameter vector order and type codes.**  dev's Doppler temperatures and
  sample thicknesses follow the densities (as on dev); thm's coherent
  backgrounds stay last.  GetParameterInfo types: 9 Doppler, 10 thickness (as
  dev), **11 cbkg (was 9)**; pyazr `_KINDS`, `bands._BAND_KINDS` (0, 1, 11),
  `AZURECalcMCMC::PARAM_*`, ParamKind (Doppler, Thickness, ThmCoherent),
  ParameterLabel, the GUI's fitting / MCMC tabs follow.
- **One Doppler parameter per `<targetInt>` line.**  thm gives every segment
  a line lists its own copy of the effect; dev made one `doppler_<n>` per
  effect.  As thm's screening and density parameters, `doppler_<n>` is now one
  per line (n = the line), moving all its copies (EData::dopplerCopies_,
  GetDopplerTargetEffect); dev's call sites use the first copy.  With one
  segment per line (dev's tests) nothing changes.
- **Nuisance priors** leave `doppler_*` and `segment_<k>_thickness` out
  (ParameterLimitsManager, as dev's name filter); analytic gradient for both
  (the THM finite-difference block skips them: a THM point has neither).
- **Shared per-point tables (3264896d)**: thm's THM vertex tables
  (`thm_jl`, `thm_rhodjl`, `thm_coul`) live in `EPoint::EDependentTables` and
  are shared by the copies like the classic rows; `thm_ps_` stays an
  immutable shared pointer and is carried by `EDependentShare` (a re-placed
  sub-point at an unchanged node keeps it); the row trimming also trims the
  THM rows.
- **Prior covariance in the API**: thm's EvaluateFilledChi2 has no parameter
  vector, so CalculateChi2RWA / CalculateChi2Physical record its chi2 for
  chiSquared.out (dev did it in its chi2 function).
- **GUI**: pairs model columns 15 = light excitation (dev), 16 = THM binding
  energy (`PairsModel::BindingEnergyColumn`, was 15), 17 = nucleus; segments
  model 24 = isTHM, 25 = Groups (`GroupsColumn`, 24 on dev).  Not built here
  (no Qt).
- **AdaptiveIntegrationGrid**: thm's level energy with dev's heavy-only
  excitation (GetHeavyExE) for the ground-transition test.
- AZURE2.cpp prints the `<levels>` line it cannot read.  NucLine does not
  call `std::ws` at eof (a newer libstdc++'s sentry sets failbit there).
- Tests: tests/convolution_unify keeps thm's grid tolerance for target
  integration and takes dev's 1e-6 for the bare calculation (CODATA 2018);
  tests/pyazr/light_excitation_test adapted to the layouts above;
  tests/pyazr/thm_spectator_window_test mirrors the CODATA 2018 uconv and
  hbar c.  The reference tables of tests/reference/thm_*_reference.py and
  tests/thm_coherent_background/check.sh still use the CODATA 2006 constants
  (they pass at their tolerances; not regenerated).
- adaptive_grid_reference_test still fails 1 of 14 checks (thm's tail bias,
  see the 2d97b1f0 merge below).

Re-pinned (each equal to a measured run of this merge):
- tests/natural_neon 3346.1 -> **3342.51** (segment 6 247.1 -> 245.029,
  3 433.353 -> 433.279, 5 27.9777 -> 27.9437, 7 947.456 -> 946.069; dev
  3342.46, segment 6 245.09).  The thm-only difference (76888ad, section A)
  is gone with dev's per-point conversion and 5-sigma kernels; the rest is
  thm's grid (segments 4-6 by up to 4e-4).  The same value the 2d97b1f merge
  measured and proposed; README updated.
- tests/shared_point_tables (new from dev) expected/chi2_fit 109.136 (dev)
  -> **109.250542525**: thm's adaptive grid.  On its 1/16 natural-neon project
  the calculation differs from dev only in the folded segments 3-6, by 1e-4
  to 2e-4 (131.190625 against 131.196); five MIGRAD calls take that to 1e-3.
  Its 1- and 4-thread outputs are identical (check c).
- Test mirrors of changed formats/constants (no pinned value moved):
  tests/thm_coherent_background (closed form: CODATA 2018 hbar c and u; it
  was off by 4e-6), tests/pyazr/thm_spectator_window_test (the same),
  tests/fixed_param_par (dev's `#columns` header of covariance.dat,
  b7e77292), tests/pyazr/light_excitation_test (the column layouts above).

Tests (this machine, GCC 11 headless):
- `ctest -LE slow`: 91 of 92 pass; adaptive_grid_reference fails 1 of 14
  checks as before (thm's tail bias).
- `tests/run_tests.sh`: 122 of 126 at first; the four failures are the
  re-pins and mirrors above (natural_neon, shared_point_tables,
  thm_coherent_background, fixed_param_par), each passing after it.
- `ctest -L slow` (one at a time): 8 of 8 pass.

## Merge of origin/dev 2d97b1f0 (November 2026)

dev between 0919dcb and 2d97b1f (208 commits): the audit fixes (items 3-32 and
tier 3: shifted sub-points at the parent's angle and moved rigidly with it,
EC final states that must be listed levels, the M1 EC spin/orbital terms and
EC signature v4, 5-sigma Gaussian kernels, the reaction rate from ln sigma
splines, lab-frame SUM terms, RATIO Jacobians, ratio-to-Rutherford lab angles,
UPOS spin and mixing fixes, physical-parameter objectives with norms and
penalties, nuisance priors by full index, normalizations.out rows, the
exact thin-target range integral, per-point lab -> c.m. conversion of a
shared `<targetInt>` line, Bohr straggling, thick targets as density 0,
kinematic-limit refusals, table-edge stencils, session-owned potentials,
exact hybrid derivatives), new features (data covariance matrices, relativistic
kinematics, angular-coefficient data isDiff 9, fitted target densities
`densityfit`/`densitygroup`, sample transmission/yield, total capture rate,
`--max-calls`, AME2020 masses, EC convection current, scattering lengths in
the segment table, the R-matrix route fixes), the speed-ups (one parallel
point pass for pyazr, forward-only Jacobian shift columns, parsing once per
session, setup scaling), shift penalties from zero (0dd28e5), exit status 1,
docs, skills, CI (artifact retention, AppImage).  All of it is in, with these
decisions:

- **Sub-point grid: thm's.**  dev fixed the narrow-resonance tail error
  (d216517: graded zones, step 10 % of the distance, moving with E_R); thm had
  fixed it in 02ae97b (geometric tails, quantized anchors).  One definition is
  kept: thm's AdaptiveIntegrationGrid (section A), as the merge rules require.
  dev's new tests/reference/adaptive_grid_reference_test runs against it and
  **fails 1 of 14 checks**: "eV level crossing the window edge" swings by
  1.685e-4 of the area (tolerance 1e-4).  Every accuracy check passes (largest
  error 1.70e-4 against 2e-4; dev's grid 7.2e-5) and the swing with the level
  inside the window is 2.2e-7: the failing swing is thm's constant tail bias
  (~1.7e-4 of the area, r^2 Gamma/(2 pi core), see TailFactor) leaving the
  window with the level, a smooth function of position (largest jump between
  neighbouring positions 9.4e-6), not the oscillation the test was written
  for.  Neither the test nor its tolerance was changed.  To decide: keep thm's
  grid and relax/split that check, or tighten thm's tailRatio (which moves
  every THM fold), or take dev's grid (which moves THM outputs).
- **Gaussian kernels to 5 sigma for every segment (dev 7d24c88).**
  TargetEffect::thmConvolutionRange stays (also 5).  This is the only change
  to THM outputs (direct segments of THM projects), see below.
- **Shared `<targetInt>` lines.**  dev now converts a line's lab energies with
  each point's own pair (170ab0c); thm's per-segment copies (76888ad) are kept
  for their bookkeeping but no conversion happens on them any more (the
  in-place ConvertSigmaToCM is gone with dev's change).  dev's density fit
  follows thm's screening scheme: `density_<n>` is the line's (n = its place
  among active lines), shared by every copy (EData::SetUpDensityParams walks
  the copies by line; GetDensityTargetEffect like GetScreeningTargetEffect);
  **differs from dev** only for a densityfit line listing no existing segment
  (no parameter).  Five dev call sites that read `GetTargetEffect(effect
  number)` (MCMC screening/density priors, labels, LM penalty rows, API
  records) use the line's first copy instead.
- **Parameter layout / type codes.**  ..., screening, **density (dev, type
  8)**, THM cbkg last, **type 9 (was 8)**.  ParamKind Density before
  ThmCoherent; MCMC PARAM_DENSITY 8, PARAM_CBKG 9; pyazr `_KINDS` 8 density,
  9 cbkg; `bands._BAND_KINDS = (0, 1, 9)`.  GetParameterInfo pushed density
  after cbkg (auto-merge) with both type 8 -- fixed.
- **Shift penalties from zero (dev 0dd28e5)**, with thm's explicit
  `segment_N_energy_shift prior_centre c` row still honoured: the constant
  term's centre is 0 unless such a row moves it (ESegment's nominal shift now
  starts at 0, not at the line's value; AddShiftPenalty, ForEachShiftPenalty,
  shifts.out, MCMC priors and pyazr `penalties()` / `nominal_shift` follow).
  sqrt(E) and shift-function terms: from 0 (thm has no centre rows for them).
  save_fit writes a shift centre row only for an explicit non-zero centre.
  Norm priors keep thm's prior_centre (`nominal_norm`).
- **Same bug, two fixes, one kept:** nuisance priors by full index (thm
  953486f's `NuisancePrior(f)`; dev's ce4462b loops dropped, `density_` added
  to the skipped names); param.fit `#parametrization` (one line, dev's
  comment); `--help` keeps thm's complete list (`--use-api` included); exit
  status 1 (dev; thm returned -1); EC just above a Coulomb threshold (dev's
  interpolation 3e9b529, with its tests/ec_threshold; thm's mirrored energy
  dropped); GUI duplicate segment lines (dev ffc86b3, thm's `fromFile` flag
  dropped; the duplicate test now also compares isTHM); MCMC tab by name
  (thm, with dev's density category).
- **AMatrixFunc adjoint:** dev's `integratedBar` refactor takes thm's
  (1 + delta_12) of an identical entrance pair for the reported
  angle-integrated cross section and 1 for the angular-coefficient sigma_int
  (the forward A_0 is unscaled).
- **ShftFunc:** dev's relativistic eta_b/rho_b (Kinematics) and its table-edge
  rule (forceExact_ over the whole stencil) on thm's structure (U ratio with
  the linear continuation to S(0)); EnergyDerivative keeps the caller's step
  (ChannelFunc's rule for closed channels).
- **ESegment:** dev's rigid sub-point shift (f56e02a); a THM segment's shift
  stays unfloored.
- **EData:** thm's THM floor (minSubEnergy) and 5-sigma THM fold window on
  dev's per-point GridExtent/SigmaCM; dev removed the zero-thickness branch
  (thick target = density 0), so thm's copy of it went too.
- **New refusals (no silent fallbacks):** isDiff 19 (THM + angular
  coefficients) on a data segment, isDiff 9/19 on a test segment, a
  covariance block on a THM segment (its chi^2 is the profiled experiment's),
  and a THM segment with `relativistic=1` (the THM vertex/line shape/HOES
  kinematics are non-relativistic).
- **AZURECalc:** PrecalculatePoints is dev's CalculatePointsInParallel (failed
  point NaN) restricted to THM segments for thmOnly; the covariance
  correction after thm's thmOnly return; the shift fast path skips THM
  experiments and covariance groups.
- **AZUREAPI:** thm's FillFromFullRWA/EvaluateFilledChi2 structure with dev's
  parallel point pass and covariance correction; CalculateChi2Physical returns
  +inf for a rejected transform and adds the same Park/abundance/
  scattering-length terms as CalculateChi2RWA (new FilledModelPenalties; thm
  added none); calculate_residuals_rwa whitens covariance rows; the residual
  Jacobian takes dev's FillShiftColumns (with `skip`) unless the project has
  THM segments, which keep thm's forward differences of EvaluateFilledChi2
  (THM rows profiled), now also for density columns.  RefreezeCovariance
  (dev) added; every pybind ConfigScope carries the session's potentials.
- pyazr 2.9.0 (dev's number; thm's public names were already in 2.8.0).
- Tests: expected files -- see below.  tests/pyazr/CMakeLists.txt: dev's new
  tests registered before thm's SKIP_RETURN_CODE block (still last).
  park_formalism and convolution_unify check.sh: thm's versions (dev's 2b
  check is a subset of thm's; both hold the grid tolerances with thm's
  reason).
- GUI (not compiled here: no Qt/Qwt): both Runtime Options boxes (THM,
  relativistic), THM tick in the data dialog with dev's l/s widgets moved one
  column, density and cbkg rows in FittingTab in engine order, MassTable /
  MassCheckDialog sources, round-trip writer for component angles/scalings.
- Docs: experimental_effects describes thm's grid (and that dev's graded tails
  are not used here), 5-sigma kernels, the exact range integral; chi_squared
  the shift centre rule; pyazr the kind table and bands; command_line the
  exit status; contributing is dev's page plus thm's test-infrastructure
  sections.  Sphinx is not installed here: the docs were not built.

Output differences, measured (OMP_NUM_THREADS 1-2, Coulomb table on; merged
= this tree, thm = 30d3bd9 built in thm-ref, dev = 2d97b1f in dev-ref; 23
regression projects, 8 THM examples):

- THM regression projects (17O, 18O_p_a_thm, 6Li_d, 7Li_p_a) and the THM
  examples n15_lacognata2007, o17_guardo, o18_lacognata2008/2010: every
  output file byte for byte as thm.  f19_pag_thm: identical but for
  normalizations.out (dev cebec4f: no malformed row per fixed norm).
  c12c12_tumino2018 111.613 -> 111.661, li6_pizzone2011 157.036 -> 157.041,
  li7_tumino2006 439.148 -> 439.149: only the direct (non-THM) segments move
  (c12c12 segments 5-8; THM segments 1-4 identical); a merged build with
  convolutionRange = 3 gives thm's files byte for byte on all four, so the
  whole difference is dev's 5-sigma kernel.
- Classic against thm: identical on 15 projects; differences:
  13N segments 9/10 1841.85 -> 1875.78, 684.258 -> 676.017 (dev f389d46,
  photon Jacobian of lab-frame capture; dev's own values exactly), total
  336454 -> 336480; target_effect_ranges 3838.66 -> 3837.11 (the merged build
  with convolutionRange = 3 gives thm's files byte for byte); thick_target_
  resonance 30.8444 -> 30.9094 (dev 3207ae1, exact range integral: the merged
  build with the first-order loss gives thm's files byte for byte);
  hybrid_potential segment 3 2239.69 -> 2238.95 and the formal parameters in
  parameters.out by up to 2.7e-3 (dev 57e31bb, hybrid S on a grid node with an
  O(dr^4) derivative: parameters.out equals dev's); channel_space
  normalizations.out (cebec4f); intEC.dat in the 17-digit format (dev
  479153b, signature v5).
- Classic against dev (19 classic projects): identical on 9; the rest are section A
  (13N segments 1/2 814.405 -> 814.349 and 39.7064 -> 39.7053, beam_profile_
  kernel 40.4396 -> 40.4393, hybrid_potential 2.46215e6 -> 2.46612e6,
  7Li_p_ay-type parameters.out in the 7th digit, 13N_capture_ay 6e-9,
  beta_delayed_alpha 9e-9, shift_subpoint_angle 4e-10, intEC.dat node counts),
  plus thick_target_resonance 30.7201 -> 30.9094 and target_effect_ranges
  3839.73 -> 3837.11 (thm's grid under dev's kernel and thickness) and
  natural_neon 3342.46 -> 3342.51 (segments 4 and 6 by 4e-4 and 2.5e-4; the
  thm pin's segment 6, 247.1, was the 76888ad difference, gone with dev's
  per-point conversion).
- Re-pinned (rows explained above, each equal to a measured run):
  tests/13N (segments 9/10 and the total: dev's values), tests/thick_target_
  resonance (30.9094, README updated).  **Not re-pinned:** tests/natural_neon
  keeps thm's pin and fails (3342.51 against 3346.1, 1.07e-3; segment 6
  245.029 against 247.1): the move is dev's 5-sigma kernel (-0.1 % on dev) and
  per-point conversion, but it was not isolated by experiment here; proposed
  new pin = the merged run above.  target_effect_ranges stays within its
  tolerance of thm's pin (4e-4) and was not touched.

Tests (this machine: GCC, no Qt; another job held ~12 GB, so heavy runs were
capped at 4.5 GB of address space):

- Build (GCC 11, Release, BUILD_GUI=OFF, pyazr on): no warnings.
- ctest -LE slow (86 tests, -j 2): 85 passed, 1 failed:
  adaptive_grid_reference (the window-edge swing above).  The pyazr tests
  ran their engine parts (pyazr imports, 2.9.0; e.g. target_density,
  thm_adjoint, thm_park, park_gradient, scattering_length print ok).
  The slow label was not run.
- tests/run_tests.sh (all regression projects and check.sh, every engine
  run under 4.5 GB): 118 passed, 3 failed -- natural_neon (pin, above),
  convolution_unify (bad_alloc under the cap; it needs more than 3 GB with
  either binary, as at the previous merges; not run uncapped here because
  of the other job), cli_help (dev's `--ec-convection-current` was accepted
  but not printed, and the combined help/docs row for the three positive
  defaults did not name `--grid-replacement` on its own line: fixed, one
  line per flag in --help and in command_line.rst; cli_help and the
  consistency test then pass).  Every THM check.sh passes, as do dev's new
  ones (covariance, relativistic, transmission_yield, shift_penalty_zero,
  target_density, ec_threshold, convolution_window, ...).
- THM check.sh work files were not compared file by file against thm (only
  pass/fail); the regression projects and the examples were.
- Not checked: the GUI (Qt/Qwt are not installed: its files were merged by
  reading only), Sphinx, macOS/Windows, the slow pyazr tests, the 12C+12C
  best fits outside the repository.

## Merge of origin/dev 0919dcb (October 2026)

dev between 8673c5c and 0919dcb: R. deBoer's Screening sub-tab of the
Fitting tab (334f104: `screening_<n>` rows, synced both ways with the
Experimental Effects tab) and tests/screening_6li (657b14b, the 6Li(p,a)
electron-screening benchmark on Cruz et al. 2005); and the Park speed-up:
Park rows of the residual Jacobian by free parameter (fc88274:
`ChannelSpace::CurrentLevelDerivatives` takes the formalism, `LevelAdjoint`
and `GroupRows` serve Park), the J > 0 wall in `Chi2On` (3ea6020), Park's
gradient from the residual Jacobian (f61690a: `GradientFromJacobian` adds the
wall's gradient), d2S/dE2 only when a gradient asks for it and Park's dS/dE
table read from the levels (ecbca5e: `CalcShiftFunctions(config, false)` for
a chi-squared alone, `ALevel::ShiftSecondDerivativesValidAt`;
`BuildShiftDerivTable` reads `GetShiftDerivative` under Park), the shift-
function derivatives for all levels in parallel and without the Coulomb memo,
shared through an exact-argument memo (45cc5cf), docs (0919dcb).  All of it
is in, with these decisions:

- **CNuc::CalcShiftFunctions: dev's structure, thm's functions.**  S is
  thm's `ChannelFunc::Shift` (continuous through threshold), recomputed only
  for a level that moved, as before.  dev's deferred, parallel derivative
  pass is kept (channel-spin sharing, `g_shiftDerivMemo`, the lazy S''), but
  S' there is thm's `ChannelFunc::ShiftDerivative` -- the dS/dE of every
  other caller on thm, `BuildShiftDerivTable` under Brune included (row
  ec416de) -- and S'' the branch the level is on (ShftFunc below threshold,
  `PEShift_d2E` above), as thm had it.  Like dev, neither uses the Coulomb
  memo: `ChannelFunc::SetUseGlobalCache(false)` (new; it switches the
  CoulFunc it makes).  So under Park dev's `BuildShiftDerivTable` shortcut
  reads the same function at the same energy that it would have
  differentiated, without the memo.
- **Brune (and every non-Park run) is untouched**: no derivative is taken
  there, S is computed as before, and `BuildShiftDerivTable` still
  differentiates with `ChannelFunc` and the memo.  Every Brune output below
  is byte for byte 031793d's.
- **Park round-off (45cc5cf).**  Park's S' no longer goes through
  `CoulFuncCache` (which matches energies within 1e-12 MeV and depends on the
  run's history), so Park values move in the last digits.  Measured: of
  the 18 THM check.sh scripts only tests/thm_park writes different numbers,
  four files of its f19_pag_thm runs (thm_experiments.out b0 2.1e-11,
  AZUREOut_aa=1_R=6.out 4.8e-11 relative: the 11th printed digit); its
  Park = Brune comparison for f19_pag_thm reads max rel 2.144e-11 (was
  0.000e+00) and the file round trip 8.267e-08 (was 8.262e-08), against
  tolerances of 1e-8 and 1e-5.  Its fits, MCMC runs and every other Park
  output (17O, 18O_p_a_thm with each feature, the on-shell examples,
  c12c12_tumino2018 under Park) are byte for byte.  No pin or tolerance
  changed.  tests/pyazr/thm_park_test.py prints the same numbers with both
  builds.
- **Park + THM gradient.**  `GradientFromJacobian` no longer declines Park,
  but it still declines every project with a THM segment:
  `ComputeResidualJacobian` (called with skipTHM = false by
  `AZURECalc::ResidualJacobian`) refuses THM segments, so `Gradient` takes
  the adjoint for the classic segments and the THM-only differences
  (Chi2On thmOnly, which adds no wall: the wall's gradient is the adjoint
  path's, AddParkPenaltyGradient), as before.  pyazr's
  `residual_jacobian`/`chi2_and_grad` (skipTHM = true, the THM rows from
  `ComputeTHMRows`) now get Park's classic rows from `LevelAdjoint`/
  `CurrentLevelDerivatives` (mode 2) instead of the bilinear fallback.
  Checked with `AZURE_GRAD_CHECK=exit` (first MIGRAD gradient, every free
  component against central differences, one thread), Park and Brune, both
  builds: tests/17O (8 widths; with the 4 energies free too) and
  examples/f19_pag_thm (25 widths; with 11 energies free too).  Widths agree
  with the differences to 8e-7 (largest |a - fd| 4e-8 of the largest
  component) in both modes and builds.  Two f19 energies (J3 L4, J6 L1) read
  1.8e-2 (Park) and 1.6e-1 (Brune) off at the default step (1e-6 relative)
  with both builds; at 1e-8 they agree to 1.6e-9 and 2.1e-8 (Park), 4.7e-8 and
  5.3e-7 (Brune): truncation of the difference near narrow levels, not the
  gradient.  The analytic gradients of the two builds are equal bit for bit
  under Brune and for 17O under Park; for f19 under Park the first gradient
  comes after the seed's 610 evaluations, which the 1e-11 round-off above
  steers to a slightly different point (1e-5 in the gradient there, the
  differences moving with it).  One more, in both builds: under Park the
  energy of tests/17O's first 1/2- level reads 1.1e-4 off the differences
  at either step (Brune: 4e-8), within dev's Park tolerance (1e-4 of the
  largest component in tests/channel_space and park_formalism: S' and S''
  are numerical derivatives); it is 031793d's value bit for bit.
- `Chi2On`: dev's wall (3ea6020) was already on thm (dcb8064); thm's
  `thmOnly` early return stays before it, comments merged.
- **GUI**: the Screening sub-tab sits between Energy Shifts and thm's
  THM Coherent Background tab; parameter order screening, then cbkg, as the
  engine (`appendScreeningParameters` before `appendCoherentParameters` in
  both `populateFromCurrentGUIState` and `loadSettings`); screening rows have
  no prior-centre column (norms and shifts only).  `getTargetIntTab` was on
  both sides (kept once).  TargetIntModel's columns 35-38 are the same on
  thm.
- tests/screening_6li/check.sh takes thm's time guard (tests/lib/guard.sh).
- docs/source/developer/level_matrix.rst: S' is ChannelFunc's on this
  branch; a THM project keeps the adjoint gradient in every formalism.

Verified after the merge (one heavy process at a time, OMP_NUM_THREADS=1,
engines under a 1.6 GB address-space cap; reference: 031793d built from the
same sources, tag thm-premerge-2026-10-22c; dev: 0919dcb's own build):

- THM: every THM example and THM test project (13 projects, 73 output files,
  86 with scatteringLengths.out) byte for byte as 031793d; the 12C+12C best
  fits too (33 files, 36): plain_6.5 92.2140340427, plain_8.5 70.7968723529,
  NCxR_8.5 100.894445567 (totals with priors), identical.  The 18 THM
  check.sh scripts pass with both binaries; of their 5220 work files 78
  differ: run.azr/thm_experiments.out with the temporary directory names,
  logs in timings only, and the four tests/thm_park files above.
- Classic (21 classic test projects, natural_neon apart): 113 of 113 files
  byte for byte as 031793d.  Against dev 0919dcb: the same 29 files as at
  the previous merge (section A), and dev 0919dcb gives dev 8673c5c's files
  on all of them.
- Suites: build/ ctest 61/61 (pyazr_thm_adjoint, pyazr_thm_park,
  pyazr_park_gradient, pyazr_thm_profiled included), build-gui ctest 82/82
  (dev's screening_fitting included), tests/run_tests.sh
  under 2.4 GB 67 passed, 0 failed (convolution_unify, natural_target and
  natural_neon left out: memory); tests/natural_target/check.sh passes on
  its own under 3.0 GB.  dev's tests/screening_6li (new), channel_space,
  park_formalism (its new wall-gradient check: both paths), 13N,
  scattering_length and electron_screening check.sh pass on thm.
- Speed (one thread, same machine): CLI evaluation f19_pag_thm 8.3 s
  (031793d 8.3-8.4 s), c12c12_tumino2018 45.9 s (46.0 s); the CLI's first
  gradient on f19 1.02 s (1.02-1.04 s); pyazr chi2_and_grad (HOES adjoint)
  f19 1.19-1.23 s (1.20-1.21 s), c12c12 (230 parameters) 1.76-1.77 s
  (1.78-1.79 s), same chi^2 and gradient sums.  THM takes no Park speed-up
  (none intended: THM stays on the level route).

## Merge of origin/dev 8673c5c (October 2026)

dev between 7e6dba8 and 8673c5c: R. deBoer's scattering lengths (dfd15a4,
448244d, 79d06ae: `scatteringLengths.out` after every calculation, the
`<scatteringLengths>` block of measured values, fittable with the adjoint),
electron screening as an experimental effect (ff356ab: `screening Ue vary
sigma` on a `<targetInt>` line, parameter `screening_<n>`), the segment
dialogs' compact observable list (0c0f8eb), tests/shift_fn, the segfault of
fits with external capture (2dcb846; thm had it since the merge of 7e6dba8);
and the speed-up: the channel-space level-matrix route (e3ff2a9, 487e194,
3ae98f2: formal and pencil routes, near-pole deflation, chosen per J-group and
energy; `AZURE_LEVEL_SPACE=1` forces the level route, `AMatrixFunc::
SetLevelSpaceOnly(true)` one object), differential sums by Legendre order
once per T-matrix (1e7eced) and their cotangents (2bc807f), the deferred
level-matrix adjoint and the Jacobian rows by free parameter (8546dfa,
a22a2a1, 095bcc2), `AZURE_PROFILE` timers (c84ba60), the ALevel stamp
(875f508), tests/channel_space and tests/13N/check.sh.  All of it is in, with
these decisions:

- **THM stays on the level-space route.**  `THMMatrixFunc` is switched to it
  in its constructor (`SetLevelSpaceOnly(true)`; `Reset` keeps it), so every
  THM point -- the forward HOES amplitude (`GetAMatrixElement` over the active
  levels) and the HOES adjoint (`GetAMatrix`, `LevelMatrixABarAdjoint`) --
  is solved by filling and LU-factoring the level matrix, as before; no THM
  object ever holds a channel-route group, so the adjoint always has the
  factorization it reads.  `THMMatrixFunc` forms its observable from its own
  amplitudes and never calls `GenMatrixFunc::CalculateCrossSection`, so the
  Legendre-order sums (`rt_wl_`), the Coulomb-interference table and the
  per-order T-cotangents (`tbar_q_`, `PointTBar`) are not on any THM path;
  `ThmPointModelGradient` and `ComputeTHMRows` use `LevelMatrixABarAdjoint`
  -> `LevelMatrixAdjoint` directly (structural, not the deferred `lmK/lmN`
  path of a scalar pass).
- **Classic segments of a THM project take dev's routes** (channel space where
  it is chosen, the Legendre-order sums, the grouped and deferred gradient):
  every output of every THM example and test, and of the 12C+12C best fits,
  is byte for byte what it was (below), so nothing forces them back.
- `LevelMatrixAdjoint(..., lambda, levelOf)`: dev's kernel argument `lambda`
  first, then thm's `levelOf` (active level -> level number).  Without
  `levelOf` the level data come from dev's per-parameter-set cache
  (`ChannelSpace::CurrentLevelArrays`); a THM block with levels left out of
  the R matrix gathers them over its active levels (`lma_active_`), the same
  arithmetic as before.
- **Parameter layout and type codes.**  Levels, abundances, norms, shifts,
  sqrt(E), shift-function parameters, **screening potentials (dev)**, THM
  cbkg (still last).  `GetParameterInfo` type codes: 7 = screening (dev),
  **8 = cbkg (was 7 on thm)**.  `ParamKind` ... Abundance, Screening,
  ThmCoherent.  MCMC kinds: PARAM_SCREENING 7 (dev), PARAM_CBKG 8 (internal;
  was 7).  pyazr `_KINDS` {7: screening, 8: cbkg}, `bands._BAND_KINDS =
  (0, 1, 8)`.  thm's by-name MCMC classification knows `screening_<n>` (its
  prior from the `<targetInt>` line, as dev's BuildAutoPriors); the GUI's
  `MCMCTab::categoryOfName` too.  The nuisance priors by full index skip
  `screening_*` as they skip norms, shifts and abundances (dev skipped it by
  substring).
- **Screening parameters belong to `<targetInt>` lines.**  thm gives every
  segment listed on a line its own copy of the effect (76888ad), so effect
  numbers and lines differ there.  A screening parameter is the line's
  (`TargetEffect::GetLineNumber`, the place among the active lines, dev's
  numbering of `screening_<n>` and of `GetParameterInfo` field 6), shared by
  every copy: `EData::GetScreeningTargetEffect(k)` is the first copy, U_e is
  set on all of them, the prior is counted once.  **Differs from dev** only
  for a screening line that lists no existing segment: no copy, so no
  parameter (dev makes one that moves nothing).
- A THM segment carries no screening factor: the factor is applied in
  `GenMatrixFunc::CalculateCrossSection`, which `THMMatrixFunc` does not use,
  so a screening line on a THM segment moves nothing there (zero derivative).
- **EPoint at threshold.**  dev's exact-threshold branch (|E| <= 1e-10 MeV:
  S_l(0) from `ScatteringLength::ZeroEnergyShift`, P = 0), which the
  scattering lengths evaluate, comes first; thm's `ChannelFunc` branch (closed
  channels, the Coulomb threshold) keeps every other energy at or near
  threshold.
- pyazr: `CalculateChi2RWA` adds the scattering-length chi^2 (dev) on thm's
  `FillFromFullRWA`/`EvaluateFilledChi2`; the screening priors are not in it,
  as on dev (`penalties()["screening"]`, added by `objective`).
  `CalculateScatteringLengthsRWA`, `ScatteringLength`, `ScatteringLengthLine`
  and the AzrModel screening helpers are dev's.
- GUI: the THM tick lives in dev's new observable row of Add/Edit Data
  Segment (after the multipolarity box) and keeps its column in Add/Edit Test
  Segment; it is disabled by observable code (data: phase shift; test: phase
  shift and angular-distribution coefficients), as thm's index rules did.
  `<scatteringLengths>` is read after thm's THM content.
- Tests: dev's tests/13N, channel_space, electron_screening and
  scattering_length check.sh take thm's time guard (tests/lib/guard.sh);
  pyazr_scattering_length and pyazr_screening exit 77 without the package or
  the engine (thm's convention; they exited 0) and are registered before the
  SKIP_RETURN_CODE block.
- docs/source/developer/level_matrix.rst notes that THM is always on the
  level route.

Output: `scatteringLengths.out` (dev) is written by every calculation, THM
projects included; it is new, every other file is unchanged.

Verified after the merge (one heavy process at a time, OMP_NUM_THREADS=1,
engines under a 1.6 GB address-space cap; reference: 97e2e1a built from the
same sources, tag thm-premerge-2026-10-22b; dev: 8673c5c's own build):

- THM: every THM example and THM test project (13 projects, 73 output files)
  byte for byte as 97e2e1a; the 12C+12C best fits too (33 files): plain_6.5
  92.2140340427, plain_8.5 70.7968723529, NCxR_8.5 100.894445567 (totals with
  priors), identical.  Each run also writes the new scatteringLengths.out.
  The 18 THM check.sh scripts of 97e2e1a pass with both binaries and give the
  same 4323 work files, apart from 13 that hold the temporary directory names
  (run.azr, one thm_experiments.out); thm_park's MIGRAD and MCMC runs
  included.
- Classic (97e2e1a's 18 classic test projects, natural_neon apart): 87 of 87
  files byte for byte as 97e2e1a, although dev's channel-space route and
  Legendre-order sums are on for them.  Against dev 8673c5c: the 29 files of
  the previous merge's list, i.e. section A (thick_target_resonance,
  hybrid_potential, 13N and its EC table, beam_profile_kernel,
  target_effect_ranges, the 7Li_p_ay-type parameters.out/param.par in the 7th
  digit, 13N_capture_ay 3e-8, beta_delayed_alpha 1.3e-7, the last also with
  `AZURE_LEVEL_SPACE=1` on dev); dev 8673c5c gives exactly dev 7e6dba8's files
  on these projects.
- Suites: build/ ctest 61/61 (pyazr_thm_adjoint, pyazr_thm_park,
  pyazr_park_gradient, pyazr_thm_profiled, pyazr_scattering_length,
  pyazr_screening included), build-gui ctest 81/81, tests/run_tests.sh under
  2.4 GB 66 passed, 0 failed (convolution_unify, natural_target and
  natural_neon left out: memory; not run this time).  dev's tests/13N,
  channel_space, electron_screening, scattering_length, shift_fn, ec_shift
  and angle_mapping check.sh pass on thm.
- Speed (one thread, the same as 97e2e1a): CLI evaluation f19_pag_thm 8 s,
  c12c12_tumino2018 46 s; pyazr chi2_and_grad (HOES adjoint) f19 1.20 ->
  1.21 s, c12c12 (230 parameters) 1.83 -> 1.79 s, same chi^2 and gradient.

## Merge of origin/dev 7e6dba8 (October 2026)

dev between 64ed575 and 7e6dba8: external-capture energy tables (687105f,
ECTable, `--no-ec-table`; shifted capture points exact), the GUI's engine
numbering of `<parameterSettings>` level entries (be817bf, e36d485), and the
speed-up: point mapping across segments (ccdd364, 92d837e), interference
products once per T-matrix (deec7fc), the rank-one level-matrix cotangent
and flat structural loops of the adjoint (e8b1c08), the energy-grouped
analytic gradient and residual Jacobian (861d506), tests/angle_mapping and
tests/ec_shift.  All of it is in, with these decisions:

- **THM and point mapping across segments.**  A THM (HOES) segment is
  mapped within itself only, and no point of another segment (THM or
  classic) is mapped onto one of its points (`EData::MapData`: THM segments
  are `ownOnly`, like a shifted capture segment).  Extending
  `SharesCalculation` instead would have to compare every THM option --
  experiment, observable and isDiff, vertex model and the DW vertex,
  kinematics, the ps nodes and theta windows, distortion, line shape, cbkg,
  the THM binding energy of the pair, R(E), weights, the 5-sigma fold -- and
  any option added later would silently map unequal points.  Two code paths
  also assume a THM point's host is in a THM segment: the THM-only chi^2 of
  `AZURECalc::Gradient` (Chi2On thmOnly; PrecalculatePoints hands out the
  unmapped THM points only, so a host in a classic segment would not be
  recalculated) and `ComputeTHMRows` (it evaluates THM segments alone, mapped
  points after their hosts).  Within one THM segment every option is shared,
  so mapping there is what thm already did (one-angle THM segments with
  identical settings are not mapped onto each other; not needed so far).
  Classic segments of a THM project are mapped across segments as on dev.
- **Grouped gradient.**  `BuildGradItems` leaves THM segments out
  (`skipTHM`): `AccumulateEGammaGradient` always (their gradient is the HOES
  adjoint or the THM-only differences, as before), `ComputeResidualJacobian`
  refuses a THM segment unless `skipTHM` (the LM fallback, as before), and
  `ComputeModelGradients` honours `skipTHM`.  No group mixes THM and classic
  points, since MapData never maps across them.  `ThmPointModelGradient`,
  `ComputeTHMRows`, the cbkg columns, profiled norms and THM penalty rows are
  unchanged.
- **Adjoint.**  dev's `PointTBar`/`TBarAdjoint`/`BilinearAdjoint` and the flat
  `LevelMatrixAdjoint(Mbar)` are kept, with thm's identical-entrance
  (1 + delta_12) factor in `PointTBar`.  `LevelMatrixAdjoint` takes thm's
  `levelOf` (active level -> level number), and the HOES adjoint enters
  through `LevelMatrixABarAdjoint` (the dense inverse lemma, then
  `LevelMatrixAdjoint`).  The Park branch is dev's, unchanged; the Brune
  terms use the stored ShftFunc/ChannelFunc shift values and thm's dS/dE
  table, as before; gammaIsRWA acts in TransformIn, before the adjoint.
- **Interference-product cache** (deec7fc): only `GenMatrixFunc::
  CalculateCrossSection` reads it; `THMMatrixFunc` forms the HOES observable
  from its own amplitudes and never goes through it, and the cache is dropped
  at every `ClearMatrices`/T-matrix store, so no THM point reads products of
  another observable.
- **CopyEDependentValuesFrom** still copies the THM tables; with THM
  segments mapped within themselves a THM point's host is in its segment.
- **EC amplitudes.**  thm's closed-entrance-channel guard (zero amplitude
  where P = 0) moved into dev's `ECPathwayAmplitude`, so it also covers
  shifted and tabulated amplitudes.
- **GUI**: dev's engine numbering of level entries replaces thm's
  `engineLevelOrder`/legacy reader (see B, f3210c7); `engineNumbering` skips a
  level without channel lines, as thm's order did; thm's units, round-trip
  level writer, prior centres and THM members are kept.  dev's
  `getPairsTab`/`getSegmentsTab` replace thm's copies.
- `--help` prints `--ec-table` (dev parses it; tests/cli_help).
  tests/angle_mapping and tests/ec_shift take thm's time guard.

Verified after the merge (one heavy process at a time, OMP_NUM_THREADS=1,
engines under a 1.6 GB address-space cap; reference: c5c747d built from the
same sources):

- THM: every THM example and THM test project (13 projects, 73 output files)
  byte for byte as c5c747d; the ¹²C+¹²C best fits too (33 files): plain_6.5
  92.2140340427, plain_8.5 70.7968723529, NCxR_8.5 100.894445567 (totals with
  priors), identical.  The 18 THM check.sh scripts of c5c747d pass with both
  binaries and give the same 4838 work files, apart from temporary paths in
  logs and run files and the analytic gradient's last bits (dev's flat
  structural loops): thm_band's band 2.7e-11, thm_park's two MIGRAD fits end
  with the same chi^2 to 12 digits and models within 7e-10 (a flat
  direction: a width 1e-4 apart).  The first gradient of that fit agrees with
  c5c747d's to ~1e-15 relative.
- Classic (c5c747d's 17 classic test projects, natural_neon apart): 71 of
  78 files byte for byte; the other 7 are 13N, 13N_capture_ay and
  beam_profile_kernel, dev's EC tables (687105f: 13N_capture_ay 1e-5 in one
  small capture value, chi^2 3.73951 -> 3.7395); with `--no-ec-table` they
  are byte for byte too.  Against dev 7e6dba8 the differences are those of
  section A only (the previous merge's list).
- Suites: build/ ctest 59/59 (pyazr_thm_adjoint, thm_park, park_gradient
  included), build-gui ctest 77/77, tests/run_tests.sh under 2.4 GB 62
  passed, 3 failed: cli_help (fixed above, rerun: pass), natural_target and
  convolution_unify (out of memory, as with dev's binary).  natural_target
  passes under 3.0 GB; convolution_unify's first model needs more than
  3.0 GB (known).  natural_neon (project and check.sh) passes under 2.4 GB.
  tests/angle_mapping/check.sh and tests/ec_shift/check.sh pass on thm.
- Speed (one thread): CLI evaluation f19_pag_thm 8.2 s, c12c12_tumino2018
  45.9 s, both as before; pyazr chi2_and_grad (HOES adjoint) f19 1.78 ->
  1.29 s, c12c12 (230 parameters) 2.81 -> 1.99 s, same values; the CLI's
  first gradient on f19 1.28 -> 1.02 s.

## Bugs from origin/dev fixed on thm (October 2026)

Found after the merge of 6f3228e; each came with dev's code, so `dev` has it
too and can take the commit (cherry-pick, or the same change by hand where
the THM lines differ).  Verified against c4890ab: the 133 output files of
every tests/ and examples/ project byte for byte (run logs differ in timings
and progress bars only); c4890ab's 29 check.sh scripts give the same 2409
files, apart from temporary names and the MCMC runs of tests/thm_park, which
are the cases fixed: its 18O THM project (two profiled norms) printed
"parameter classification produced 12 entries but there are 10 parameters"
and now gets the 0.5 keV energy spread (other samples), and its J < 0 start
is now refused (the check expects that since baa6642).

| commit | bug (on dev too) | fix | test |
|---|---|---|---|
| 53154cb | Every CLI prompt loop (menu, MINOS variance, reaction-rate pairs and temperatures, MCMC walkers, steps, threads, RWA) printed its prompt forever at end of input: one MCMC run wrote 10 GB | `src/AZURE2.cpp`: `readAnswer` stops with "ERROR: the input ended (end of file) at the prompt ..." and exit 1; prompts whose blank answer is a default take it, as before; the temperature-file prompt stops too (it returned with no file) | tests/cli_eof (c4890ab: 26 failures) |
| fcf3a24 | MCMC `BuildAutoPriors` re-walked the levels and segments: a free THM norm (profiled, not a parameter) got an entry, a coherent background none, so the lengths differed ("12 entries but there are 10 parameters"), and the automatic priors and the level-energy spread were dropped (walkers 1 MeV apart).  **On dev** the walk is the same but cannot mismatch (no THM); the sqrt(E) block is classified correctly in both.  `AZURELabel::Parameter` had the same walk | kinds from the parameter names; `PARAM_SHIFT_SQRT` = 4, `PARAM_CBKG` = 5 (GetParameterInfo's codes); cbkg: no automatic prior, a user prior kept; labels skip profiled norms and name cbkg | tests/mcmc_classification (c4890ab fails the two THM runs) |
| eb7deb6 | MCMC under `--use-park` from a point with J <= 0: no word, acceptance 0, the chain is the start | the start ensemble is evaluated first; levels with J <= 0 at the start are named with their J; with no walker of finite posterior the run is refused ("MCMC not started", no samples.mcmc).  `RunMCMCSampling` returns bool; an MCMC failure exits the CLI non-zero (it exited 0); the GUI gets the report through `SetGUIMessageCallback` | tests/mcmc_park_start (c4890ab: 11 failures) |
| dcb8064 | `AZURECalc::Chi2Value` lacked the Park J > 0 penalty that `operator()` adds: the finite-difference fallback of `Gradient()` (all of it when the adjoint bails) differentiated another function beyond the wall, and LM / GSL-LM took the penalty-free value as their cost | penalty added (not to the THM-only part, which adds to an analytic gradient that has it).  Dev's `Chi2Value` has no `thmOnly`: add it unconditionally under Park | tests/reference/park_chi2_fallback (ctest; c4890ab: 4 failures) |
| 90cfb62 | pyazr_sqrt_shift, pyazr_park_gradient, pyazr_sav_readers were registered after the SKIP_RETURN_CODE 77 block and exited 0 (or crashed) without numpy or the engine.  **On dev** there is no 77 convention at all (the block is thm's) | block moved to the end of tests/pyazr/CMakeLists.txt; the three exit 77 | numpy / `_azure2` hidden: 77 in all six cases |
| 8d6f382 | `AZURECalc::WriteParameters` wrote `param.fit` (the 100-evaluation snapshot of a fit) without the `#parametrization` line that 8cdc841 gave `param.par`/`param.sav`: a `--use-park` run that read it back took its Park amplitudes for Brune's and converted them (the p+p Park fit of tests/park_formalism resumed at chi2 1696.44 instead of 103.129).  **On dev**: same WriteParameters (6f3228e) | the line written as WriteUserParameters writes it | tests/park_formalism section 2b (7f8037e fails its 4 checks) |
| 9490640 | `--help` printed `--no-long-wavelenth` (the parser takes `--no-long-wavelength`) and lacked `--help`, `--covariance-band`, `--scale-covariance`, `--use-api`.  **On dev**: the same help text and parser (6f3228e) | help lists exactly the parser's flags | tests/cli_help (help, parser literals, command_line.rst) |
| 482eabe | GUI MCMC tab: kinds by substring, so `segment_N_energy_shift_sqrt` showed the energy shift's automatic prior (the sampler overrode it).  **On dev**: same code | `MCMCTab::categoryOfName`, by the name's form as BuildAutoPriors | tests/gui/mcmc_priors |
| 122bd12 | GUI MCMC tab, Load RWA Parameters: `param.par` read by row position, so the `#parametrization` line became the first parameter's value and every value moved a row.  **On dev**: same reader | `ReadUserParameters` (by name) + `ReconcileBasis` | tests/gui/mcmc_priors |

Results that move: none for a converging run.  Behaviour: a CLI fed short
input now exits 1 instead of looping; a refused or failed MCMC run exits
non-zero; under Park, LM's cost includes the J penalty (zero where J > 0).

## Merging thm into dev later

1. `git fetch`; merge the newest origin/dev into thm first (as 3628a87 did),
   not the other way round, and get CI green on thm.
2. Re-read sections A-D: each row is a deliberate difference that dev takes
   with the merge.  Re-pin dev's tests/*/expected that A moves (the thm pins
   are already moved; a new dev test may need it).
3. Parameter layout and type codes: classic as dev; check that a newer dev
   did not add parameters after the screening block (the cbkg block must
   stay last) or reuse type code 8 (cbkg; 5, 6 and 7 are dev's shift_fn,
   abundance and screening); the pair record keeps dev's fields first
   (binding energy last, field 19).  THMMatrixFunc must stay on the
   level-space route (SetLevelSpaceOnly, see the merge of 8673c5c).
4. Engine seams to re-check: EData.cpp/EDataThm.cpp hooks, EData::MapData
   (THM segments map within themselves only), the energy groups of
   AZUREGrad.cpp (THM segments left out), AMatrixFunc::LevelMatrixAdjoint
   (levelOf) and LevelMatrixABarAdjoint (the HOES adjoint), EPoint THM
   tables (CalcEDependentValues, CopyEDependentValuesFrom, mapped points),
   ReplaceSubPointGrids and THM sub-points, AZURECalc (Chi2On thmOnly,
   Gradient)/AZUREGrad/AZUREAPI THM paths, CNuc TransformIn (gammaIsRWA,
   Park), ALevel per-channel vectors (gamma_is_rwa_), ChannelFunc and
   ShftFunc, CNuc::CalcShiftFunctions (ChannelFunc's S and S' inside dev's
   parallel derivative pass, no Coulomb memo there; see the merge of
   0919dcb), the Coulomb-table keys of the THM memos, Config.h flags;
   AZURECalc::GradientFromJacobian must keep declining THM projects (or
   learn their rows) now that it takes Park.
5. GUI: SegmentsDataData/SegmentsTestData fields (isTHM last before any new
   dev field), level lines 32/33, parameterSettings rows and prior_centre,
   the Runtime Options box; open + save of a classic project must stay byte
   for byte (tests/gui/thm_opt_in, parameter_settings, level_precision).
6. pyazr: `_KINDS`, `_NFIELDS` (Parameter 16, Pair 20), `penalties()`, the
   sav readers (by name), save_fit (both verifications), mode-dependent
   vector lengths.
7. Tests and CI: tests/pyazr/CMakeLists.txt (AZURE2_BIN, the SKIP 77 block
   last), tests/reference, tests/gui, .github/workflows/build.yml (one Python
   pinning block); run build/ and build-gui ctest and tests/run_tests.sh under
   the memory cap (convolution_unify, natural_target and natural_neon need
   more than 2.4 GB of address space with dev's binary too).
9. A test toy model must be one compound nucleus per level (dev's rule).
8. THM + `--use-park` is allowed (see "THM under Park"); dev has neither
   THM nor gammaIsRWA, so only the adaptive-grid row there touches dev's
   classic Park projects.

## THM under Park (October 2026)

The refusal of THM segments under `--use-park` (EData.cpp, data and test
segments) is lifted.  Park's level matrix is Brune's with diag(sqrt J) on both
sides and the HOES amplitude is bilinear in the amplitudes with level-diagonal
factors only (vertex boundary, line shape), so it is invariant; the coherent
background, kinematic factors, R(E), the DW vertex and the windows do not see
the amplitudes.  Three readers of the amplitudes did not know the
parametrization and are fixed:

| change | before | effect |
|---|---|---|
| ThmLevelWidth (line-shape pole): Gamma = sum 2 P gamma^2 under Park; the parametrization is part of its memo key | Brune's 2 P gamma^2 / (1 + sum gamma^2 dS/dE) of Park's amplitudes: Gamma low by J | 18O_p_a_thm + lineshape: 2.6 % in the model |
| gammaIsRWA under Park: Brune's amplitude in the file (row above) | read as Park's | 18O_p_a_thm (all channels amplitudes): chi2 2228.94 vs 2630.55 |
| AdaptiveIntegrationGrid, fit-time anchors (EPoint::RefreshSubPointGrid): GridConfig::parkAmplitudes gives the widths of the same Brune parameters (particle sum / (J + sum_open gamma^2 dS/dE), radiative / J); the parametrization is part of the anchor memo key | Park's amplitudes divided by Brune's 1 + sum: a narrower width, other quantized anchors, another sub-point grid | f19_pag_thm 1e-5, li7_tumino2006 4e-4, tests/17O 2e-5 in the model.  **Differs from dev** for a classic Park project with target effects whose levels have J well below 1 (dev's grid is the "before" column) |

Brune and classic results are unchanged against e591d33: the 133 output files
of every tests/ and examples/ project byte for byte, and the logs of the 28
other check.sh scripts identical apart from temporary names and timings.  Under Park every THM
example and test project gives Brune's model to 1e-8 (most bit for bit;
f19_pag_thm 4e-10, levels with J down to 0.02); tests/7Li_p_a and tests/6Li_d
have input widths beyond Brune's reach ("Denominator less than zero"), so
their .azr is a different model in the two modes and they are compared
through converted parameter files (8e-8).  tests/thm_park/check.sh,
tests/pyazr/thm_park_test.py (THM Jacobian under Park: finite differences and
the chain rule through J), tests/gui/thm_opt_in (the Park and THM boxes
together) and thm_workspace (GUI diagnostics in both modes).

