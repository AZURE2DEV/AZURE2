# thm -> dev: differences for classic projects

For the merge of `thm` into `dev`. "Classic" means a project without a `<thm>`
block and without a THM segment (isDiff >= 10). Everything listed here changes
something for such a project. The THM machinery itself (HOES cross section,
`<thm>` keys, experiments, R(E), the DW vertex, the ps and angle windows,
cbackground=, thm_experiments.out) is new and does nothing unless it is used.

`origin/dev` 6f3228e is merged into `thm` (3628a87, October 2026; see
"Merge of origin/dev 6f3228e" below), so all of dev is already in `thm`.
The dev changes merged by 18c49ea, f56c150, da72def, 497548a, e191991 and
3628a87 are not differences: the Wigner-limit bound, the elastic identical-pair factor
(7c34992), the adaptive-grid convergence work, target-integration anchoring,
the beam profile, the LM energy-shift Jacobian, CR line endings and the EC
cache guard.

## A. Results that move

| commit | change | reason | re-pinned (tests/*/expected) |
|---|---|---|---|
| b6cc41b | ShftFunc below threshold from the U(a,b,z) recurrence with an exponent carry; exact S_l(0) at E = 0; linear interpolation for eta > 300; an open channel with eta > 100 is closed (P = 0, S = 2 S(0) - S(-e)) | NaN and underflow near threshold; ~1e-9 noise in S became 0.03-0.5 % in dS/dE (Brune transform) | 22c3d03: 13N, 13N_capture_ay, 7Li_p_ay, beam_profile_kernel, data_file_comments, hybrid_potential, polarization_product, target_effect_ranges, thick_target_resonance (all < 1e-3) |
| 12fe484 | ChannelFunc: one helper for S, P and dS/dE at a channel energy (e < 0, e = 0, eta > 100 closed, Coulomb), used by the CNuc transforms, the boundary conditions, AZUREGrad, the adaptive grid, EPoint and ECIntegral; GSL U trusted only for eta <= 50 or z >= 0.1 | a level at or near threshold aborted the CLI or gave NaN; dS/dE was 0 for a level exactly at threshold | none; a width that cannot be converted (P = 0) becomes gamma = 0 with a WARNING |
| ec416de | dS/dE central difference with starting step \|e\|/4 within [1 eV, 1 keV] (was 1 eV) | round-off dominated: 1.6e-3 in chi2 between compilers | 7Li_p_ay, beam_profile_kernel, data_file_comments, target_effect_ranges; 15bb0a2: 13N, hybrid_potential (< 7e-5) |
| 02ae97b | adaptive grid: geometric tails beyond the resonance core; anchors quantised (width to 1.25^n, energy to width/4); only the finest of overlapping lattices kept; a pure Gaussian convolution whose `<targetInt>` line omits the multiplier uses a 2-width core (an explicit value, and target integration, keep 20) | narrow resonances folded 16-38 % high; the union of lattices tripled the sub-points | thick_target_resonance 19.0249 -> 30.8444, hybrid_potential (+1.6e-3) |
| 02ae97b, 528b672 | sub-point grid refresh: a folded point rebuilds its sub-points when a fit moves a narrow level's anchor (Brune; 528b672 formal parameters too). Not for components, beam profiles, energy-shifted segments or EC models | the grid was built once at the input parameters | fits only; no fixed-parameter pin |
| df8c2c8 | (1 + delta_12) = 2 for every cross section out of an identical entrance pair (reactions, capture, differential); none for an identical exit pair; the adjoint gradient follows; a rate run prints a note | dev applied it to elastic only: 12C(12C,a), d(d,p), d(d,g) were a factor 2 low | none (new tests/identical_entrance_reaction) |
| 76888ad | a `<targetInt>` line listing several segments gives each segment its own copy, so sigma lab -> c.m. is converted once | the shared object was converted once per listed segment | target_effect_ranges 3861.5 -> 3838.66. dev fixed the same bug with a once-only guard (f0444ab, pin 3838.73); both are kept, and the 7e-5 between the pins is the shift-function rows above. A line shared by segments of different entrance pairs still differs: thm converts each copy with its own pair's factor, dev uses the first segment's |
| 79358bb | EPoint::CalcEDependentValues clears its L_o, P and phase tables before refilling | a second call appended copies; lookups read the stale first copy | none. dev fixed the same doubling in RecalcEDependentValues (797c2ea, rows emptied in place); since the merge both functions use dev's ClearEDependentRows, so no output difference remains |
| 953486f | nuisance priors looked up by full Minuit index in chi2, gradient and LM setup | a fixed parameter before the prior's moved the prior to another parameter (also on dev) | none (new tests/nuisance_prior) |
| 178ae16 | AZURECalc::FixedMask from the minimizer: LM, GSL-LM, the band covariance and the FD gradient honour the fixed flags | --use-lm / --use-gsl-lm fitted fixed parameters | none (new tests/fixed_param_par); covariance.dat has one row per free parameter |

## B. Behaviour and format (no physics change)

- a0cd873: the CLI exits with -1 when a run fails (fill, initialize or transform abort); it always returned 0.
- d284224: a malformed `<potential>` block makes ReadConfigFile return -2; the binding refuses the file instead of running with defaults.
- 953486f: Calculate prints "Total Chi-Squared:" (the total with priors) to 12 significant digits, not 6.
- a2b8c81: "Denominator less than zero while transforming" names the channel and its saturation width.
- c465277: GslQuiet replaces the process-wide GSL handler swaps (WhitFunc, ShftFunc), which were not thread-safe. Values unchanged.
- e3f43df: segment norms and energy shifts are recognised by their name form, not a substring (no classic name contains the words).
- GUI:
  - 5832bee, ee092f9, 2887d99: level lines have 33 fields (32: THM binding energy, 33: gammaIsRWA; both 0 for classic). Dev's engine ignores trailing tokens.
  - d2d30a9, 2887d99, 3e7bd97: doubles are saved with round-trip precision (was 6 digits) and fixed-width columns keep a space; opening no longer rounds free level parameters. A GUI-saved classic project can therefore give a different chi2 than on dev (7Li_p_a: 2137.82733857 vs 2137.82673527 after open + save on dev).
  - f3210c7: `<parameterSettings>` level rows use AZURE2's numbering; legacy rows are recognised by value.
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
- c8d8db5: the GUI keeps a segment line's `sqrtshift` block on save.  **Differs
  from dev**, whose GUI drops it (the tokens never reached the model).

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

## Merging thm into dev later

1. `git fetch`; merge the newest origin/dev into thm first (as 3628a87 did),
   not the other way round, and get CI green on thm.
2. Re-read sections A-D: each row is a deliberate difference that dev takes
   with the merge.  Re-pin dev's tests/*/expected that A moves (the thm pins
   are already moved; a new dev test may need it).
3. Parameter layout and type codes: classic as dev; check that a newer dev
   did not add parameters after the sqrt(E) block (the cbkg block must stay
   last) or reuse type code 5.
4. Engine seams to re-check: EData.cpp/EDataThm.cpp hooks, EPoint THM
   tables and RefreshSubPointGrid, AZURECalc/AZUREGrad/AZUREAPI THM paths,
   CNuc TransformIn (gammaIsRWA, Park), ChannelFunc users, Config.h flags.
5. GUI: SegmentsDataData/SegmentsTestData fields (isTHM last before any new
   dev field), level lines 32/33, parameterSettings rows and prior_centre,
   the Runtime Options box; open + save of a classic project must stay byte
   for byte (tests/gui/thm_opt_in, parameter_settings, level_precision).
6. pyazr: `_KINDS`, `_NFIELDS`, `penalties()`, the sav readers (by name).
7. Tests and CI: tests/pyazr/CMakeLists.txt (AZURE2_BIN), tests/reference,
   tests/gui, .github/workflows/build.yml; run build/ and build-gui ctest and
   tests/run_tests.sh under the memory cap.
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

