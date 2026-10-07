# thm -> dev: differences for classic projects

For the merge of `thm` into `dev`. "Classic" means a project without a `<thm>`
block and without a THM segment (isDiff >= 10). Everything listed here changes
something for such a project. The THM machinery itself (HOES cross section,
`<thm>` keys, experiments, R(E), the DW vertex, the ps and angle windows,
cbackground=, thm_experiments.out) is new and does nothing unless it is used.

`origin/dev` (d801ba5) is the merge base, so all of dev is already in `thm`.
The dev changes merged by 18c49ea, f56c150, da72def, 497548a and e191991 are
not differences: the Wigner-limit bound, the elastic identical-pair factor
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
| 76888ad | a `<targetInt>` line listing several segments gives each segment its own copy, so sigma lab -> c.m. is converted once | the shared object was converted once per listed segment | target_effect_ranges 3861.5 -> 3838.66 |
| 79358bb | EPoint::CalcEDependentValues clears its L_o, P and phase tables before refilling | a second call appended copies; lookups read the stale first copy | none |
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

## To check at the merge

- 79358bb: whether dev ever recomputes a classic sub-point twice (energy-shifted folded segments); if not, the row moves to "no effect".
- 02ae97b: the GUI's Add Experimental Effect dialog writes the multiplier (default 20), so GUI-made projects keep the 20-width core.
- 12fe484: the narrower GSL U trust region could move S slightly for light pairs at large eta that no test covers.
