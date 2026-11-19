# Qt GUI parity checklist

This list records every user-visible capability of the Qt setup utility AZURESetup as it stands on `feature/upgrade`
at 4bf7aa06, before the Qt GUI is removed from the upgrade branch (decision 3 of HANDOFF 2026-11-07). The Phase 9 web
UI is checked against it, item by item: a box is ticked when the web UI does the same thing, or when the item is
dropped on purpose and the reason is written next to it. Each item names the source file it comes from. The Qt code
stays on `dev` and `qt5` for reference; any file can be read with `git show 4bf7aa06:gui/src/<file>` or
`git show 4bf7aa06:gui/include/<file>`. Files under `tests/gui/` are named by their basename.

The Qt program came back on 2026-11-18 on the facade `azure2::front` (`docs/upgrade/qt_gui.md`). Its checks are
the lines marked *Qt (2026-11-18): ...* under an item: what `gui_window` (`tests/gui/gui_window_test.cpp`, label
`gui`) or the screenshots of `work/qt-shots/` showed of the Qt program; a box stays as it was (it records the web
page). An item without such a line was not checked in the Qt program.

Twenty items are ticked as dropped by the PI's decision of 2026-11-13: Qt-only behaviour, the grammar of format
1 and 2 that the converter handles, or what a format-3 property replaces; the reason is written under each. The
navigation of 2026-11-19 (`web_navigation.md`) ticks two (the logo's About, File › About) and notes under the
items of sections 2, 3, 4, 5 and 10 where the page now offers them.

## 1. The .azr file

The GUI rebuilds the whole file from its tabs on every save (`AZURESetup::writeProject`). Blocks are read in the
order below; `<thm>`, `<scatteringLengths>` and `<externalCapture>` get separate passes from the top of the file,
because the engine also searches the whole file for them.

### 1.1 Blocks read and written

- [ ] `<config>`: 11 fixed lines: `true|false` (A-matrix), output directory, checks directory, then eight checks
      (compound nucleus, boundary conditions, data, Lo-matrix/penetrability, Legendre, Coulomb amplitudes, reaction
      pathways, angular distributions), each `none|screen|file`. Text after `#` on the directory lines is ignored;
      each line is written back with the GUI's own `#` comment, padded to 100 characters (`AZURESetup.cpp`)
- [ ] `<config>` optional `key=value` lines: `includeAllJPi` and `maxOrbitalMomentum` (always written);
      `abundanceSumSigma` (only when > 0); `relativistic=1` (only when on); `relativisticReferencePair` (only when on
      and > 0; no editor, carried through); `priorCovariance=<path>` (only when set) (`AZURESetup.cpp`)
- [ ] An empty output or checks directory is set to the project directory on save, with a message box
      (`AZURESetup.cpp`)
- [ ] `<potential>`: `useHybridPotential`, `useAdaptiveGrid` (not in format 3: removed 2026-11-12, every fold is on the
      shared mesh; the web UI has no control for it), then the default shape (`potentialType` 0 Woods-Saxon /
      1 Gaussian, `V0`, `R`, `a` or `r0`), then one `pair=<key>` section per pair with its own setting. The master
      switch is forced on when any pair is on. Optional on read (`NuclearPotentialTab.cpp`, `AZURESetup.cpp`)
- [ ] `<levels>`: one line per channel; the pairs have no block of their own and are rebuilt from these lines.
      Columns: 1 J, 2 parity, 3 E_x (MeV), 4 fix-energy flag, 5 `1`, 6 pair key, 7 2s, 8 2l, 9 level id, 10 active
      flag, 11 fix-width flag, 12 width (eV, or ANC in fm^-1/2, or amplitude in MeV^1/2), 13-16 light J, light
      parity, heavy J, heavy parity, 17 heavy E_x (MeV), 18-19 light and heavy atomic mass (u), 20-21 light and heavy
      Z, 22 separation energy of pair 1 (MeV), 23 separation energy of this pair (MeV), 24-26 placeholders, 27 pair
      type (0 particle, 10 gamma, 20 beta), 28 channel radius (fm), 29-30 light and heavy g factor, 31 external-capture
      multipolarity mask (bit 0 E1, bit 1 M1, bit 2 E2) (`LevelsTab.cpp`)
- [ ] `<levels>` trailing columns, told apart by count as the engine's reader (io/ LevelRows.h) does: one extra = light-particle E_x;
      two = THM binding energy B (MeV) and the width-is-amplitude flag; three = both, light E_x last. The GUI always
      writes 33 columns, and a 34th (light E_x, MeV) only on lines of a pair whose light E_x is nonzero
      (`LevelsTab.cpp`)
- [ ] A pair is identified by its key (column 6): the first line with a key defines it; pairs alike in everything
      but B stay distinct (`LevelsTab.cpp`, `PairsTab.cpp`)
- [ ] A level belongs to the compound nucleus of its channels' pairs; levels of several nuclei are written as one
      group per nucleus, each sorted by J, parity, E; the level id column is renumbered and a blank line ends each
      level (`LevelsTab.cpp`)
- [ ] The maximum l, gamma multipolarity and multipolarities per decay are inferred from the channels read
      (`LevelsTab.cpp`)
- [ ] `<abundances>` (optional, between `</levels>` and `<segmentsData>`): `[<sample>:] Z A value vary sigma` and
      `[<sample>:] sum <sigma>`; written only when several nuclei or a non-default entry exist, with a `# <nucleus>`
      comment per line (`PairsTab.cpp`)
- [ ] `<segmentsData>`: active, entrance key, exit key (-1 total capture), E_low, E_high (MeV, lab), θ_low, θ_high
      (deg), isDiff (+10 for THM), phase J and l (isDiff 2 only), norm, vary-norm, norm error (%), energy shift,
      shift error (MeV), vary-shift, data file (double-quoted when it has a space), composite flag
      `[op [-1] n {entrance exit angle [scaling]}]`, UPOS flag `[L J_f δ]`, keyword tail (`SegmentsTab.cpp`)
- [ ] Composite components read in three layouts (2, 3 or 4 tokens; `-1` marks the 4-token scaling form); the
      4-token form is written only when a component has a scaling; a ratio writes one component (`SegmentsTab.cpp`)
- [ ] `<segmentsTest>`: active, entrance, exit, E_low, E_high, E_step (MeV, lab), θ_low, θ_high, θ_step (deg),
      isDiff (+10 THM), phase J l (isDiff 2) or maximum order (isDiff 3), composite block, keyword tail
      (`SegmentsTab.cpp`)
- [ ] `<scatteringLengths>` (optional): `active pairKey s value sigma [imValue imSigma]` (fm; s = -1 coherent). Read
      wherever it sits; written after `</segmentsTest>` only when rows exist, with one header comment
      (`ScatteringLengthsModel.cpp`, `AZURESetup.cpp`)
- [ ] `<targetInt>`: active, `"segment list"`, integration points, convolution flag, σ (number or expression),
      target-integration flag, density (atoms/cm²), `"stopping equation"`, n parameters, parameters, attenuation
      flag, n, Q coefficients, energy-dependent convolution flag, `"σ(E) equation"`, n, coefficients, straggling flag,
      straggling coefficient, width multiplier, points per width (`TargetIntTab.cpp`)
- [ ] `<targetInt>` optional tail, written only when used: `sigmalab` (the `"ranges" blend tolerance` tail has no
      format-3 counterpart: deleted with decision D3, nothing to build),
      `screening U_e vary sigma`, `densityfit vary sigma [max]`, `densitygroup <label>`, `doppler T [vary sigma
      [Tmax]]`, `beamprofile n {ξ ω α w} σ_det truncation detailedBalance`, `udr "<file>" L burst channel centred`
      (`TargetIntTab.cpp`)
- [ ] `<thm>` (optional): kept as a list of raw lines and written verbatim after `</targetInt>` unless the THM
      workspace changes it; a project without it is written without it (`AZURESetupThm.cpp`)
- [ ] `<thm>` keys the GUI edits: `entranceL`, `vertex` (`real` read as `perlevel`), `kinematics`,
      `coulombIntegral`, `spectatorEnergy`, `spectatorEnergy[<pair>]`, `weight[<data line>]`,
      `weightTest[<test line>]`, `experiment[<name>] ...` lines (`ThmSettings.h`)
- [ ] `experiment[...]` keys the GUI edits: `segments=`, `background=`, `theta=`, `beam=`, `target=`, `spectator=`,
      `Ebeam=`, `lineshape=on`, `ps=`, `psNodes=`, `spectatorAngles=`, `spectatorAngleNodes=`, `distortion=`,
      `opticalAA=`, `opticalSF=`, `spectatorAngle=`, `distortionRef=`, `distortionRatio=`, `boundState=`,
      `vertexModel=`, `cbackground=`; unknown keys kept verbatim (`ThmSettings.h`)
- [ ] No form-factor block exists in the file: the THM form factor M_l(E) is computed and plotted only on the
      Diagnostics page (`ThmDiagnosticsPage.cpp`)
- [ ] `<parameterSettings>`: comment header, `# level numbering: engine`, then
      `name value lower upper error fitError nuisance category minuitIndex`, then optional
      `segment_N_norm|segment_N_energy_shift prior_centre <value>` rows (`FittingTab.cpp`)
- [x] `<parameterSettings>` read: older 12- and 14-token rows (names with spaces); an unmarked block is read by
      engine numbering or by GUI row, whichever matches more current values; rows for segments beyond the table and
      for no active level are dropped (`FittingTab.cpp`)
      *Dropped (the PI, 2026-11-13): old 12/14-token rows: the converter's job.*
- [x] `<lastRun>`: run mask (calculation mode, minimizer, formalism and option bits), `"param file"`,
      `"integrals file"`, rate entrance and exit keys, `0|1 "temperature file"`, T_min T_max T_step (GK; `-1.` when
      blank) (`AZURESetup.cpp`)
      *Dropped (the PI, 2026-11-13): <lastRun> run mask: the job request carries its settings; the page remembers the
      last one.*
- [x] `<mcmc>` (MCMC builds): `nwalkers`, `nsteps`, `chainspread`, `energyspread`, `nthreads`, `usereducedwidths`,
      `<parameters>` rows `name value priorMean priorStd useGaussian` (`MCMCTab.cpp`)
      *Dropped (the PI, 2026-11-13): <mcmc> of numcmc: numcmc goes (decision 3b); sampler settings belong to the job.*
- [x] `<externalCapture>` (legacy, read only): its multipolarity masks are moved onto the pairs (`PairsTab.cpp`)
      *Dropped (the PI, 2026-11-13): legacy block: the converter's job.*

### 1.2 Segment keyword tails

- [x] Whitelist of tail keywords carried verbatim (case-insensitive match): `sqrtshift`, `sample`, `natural`,
      `shiftgroup`, `shiftfngroup`, `normgroup`, `covariance`, `shiftfn`, `total`, `absolute`, `spin`,
      `transmission`, `yield`, `thicknessgroup`, `mixture` (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): segment keyword tails: format 1/2 only; format 3 has fields.*
- [x] Block lengths: `sqrtshift` 4 tokens, `shiftfn "<f>" N ...` 3+3N, retired `mixture` 3+2n, `transmission` and
      `yield` 2 or 4 (fit form), every other keyword 2 (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): segment keyword tails: format 1/2 only; format 3 has fields.*
- [x] Edited by dialogs: `sample`/`natural`, `covariance`, `spin`, the four group labels, `transmission`/`yield`;
      each edit rewrites its own block only and a new block goes before the `sample` block (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): segment keyword tails: format 1/2 only; format 3 has fields.*
- [x] Not edited, carried only: `sqrtshift`, `shiftfn`, `total`, `absolute`, `mixture` (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): segment keyword tails: format 1/2 only; format 3 has fields.*
- [x] A quoted path with spaces in a `covariance` block stays one token (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): segment keyword tails: format 1/2 only; format 3 has fields.*

### 1.3 Round-trip guarantees and known differences

- [ ] Open + save without edits reproduces the file; checked for segments, `<thm>`, scattering lengths, THM
      opt-in, segment renumbering and the Experimental Effects tail (`segments_roundtrip_test.cpp`,
      `thm_block_test.cpp`, `scattering_lengths_test.cpp`, `thm_opt_in_test.cpp`)
      *Qt (2026-11-18): verified: T1, T2, T3, 13N saved byte for byte, and again after every record was written back through its form (81 records)*
- [ ] Inactive lines (levels, segments, effects, scattering lengths) are kept with their flags (`*Model.cpp`)
- [x] Unknown segment tails are kept only if they start with a whitelisted keyword (`SegmentKeywords.h`)
      *Dropped (the PI, 2026-11-13): format-1/2 write-back artefacts: no rewrite of untouched lines in format 3.*
- [ ] Numbers are written with the fewest digits (6 or more) that read back to the same double; fixed-width columns
      get a separating space when full (`RoundTripNumber.h`)
      *Qt (2026-11-18): verified for <levels>: a value not edited keeps its token as read; a new value is written in the shortest form that reads back (LevelsFile.cpp); one width edited changes one row*
- [ ] Comments: kept in `<thm>` (inline comments of unchanged keys too); dropped in `<scatteringLengths>`,
      `<abundances>` and `<parameterSettings>`, which write their own header comments (`AZURESetupThm.cpp`,
      `ThmSettings.h`)
- [x] Known write-back difference: older files hold `1 0 0` in `<levels>` columns 24-26; the GUI writes `0 0 0.0`
      (`gui/CLAUDE.md`, `LevelsTab.cpp`)
      *Dropped (the PI, 2026-11-13): format-1/2 write-back artefacts: no rewrite of untouched lines in format 3.*
- [x] Known write-back difference: levels re-sorted by nucleus, J, parity, E and the level id renumbered
      (`LevelsTab.cpp`)
      *Dropped (the PI, 2026-11-13): format-1/2 write-back artefacts: no rewrite of untouched lines in format 3.*
- [ ] Known write-back difference: pyazr-written files format numbers differently. `gui/CLAUDE.md` still says the
      GUI saves 6 significant digits; that predates `roundTripNumber` (see section 15) (`gui/CLAUDE.md`)
- [ ] A saved width conversion snapshot (`GuiStateSnapshot`) writes the current state to a temporary `.azr` without
      changing the open file name, title or working directory (`AZURESetup.cpp`)

## 2. Main window and File menu

- [ ] Window title `AZURE2 -- <path>`, `-- untitled`, `-- Running`, `-- Running MCMC` (`AZURESetup.cpp`)
- [ ] Tab order: Particle Pairs, Levels and Channels, Segments, Experimental Effects, Nuclear Potential (only with
      the hybrid method on), Fitting Settings, Calculate, MCMC (MCMC builds), Plot (Qwt builds) (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): the first six tabs in this order and with these names (Fitting
      Settings and Calculate are one tab, Calculate and Fit; MCMC is its sampler; no Nuclear Potential tab), then Rate
      and ENDF, Physics and THM, Covariances and Project.*
      *Qt (2026-11-18): verified (screenshots): Particle Pairs, Levels and Channels, Segments, Experimental Effects, Nuclear Potential (always shown), THM (new), Fitting Settings, Calculate; MCMC is a command of Calculate; no Plot tab (no Qwt)*
- [x] Logo button at the right of the tab bar opens About (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): the logo at the left of the top bar opens About.*
- [x] File > About AZURE2...: version `v<version> (<git describe>)` (`AboutAZURE2Dialog.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): File › About AZURE2… and Help › About AZURE2…: the version of the
      program and of pyazr, the reference, the documentation.*

- [ ] File > New Project (Ctrl+N): resets every tab, A-matrix, THM off (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): File › New project from this one… (a new file with the pairs and
      levels, no data; POST /api/project/new); an empty project is `pyazr ui new.azr`.*
- [ ] File > Open... (Ctrl+O); a read error resets and shows "Can't Access File" (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): File › Open… (a .azr of the folder, POST /api/project/open); a
      refusal is shown in the dialog. No Ctrl+O (the browser's).*

      *Qt (2026-11-18): verified: an old file is refused with the converter's command (and offered to the converter in the dialog path, not tested)*
- [ ] File > Open Recent: five most recent files, stored in QSettings, plus Clear (`AZURESetup.cpp`)
- [ ] File > Save (Ctrl+S; Save As when untitled), Save As... (Shift+Ctrl+S), Quit (Ctrl+Q) (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): File › Save (Ctrl+S) and Save as… (a file of the folder); Quit is
      closing the browser tab and Ctrl-C on the server.*
- [ ] Opening sets the working directory to the project's directory (`AZURESetup.cpp`)
- [ ] First non-`--` command-line argument is opened as the project (`main.cpp`)
- [x] After open: a non-modal warning when the last run used GSL Coulomb functions, with "Use the default
      functions" (clears the bit) and "Keep GSL" (`AZURESetup.cpp`)
      *Dropped (the PI, 2026-11-13): the GSL bit lived in <lastRun>; use_gsl_coul is a session option shown with the
      others.*
- [ ] Documentation > Show Documentation For Current Tab (Ctrl+D): HTML help by tab index 0-5 (Pairs, Levels,
      Segments, Effects, Calculate, Plot). The indices predate the Nuclear Potential and Fitting tabs, so on
      Fitting it shows the Calculate text and on Calculate the Plot text (`AZURESetup.cpp`, `InTabDocs.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): Help › Documentation for this tab (F1): the user guide's section of
      the open tab (`web_interface.rst`, "The tabs").*
- [ ] Documentation > Open AZURE Website... opens https://azure.nd.edu (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): Help › Documentation opens the documentation site.*
- [x] Status bar message after a segment renumbering removes effect lines or THM experiments (`AZURESetup.cpp`)
      *Dropped (the PI, 2026-11-13): no renumbering: records are named.*

## 3. Configure menu

- [ ] Formalism > A-Matrix / R-Matrix (exclusive); R-Matrix refused while "Ignore external width" is on
      (`AZURESetup.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): the Configure menu's Run settings… opens the run settings
      (formalism, checks, directories, runtime options) of the Project tab.*
- [ ] Checks...: eight drop-downs None / Screen / File: Compound Nucleus, Boundary Conditions, Data, L-Matrix/Phases/
      Penetrabilities, Legendre Polynomials, Coulomb Amplitudes, Reaction Pathways, Angular Distributions
      (`EditChecksDialog.cpp`)
- [ ] Directories...: Output Directory and Checks Directory with Choose... (`EditDirsDialog.cpp`)
- [ ] Runtime Options...: Use GSL Coulomb functions; Use Brune formalism; Use Park parametrization (implies Brune);
      Ignore external width if internal width is zeroed (switches to A-matrix with a message); Use RMC capture
      formalism (neutron capture only); Do not perform parameter transformations; Use Hybrid Coulomb method (shows
      the Nuclear Potential tab); Use adaptive integration grid for target effects; Use Trojan Horse Method
      (`EditOptionsDialog.cpp`, `AZURESetup.cpp`)
- [ ] Runtime Options > Kinematics box: Relativistic kinematics, noted as saved in the project (`EditOptionsDialog.cpp`)
- [ ] Runtime Options > Prior Parameter Covariance: file field + Choose..., relative to the project directory
      (`EditOptionsDialog.cpp`)
- [ ] Switching THM off on a project with THM content warns, lists the content, and only hides the controls
      (`AZURESetupThm.cpp`)
- [ ] THM opt-in is turned on at open when the project has a `<thm>` block, a THM segment, a binding energy or an
      amplitude width (`AZURESetupThm.cpp`)
- [ ] THM Workspace... (visible only with THM on); refused with the engine's message when `<thm>` has a bad line
      (`AZURESetupThm.cpp`)
- [ ] Build Coulomb Table...: asks first (about 5 min, 88 MB), builds on all threads with a progress dialog and
      Cancel, reports the path and time; later runs use the table (`AZURESetup.cpp`)

## 4. Particle Pairs tab

- [ ] Table of pairs numbered by key; shown columns: Light Spin, Light Particle, Heavy Spin, Heavy Particle, Heavy
      Excitation, Separation Energy, Channel Radius, THM Binding; Compound Nucleus only with several nuclei. The
      hidden column 15 is Light Excitation since the dev merge (comment says "THM binding energy"), so B shows and
      the light E_x does not: probably a regression (`PairsTab.cpp`, `PairsModel.h`)
- [ ] `+` opens Add a Particle Pair; `-` deletes the selected pair; double-click edits (`PairsTab.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): Add pair and Remove pair at the table, the selected pair's numbers
      beside it.*
- [ ] Pair dialog fields: Particle Pair Type (Particle, Particle / Particle, Gamma / Beta Decay); light and heavy J
      (validated), parity (- / +), Z, M (atomic mass, u); Heavy Excitation [MeV]; Light Excitation [MeV] (particle
      pairs only); Separation Energy [MeV]; Channel Radius [fm] (`AddPairDialog.cpp`)
- [ ] Gamma pair: light fields fixed (J 1, +, Z 0, M 0), S 0, radius 0; External Capture Multipolarities E1 and E2
      shown (M1 hidden) (`AddPairDialog.cpp`)
- [ ] Beta pair: light J 1/2, +, M 0.0005; heavy excitation fixed at 0 (`AddPairDialog.cpp`)
- [ ] Isotope fields ("p", "a", "7Li", "Li7") fill Z and the AME2020 atomic mass when edited (`AddPairDialog.cpp`)
- [ ] "Fill masses and S from AME2020": both masses and S = m1 + m2 - M; heavy mass only for gamma/beta
      (`AddPairDialog.cpp`)
- [ ] Live mass check line: table S, "consistent" within 1 keV, orange warning beyond; mass offset flagged above
      1e-3 u (`AddPairDialog.cpp`, `MassTable.h`)
- [ ] Refusals: first pair must be particle-particle; duplicate pair; deleting pair 1 when pair 2 is a gamma pair
      (`PairsTab.cpp`)
- [ ] Separation-energy cell: tooltip against AME2020, orange text when off by more than 1 keV (`PairsModel.cpp`)
- [ ] "Masses from AME2020..." dialog: per pair Fill tick, light and heavy mass (old → table), S, table S,
      S - table (keV); option to also set S; "Fill ticked pairs" with a count (`MassCheckDialog.cpp`)
- [ ] Pair labels: `²⁷Al(p,γ)` style with element symbols, `*(E)` mark for an excited light particle,
      `[E_x MeV]` (`PairsModel.cpp`)
- [ ] Editing or adding a pair regenerates every level's channels, keeping widths and fix flags of channels that
      survive; deleting a pair removes its channels and shifts higher keys (`LevelsTab.cpp`)
- [x] Deleting a pair moves per-pair nuclear-potential settings down with the keys (`NuclearPotentialTab.cpp`)
      *Dropped (the PI, 2026-11-13): potential records name their pair; no keys to shift.*
- [ ] Isotopic abundances box (only with more than one compound nucleus): Compound Nucleus (target pair),
      Abundance, Fit?, Prior Sigma; parameter `abundance[_<sample>]_<Z>_<A>` (`PairsTab.cpp`, `AbundancesModel.cpp`)
- [ ] Samples: combo, "+ sample" (one-word name, no `:` or `#`), "- sample" (confirm), per-sample sum-to-one σ
      ("default" = use the default) (`PairsTab.cpp`)
- [ ] Default sum-to-one prior σ spin box (0 = off), 4 decimals (`PairsTab.cpp`)

## 5. Levels and Channels tab

- [ ] Levels table (sortable except the first two columns): Include?, Fix?, Level Spin (J^π with ⁺/⁻), Energy
      [MeV]; Compound Nucleus column and a nucleus filter only with several nuclei (`LevelsTab.cpp`,
      `LevelsModel.cpp`, `LevelsHeaderView.h`)
      *Qt (2026-11-18): shown in the screenshot of T2*
- [ ] `+` Add a Level, `-` delete, double-click Edit a Level (`LevelsTab.cpp`)
      *Web page (2026-11-19, `web_navigation.md`): Add level and Remove level at the table (Remove level new: the page
      had no control), the selected level and its channels beside it.*
- [ ] Level dialog: Compound Nucleus (several nuclei; fixed on edit), Excitation Energy [MeV], Spin-Parity
      (J^π) list of only the J^π that give channels; Accept disabled until the energy is a number
      (`AddLevelDialog.cpp`)
- [ ] With no pairs the J^π list falls back to 0⁺ ... 10⁻ with a note; an edited level keeps its own J^π marked
      "(no channels)" (`AddLevelDialog.cpp`)
- [ ] Duplicate level refused (`LevelsTab.cpp`)
- [ ] Channel Configuration: Maximum Orbital Momentum (0-10, default 2), Maximum Gamma Multipolarity (1-10),
      Maximum Gamma Multipolarities Per Decay (1-10) (`LevelsTab.cpp`)
      *Qt (2026-11-18): shown in the screenshot of T2*
- [ ] "Hard-sphere scattering in every J^π (levels or not)" (default on, saved as `includeAllJPi`) (`LevelsTab.cpp`)
- [ ] Missing-J^π line: all covered / hard-sphere only list / red "not in the calculation" list (`LevelsTab.cpp`)
      *Qt (2026-11-18): shown in the screenshot of T2 (the red list)*
- [ ] Channel generation rules: particle (l, s) with parity and l ≤ l_max, identical particles l + s even; beta
      Fermi and Gamow-Teller; gamma E/M multipoles with integer/half-integer matching (`LevelsTab.cpp`)
- [ ] Channels In Selected Level table: Fix?, Channel Pair, s, l (`LevelsTab.cpp`, `ChannelsModel.cpp`)
- [ ] Channel Details text: level, pair key, s and l or radiation type (E1, M1, ...), internal or internal +
      external capture, Fermi/Gamow-Teller, particle spins, Z, M, g, excitation, S, radius (`LevelsTab.cpp`)
- [ ] Width field label and unit by channel kind: Partial Width (eV); ANC (fm^-1/2) below threshold; Mu (nm) and
      Q (b) for M1 / E2 gamma to the pair's own level; B for beta; Width (MeV^1/2) for an amplitude
      (`ChannelDetails.cpp`)
- [ ] Wigner Limit button: Γ_W = 2 P γ²_W in eV, or √γ²_W in MeV^1/2 for amplitudes; "N/A" for capture, beta and
      bound channels (`LevelsTab.cpp`)
      *Qt (2026-11-18): now physics::single_particle_width through the facade; not run in a test*
- [ ] Fix All Widths, Fix All Energies (`LevelsTab.cpp`)
- [ ] NDS button: picks the nucleus (several nuclei), opens Import Levels from IAEA Database; reports added and
      skipped duplicates (`LevelsTab.cpp`)
- [ ] IAEA dialog: Nuclide field, Query IAEA Database (network, `nds.iaea.org/relnsd/v1/data?fields=levels`), table
      Select / Energy (MeV) / J^π (IAEA) / J choice / parity choice, Add Selected Levels (`IAEALevelDialog.cpp`,
      `IAEALevelData.cpp`)
- [ ] Export LaTeX: table of levels (sorted by E) × channel types with widths, amplitude widths marked; Copy to
      Clipboard or Save to File (.tex) (`LevelsTab.cpp`)

## 6. Segments tab

### 6.1 Segments From Data

- [ ] Table: active check, Reaction (pair labels, composite sum/ratio, total capture, sample mark, UNDEFINED in red),
      Energy Range, Angle Range, Data Type, Data File, Data Norm. (red when varied, grey when tied), Vary Norm.?,
      Energy Shift [MeV], Groups (`SegmentsTab.cpp`, `SegmentsDataModel.cpp`)
- [ ] Tooltip on a tied norm or shift: "Not used: segment k takes the ... of segment m" (`SegmentsDataModel.cpp`)
- [ ] Measured scattering lengths appear as rows after the segments in the same table (`SegmentsTab.cpp`)
      *Qt (2026-11-18): model data are the Model Data sub-tab of Segments (records of <modelData>)*
- [ ] Buttons: `+`, `-`, up, down, Check All, Uncheck All, EXFOR; Filter by entrance and exit pair
      (`SegmentsTab.cpp`)
- [x] Moving or deleting a line renumbers what names segments: effect segment lists (line dropped when empty), THM
      `segments=` and `weight[k]`, `segment_N_*` settings and prior centres (`SegmentsTab.cpp`, `AZURESetup.cpp`)
      *Dropped (the PI, 2026-11-13): no renumbering: records are named.*
- [ ] Add/Edit a Segment From Data: Entrance Pair and Exit Pair drop-downs (by name; exit lists reactions and
      "Total capture (all γ pairs)") (`AddSegDataDialog.cpp`, `PairKeyCombo.cpp`)
- [ ] Lab Energy [MeV] low/high; Lab Angle [degrees] low/high, retitled C.M. Angle for C.M. data and A_y
      (`AddSegDataDialog.cpp`)
- [ ] Data Type: Angle Integrated, Differential (+ Frame Lab / C.M.), Phase Shift (J, l, s), Analyzing Power,
      Angular Coefficients, Scattering Length; isDiff codes 0-7, 9; unknown codes shown as "Observable n"
      (`AddSegDataDialog.cpp`)
- [ ] Multipolarity All / E1 only / E2 only for angle-integrated capture (codes 0/5/6) (`AddSegDataDialog.cpp`)
- [ ] Total capture allows only angle-integrated data (`AddSegDataDialog.cpp`)
- [ ] THM (modified R-matrix) tick, THM projects only, disabled for phase shifts (`AddSegDataDialog.cpp`)
- [ ] Data Norm., Vary Norm? (disabled for A_y), Norm. Error [%] (shown when varied) (`AddSegDataDialog.cpp`)
- [ ] Energy Shift [MeV], Energy Shift Error [MeV], Vary Energy Shift? (`AddSegDataDialog.cpp`)
- [ ] Data File with Choose...; Covariance File with Choose... (`AddSegDataDialog.cpp`)
- [ ] Groups: Normalization, Energy Shift, Shift Function drop-downs listing the project's labels; no blank or
      quote; energy-shift and shift-function groups exclusive (`AddSegDataDialog.cpp`)
- [ ] Sample (abundance weights): none / default sample / named samples / "(not in the Pairs tab)"; disabled for
      ratios, phase shifts, A_y (`AddSegDataDialog.cpp`)
- [ ] Sample row (neutral entrance, angle-integrated): None / Transmission / Self-shielded yield, n [atoms/b] > 0,
      Fit, Prior sigma, thickness Group; invalid block refused with a warning (`SampleRow.h`, `SegmentsTab.cpp`)
- [ ] UPOS: Secondary Decay L, Final State J (Ic), Mixing Ratio (δ) (`AddSegDataDialog.cpp`)
- [ ] Advanced Settings: Sum or Ratio, components (entrance, exit) with Add/Remove, Scaling (sum), fixed denominator
      Angle [degrees] (ratio of differential data) (`AddSegDataDialog.cpp`)
- [ ] Scattering Length mode: entrance pair only, Channel Spin (allowed s with parity, or coherent), Re a0 [fm],
      uncertainty, Im a0 [fm], uncertainty (`AddSegDataDialog.cpp`)
- [ ] Scattering-length rows flagged red with a tooltip when the pair is missing, not a particle pair, or the
      spin is not allowed; kept and written as they are (`ScatteringLengthsModel.cpp`)
- [ ] Duplicate line refused (`SegmentsTab.cpp`)

### 6.2 Segments Without Data

- [ ] Table: active, Reaction, Energy Range, Energy Step, Angle Range, Angle Step, Data Type (`SegmentsTestModel.cpp`)
- [ ] Buttons `+`, `-`, up, down; double-click edits (`SegmentsTab.cpp`)
- [ ] Dialog: pairs, Lab Energy [MeV] low/high/step, Lab Angle [degrees] low/high/step; Data Type Angle Integrated,
      Differential (Frame Lab / C.M.), Phase Shift (J, l, s), Angular Distribution Coefficients (Maximum Order),
      Analyzing Power; THM tick; Sample (weights); Sample row (kind and n, no fit); Advanced Mode (Sum/Ratio)
      (`AddSegTestDialog.cpp`)

### 6.3 EXFOR import

- [ ] Retrieve Experimental Data from EXFOR: entrance/exit pair, Target, Reaction (filled from the pairs, `EL` for
      elastic), Quantity (`SIG`, `DA`, editable), Search EXFOR (network, `nds.iaea.org/exfor/x4list`)
      (`ExforDialog.cpp`, `ExforData.cpp`)
- [ ] Results table: Author / Year, Reaction, Points, E range (MeV, lab), DatasetID, Reference; Show / Preview
      Dataset (`x4get`, CSV) with raw CSV and AZURE2-format panes (`ExforDialog.cpp`)
- [ ] Conversions: energy to MeV, cross section to b or b/sr, S factor to cross section, ratio to Rutherford to
      cross section in the angle column's frame; c.m. energy used when no lab energy (`ExforData.cpp`)
- [ ] Save to data folder... as `<DatasetID>_<reaction>.dat` (E, angle, σ, error) (`ExforDialog.cpp`)

## 7. Experimental Effects tab

- [ ] Table of effect lines with `+`, `-`, double-click edit; columns include Segment List, Int. Points, flags and
      parameters of each effect, Effects summary (`TargetIntTab.cpp`, `TargetIntModel.cpp`)
      *Qt (2026-11-18): targets, beams and detectors are records listed and edited by their form (screenshot of T2); a record broken in its form is named by the engine in the diagnostics panel (verified on a data set: E_UNKNOWN_UNIT)*
- [x] Segments List ("5, 8, 13", "4-12", "3,6,7-14"); Number of Integration Points (`AddTargetIntDialog.cpp`).
      Not to be built: Apply in Energy Ranges, Blend Width, Auto Tolerance (deleted with decision D3: every point of
      a folded data set is folded)
      *Dropped (the PI, 2026-11-13): no renumbering: records are named.*
- [ ] Adaptive grid: Resonance Width Multiplier, Points Per Width (`AddTargetIntDialog.cpp`)
- [ ] Gaussian convolution: σ [MeV, lab] as a number or expression of lab energy x; parameters of σ(E) table;
      lab or legacy c.m. reading; Convert to lab form (`AddTargetIntDialog.cpp`)
- [ ] Target integration: Active Density [atoms/cm²]; Effective Stopping Cross Section [MeV cm²/atoms] equation
      and parameter table; Target Element or Compound Formula + Active Element; Fetch from ERYA (built-in SRIM
      tables; message when not built) (`AddTargetIntDialog.cpp`)
- [ ] Calculate ΔE: Energy [keV] → ΔE [keV] from density and stopping parameters (`AddTargetIntDialog.cpp`)
- [ ] Energy straggling with Coefficient (0 = Bohr) (`AddTargetIntDialog.cpp`)
- [ ] Fitted Target Density: Density Is a Fit Parameter, Fit Density, Prior Sigma [atoms/cm² or %], Grid Maximum,
      Density Group; refused for density 0 or with beam profile / UDR (`AddTargetIntDialog.cpp`)
- [ ] Attenuation Coefficients: number and values (`AddTargetIntDialog.cpp`)
- [ ] Beam Profile: skewed-Gaussian components (Location ξ, Scale ω [MeV], Skewness α, Weight), Detector
      Resolution Sigma [MeV], Truncation [s.d.], detailed-balance weight; window from data columns 5-6
      (`AddTargetIntDialog.cpp`)
- [ ] SAMMY UDR: file with Browse..., Flight Path [m], Burst FWHM [ns], Channel Width [ns] (negative = per decade),
      re-centre (`AddTargetIntDialog.cpp`)
- [ ] Electron Screening: U_e [eV], Adiabatic Limit button (from electron binding energies), Fit U_e, Prior
      Sigma [eV] (`AddTargetIntDialog.cpp`)
- [ ] Free-Gas Doppler Broadening: Effective Temperature T [K] > 0, Temperature Is a Fit Parameter, Fit T, Prior
      Sigma [K], Grid Maximum [K] (> T) (`AddTargetIntDialog.cpp`)

## 8. Nuclear Potential tab (hybrid method on)

- [ ] Applies to: Default (every pair without its own) or Pair k with state on / on, inherited / off; list follows
      the pairs as they load (`NuclearPotentialTab.cpp`)
- [ ] Use hybrid potential here; Follow Default (drop a pair's own setting) (`NuclearPotentialTab.cpp`)
- [ ] Potential Type Woods-Saxon (V₀ MeV, R fm, a fm) or Gaussian (V₀ MeV, r₀ fm); numbers > 0 required
      (`NuclearPotentialTab.cpp`)
- [ ] Apply Settings, Reset to Default; summary line of which pairs use the potential (`NuclearPotentialTab.cpp`)

## 9. Fitting Settings tab

- [ ] Sub-tabs: Level Parameters, Abundances, Normalization, Energy Shifts, Screening, Target Densities, Doppler
      Temperatures, Sample Thicknesses, THM Coherent Background (THM on) (`FittingTab.cpp`, `FittingTabThm.cpp`)
- [ ] Columns: Parameter, Value, Lower Limit, Upper Limit, Error, Fit Error (read-only), Use as Nuisance; Prior
      Centre on Normalization and Energy Shifts (`FittingTab.cpp`)
- [ ] Only free parameters are listed; names `Level N Energy (MeV)`, `Level N Channel M Width (eV|MeV^(1/2))`,
      `segment_N_norm`, `segment_N_energy_shift`, abundance, `screening_n`, `density_n`, `doppler_n`,
      `segment_<key>_thickness` (`FittingTab.cpp`)
- [ ] Level numbering follows the engine (active levels, J-group order) (`FittingTab.cpp`, `LevelsTab.cpp`)
- [ ] Group members shown greyed, not editable, "tied to segment k (<keyword> <label>)" (`FittingTab.cpp`)
- [ ] Edits go back to the source tab: value, error and vary flag of norms, shifts, abundances, screening,
      density (whole group), Doppler T, thickness, coherent background; non-numbers refused (`FittingTab.cpp`)
- [ ] Refresh from Current; Clear Limits (all limits to 0, count reported) (`FittingTab.cpp`)
- [ ] Load from .sav file: backs up the `.azr` as `.azr.bk`, converts amplitudes to physical values against the
      current GUI state, NaN set to 0 with a warning, verification report (`FittingTab.cpp`, `AZURESetup.cpp`)

## 10. Calculate tab

- [ ] Calculation Type: Calculate Segments From Data, Fit Segments From Data, Calculate Segments Without Data,
      Perform MINOS Error Analysis, Calculate Reaction Rate (`RunTab.cpp`)
      *Qt (2026-11-18): verified: the engine's commands under the old titles; calculate and fit run from the tab*
- [ ] Minimizer (fit only): Minuit2, Minuit2 (analytic grad.), Levenberg-Marquardt, GSL Trust-Region, NLopt SBPLX,
      COBYLA, BOBYQA, NEWUOA, PRAXIS, Nelder-Mead (NLopt builds); MINOS forces Minuit2 (`RunTab.cpp`)
      *Qt (2026-11-18): the engine's fitters (front::Fitters) in the Fitting Settings tab; a call limit verified*
- [ ] Chi-Sq Variance (MINOS only, default 1.0) (`RunTab.cpp`)
- [ ] Uncertainty (band) for fit, extrapolation and MINOS; Scale covariance (fit/MINOS with band); Set Wigner
      limits (`RunTab.cpp`)
- [ ] Parameters File: Create New / Use: <file> with Choose... (`RunTab.cpp`)
- [x] External Capture Integrals File: Create New / Use: <file> (hidden for rates) (`RunTab.cpp`)
      *Dropped (the PI, 2026-11-13): the session caches intEC.dat by its signature.*
- [ ] Reaction Rate Parameters: Entrance Pair, Exit Pair (drop-downs), Minimum / Maximum Temperature (GK),
      Temperature Step (GK), or Use Temperature File (`RunTab.cpp`)
- [ ] Save and Run (Alt+R) saves first, then validates: no active segments; total capture without capture pair;
      undefined pair key; elastic rate; missing output/checks directory; missing parameter, integrals or temperature
      file; angular-distribution coefficients with convolution or target integration (`AZURESetup.cpp`)
- [ ] Stop AZURE2 sets the engine's stop flag; Run disabled while running; log scrolled to the end (`AZURESetup.cpp`,
      `AZUREMainThread.h`)
      *Web page (2026-11-19, `web_navigation.md`): Stop in the job's panel and in the status line on every tab.*

      *Qt (2026-11-18): verified: a fit stopped from the window ends at its best point, fitStatus.out says stopped; run controls disabled while running*
- [ ] Runtime log pane, `\r` lines overwrite the last line, mouse filtered while running (`FilteredTextEdit.h`)
- [ ] Files the run writes are the engine's in the output directory (`param.sav`, `parameters.out`,
      `normalizations.out`, `param.errors`, `covariance_matrix.out`, `reactionrates.out`, `ECint.dat`, `AZUREOut_*`)
      (`InTabDocs.cpp`)
      *Qt (2026-11-18): verified: the files of calculate from the window equal AZURE2 calculate's, 9 of 9 (T1)*

## 11. MCMC tab (MCMC builds)

- [ ] Sub-tabs Parameters, Sampling Settings, Results; four "?" help buttons (`MCMCTab.cpp`)
- [ ] Parameters table: Parameter, Current Value, Prior Mean, Prior Std, Use Gaussian Prior, Category
      (`MCMCTab.cpp`)
- [ ] Load Physical Parameters; Load RWA Parameters (reads `output/param.par` by name) (`MCMCTab.cpp`)
- [ ] Automatic priors for norms, energy shifts, sqrt(E) shifts, abundances, with prior centres honoured
      (`MCMCTab.cpp`)
- [ ] Number of Walkers 2-1000 (50), Number of Steps 100-10⁶ (10000), Initial Chain Spread 0.001-50 % (5 %),
      Level Energy Spread 0.001-1 keV, Number of Threads 1-N (min(4, N)) (`MCMCTab.cpp`)
- [ ] Run MCMC Sampling (saves, validates as Calculate), Stop Sampling, Fresh Start (deletes
      `output/samples.mcmc`), Use Reduced Widths (`MCMCTab.cpp`, `AZURESetup.cpp`, `AZUREMCMCThread.cpp`)
- [ ] Progress bar, status, current iteration, log probability, log likelihood, log prior; Sampling Log
      (`MCMCTab.cpp`)
- [ ] Results: Statistics (Mean, Std Dev, 2.5 %, 50 %, 97.5 %), Chain Data, Refresh Results from File, warning
      row when stopped early (`MCMCTab.cpp`)

## 12. Plot tab (Qwt builds)

- [ ] Segment lists From Data and Without Data, each filtered by In and Out pair; multi-select; Draw, Clear
      (`PlotTab.cpp`)
- [x] Reads `output/AZUREOut_aa=<in>_R=<out>.out` or `_TOTAL_CAPTURE.out` (data), `.extrap` (no data), and the
      `.band` sibling; block within the file by segment order (`PlotTab.cpp`, `AZUREPlot.cpp`)
      *Dropped (the PI, 2026-11-13): no AZUREOut files read: arrays from the API.*
- [ ] X axis: CoM Energy, Excitation Energy, CoM Angle; log scale (`PlotTab.cpp`)
- [ ] Y axis: Cross Section or S-Factor; log scale; title per quantity: Cross Section [b] / [b/sr], S-Factor
      [MeV b] / [MeV b/sr], Phase Shift [degrees], Analyzing Power, Angular Distribution Coefficient
      (`AZUREPlot.cpp`)
- [ ] Analyzing power in the selection forces a linear axis and cross section (`PlotTab.cpp`)
- [ ] Display: Grid, Legend, Levels (markers at E_x or E_c.m. per entrance pair, not on angle axis), Uncertainty
      (1σ band) (`PlotTab.cpp`, `AZUREPlot.cpp`)
- [ ] Curve list with Customize... / double-click: label, data and calculation colour, marker (11 shapes), marker
      size, line width (`PlotTab.cpp`)
- [ ] THM segments labelled "(THM)"; label from the data file name (`PlotTab.cpp`)
- [ ] Zoom by rubber band, right-click zoom out, middle-button pan, coordinate tracker (`AZUREPlot.cpp`)
- [ ] Export... (PDF, SVG, PS, image formats); Print... (landscape) (`AZUREPlot.cpp`)

## 13. THM Workspace (Configure > THM Workspace...)

- [ ] Dialog with pages Model, Experiments, Channels, Diagnostics (engine-API builds); Help opens the local or online
      user guide; Accept validates every page and applies; pages disabled when no THM segment exists
      (`ThmWorkspace.cpp`)
- [ ] Model: Entrance partial waves (incoherent / coherent), Vertex boundary (constant / perlevel / onshell),
      Kinematic factors (lacognata / triple / kf3body / lambda32), External Coulomb term, Spectator energy (MeV)
      (`ThmModelPage.cpp`)
- [ ] Model: Spectator energy per pair table; Energy-dependent weight per segment table (Data/Test, Segment, weight
      file E_cm MeV, w) with file chooser; duplicates and bad weight files refused (`ThmModelPage.cpp`,
      `ThmSettings.h`)
- [ ] Experiments: list (Name, Segments, Background, Three-body reaction), Add, Remove; Name letters, digits,
      `_-.+` (`ThmExperimentsPage.cpp`)
- [ ] Experiments: segment picker of THM data segments with a free norm in no other experiment; Background none /
      const / linear / quadratic; Exit angle all or c.m. window (`ThmExperimentsPage.cpp`)
- [ ] Three-body reaction: Beam, Target, Spectator (AME2020 nuclide or Z,A,mass), Beam energy (lab, MeV); derived
      B(x+s), E_qf, ζ at the data ends; Coulomb line shape (needs Brune) (`ThmExperimentsPage.cpp`)
- [ ] Spectator acceptance: point / Hulthén (custom a, b in fm⁻¹) / Gaussian (FWHM, MeV/c) / table; p_min, p_max
      (MeV/c); nodes; ⟨T_s⟩; Directions one / lab or c.m. window / lab or c.m. table with nodes
      (`ThmExperimentsPage.cpp`)
- [ ] Distortion: None / Coulomb / Optical / Table; Vertex plane or distorted wave; spectator angle quasi-free /
      lab / c.m.; E_ref (auto); bound state Whittaker / Yukawa with r_min (fm); optical a+A and s+F plane / Coulomb /
      global / Woods-Saxon; ratio DWBA/PWBA or DWBA; R(E) at the data ends (`ThmExperimentsPage.cpp`)
- [ ] Coherent background table: J^π, Exit, Channels (all / one), Form (const / linear), Re/Im c0, c1 with fixed
      ticks; parameters appear in the Fitting tab (`ThmExperimentsPage.cpp`, `FittingTabThm.cpp`)
- [ ] Channels: binding energy B (MeV) per pair with its THM segments; width-is-amplitude flag per channel (Level,
      Pair, l, s, Width, Amplitude); written to `<levels>` columns 32-33 (`ThmChannelsPage.cpp`)
- [ ] Diagnostics: Segment, Compute (temporary copy, off the GUI thread, busy bar), J^π group; plots of entrance
      vertex |M_l|² with nodes, HOES vs on-shell, line shape |N_C|² and ζ(E), weight w(E), spectator window,
      distortion R(E), angular distribution at E (`ThmDiagnosticsPage.cpp`, `ThmDiagnostics.cpp`)

## 14. Cross-cutting

- [ ] Runs on a worker thread with its own copy of the Config (every field the GUI sets, including hard-sphere,
      abundance prior, relativistic, prior covariance); output streamed to the log (`AZUREMainThread.h`)
      *Qt (2026-11-18): verified: front::Run in a QThread; the window's event loop turned during the run; the run reads the saved file*
- [x] File > Save before every run; runs use the saved file (`AZURESetup.cpp`)
      *Dropped (the PI, 2026-11-13): jobs read the working text; the user's file is written only by Save.*
- [ ] Pair keys are never typed: drop-downs list "k: label [E MeV]"; a key beyond the list kept as
      "(no such pair)" (`PairKeyCombo.cpp`)
- [ ] Azure theme (Fusion style, palette, watermark); `AZURE2_THEME=system` turns it off (`AzureTheme.cpp`)
- [ ] File dialogs: Qt's own on macOS, native elsewhere; `AZURE2_NATIVE_DIALOGS=1|0` forces (`FileDialogPolicy.cpp`)
- [ ] `AZURE2_SCREENSHOT_DIR=<dir>`: open the project, save each tab as `tab<n>_<name>.png` after 1.5 s, quit
      (`main.cpp`, `AZURESetup.cpp`)
      *Qt (2026-11-18): the grab of every tab verified (saveTabScreenshots, AZURE2_QT_SHOTS in the test); the variable itself not run*
- [ ] High-DPI scaling; macOS font substitutions (`main.cpp`)
- [ ] Network fetches: IAEA levels (`nds.iaea.org/relnsd`), EXFOR (`nds.iaea.org/exfor`), website and THM help
      opened in the browser (`IAEALevelData.cpp`, `ExforData.cpp`, `AZURESetup.cpp`, `ThmWorkspace.cpp`)
- [ ] Embedded AME2020 table, the same file pyazr reads (`MassTable.cpp`)
- [ ] Per-tab help text (HTML) for Pairs, Levels, Segments, Experimental Effects, Calculate, Plot; Fitting and MCMC
      have their own text (`InTabDocs.cpp`, `FittingTab.cpp`, `MCMCTab.cpp`)

## 15. Tests that pin GUI behaviour

- `density_fitting_test.cpp`: `density_<n>` rows, group sharing, % prior kept, Minuit index, values written back.
- `doppler_fitting_test.cpp`: `doppler_<n>` rows only for the fit form; edits both ways; values written back.
- `exfor_rutherford_test.cpp`: ratio-to-Rutherford conversion in the frame of the angle column (canned CSV).
- `file_dialog_policy_test.cpp`: Qt file dialogs on macOS, native elsewhere; env override.
- `gsl_coulomb_warning_test.cpp`: GSL-Coulomb warning on open, its two buttons, the bit written or kept.
- `level_dialog_test.cpp`: level dialog lists only J^π with channels; fallback list; own J^π kept on edit.
- `level_precision_test.cpp`: save keeps every `<levels>` number exact; χ² unchanged; snapshot leaves project path.
- `levels_compat_check.cpp` (not a ctest): write-back of real projects, J^π list, channel regeneration.
- `light_excitation_test.cpp`: light E_x column, `*(E)` label, identical pairs, 33/34-column layout.
- `main_thread_config_test.cpp`: the worker's Config copy carries the hard-sphere option.
- `mass_table_test.cpp`: AME2020 values, isotope fields, Fill button, S tooltip/colour, masses dialog.
- `mcmc_priors_test.cpp`: MCMC prior kinds by name (norm, shift, sqrt shift, cbkg); prior centre; param.par by name.
- `natural_target_test.cpp`: several nuclei, abundance block round trip, levels per nucleus, sample tail, abundance
  parameters, full-precision levels.
- `nuclear_potential_tab_test.cpp`: per-pair selector follows pairs read after `<potential>`.
- `pair_identity_test.cpp`: pairs identified by key; twins differing only in B or active flag kept.
- `pair_key_combo_test.cpp`: pair drop-downs show names, return keys, keep unknown keys.
- `parameter_settings_numbering_test.cpp`: `Level N` rows follow engine numbering; old blocks read either way.
- `parameter_settings_test.cpp`: `<parameterSettings>` round trip with engine level/channel numbering.
- `runtime_options_test.cpp`: relativistic option and prior covariance in Runtime Options, read and written.
- `scattering_lengths_test.cpp`: allowed spins, block round trip anywhere in the file, flagged rows, editors.
- `screening_fitting_test.cpp`: `screening_<n>` rows, edits both ways, values written back.
- `segment_groups_test.cpp`: normgroup/shiftgroup/shiftfngroup labels, dialog edits only its block, tied rows.
- `segment_observable_test.cpp`: observable + frame + multipolarity map to isDiff codes and back unchanged.
- `segments_move_test.cpp`: moving a data segment keeps every field (UPOS included).
- `segments_renumber_test.cpp`: move/add/delete renumbers THM, effects and settings; saves byte for byte.
- `segments_roundtrip_test.cpp`: keyword tails (sqrtshift, sample, groups, covariance, shiftfn) verbatim.
- `target_int_tab_test.cpp`: Experimental Effects optional tail (ranges, beam profile, etc.) read and written.
- `thickness_fitting_test.cpp`: transmission/yield blocks, thickness groups, `segment_<key>_thickness` rows.
- `thm_block_test.cpp`: `<thm>` carried verbatim; no block added to a project without one.
- `thm_opt_in_test.cpp`: THM opt-in from content, warning on switching off, byte-identical saves, Park + THM.
- `thm_options_dialog_test.cpp`: Model page and ThmSettings parse/compose; engine accepts the GUI's block.
- `thm_workspace_test.cpp`: experiments, channels, line shape, spectator window, distortion, exit angle,
  Diagnostics against the engine; untouched workspace saves byte for byte.
