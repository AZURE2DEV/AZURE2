# Skills migration: old section → new location

Phase 0 rewrote `.claude/skills/` into the ten-skill set of `RULES.md` §3 plus a trimmed `rmfit-campaign`. This
table lets a reviewer check that nothing was dropped: every section of every old skill (line ranges refer to the old
files on `feature/upgrade` before the rewrite, e.g. `git show ff75fc28:.claude/skills/azure2-eval/SKILL.md`) has a row
naming where its content now lives, or why it was dropped. "Dropped" is used only for history with no effect on
current or older builds (directory provenance, job IDs, hashes, dated narratives, results of particular campaigns).
Bug fixes that change results for users of older builds are kept as version notes
(`builds before YYYY-MM-DD (expires: 2027-02-01)`). Content that describes how to change the code went to
`CLAUDE.md`, `tests/CLAUDE.md` or `gui/CLAUDE.md`; a few behaviour facts missing from the Sphinx docs were added to
`docs/source/reference/command_line.rst`.

## Old skill → new owner

| Old skill (lines) | New owner(s) |
|---|---|
| azure2-eval (3,721) | split by item: azure2-setup, azure2-calculate, azure2-fit, azure2-check, azure2-uncertainty, azure2-plot, azure2-targets, azure2-exchange, azure2-thm; developer notes to `tests/CLAUDE.md` |
| azr-project-builder (172) | azure2-setup |
| nds-explorer (321) | azure2-setup (`references/nds-exfor.md`) |
| r-matrix-analysis (367 + references) | §1–5, 9, 11 and key papers: azure2-setup; §6–8, 10 and formalism: azure2-fit |
| azure2-legacy-params (117 + script) | azure2-fit (`references/legacy-params.md`, `scripts/convert_sav_names.py`) |
| azure2-drop-dataset (159) | azure2-fit (`references/drop-dataset.md`) |
| data-covariance (154) | azure2-uncertainty (`references/data-covariance.md`) |
| azure2-gui-review (203) | lines 1–152: azure2-plot; 153–203: `gui/CLAUDE.md` |
| ferdinand-gnds (157) | azure2-exchange |
| azure2-thm (458 + references) | azure2-thm (same name, aligned with `docs/upgrade/THM.md`) |
| resonance-search (310 + scripts) | azure2-resonance-search (renamed; scripts moved) |
| rmfit-campaign (2,280) | rmfit-campaign (procedure); lines 62–161: azure2-fit; 306–348: azure2-resonance-search; general rules from its log: the owner skills |

Paths in the "New location" column are relative to `.claude/skills/` unless they start with `docs/`, `tests/`,
`gui/` or name a `CLAUDE.md`. Rows marked "(spill)" were moved by one writer to the skill that owns them.

## azure2-setup

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 33-44 | Particle pair masses must be precise atomic masses | azure2-setup/SKILL.md#Masses and separation energies |
| azure2-eval | 45-71 | Several compound nuclei in one .azr (natural targets) | azure2-setup/references/data-definition.md#3; behaviour docs/source/reference/data_formats.rst, docs/source/examples/natural_neon.rst |
| azure2-eval | 86-97 | Kinematics is non-relativistic unless relativistic=1 | azure2-setup/SKILL.md#Kinematics |
| azure2-eval | 98-127 | A J^pi is in the calculation only if an active level carries it (includeAllJPi) | azure2-setup/references/azr-levels.md#6 (+ version note for the zero-width singular level); SKILL.md#Levels and channels |
| azure2-eval | 163-179 | Input is LAB frame; add_extrapolation lab default; run from .azr dir; output/ must exist | azure2-setup/SKILL.md#Frames, units and directories (anecdote "caught this firsthand" dropped: narrative) |
| azure2-eval | 323-339 | An external-capture final state must be a level of the project | azure2-setup/references/azr-levels.md#7 (history as version note, builds before 2026-10-26) |
| azure2-eval | 406-414 | pyazr.nds ratio-to-Rutherford data with a lab angle column | azure2-setup/references/nds-exfor.md#4 (version note, builds before 2026-10-27) |
| azure2-eval | 664-679 | Old segment-line layouts are read right; test lines skip data-only keywords | azure2-setup/references/azr-segments.md#2 (Old layouts) and #3 |
| azure2-eval | 680-686 | An external-capture final state matching two levels is refused | azure2-setup/references/azr-levels.md#7 |
| azure2-eval | 702-717 | A lab angle the forward-branch kinematics cannot reach is refused | azure2-setup/SKILL.md#Frames, units and directories (version note builds before 2026-10-27) |
| azure2-eval | 766-784 | relativistic=1: lab angles past the relativistic maximum are refused | azure2-setup/SKILL.md#Kinematics (version note builds before 2026-10-29; personal path of the kinematics calculator dropped) |
| azure2-eval | 785-809 | Masses and separation energies are checked at load; filling from AME2020 | azure2-setup/SKILL.md#Masses and separation energies; behaviour docs/source/user_guide/particle_pairs.rst |
| azure2-eval | 994-1000 | Excited light particle: level-line column 32 | azure2-setup/SKILL.md#Levels and channels; azr-levels.md#2 |
| azure2-eval | 1726-1732 | The <levels> gamma field is NOT a reduced-width amplitude | azure2-setup/references/azr-levels.md#4 |
| azure2-eval | 2606-2614 | .azr file anatomy: <config>, A-matrix vs R-matrix | azure2-setup/references/azr-segments.md#1 |
| azure2-eval | 2615-2627 | <levels> token map, J stored as J, 2S 2L doubled | azure2-setup/references/azr-levels.md#2 (12C+alpha anecdote reduced to the 3- -> 3/2- trap) |
| azure2-eval | 2628-2638 | <segmentsData> layout; fixed-width formatting for the GUI | azure2-setup/references/azr-segments.md#2 |
| azure2-eval | 2639-2660 | Read the tail tokens; sqrtshift | azure2-setup/references/azr-segments.md#2 Keyword blocks; behaviour docs/source/reference/data_formats.rst |
| azure2-eval | 2661-2689 | sample 1 / sample <name> (natural targets, samples) | azure2-setup/references/data-definition.md#3 |
| azure2-eval | 2690-2712 | shiftgroup, shiftfngroup | azure2-setup/references/azr-segments.md#2 Keyword blocks (rmfit shard note: dated history dropped, rmfit-campaign owns sharding) |
| azure2-eval | 2713-2723 | normgroup | azure2-setup/references/azr-segments.md#2 Keyword blocks |
| azure2-eval | 2724-2761 | covariance <file> (rows, relvector freeze, refreeze, Peelle) | azure2-setup/references/azr-segments.md#2 (Covariance files paragraph); building files: azure2-uncertainty; behaviour docs/source/user_guide/chi_squared.rst |
| azure2-eval | 2762-2784 | shiftfn; .sav never frees or fixes a shift | azure2-setup/references/azr-segments.md#2 Keyword blocks (history of the penalty change dropped; penalty from zero owned by azure2-fit) |
| azure2-eval | 2785-2798 | isAdvanced composite; isUPOS; reading a tail | azure2-setup/references/azr-segments.md#2 Composite and UPOS tails |
| azure2-eval | 2799-2808 | An extrapolation cannot reproduce a UPOS segment; tail audit one-liner | azure2-setup/references/azr-segments.md#2 |
| azure2-eval | 2809-2828 | <segmentsTest>; <targetInt> key matching; mirrored test segments; <parameterSettings>; <mcmc>; param.sav key collision | azure2-setup/references/azr-segments.md#1, #4 (path of the campaign's mirroring script dropped: personal campaign directory) |
| azure2-eval | 2823-2831 | Data file format (.dat) | azure2-setup/references/azr-segments.md#5 |
| azure2-eval | 2993-3006 | A total must be summed over every channel it measures | azure2-setup/references/data-definition.md#1 |
| azure2-eval | 3007-3040 | Two segment lines that share a data file are NOT summed | azure2-setup/references/data-definition.md#2 (archive directory names, AFS path and commit hash dropped: provenance) |
| azure2-eval | 3041-3057 | Neutron totals: prefer total 1; stored output not reproducible | azure2-setup/references/data-definition.md#2 |
| azr-project-builder | 1-24 | Building an AZURE2 project from nothing; no <pairs> section; azrbuild.py | azure2-setup/references/azr-levels.md#1 (azrbuild.py named as external helper, its tree path dropped) |
| azr-project-builder | 25-61 | The 31 fields of a level line; THM fields 32, 33 | azure2-setup/references/azr-levels.md#2 (m1 m2 corrected to ATOMIC masses per docs/source/user_guide/particle_pairs.rst; the old "nuclear masses" row contradicted the docs and the mass check) |
| azr-project-builder | 62-105 | Rules the code enforces silently | azure2-setup/references/azr-levels.md#3 |
| azr-project-builder | 106-117 | Channel enumeration | azure2-setup/references/azr-levels.md#5 |
| azr-project-builder | 118-138 | Minimum viable project (azrbuild example) | azure2-setup/references/azr-levels.md#1 (API names listed; full example dropped: the module is outside the repo) |
| azr-project-builder | 139-154 | Verifying a newly built file | azure2-setup/references/azr-levels.md#8 |
| azr-project-builder | 155-165 | What this does not cover; channel radius start value | azure2-setup/references/azr-levels.md#9; SKILL.md Procedure step 3 |
| azr-project-builder | 166-172 | See also | azure2-setup/SKILL.md#References |
| nds-explorer | 1-41 | Intro; the two services at a glance | azure2-setup/references/nds-exfor.md (intro, #1) |
| nds-explorer | 42-74 | Searching EXFOR | azure2-setup/references/nds-exfor.md#2 |
| nds-explorer | 75-122 | Fetching a dataset; RTH; units | azure2-setup/references/nds-exfor.md#3, #4 (RTH history as version note builds before 2026-10-27; campaign directory names dropped) |
| nds-explorer | 123-171 | Converting to AZURE2 | azure2-setup/references/nds-exfor.md#4 |
| nds-explorer | 172-188 | Level structure, gammas, ground state | azure2-setup/references/nds-exfor.md#5 |
| nds-explorer | 189-208 | Getting the paper behind a dataset | azure2-setup/references/nds-exfor.md#6 |
| nds-explorer | 209-231 | Frames: EXFOR's target is not always AZURE2's target | azure2-setup/references/nds-exfor.md#4 (light/heavy; nonexistent nuclear_mass() replaced by pyazr.masses.atomic_mass) |
| nds-explorer | 232-262 | Two conversion bugs fixed (S-factor unit; za() light projectiles); cross-check lesson | azure2-setup/references/nds-exfor.md#3 (unit rule), #4 (charges, cross-check); dataset IDs of the bug history dropped |
| nds-explorer | 263-310 | Common failure modes (Cloudflare, local master files, x4list HTML, frames, invented errors, stale EC cache) | azure2-setup/references/nds-exfor.md#7, #3, #4 (CRC cluster path of the mirror and dates dropped; external diagnose() helper reduced to the checks it performs) |
| nds-explorer | 311-321 | See also | azure2-setup/references/nds-exfor.md intro |
| r-matrix-analysis | 1-22 | Purpose (theory counterpart) | azure2-setup/references/rmatrix-picture.md intro |
| r-matrix-analysis | 23-64 | 1. What R-matrix theory is | azure2-setup/references/rmatrix-picture.md#1 |
| r-matrix-analysis | 65-103 | 2. Vocabulary | azure2-setup/references/rmatrix-picture.md#2 |
| r-matrix-analysis | 104-135 | 3. Resonances to levels | azure2-setup/references/rmatrix-picture.md#3 |
| r-matrix-analysis | 136-171 | 4. Selecting levels | azure2-setup/references/rmatrix-picture.md#4 |
| r-matrix-analysis | 172-198 | 5. Capture vs particle channels | azure2-setup/references/rmatrix-picture.md#5 |
| r-matrix-analysis | 282-299 | 9. Angular distributions | azure2-setup/references/rmatrix-picture.md#6 |
| r-matrix-analysis | 334-357 | 11. Connecting to AZURE2 inputs | azure2-setup/references/rmatrix-picture.md#7 |
| r-matrix-analysis | 358-367 | Reference files | azure2-setup/SKILL.md#References (formalism-and-fitting.md is azure2-fit's) |
| r-matrix-analysis | references/key-papers.md | Key papers digest | azure2-setup/references/key-papers.md (copied; skill names updated) |

## azure2-calculate

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 6-15 | Running AZURE2 evaluations (CLI vs pyazr) | azure2-calculate/SKILL.md#two-ways-to-run-the-engine |
| azure2-eval | 16-30 | pyazr does not run the binary; reference docs | azure2-calculate/SKILL.md#two-ways-to-run-the-engine |
| azure2-eval | 139-151 | R-matrix route and exit status | azure2-calculate/references/physics-options.md#route; exit status in SKILL.md Procedure A step 6; docs/source/reference/command_line.rst |
| azure2-eval | 152-162 | Composite segments (RATIO 0, Scaling digits, old component lists) | azure2-calculate/references/observables.md#composite-segments; older builds in azure2-calculate/references/older-builds.md |
| azure2-eval | 180-201 | STOP before any extrapolation (<targetInt> applies to test segments) | azure2-calculate/SKILL.md#extrapolation step 1 (directory provenance dropped) |
| azure2-eval | 202-207 | Delete intEC.extrap when the grid changes | azure2-calculate/SKILL.md#extrapolation step 2 and #caches |
| azure2-eval | 208-214 | intEC.dat holds every digit (signature v5) | azure2-calculate/references/physics-options.md#cache-signatures (version note) |
| azure2-eval | 215-235 | Several sessions can be open at once (EC cache sharing in old builds) | azure2-calculate/SKILL.md#procedure-b-pyazr (version note) |
| azure2-eval | 236-248 | CLI does not read Runtime Options; three parametrizations; RMC; pyazr defaults vs archive flags | azure2-calculate/SKILL.md Procedure A step 1 and Procedure B first bullet; bases -> azure2-exchange |
| azure2-eval | 249-260 | A pyazr session can rewrite a live fit's cache | azure2-calculate/SKILL.md#procedure-b-pyazr |
| azure2-eval | 261-296 | --no-long-wavelength / EC operator / M1 EC bugs | azure2-calculate/references/physics-options.md#external-capture-operator; azure2-calculate/references/older-builds.md#external-capture-m1 |
| azure2-eval | 352-371 | Reaction rates (mode 5) were low | azure2-calculate/SKILL.md#reaction-rates; azure2-calculate/references/older-builds.md#rates |
| azure2-eval | 388-405 | Lab-frame SUM composites Jacobians | azure2-calculate/references/observables.md#composite-segments; older-builds.md#composite-segments-and-frames |
| azure2-eval | 415-428 | Capture analyzing power coupling order | azure2-calculate/references/older-builds.md#capture-analyzing-power-and-upos; spin-1/2 rule in observables.md |
| azure2-eval | 429-437 | UPOS secondary-gamma spin from decay index | azure2-calculate/references/older-builds.md#capture-analyzing-power-and-upos |
| azure2-eval | 452-464 | A param.sav in an extrapolation moved test grids | azure2-calculate/SKILL.md#extrapolation step 3; older-builds.md |
| azure2-eval | 495-511 | RATIO denominators Jacobian | azure2-calculate/references/observables.md#composite-segments; older-builds.md |
| azure2-eval | 512-524 | Hybrid-potential Coulomb derivative | azure2-calculate/references/physics-options.md#hybrid-nuclear-potential (version note) |
| azure2-eval | 539-546 | normalizations.out malformed rows | azure2-calculate/references/outputs.md#normalizationsout |
| azure2-eval | 603-614 | Phase-shift segment is one (J,l,s) wave | azure2-calculate/references/observables.md#phase-shifts; older-builds.md |
| azure2-eval | 615-627 | UPOS mixing at every order; l > 10 flagged | azure2-calculate/references/observables.md#upos-secondary-γ-rays-and-high-partial-waves; older-builds.md |
| azure2-eval | 628-638 | relativistic=1 in intEC signature | azure2-calculate/references/physics-options.md#relativistic-kinematics |
| azure2-eval | 639-654 | E1-M1 interference of EC verified (m-scheme truth traps) | behaviour: docs/source/theory/external_capture_operator.rst + tests/ec_e1m1_phase; truth-building traps -> spill-DEV.md |
| azure2-eval | 655-663 | Hybrid potential bends scattering functions only | azure2-calculate/references/physics-options.md#hybrid-nuclear-potential; docs/source/user_guide/pyazr.rst |
| azure2-eval | 687-701 | pyazr sessions own hybrid potential; write_output_files penalties; nan in parameterSettings | physics-options.md#hybrid-nuclear-potential; outputs.md#writing-the-files-from-pyazr; nan-limit behaviour: docs (dropped from skill, fixed behaviour, no user action) |
| azure2-eval | 718-738 | Total capture exit -1 with data type; component sub-point angle; Wigner limit particle-only; shift stencil | observables.md#total-capture; older-builds.md; Wigner-limit/shift-stencil remarks: dropped (internal fixes, no user action) |
| azure2-eval | 739-748 | E_x column includes E* | observables.md#excitation-energy-column |
| azure2-eval | 835-859 | Workflow A — interactive CLI, prompt order, modes table | azure2-calculate/SKILL.md#procedure-a-the-cli |
| azure2-eval | 896-904 | Total radiative capture rate | azure2-calculate/SKILL.md#reaction-rates |
| azure2-eval | 1033-1041 | Write the flags out (zsh) | azure2-calculate/SKILL.md Procedure A step 2 |
| azure2-eval | 1042-1058 | Prompt ordering silently ignores parameter file; y/n validation; --covariance-band | azure2-calculate/SKILL.md Procedure A steps 4-5; docs/source/reference/command_line.rst |
| azure2-eval | 1006-1032 (spill from azure2-check) | parameter-file prompt; printf recipes modes 1/2/3 | azure2-calculate/SKILL.md Procedure A step 4 |
| azure2-eval | 1059-1104 | Workflow B — pyazr (install, sessions, parse once, copies) | azure2-calculate/SKILL.md#two-ways-to-run-the-engine and #procedure-b-pyazr |
| azure2-eval | 1105-1147 | Two parameter conventions; transform_rwa length; short vectors | azure2-calculate/SKILL.md#two-parameter-vectors |
| azure2-eval | 1218-1231 | Data mode vs extrapolation mode; segment indexing | azure2-calculate/SKILL.md#data-mode-and-extrapolation-mode |
| azure2-eval | 1320-1337 | Writing the run's output files | azure2-calculate/SKILL.md#observables; outputs.md#writing-the-files-from-pyazr |
| azure2-eval | 1355-1369 | Calculating observables | azure2-calculate/SKILL.md#observables |
| azure2-eval | 1370-1391 | Angular-distribution coefficients as data (observable 9) | observables.md#legendre-coefficients-as-data-observable-9 |
| azure2-eval | 1392-1467 | Analyzing power (observable 7) | observables.md#analyzing-power-a_y-observable-7; docs/source/theory/polarization_*.rst |
| azure2-eval | 2056-2080 | Defining extrapolation grids; isDiff codes | azure2-calculate/SKILL.md#extrapolation step 4 |
| azure2-eval | 2148-2173 | intEC.extrap must be deleted before every extrapolation run | azure2-calculate/SKILL.md#extrapolation step 2 |
| azure2-eval | 2295-2323 | Decomposing a cross section; capture-specific | azure2-calculate/SKILL.md#decomposing-a-cross-section |
| azure2-eval | 2324-2344 | S-factor: never re-derive the conversion | azure2-calculate/SKILL.md#s-factor |
| azure2-eval | 2832-2857 | Output files; norm multiplies data; n* | azure2-calculate/references/outputs.md |
| azure2-eval | 2858-2878 | param.fit autosave semantics and restart | semantics: outputs.md#paramfit; restart recipe -> spill-azure2-fit.md |
| azure2-eval | 2879-2888 | missing checks/ directory; normalizations.out; MINOS; rates; caches | azure2-calculate/SKILL.md Procedure A step 3; outputs.md table |
| azure2-eval | 3103-3114 | Examples shipped with pyazr (table) | azure2-calculate/references/examples.md |
| azure2-eval | 3115-3118 | <targetInt> applied to extrapolation segments; folded shape via dummy data | azure2-calculate/SKILL.md#extrapolation step 1 |
| azure2-eval | 3120-3128 | add_level 16-digit energies misread | spill-azure2-fit.md |
| azure2-eval | 3194-3206 | Closed neutral channels closed form | physics-options.md#coulomb-functions (perf numbers dropped: performance history) |
| azure2-eval | 3207-3222 | Universal Coulomb table; --gsl-coul wrong | SKILL.md#caches; physics-options.md#coulomb-functions; docs/source/developer/coulomb_table.rst |
| azure2-eval | 822-824 (spill from azure2-targets) | --gsl-coul speed/accuracy tradeoff | physics-options.md#coulomb-functions |
| azure2-eval | 3279-3294 | EC energy table; shifted capture segments; segfault range | physics-options.md#external-capture-energy-table; docs/source/developer/ec_table.rst; segfault-range note dropped (a dev-only window fixed the next day) |
| azure2-eval | 3295-3318 | S-wave scattering lengths | observables.md#s-wave-scattering-lengths; docs/source/theory/scattering_length.rst |
| azure2-eval | 3319-3323 | pyazr reload of intEC.dat gives wrong chi2 | physics-options.md#cache-signatures |
| azure2-eval | 3455-3469 | isDiff 8 removed; enter P as inverse A_y | observables.md#analyzing-power (last bullets); archive file list dropped (provenance) |
| azure2-eval | 3634-3647 | Adaptive cross-section tables | azure2-calculate/SKILL.md#adaptive-cross-section-tables |
| azure2-eval | 885-895 | mode-5 reaction-rate grid (spill from fit) | azure2-calculate/SKILL.md#Reaction rates |

## azure2-fit

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 72-85 | A width of exactly 0 is a fixed parameter | azure2-fit/SKILL.md#2. Free and fixed |
| azure2-eval | 297-307 | Energy-shift penalties are measured from zero | azure2-fit/references/normalizations-and-shifts.md#The two homes (version note); SKILL.md#3 |
| azure2-eval | 372-387 | Physical-parameter objectives ignored the normalizations | azure2-fit/references/pyazr-fitting.md#χ², residuals and the objective (version note) |
| azure2-eval | 479-494 | MIGRAD nuisance priors on level parameters hit the wrong parameter | azure2-fit/references/running-fits.md#Nuisance priors on level parameters (version note) |
| azure2-eval | 547-560 | pyazr penalties() shiftgroup follower; teichmann_wigner 1.5 | normalizations-and-shifts.md#The two homes (shiftgroup); pyazr-fitting.md#A subset fit (Wigner units) |
| azure2-eval | 860-884 | Mode 2's MIGRAD has no call limit; plateau not convergence | azure2-fit/references/running-fits.md#Convergence; SKILL.md#5; docs/source/user_guide/fitting.rst |
| azure2-eval | 885-895 | Mode 5 tabulates sigma(E) on a resonance-aware grid | spilled to S/spill-azure2-calculate.md (reaction rates) |
| azure2-eval | 1232-1262 | Inspecting a model; m.parameters RuntimeError; LevelKey | pyazr-fitting.md#Inspecting the model; parameters-and-sav.md#Which physical level |
| azure2-eval | 1263-1319 | χ², residuals, objective; threads; fit against objective; penalty rows | pyazr-fitting.md#χ², residuals and the objective; SKILL.md#3 |
| azure2-eval | 1338-1354 | Free-vector helpers | pyazr-fitting.md#Free-vector helpers |
| azure2-eval | 1468-1512 | Editing the model: AzrModel add/remove/deactivate; J-group channel set; existing pairs; gamma seed | editing-levels.md#Adding, removing, deactivating; SKILL.md#6 |
| azure2-eval | 1513-1528 | Deactivating renumbers names; how names are assigned (FillMnParams) | parameters-and-sav.md#How the names are assigned; #Remapping after a level edit |
| azure2-eval | 1529-1539 | Read .sav by NAME; parametrization tag row | parameters-and-sav.md#Reading a parameter file (dated commit window replaced by description of the three file forms) |
| azure2-eval | 1540-1553 | load_sav, not slicing full_rwa_from_sav | parameters-and-sav.md#Loading a parameter file in pyazr |
| azure2-eval | 1554-1567 | pyazr names vs param.sav names; plain copy unsafe after edit | parameters-and-sav.md#Which physical level; #Remapping after a level edit |
| azure2-eval | 1568-1606 | Removing a level: file vs runtime table; other .azr edits; total-capture observable; analyzing-power | editing-levels.md#Adding, removing, deactivating; #Other .azr edits |
| azure2-eval | 1623-1631 | Adding/removing segments invalidates EC cache | editing-levels.md#Other .azr edits |
| azure2-eval | 1632-1655 | Removing a segmentsData row breaks param.sav | parameters-and-sav.md#Segment keys |
| azure2-eval | 1656-1671 | save_fit snapshots and verifies; beta sign bug | editing-levels.md#save_fit (version note) |
| azure2-eval | 1672-1725 | Fixed parameters fixed in rwa space; EC integrals; ignore_externals; sub-point grid re-placement; graded tail | editing-levels.md#save_fit (two version notes); docs/source/user_guide/experimental_effects.rst |
| azure2-eval | 1733-1738 | A Parameter's pair is the engine number | parameters-and-sav.md#Which physical level |
| azure2-eval | 1739-1755 | Fitted norms/shifts only in .sav; apply_fit; save_fit_to_azr.py | normalizations-and-shifts.md#Consequences; editing-levels.md#apply_fit and the example script |
| azure2-eval | 1756-1792 | Background pole restricted to the channels that need it (16O case); fix one width; renumbering | background-poles.md#Restricting a pole |
| azure2-eval | 1854-1927 | Normalizations live in two places + consequences + audit script | normalizations-and-shifts.md#The two homes; #Consequences |
| azure2-eval | 1928-2001 | Fitting a data set's energy scale (scan, form, seed, checks, one calibration per experiment) | normalizations-and-shifts.md#Fitting a data set's energy calibration |
| azure2-eval | 2002-2016 | What a snapshot cannot carry: norms, check dump keywords | editing-levels.md#What a snapshot still cannot carry |
| azure2-eval | 2017-2034 | Inactive levels lines corrupted bakes; verify bake with blank file | editing-levels.md#Adding... (Inactive lines, version note); #What a snapshot still cannot carry |
| azure2-eval | 2035-2040 | m.datasets[i] vs active_datasets | pyazr-fitting.md#χ², residuals and the objective; SKILL.md#Pitfalls |
| azure2-eval | 2041-2055 | Verify tolerance hides fragility | editing-levels.md#save_fit |
| azure2-eval | 2081-2141 | Fitting from Python: pyazr does not fit; objective; Wigner; least_squares subset; iminuit; rwa space; MCMC pattern | pyazr-fitting.md (intro, #A subset fit); SKILL.md#5 |
| azure2-eval | 2142-2147 | Three traps header (BBN campaign) | dropped: provenance header only; traps placed individually |
| azure2-eval | 2174-2206 | Zero Jacobian column freezes least_squares; zero amplitudes never move (n1 ANC) | pyazr-fitting.md#Traps; SKILL.md#2 |
| azure2-eval | 2207-2214 | residual_jacobian raises | pyazr-fitting.md#Traps |
| azure2-eval | 2215-2230 | Fitting a model built from hand-chosen seeds (stages) | pyazr-fitting.md#A model built from hand-chosen seeds; SKILL.md#4 (external Fitter tool path dropped: not in repo) |
| azure2-eval | 2231-2294 | Interfering channels stuck on wrong branch; guard vs runaway | interference.md#Wrong branches; #Guard or runaway; SKILL.md#4 |
| azure2-eval | 2889-2927 | Long fits: memory, leak, kill/restart from param.fit | running-fits.md#Long fits on a batch cluster (version note; site-specific qstat/skill references generalised) |
| azure2-eval | 3223-3240 | GUI-saved limits and nuisance priors on other levels | parameters-and-sav.md#Limits and priors written by the GUI (version note); docs/source/user_guide/fitting.rst; archive scan counts dropped (archive-specific) |
| azure2-eval | 3324-3335 | Nuisance priors on ANCs/gamma widths wrong penalties | running-fits.md#Nuisance priors (superseded cause merged with 479-494; test recipe kept) |
| azure2-eval | 3422-3432 | Threaded fits with energy shifts could crash | running-fits.md#Long fits (version note) |
| azure2-eval | 3589-3618 | Level-merge tolerance 1 keV | parameters-and-sav.md#Two levels of one J^π closer than 1 keV (updated to current behaviour per docs/source/user_guide/levels_channels.rst; version note); local patch/binary paths dropped (machine-specific history) |
| azure2-eval | 3648-3655 | Correction: fixed energy shifts in .sav are applied | normalizations-and-shifts.md#Consequences; SKILL.md#2 |
| azure2-eval | 3656-3687 | --use-gradient G2 seeding, same-minimizer rule, diagnostics, parallel chi2 | running-fits.md#The analytic gradient; diagnostics also added to docs/source/reference/command_line.rst |
| azure2-eval | 3688-3721 | --use-gradient on penalized norms, A_y with target integration, many shifts; FD step; speed tip | running-fits.md#The analytic gradient (version note; internal fix list (1)-(7) dropped as code history; points-per-width tip kept) |
| rmfit-campaign | 62-161 | First fits: starting parameters; rmfit init; objective options; policy; cost check | azure2-fit/references/first-fits.md; SKILL.md#1 |
| r-matrix-analysis | 199-232 | §6 Interference (same-J, different-J, ghost) | interference.md |
| r-matrix-analysis | 233-254 | §7 Background poles | background-poles.md#What a background pole means |
| r-matrix-analysis | 255-281 | §8 Boundary conditions and Brune | interference.md#Formal, observed, Brune and Park parameters |
| r-matrix-analysis | 300-333 | §10 Common pitfalls | interference.md#Pitfalls; last bullet (per-point judging) → see note: judging per point belongs to azure2-check, kept here as a one-line pointer and spilled in full |
| r-matrix-analysis | references/formalism-and-fitting.md | all | azure2-fit/references/formalism.md (copied; cross-references updated) |
| azure2-legacy-params | SKILL.md all | old-format parameter files | azure2-fit/references/legacy-params.md; SKILL.md#1 |
| azure2-legacy-params | convert_sav_names.py | converter | azure2-fit/scripts/convert_sav_names.py (moved; docstring hash/date removed) |
| azure2-drop-dataset | SKILL.md all | dropping a data set | azure2-fit/references/drop-dataset.md; SKILL.md#7 |
| rmfit-campaign | 682-695 | MINUIT initial steps; narrow levels drift (spill) | azure2-fit/SKILL.md#Pitfalls, short |
| rmfit-campaign | 1101-1130 | inert background poles: levelFix vs channelFix; zero widths never move (spill) | azure2-fit/references/background-poles.md#Restricting a pole to the channels that need it |
| rmfit-campaign | 1191-1203, 2218-2231 | seed dead channels in γ_W units; tiny physical seeds; Denominator warning (spill) | azure2-fit/references/background-poles.md#Restricting a pole to the channels that need it |
| rmfit-campaign | 1547-1562, 1980-1983 | threshold guard on free energies; Brune caps an l=0 width (spill) | azure2-fit/SKILL.md#Pitfalls, short |
| azure2-eval | 3122-3128 | add_level 16-digit energy trap (spill from calculate) | azure2-fit/references/editing-levels.md#Adding, removing, deactivating |
| azure2-eval | 2858-2878 (restart part) | restarting a killed fit from param.fit (spill from calculate) | azure2-fit/references/running-fits.md#Long fits on a batch cluster |
| azure2-legacy-params | 1-41 | The trap | azure2-fit/references/legacy-params.md#The two naming schemes |
| azure2-legacy-params | 42-54 | Detect | azure2-fit/references/legacy-params.md#Detect |
| azure2-legacy-params | 55-71 | Convert | azure2-fit/references/legacy-params.md#Convert |
| azure2-legacy-params | 72-111 | Verify (all three, in order) | azure2-fit/references/legacy-params.md#Verify (all three, in order) |
| azure2-legacy-params | 112-117 | Related | azure2-fit/references/legacy-params.md#Related |
| azure2-drop-dataset | 1-31 | The rule | azure2-fit/references/drop-dataset.md#The rule: deactivate, never delete |
| azure2-drop-dataset | 32-63 | The drop did not reach the model being run | azure2-fit/references/drop-dataset.md#The drop must reach the model being run |
| azure2-drop-dataset | 64-99 | Doing it | azure2-fit/references/drop-dataset.md#Doing it |
| azure2-drop-dataset | 100-134 | Verify: the difference must be exactly the dropped segment | azure2-fit/references/drop-dataset.md#Verify: the difference is exactly the dropped segment |
| azure2-drop-dataset | 135-146 | Running the refit afterwards | azure2-fit/references/drop-dataset.md#The refit afterwards |
| azure2-drop-dataset | 147-153 | Record it | azure2-fit/references/drop-dataset.md#Record it |
| azure2-drop-dataset | 154-159 | Bringing a data set back | azure2-fit/references/drop-dataset.md#Bringing a data set back |

## azure2-check

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 128-133 | A `nan` point / Total-Chi-Squared: nan means a calculation FAILED | azure2-check/SKILL.md#Failures that look like results (version note builds before 2026-11-01) |
| azure2-eval | 1001-1005 | CODATA 2018 constants | azure2-check/SKILL.md#Builds, constants and logs (version note builds before 2026-11-05) |
| azure2-eval | 1006-1032 | external parameter file prompt; printf recipes for modes 1/2/3 | spilled to azure2-calculate (S/spill-azure2-calculate.md) |
| azure2-eval | 1607-1622 | Never fit a derived/extrapolated point alongside the data | azure2-check/SKILL.md#Data and their definitions (archive file names dropped as provenance) |
| azure2-eval | 1793-1853 | A background pole can run away to saturation | azure2-check/SKILL.md#Saturated background poles + references/background-pole-saturation.md (directory paths and hash dropped) |
| azure2-eval | 2345-2354 | Evaluation recipes; Is a level needed? frozen Δχ² | azure2-check/SKILL.md#Is a level needed? Is a new one warranted? |
| azure2-eval | 2355-2372 | Changing the channel radius; second set_channel_radius undid the first | azure2-check/SKILL.md#Channel-radius sensitivity (version note builds before 2026-11-02) |
| azure2-eval | 2384-2388 | Is a new level/channel warranted? | azure2-check/SKILL.md#Is a level needed? Is a new one warranted? |
| azure2-eval | 2389-2409 | Dimensionless widths (dimensionless_widths, WidthTable) | azure2-check/SKILL.md#Wigner and Weisskopf limits; behaviour also in docs/source/user_guide/pyazr.rst "Dimensionless widths" |
| azure2-eval | 2410-2437 | GUI/width form vs API form of the Wigner limit; γ-ray strengths | azure2-check/references/widths.md |
| azure2-eval | 2928-2940 | Verifying a run | azure2-check/SKILL.md#The checklist item 1 |
| azure2-eval | 2941-2958 | Data points too close to a channel threshold explode χ² | azure2-check/SKILL.md#Failures that look like results (links docs/source/reference/troubleshooting.rst) |
| azure2-eval | 2959-2967 | Before submitting a fit job after a structural edit, run a plain calculate | azure2-check/SKILL.md#Failures that look like results |
| azure2-eval | 2968-2978 | A starting point far above the baseline | azure2-check/SKILL.md#Failures that look like results |
| azure2-eval | 2979-2986 | Tolerance for "the numbers don't match" | azure2-check/SKILL.md#Tolerance for "the numbers do not match" |
| azure2-eval | 2987-2992 | Checking what the data and the figures mean (intro) | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 3058-3062 | A large, angle-dependent normalization is a symptom | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 3063-3068 | Find outliers by comparing data sets | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 3069-3075 | Rebuild combined data files with a script, diff row by row | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 3094-3102 | Legendre coefficients are a data-consistency tool | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 3336-3384 | Which experimental-effect kernel; target integration looks wrong for a thin target | spilled to azure2-targets (S/spill-azure2-targets.md) |
| azure2-eval | 3385-3391 | Run logs: torn "Segment #" progress bars | azure2-check/SKILL.md#Builds, constants and logs (version note builds before 2026-11-06) |
| azure2-eval | 3584-3588 | include/Constants.h defines pi = 3.141592650 | azure2-check/SKILL.md#Builds, constants and logs (merged with CODATA; now a version note: pi exact since 2026-11-05) |
| (CLAUDE.md) | Rules: Frames and units | "full list of conventions is in azure2-eval Golden rules" | azure2-check/SKILL.md#Frames and units (new compact table, sourced from docs/source/reference/*.rst) |
| rmfit-campaign | 1265-1267 | secondary-γ detectors from a J=1/2 state must agree (spill) | azure2-check/SKILL.md#Data and their definitions |
| rmfit-campaign | 1279-1293 | one bad point can hold a converged fit hostage (spill) | azure2-check/SKILL.md#Data and their definitions |
| rmfit-campaign | 1834-1836 | data handed over between analyses: compare by value (spill) | azure2-check/SKILL.md#Data and their definitions |
| azure2-eval | 639-654 (trap 3) | ANC must be physical for the channel radius (spill from calculate) | azure2-check/SKILL.md#Failures that look like results |
| r-matrix-analysis | 316-333 (last bullet of §10) | judge a fit point by point (spill from fit) | azure2-check/SKILL.md#Data and their definitions |

## azure2-uncertainty

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 749-765 | pyazr transform_rwa in a sampler no longer grows memory | azure2-uncertainty/SKILL.md#4. Sampling on the objective (version note builds before 2026-10-28; diagnosis cache_stats) |
| azure2-eval | 981-993 | Prior covariance (SAMMY PUP / sequential Bayesian update) | azure2-uncertainty/SKILL.md#6. Prior parameter covariance (behaviour in docs/source/user_guide/chi_squared.rst "Prior parameter covariance") |
| azure2-eval | 2373-2383 | Sampling the radius in a likelihood leaked memory | azure2-uncertainty/SKILL.md#4. Sampling on the objective (version note builds before 2026-11-04) |
| data-covariance | 1-4 | frontmatter | azure2-uncertainty/SKILL.md description (data covariance triggers) |
| data-covariance | 6-16 | Standing rule: no data set without a covariance file; pointers | azure2-uncertainty/SKILL.md#1. Data covariance + references/data-covariance.md intro (DeBoer quote dropped: quote) |
| data-covariance | 18-30 | 1. Inventory first | azure2-uncertainty/references/data-covariance.md#1 |
| data-covariance | 32-61 | 2. What goes into C, and what stays out (table) | azure2-uncertainty/references/data-covariance.md#2 (summary in SKILL.md#1) |
| data-covariance | 63-84 | 3. Our own data: recipe | azure2-uncertainty/references/data-covariance.md#3 (worked-example campaign paths dropped: outside repo, personal provenance) |
| data-covariance | 86-109 | 4. Literature data: recipe | azure2-uncertainty/references/data-covariance.md#4 (named-paper examples generalised) |
| data-covariance | 111-121 | 5. Checks before the fit | azure2-uncertainty/references/data-covariance.md#5 |
| data-covariance | 123-154 | 6. Traps seen so far | azure2-uncertainty/references/data-covariance.md#6 (campaign paths and dates dropped; numbers kept) |
| (new) | — | Hessian/LM covariance, MINOS, profile likelihood, bands (pyazr.bands), rate propagation, model averaging (pyazr.modelavg) | azure2-uncertainty/SKILL.md#2,3,5,7 (grounded in docs/source/reference/output_files.rst, command_line.rst, pyazr/bands.py, pyazr/modelavg.py) |
| rmfit-campaign | 2155-2217 | relvector covariance + free norm biased; t0 method (spill) | azure2-uncertainty/references/data-covariance.md#6. Traps |

## azure2-plot

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 825-834 | Showing a result in the GUI: recalculate first | azure2-plot/SKILL.md#showing-a-result-in-the-gui |
| azure2-eval | 3076-3093 | Check every figure against the fit column; UPOS, coarse grids, threshold spikes; data x norm at shift | azure2-plot/SKILL.md#rules-for-every-figure |
| rmfit-campaign | 1618-1624 (spill) | Figures of UPOS data need data-mode curves | azure2-plot/SKILL.md#rules-for-every-figure |
| azure2-gui-review | 1-5 | frontmatter | azure2-plot/SKILL.md frontmatter (merged into azure2-plot description) |
| azure2-gui-review | 6-35 | The rule that keeps getting broken (GUI plots output/) | azure2-plot/SKILL.md#showing-a-result-in-the-gui steps 1-5 |
| azure2-gui-review | 36-68 | Doing it (close by pid, calculate, confirm, open) | azure2-plot/SKILL.md#showing-a-result-in-the-gui code block (~/bin/AZURE2 -> $AZURE2_BINARY) |
| azure2-gui-review | 69-87 | Open it on the screen the person is looking at | azure2-plot/SKILL.md#open-it-on-the-screen-the-person-is-looking-at |
| azure2-gui-review | 88-121 | Things that make the picture wrong | azure2-plot/SKILL.md#things-that-make-the-gui-picture-wrong-without-a-warning |
| azure2-gui-review | 122-136 | Budget the time | azure2-plot/SKILL.md#budget-the-time (qsub -> batch system) |
| azure2-gui-review | 137-152 | Before handing the window over; segment keys | azure2-plot/SKILL.md#before-handing-the-window-over |
| (new) | - | figure recipes: pulls, angular distributions, S with bands and decomposition, rate ratios | azure2-plot/references/recipes.md (written from pyazr API: calculate_rwa, segment_norms, extrapolation_bands, angular_distribution) |

## azure2-targets

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 308-322 | An energy shift moved sub-points to the wrong angle | azure2-targets/SKILL.md#9 (version note); behaviour: docs/source/user_guide/experimental_effects.rst |
| azure2-eval | 340-351 | Gaussian kernels now reach 5 sigma | azure2-targets/SKILL.md#4 (rule + version note); docs/source/user_guide/experimental_effects.rst "Gaussian Beam Spread" |
| azure2-eval | 438-451 | Sub-points moved by the shift at their own energy | azure2-targets/SKILL.md#9 |
| azure2-eval | 525-538 | Thick target with beam spread zeroed first point | azure2-targets/SKILL.md#1 (definition), #3 (ε ≤ 0 refused), #11 (version note) |
| azure2-eval | 561-578 | Straggling coefficient is c.m.; 0 means Bohr | azure2-targets/SKILL.md#2 |
| azure2-eval | 579-590 | Thin-target energy loss is the exact range integral | azure2-targets/SKILL.md#1, #3, #11 |
| azure2-eval | 591-602 | targetInt line shared by several entrance pairs converts per point | azure2-targets/SKILL.md#2 |
| azure2-eval | 810-821 | Attenuation (Q) coefficients act on the lab distribution | azure2-targets/SKILL.md#10 |
| azure2-eval | 822-824 | --gsl-coul speed/accuracy tradeoff | spilled to spill-azure2-calculate.md |
| azure2-eval | 905-934 | Sample transmission and self-shielded yield (+ builds before ignore) | azure2-targets/SKILL.md#6 |
| azure2-eval | 935-960 | Free-gas Doppler broadening | azure2-targets/SKILL.md#1 table, #6; sub-point count example in #4 |
| azure2-eval | 961-963 | Fitted Doppler T ... (intro) | azure2-targets/SKILL.md#7 |
| azure2-eval | 964-971 | Fitted Doppler T | azure2-targets/SKILL.md#7 table |
| azure2-eval | 972-978 | Fitted thickness | azure2-targets/SKILL.md#7 table |
| azure2-eval | 979-980 | rmfit: SingleSession only | azure2-targets/SKILL.md#7 |
| azure2-eval | 3130-3149 | Beam-profile experimental effect | azure2-targets/references/udr-and-beam-profile.md#Beam profile; GUI digits note in azure2-targets/SKILL.md#11 |
| azure2-eval | 3150-3157 | Convolution coefficients bound ascending | azure2-targets/SKILL.md#2 |
| azure2-eval | 3158-3179 | Narrow resonance inside a target integration mis-integrated | azure2-targets/SKILL.md#4 (version note, resonanceWidthMultiplier, .so rebuild rule); archive scope list dropped: archive-specific counts |
| azure2-eval | 3180-3193 | Target sub-points recomputed Coulomb functions | azure2-targets/SKILL.md#11 (version note + timing recipe) |
| azure2-eval | 3241-3258 | Electron screening | azure2-targets/SKILL.md#7 table and bullets; #8 (also folds test segments); docs "Electron Screening" |
| azure2-eval | 3259-3278 | Fitted target density (+ first evaluation at shift 0) | azure2-targets/SKILL.md#7 |
| azure2-eval | 3392-3406 | Target-integration windows: outside a point is sigma(E_bar) | azure2-targets/SKILL.md#5 |
| azure2-eval | 3407-3421 | Automatic windows ("auto") | azure2-targets/SKILL.md#5 |
| azure2-eval | 3433-3454 | Compound targets: stopping and density count the same thing | azure2-targets/SKILL.md#3 |
| azure2-eval | 3470-3490 | One Gaussian beam spread; sigmalab | azure2-targets/SKILL.md#2; docs "Gaussian Beam Spread" (normalisation numbers) |
| azure2-eval | 3491-3500 | Energy-dependent convolution, legacy convention | azure2-targets/SKILL.md#2 |
| azure2-eval | 3501-3536 | SAMMY UDR | azure2-targets/references/udr-and-beam-profile.md#UDR; azure2-targets/SKILL.md#6; docs "User-defined numerical resolution function"; branch/merge history dropped |
| azure2-eval | 3537-3558 | Channel-threshold blow-up caught by sub-point grids | azure2-targets/SKILL.md#4 (version note + AZURE_DUMP_NONFINITE); env var added to docs/source/reference/command_line.rst; 1e100 to Minuit2 in docs |
| azure2-eval | 3619-3633 | targetInt line shared by several segments convolved too narrow | azure2-targets/SKILL.md#11 |
| azure2-eval | 3336-3384 | Which experimental-effect kernel; target integration is not a bug; numpy check; grid-collapse regression; constant-loss foils (spilled by azure2-check writer) | azure2-targets/SKILL.md#1 "Which kernel"; references/target-integration-check.md (recipe + version note) |
| rmfit-campaign | 707-743 | narrow resonance in <targetInt>: grid-setting cost table (spill) | azure2-targets/SKILL.md#4. Narrow resonances and sub-point grids |

## azure2-thm

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-thm | 1-4 | frontmatter description | azure2-thm/SKILL.md frontmatter (shortened to 564 chars; old-skill names replaced) |
| azure2-thm | 6-25 | intro, what THM data can and cannot give | azure2-thm/SKILL.md#intro (pointers now to references/thm-block.md, azure2-setup, azure2-fit, docs/upgrade/THM.md); new "Rules that hold throughout" block from docs/upgrade/THM.md |
| azure2-thm | 26-44 | When to use, Prerequisites | azure2-thm/SKILL.md#When to use ("thm branch" generalised to "a branch that carries the THM code") |
| azure2-thm | 45-117 | Stage 1 gather and vet; pitfalls met | azure2-thm/SKILL.md#Stage 1 (nds-explorer -> azure2-setup) |
| azure2-thm | 118-175 | Stage 2 build the project | azure2-thm/SKILL.md#Stage 2 (resolution recorded in c.m. per THM.md; lab σ for present <targetInt>) |
| azure2-thm | 176-205 | Stage 3 baseline model table | azure2-thm/SKILL.md#Stage 3 ("pre-Sep-2026" -> "former default"; kinematics legacy note; background ≥ 0) |
| azure2-thm | 206-254 | Stage 4 fit | azure2-thm/SKILL.md#Stage 4 (+ step 5: no energy shifts on THM segments) |
| azure2-thm | 255-285 | Stage 5 derive quantities, anchor formula | azure2-thm/SKILL.md#Stage 5 |
| azure2-thm | 286-349 | Stage 6 model dependence, driver, weights, robustness | azure2-thm/SKILL.md#Stage 6 |
| azure2-thm | 350-372 | Stage 7 cross-checks | azure2-thm/SKILL.md#Stage 7 (reproduce-first also stated up front) |
| azure2-thm | 373-391 | Stage 8 report | azure2-thm/SKILL.md#Stage 8 |
| azure2-thm | 392-417 | Pitfalls checklist | azure2-thm/SKILL.md#Pitfalls checklist (hash 49b6c30 dropped, test path kept; 4 THM.md items added) |
| azure2-thm | 418-443 | Pointers | azure2-thm/SKILL.md#Pointers (+ thm-block.md, THM.md rows) |
| azure2-thm | 444-459 | Memory and process discipline | azure2-thm/SKILL.md#Memory and process discipline |
| azure2-thm | references/calibration.md | all | azure2-thm/references/calibration.md ("(October 2026)" dropped) |
| azure2-thm | references/model-choices.md | all | azure2-thm/references/model-choices.md ("pre-Sep-2026", "before Oct 2026", hash df8c2c8 removed; content kept) |
| azure2-thm | references/provenance-checklist.md | all | azure2-thm/references/provenance-checklist.md (unchanged) |
| azure2-thm | references/reproduction.md | all | azure2-thm/references/reproduction.md (unchanged) |
| azure2-eval | 2518-2525 | Trojan Horse (THM) segments (intro) | azure2-thm/references/thm-block.md#intro |
| azure2-eval | 2526-2537 | A THM segment adds 10 to its observable code | azure2-thm/references/thm-block.md#THM segments |
| azure2-eval | 2538-2543 | The analytic gradient does not cover THM | azure2-thm/references/thm-block.md#THM segments |
| azure2-eval | 2544-2558 | Widths may be reduced width amplitudes | azure2-thm/references/thm-block.md#Widths given as reduced-width amplitudes |
| azure2-eval | 2559-2561 | Never convert an rwa-declared width | azure2-thm/references/thm-block.md#Widths given as reduced-width amplitudes |
| azure2-eval | 2562-2565 | Parameter.wigner_limit is gamma^2_W in MeV | azure2-thm/references/thm-block.md#Widths given as reduced-width amplitudes |
| azure2-eval | 2566-2605 | The <thm> block selects the HOES observable variant; AzrModel editors; reference cases | azure2-thm/references/thm-block.md#The <thm> block, #Editing it from Python, #Reference projects ("default since 2026-09-25" dropped: pure history, the default is stated) |

## azure2-resonance-search

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| resonance-search | 1-4 | frontmatter description | azure2-resonance-search/SKILL.md frontmatter (rewritten, 591 chars) |
| resonance-search | 6-17 | purpose; table of tools in other skills | azure2-resonance-search/SKILL.md#intro (old skill names -> rmfit-campaign, azure2-fit, azure2-check, azure2-setup) + PYTHONPATH via $RMFIT_HOME/$AZURE2_SRC |
| resonance-search | 18-19 | Worked history (13C+alpha directories) | dropped: session/directory provenance (RULES §3.5); the numbers it carried are kept in references/ |
| resonance-search | 21-39 | §0 Ground rules 1-5 | azure2-resonance-search/SKILL.md#0 Ground rules (quotes and dates removed; "narrow before broad grows into broad strength" lesson added to rule 3) |
| resonance-search | 41-81 | §1 Locate: residual dump, rank_map, ranking method, density, levelscan equivalent | azure2-resonance-search/SKILL.md#1 Locate |
| resonance-search | 83-126 | §1a data-only peak finder; blind recall table; usage | SKILL.md#1a (procedure, use) + references/locators.md#Peak finder (table, precision note, agreement example, open improvements) |
| resonance-search | 128-169 | §1b blur; scripts; validity numbers; blur_map scoring; global sign; first finding; curve grouping | SKILL.md#1b (procedure, rules) + references/locators.md#Blur validity, #Global sign table; curve grouping -> SKILL.md#1 "Curve grouping" |
| resonance-search | 171-196 | §2 Characterise a-f | azure2-resonance-search/SKILL.md#2 Characterise |
| resonance-search | 198-207 | §3 Hypotheses table | azure2-resonance-search/SKILL.md#3 Hypotheses |
| resonance-search | 209-227 | §4 Screen a-d | azure2-resonance-search/SKILL.md#4 Screen (dates dropped) |
| resonance-search | 229-246 | §5 level add with control | azure2-resonance-search/SKILL.md#5 Test (+ total vs per-channel floor, free energy check) |
| resonance-search | 248-262 | §6 Verdict | azure2-resonance-search/SKILL.md#6 Verdict (gain written as control minus add so that a positive gain is an improvement) |
| resonance-search | 264-270 | §7 After an adoption | azure2-resonance-search/SKILL.md#7 (+ propagate step) |
| resonance-search | 272-285 | §8 Traps | azure2-resonance-search/SKILL.md#8 Pitfalls ("pre-09-19 adds were fixed" -> check that an added level's energy is free) |
| resonance-search | 287-294 | §9 Costs | azure2-resonance-search/SKILL.md#9 Costs |
| resonance-search | 296-310 | §10 Refinement log | dropped as log; lifted rules: cluster by density map not interval merging (SKILL.md#1), recovery needs the chance rate at the same tolerance (SKILL.md#1a), curve grouping by detector (SKILL.md#1) |
| resonance-search | scripts/*.py (8 files) | dump_residuals, rank_map, peakfind, eval_peaks, make_broad_only, blur_compare, blur_map, blur_figure | azure2-resonance-search/scripts/ (moved; hard-coded /users/... sys.path entries replaced by $AZURE2_SRC or the script's own directory; DeBoer quotes and dates removed from docstrings; behaviour unchanged, all compile) |
| rmfit-campaign | 306-314 | §3c intro: Ex 9.4 case narrative | azure2-resonance-search/references/width-floor-recipe.md#intro (as a worked case with numbers, narrative and dates dropped) |
| rmfit-campaign | 315-319 | 1. Localise and characterise | width-floor-recipe.md#1 |
| rmfit-campaign | 320-329 | 2. Width-floor level add, every J^π | width-floor-recipe.md#2 (`setenv` -> env var; total vs per-channel floor) |
| rmfit-campaign | 330-335 | 3. Cross-checks a-c | width-floor-recipe.md#3 |
| rmfit-campaign | 336-337 | 4. Full free refit | width-floor-recipe.md#4 |
| rmfit-campaign | 338-341 | 5. Propagate (refresh_frozen recipe) | width-floor-recipe.md#5 (script path dropped, recipe kept) |
| rmfit-campaign | 342-345 | 6. Then the next region | width-floor-recipe.md#6 |
| rmfit-campaign | 347-348 | Cost; full narrative directories | width-floor-recipe.md (cost kept; directory list dropped: provenance) |

## azure2-exchange

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 134-138 | --no-brune on a Brune fit may have no standard equivalent | azure2-exchange/SKILL.md#1 "Rules for the standard basis"; docs/source/reference/troubleshooting.rst |
| azure2-eval | 465-478 | Non-Brune transform handed eigenvalues in file order | azure2-exchange/SKILL.md#1 (version note); docs/source/user_guide/levels_channels.rst |
| azure2-eval | 1148-1156 | Park's parametrization (intro) | azure2-exchange/SKILL.md#1 table; docs/source/reference/command_line.rst "Brune and Park parametrizations" |
| azure2-eval | 1157-1164 | A Park amplitude is that channel's observed width alone | azure2-exchange/SKILL.md#1 "Choose Park" |
| azure2-eval | 1165-1175 | J must stay positive | azure2-exchange/SKILL.md#1 "Rules for Park fits" |
| azure2-eval | 1176-1201 | Parameter files carry the basis (param.fit untagged; plain row) | azure2-exchange/SKILL.md#2 |
| azure2-eval | 1202-1206 | Analytic derivatives work | azure2-exchange/SKILL.md#1 |
| azure2-eval | 1207-1217 | Same physics; GUI; transform_out | azure2-exchange/SKILL.md#1; pointer to a private comparison directory dropped (personal path) |
| azure2-eval | 2438-2463 | Comparing against an ENDF/B evaluation: getting and parsing | azure2-exchange/references/endf-comparison.md |
| azure2-eval | 2464-2506 | Elastic scattering of a charged projectile (MF=6 LAW=5) | references/endf-comparison.md; σ_NI sign rule in azure2-exchange/SKILL.md#4 |
| azure2-eval | 2507-2517 | Always check the tabulated energy grid | azure2-exchange/SKILL.md#4 and references/endf-comparison.md#Energy coverage |
| ferdinand-gnds | 1-17 | Install table (paths, versions) | azure2-exchange/SKILL.md#3 ($FERDINAND_HOME, PYTHONPATH); install dates dropped |
| ferdinand-gnds | 20-39 | Running | azure2-exchange/SKILL.md#3 |
| ferdinand-gnds | 41-56 | Known problem 1: mixed masses | azure2-exchange/SKILL.md#3 problem 1 (archive file count dropped) |
| ferdinand-gnds | 58-92 | Known problem 2: settings lost; token diff script | azure2-exchange/SKILL.md#3 problem 2; script and token list in references/ferdinand-checks.md |
| ferdinand-gnds | 94-113 | Reading GNDS with FUDGE | references/ferdinand-checks.md |
| ferdinand-gnds | 115-157 | Cross-section reconstruction; no capture; Legendre fails | azure2-exchange/SKILL.md#3 problems 3-4; references/ferdinand-checks.md#Cross-section reconstruction |

## rmfit-campaign

Legend: SK = rmfit-campaign/SKILL.md; LT = references/level-tests.md; PO = references/policy-and-objective.md;
CJ = references/cluster-jobs.md (all under rmfit-campaign/). "spilled to X" = handed to the X writer.

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| rmfit-campaign | 1-20 | frontmatter, intro, PYTHONPATH | SK intro (paths -> $RMFIT_HOME, $AZURE2_ROOT) |
| rmfit-campaign | 21-41 | 1. Read the state | SK §1 |
| rmfit-campaign | 42-61 | 2. Starting a campaign: prerequisites | SK §2 |
| rmfit-campaign | 62-161 | First fits: starting parameters | azure2-fit (assigned there); link in SK §2 item 5 |
| rmfit-campaign | 162-216 | 3. The loop (table, candidates.json, levelscan, relocate, additions) | SK §3; candidates/levelscan detail LT §1-2 |
| rmfit-campaign | 217-305 | 3b. Idle-time audits | SK §3b (quote and worked-example paths dropped) |
| rmfit-campaign | 306-348 | 3c. Width-floor recipe | azure2-resonance-search (assigned there); link in SK §3c |
| rmfit-campaign | 349-381 | 4. Judging results | SK §4 |
| rmfit-campaign | 382-398 | 5. What the human decides | SK §5 |
| rmfit-campaign | 399-416 | 6. Diagnostics and benchmarks | SK §6 (oracle script path dropped) |
| rmfit-campaign | 417-476 | 7. Known limits (covariances, shifts, conv, cost) | SK §7; covariances PO §4; shifts/Jacobian PO §5 |
| rmfit-campaign | 477-537 | 8. Gotchas | SK §8; pkill/quota CJ §3; Park PO §6 |
| rmfit-campaign | 538-539 | 9. Refinement log (header) | dropped: dated log removed; durable rules lifted (rows below) |
| rmfit-campaign | 540-542 | 2026-09-05 -- pyazr `engine_level_keys` counted inactive levels: the live 13C+a  | dropped: pyazr bake/inactive-level bug, fixed; current rule (purge inactive lines) SK §2 |
| rmfit-campaign | 543-544 | 2026-09-05 -- sharding: the convolved Cierjacks segment (17 of 23 s per evaluati | dropped: sharding history; current convolution rules PO §5 |
| rmfit-campaign | 545-546 | 2026-09-06 -- flips were credited with the incumbent's own continuation gain: ad | lifted: same-budget control polish; GN vs frozen screens -> SK §3 |
| rmfit-campaign | 547-549 | 2026-09-06 -- level tests: seed only the new level's keys; window scans sharded | dropped: per-campaign results; method superseded by LT §4 |
| rmfit-campaign | 550-556 | 2026-09-07 -- recovery benchmark: an 8-flip perturbation polished to a minimum s | lifted: benchmark tools -> SK §6; donor-basin rule -> SK §3 |
| rmfit-campaign | 557-559 | 2026-09-07 -- MANA 13C(a,a): Heil's shapes at ~0.67x scale with 1% errors; three | dropped: data-specific note; rule (angle-dependent norms = evaluation question) SK §5 |
| rmfit-campaign | 560-564 | 2026-09-07 -- round 3 of the signs benchmark (multi-level screen): 17 of 404 imp | lifted: drift outcome -> SK §4 |
| rmfit-campaign | 565-573 | 2026-09-07 -- Validation A PASSED on the signs-only benchmark: 777,349 -> 771,19 | lifted: judge by rho too; multi-level GN + polish quota -> SK §3, §4 |
| rmfit-campaign | 574-575 | 2026-09-07 -- `verify_candidate` from a heredoc script fails (`spawn` needs a re | lifted: spawn needs a real __main__ file -> CJ §1 |
| rmfit-campaign | 576-581 | 2026-09-07 -- oracle on the 728,879 incumbent: the *old* incumbent's 7/2+ p5 L1  | lifted: crossover from a second basin -> SK §3 (step 6b, rules) |
| rmfit-campaign | 582-587 | 2026-09-08 -- first crossover round on 13C+a: 33 donor patterns, 8.1 h, 6 accept | lifted: sign-converge -> perturb -> crossover -> polish -> SK §3 |
| rmfit-campaign | 588-593 | 2026-09-08 -- DeBoer's guidance on underestimated error bars (BandH stays in): a | lifted: robust objectives -> PO §3 |
| rmfit-campaign | 594-599 | 2026-09-08 -- the unscaled Sivia loss is degenerate when chi2/N >> 1 (see sectio | lifted: unscaled Sivia degenerate when chi2/N >> 1 -> PO §3 |
| rmfit-campaign | 600-608 | 2026-09-08 -- objective scan on 13C+a (`9-8-26_rmfit_robust/objective_scan.md`): | lifted: objective scan rule -> PO §3 |
| rmfit-campaign | 609-615 | 2026-09-08 -- freeing normalizations is a structure change: edit a copy of the b | lifted: freeing norms is a structure change -> SK §2 |
| rmfit-campaign | 616-625 | 2026-09-09 -- PLANNED (DeBoer): residual-driven level search, `rmfit levelscan`: | lifted: finite trial width physics -> LT §2; plan-file path dropped |
| rmfit-campaign | 626-635 | 2026-09-09 -- narrow levels drift: the 13C+a 9/2+ (Gamma_n ~0.4 keV, seen only b | spilled to azure2-fit (narrow-level drift, no gradient pull); relocate SK §3 |
| rmfit-campaign | 636-642 | 2026-09-10 -- robust campaign (global-scaled Sivia, free Heil norms) finished at | lifted: perturbation size for donor basins -> SK §3 |
| rmfit-campaign | 643-648 | 2026-09-10 -- with the error bars corrected by hand (BandH x2, MANA 3% systemati | lifted: plain chi2 once errors corrected; record error-bar changes -> SK §4, §5 |
| rmfit-campaign | 649-656 | 2026-09-10 -- `rmfit levelscan` built (rmfit/levelscan.py, test_levelscan_13n.py | lifted: engine level numbering, 16-digit energies -> LT §1, §3; add/remove_level owned by azure2-fit |
| rmfit-campaign | 657-660 | 2026-09-10 -- `levelscan.review_figure` + `Campaign.review_figure`: written by ` | lifted: review figures -> SK §3 |
| rmfit-campaign | 661-664 | 2026-09-11 -- the exported `.azr` of a campaign carries the fitted `<levels>` bu | lifted: exports carry nominal norms -> SK §3; -hold_jid -> CJ §1 |
| rmfit-campaign | 665-681 | 2026-09-11/12 -- PERTURBATION SIZE, both ends measured on the same sign-converge | lifted: perturbation size; stop arm when sign distance stalls -> SK §3 |
| rmfit-campaign | 682-695 | 2026-09-11 -- RETRACTION (DeBoer challenged the claim, and he was right): an ear | spilled to azure2-fit (MINUIT initial steps vs narrow-level basin) |
| rmfit-campaign | 696-706 | 2026-09-12 -- a model variant written into a SUBDIRECTORY needs the base model's | lifted: variant subdirectory needs data/ link -> PO §7; --opt=value -> CJ §1 |
| rmfit-campaign | 707-743 | 2026-09-12 -- THE TEMPLATE SCREEN AND THE NARROW-RESONANCE INTEGRATION FIX INTER | lifted: --grid 5,20 for screening -> LT §2; cost table spilled to azure2-targets |
| rmfit-campaign | 744-760 | 2026-09-13 -- A FREE ENERGY SHIFT IN THE .azr IS NOT A FREE PARAMETER IN A CAMPA | lifted: shifts_free policy; pre-flight through same code path -> SK §2 |
| rmfit-campaign | 761-773 | 2026-09-13 (later) -- the same shift test was void a SECOND time: `run_polish` h | lifted: exactly-0 fitted value is not a result; dead-columns line -> SK §4 |
| rmfit-campaign | 774-784 | 2026-09-13 -- OPEN ISSUE, single-session path only: `SingleSession.residuals(x,  | lifted: single-session Jacobian with free shifts is slow -> PO §5 |
| rmfit-campaign | 785-800 | 2026-09-13 -- A POSITIONAL .sav IS INVALID AFTER A LEVEL ADD/REMOVE. AZURE2's pa | lifted: positional .sav invalid after add/remove; init guard -> SK §2; csh eval -> CJ §1 |
| rmfit-campaign | 801-808 | 2026-09-13 -- EVERY STRUCTURE MOVE IS JUDGED AGAINST A SAME-BUDGET CONTROL POLIS | lifted: every move vs same-budget control -> SK §3 |
| rmfit-campaign | 809-819 | 2026-09-13 -- THE REDESIGNED LEVEL SEARCH IS WIRED INTO `rmfit levelscan` with t | lifted: --differential-window, --grid, --scales, GN slots -> LT §2 |
| rmfit-campaign | 820-829 | 2026-09-13 -- REDESIGNED SCREEN, FIRST END-TO-END RESULT (13N, the removed 3/2-  | dropped: fixture validation result (numbers only) |
| rmfit-campaign | 830-846 | 2026-09-14 -- REDESIGNED SCREEN, FIRST REAL USE (13C+a at Ex 9.862 MeV, job 1441 | lifted: frozen-screen limit; pattern count is the cost -> LT §2, §3 |
| rmfit-campaign | 847-861 | 2026-09-14 -- THE TEMPLATE SCREEN IS TWO-STAGE (`template_screen(stage_scales=(1 | lifted: two-stage screen, --max-signs -> LT §2 |
| rmfit-campaign | 862-878 | 2026-09-14 -- FIRST LEVEL ACCEPTED BY `structure add` ON 13C+a AFTER THE FROZEN  | lifted: frozen template is a lower bound; jpi null assigns J^pi -> LT §1, §3 |
| rmfit-campaign | 879-902 | 2026-09-14 -- THE TEMPLATE SCREEN'S AMPLITUDES WERE WRONG BY ~100x ON 13C+a; BOT | lifted: Wigner-route amplitudes; drop-in check; grep Denominator -> LT §2, §3 |
| rmfit-campaign | 903-916 | 2026-09-14 -- FIX VALIDATED ON THE REAL CASE (login node, 6 shards, grid 5/20, 1 | lifted: global gain = existence, local = J^pi -> LT §3 |
| rmfit-campaign | 917-934 | 2026-09-14 -- GATE C' (`13C+a/9-14-26_levelscan_986_fixed`, job 1441597): the re | lifted: same; candidate on top of same-J^pi level -> LT §3 |
| rmfit-campaign | 935-953 | 2026-09-14 -- THE SCREEN'S LOCAL TOP CANDIDATE, A 9/2+ 2 keV FROM AN EXISTING 9/ | lifted: PROVISIONAL lineshape hypothesis -> LT §5 |
| rmfit-campaign | 954-977 | 2026-09-14 -- FIRST LEVEL FOUND FROM RESIDUALS ALONE (`13C+a/9-14-26_levelscan_r | lifted: dead-channel rule, --merge-kev -> LT §2, §4 |
| rmfit-campaign | 978-986 | 2026-09-14 (DeBoer) -- BROAD STRUCTURE FIRST.  "Don't get stuck investigating th | lifted: broad structure first -> LT §6 |
| rmfit-campaign | 987-998 | 2026-09-14 (evening) -- STATE AND HANDOVER.  Campaign best: 366,565, 95 levels | dropped: session state/handover narrative; queue-check rule -> SK §1 |
| rmfit-campaign | 999-1008 | 2026-09-14 (evening) -- READ THE DRIVER'S campaign.log, NOT THE SGE STDOUT FILE. | lifted: campaign.log is the truth -> SK §1, CJ §4 |
| rmfit-campaign | 1009-1022 | 2026-09-14 (night) -- THE BROAD-ONLY SEARCH PAID OFF AT ONCE (`13C+a/9-14-26_lev | lifted: GN underestimates broad levels ~3.5x -> LT §3 |
| rmfit-campaign | 1023-1037 | 2026-09-14 (night) -- THE BROAD-ONLY SEARCH FOUND THE CAMPAIGN'S LARGEST GAIN | lifted: broad same-J^pi poles share strength; not every run is a level -> LT §3, §6 |
| rmfit-campaign | 1038-1047 | 2026-09-14 (late) -- SIGN ROUNDS NO LONGER PROPOSE FLIPS OF ~ZERO AMPLITUDES (se | lifted: sign_eps exclusion -> SK §3 |
| rmfit-campaign | 1048-1052 | 2026-09-15 -- A DETACHED (nohup) SUBMITTER CANNOT qsub ON THIS CLUSTER: "job rej | lifted: nohup submitter has no token; confirm job id -> CJ §1, §2 |
| rmfit-campaign | 1053-1066 | 2026-09-15 -- SECOND BROAD LEVEL FROM THE SAME SEARCH FAMILY (`13C+a/9-15-26_lev | lifted: window delta exaggerates; tau follows effective k; levels at Wigner limit -> LT §4, §6 |
| rmfit-campaign | 1067-1081 | 2026-09-16 -- A SIGN FLIP IN A CHANNEL THE DATA DO NOT SEE IS A NO-OP, AND SIZE  | lifted: flip_floor; judge by response -> SK §3 |
| rmfit-campaign | 1082-1100 | 2026-09-16 -- WHEN SINGLE-LEVEL ADDS STOP PAYING, MAP THE CHI^2 AND AUDIT THE CA | lifted: band map + theta2 audit after two null searches -> SK §3b |
| rmfit-campaign | 1101-1130 | 2026-09-16 -- THE MODEL'S BACKGROUND POLES WERE INERT, AND THAT IS WHY THE FIT W | lifted: count free keys per level -> SK §3b; pole seeding -> PO §2; spilled to azure2-fit |
| rmfit-campaign | 1131-1141 | 2026-09-16 -- level tests accepted broad absorbers on the 13C+a high-Ex window | lifted: unphysical gate (at-bound widths) -> LT §5 |
| rmfit-campaign | 1142-1165 | 2026-09-16 -- TWO SAME-J^pi LEVELS WITHIN ~1 keV MAKE A FIT UNWRITABLE.  The 13C | lifted: same-J^pi pair < 1 keV unwritable; bake guards -> SK §3; removal screens blind -> LT §6 |
| rmfit-campaign | 1166-1175 | 2026-10-06 — Park parametrization available (`"use_park": true`) | lifted: Park mode -> PO §6 |
| rmfit-campaign | 1176-1190 | 2026-09-17 — audit the pair table before fitting (13C+α) | lifted: pair-table audit -> SK §3b item 1 |
| rmfit-campaign | 1191-1203 | 2026-09-18 — seed dead channels in units of γ_W, never in keV (13C+α) | lifted: seed in units of gamma_W -> PO §2; spilled to azure2-fit |
| rmfit-campaign | 1204-1218 | 2026-09-18 — leave-one-dataset-out: totals are not enough (13C+α) | lifted: dataset-off reading rules -> SK §4 |
| rmfit-campaign | 1219-1231 | 2026-09-18 — a per-dataset table must reproduce a number you already know | lifted: check a new diagnostic against a known total -> SK §4 |
| rmfit-campaign | 1232-1243 | 2026-09-18 — after a session restart, qsub may lose its Kerberos cache | lifted: Kerberos cache; nohup + done-file; -hold_jid -> CJ §1, §2 |
| rmfit-campaign | 1244-1267 | 2026-09-18 — new data in a channel the model never used (13C+α, (α,α₁γ)) | lifted: seeding a new channel -> PO §2; sanity after structural edit -> SK §4; isotropy check spilled to azure2-check |
| rmfit-campaign | 1268-1278 | 2026-09-18 — campaign directories go directly under the reaction directory | lifted: campaign directory directly under reaction directory -> SK §2 |
| rmfit-campaign | 1279-1293 | 2026-09-19 — one bad point can hold a converged fit hostage (13C+α, MANA) | lifted: per-segment max point chi2; never edit a read data file; frozen first -> SK §3, §3b; spilled to azure2-check |
| rmfit-campaign | 1294-1306 | 2026-09-19 — levels added by the search had FIXED energies (method gap) | lifted: print new level parameters -> SK §4; levelFix version note -> LT §4 |
| rmfit-campaign | 1307-1325 | 2026-09-19 — run_level_add now fits the new level's energy; and restarts are not | lifted: released energy window, relocated verdict -> LT §4; restart test: no change of practice (dropped) |
| rmfit-campaign | 1326-1347 | 2026-09-19 — AZURE2's 1 keV level-merge tolerance is permanent; rmfit now guards | lifted: merge_pairs/close_pairs, bake key-set message -> SK §3, §4; shared binaries -> CJ §1; one-level template on existing -> LT §3; 1 keV rule owned by azure2-fit |
| rmfit-campaign | 1348-1368 | 2026-09-20 -- `add_level_variant` NEVER SEEDED OMITTED CHANNELS, whatever `Level | lifted: RMFIT_SEED_OMITTED -> LT §4; zero-strength channel with free norm -> PO §2 |
| rmfit-campaign | 1369-1378 | 2026-09-20 — after a login-node REBOOT the Kerberos caches are gone | lifted: reboot removes caches -> CJ §2 |
| rmfit-campaign | 1379-1395 | 2026-09-20 (later) -- "EFFECTIVE PARAMETERS" IN tau WHEN CHANNELS ARE SEEDED (De | lifted: effective k -> LT §4 |
| rmfit-campaign | 1396-1412 | 2026-09-20 — a template screen can manufacture a level (13C+α, 7/2⁺ at 10.42) | lifted: templates from carrying channels, theta2-capped -> LT §2, §3 |
| rmfit-campaign | 1413-1439 | 2026-09-20 (evening) — both fixes from the entry above are installed | lifted: collapse vs narrow verdicts; min_level_width_ev -> LT §5; PO §1 |
| rmfit-campaign | 1440-1452 | 2026-09-21 — first real use of the capped templates; a level add can be a SPLIT | lifted: splits signature -> LT §5; poll, not tail -F -> SK §1 |
| rmfit-campaign | 1453-1477 | 2026-09-21 (morning) — "splits" is now a verdict; the screen's "nan" explained | lifted: splits verdict, GN-slot nan -> LT §2, §5; pkill kills own shell -> CJ §3; install order -> SK §8 |
| rmfit-campaign | 1478-1532 | 2026-09-21 (night) — J^pi reassignment study; two-channel templates; when a flag | lifted: drop-in before believing a flag; pair family; local-ranking threshold; exit code; sequential adds -> LT §2-§4; temp file + mv -> SK §8 |
| rmfit-campaign | 1533-1546 | 2026-09-22 -- "BROAD" VERDICT (DeBoer).  Report figures of the 13C+a high-Ex win | lifted: broad verdict, max_level_width_ev -> LT §5, PO §1; range scans -> LT §2 |
| rmfit-campaign | 1547-1562 | 2026-09-23 -- THRESHOLD GUARD on free level energies.  In `13C+a/9-22-26_ntot_fi | lifted: threshold_margin_mev -> PO §1; bisection recipe dropped (one-off); divergence spilled to azure2-fit |
| rmfit-campaign | 1563-1596 | 2026-09-24 — the level search never estimated a width; `--gamma auto` now does ( | lifted: --gamma auto; absorber signature -> LT §2; PO §1 |
| rmfit-campaign | 1597-1617 | 2026-09-24 (later) — stripping broad levels to "rebuild narrow" fails; rank resi | lifted: broad first; --rank density, --seg-match, --ex-range -> LT §6, §2 |
| rmfit-campaign | 1618-1624 | 2026-09-24 -- FIGURES OF UPOS DATA NEED DATA-MODE CURVES.  A <segmentsTest> extr | spilled to azure2-plot; UPOS data-mode rule also SK §3b |
| rmfit-campaign | 1625-1643 | 2026-09-25 — data-definition fixes outrank structure; measure them frozen first | lifted: data-definition fixes measured frozen first -> SK §3 |
| rmfit-campaign | 1644-1661 | 2026-09-25 — a level pinned at the threshold guard: test both sides, same budget | lifted: threshold guard both sides -> LT §6; export overwrites seed -> SK §2 |
| rmfit-campaign | 1662-1678 | 2026-09-25 — zero amplitudes never move: audit closed channels below thresholds | lifted: ANC seeding below thresholds -> PO §2 |
| rmfit-campaign | 1679-1689 | 2026-09-25 — two copies of one conversation: name one owner | lifted: one owner when two sessions share a campaign -> SK §1 |
| rmfit-campaign | 1690-1699 | 2026-09-25 — width caps by energy window (`BoundsPolicy.theta2_cap_windows`) | lifted: theta2_cap_windows -> PO §1 |
| rmfit-campaign | 1700-1710 | 2026-09-25 — copying a campaign directory: repoint the ledger's structure paths | lifted: copying a campaign directory -> PO §7 |
| rmfit-campaign | 1711-1725 | 2026-09-27 — per-level theta^2 cap for added levels; the narrow loop | lifted: new_level_theta2_cap -> PO §1; narrow loop, re-screen offsets -> LT §5, §6; pkill -> CJ §3 |
| rmfit-campaign | 1726-1739 | 2026-09-28 — peak-weighted chi² and a per-channel width bound for narrow-level t | lifted: point_weights, width_bound_windows -> PO §1, §3; narrow method -> LT §6 |
| rmfit-campaign | 1740-1755 | 2026-09-28 (later) — three traps from the run-off comparison | lifted: _r<n> suffix -> SK §8; width-bound pins -> LT §5; __main__ guard -> CJ §1 |
| rmfit-campaign | 1756-1765 | 2026-09-29 — csh job scripts: never put JSON in a shell variable | lifted: csh JSON, @file.json, dry run -> CJ §1; SK §2 |
| rmfit-campaign | 1766-1797 | 2026-09-30 — four lessons from the 12C+α campaign (`12C+a_onefile/9-8-26_rmfit_c | lifted: bg_theta2_cap choice -> SK §2; verify_rtol, in-memory override -> SK §8; candidate above data range -> LT §4 |
| rmfit-campaign | 1798-1804 | 2026-09-30 — export refuses a model with a free width at exactly 0 (guard added) | lifted: bake nudges zero widths -> SK §3 |
| rmfit-campaign | 1805-1847 | 2026-09-30 — what the 13C+α window campaigns taught about windowed fits (both se | lifted: windowed fits, pole coverage, objective switching, copies -> PO §1, §3, §7; stop signal -> LT §5; data hand-over spilled to azure2-check |
| rmfit-campaign | 1848-1874 | 2026-10-01 — CRC default Python moved to 3.12.15; its sqlite3 is broken on compu | lifted (machine paths dropped): prove interpreter and imports on the compute node -> CJ §1 |
| rmfit-campaign | 1875-1896 | 2026-10-02 — idle-time reviews are part of the process (section 3b); four lesson | lifted: free unpenalized norms; windows act by energy; sharded conv; count patterns -> SK §3b, §3; PO §1, §5 |
| rmfit-campaign | 1897-1910 | 2026-10-03 — local release: judge a flip or a new level after its neighbours hav | lifted: local release -> LT §4 |
| rmfit-campaign | 1911-1925 | 2026-10-05 — a released level that collapses never answers "is a BROAD level wan | lifted: collapse answers narrow only; floors; per-channel bounds; excess attribution -> LT §4-§6; PO §1 |
| rmfit-campaign | 1926-1929 | 2026-10-05 later: the first floor run exposed a trap -- run_level_add's control  | lifted: bounds exempt level-off; broken control -> PO §1; SK §3 |
| rmfit-campaign | 1930-1939 | 2026-10-05 later still: any script that opens an rmfit evaluator (SingleSession, | lifted: __main__ guard -> CJ §1; carry fits by key -> SK §8 |
| rmfit-campaign | 1940-1951 | 2026-10-06 — data-only window fits are a screen, not a start vector | lifted: window fit is a screen, not a start vector -> LT §4 |
| rmfit-campaign | 1952-1958 | 2026-10-07 — a candidate is judged on the whole data set, not on the data set it | lifted: whole-data-set judgement; redistributes is a flag -> SK §4 |
| rmfit-campaign | 1959-1962 | 2026-10-07 note: AZURE2 commit 692e7f8 (10-06) adds a `parametrization` first ro | lifted: never read a .sav by position -> SK §2 |
| rmfit-campaign | 1963-1967 | 2026-10-07 — the Ex 9.4 recipe is now section 3c | dropped: pointer to section 3c (moved to azure2-resonance-search) |
| rmfit-campaign | 1968-1984 | 2026-10-10 — near-degenerate BROAD same-J^pi levels are separated at write time, | lifted: bake separates broad near-degenerate pairs -> SK §3, §4 |
| rmfit-campaign | 1985-1999 | 2026-10-10 — shard layout: composite sums are split, free-shift segments get the | lifted: shard layout env vars -> PO §5; per-shard timing -> SK §6 |
| rmfit-campaign | 2000-2005 | 2026-10-11 — a `+` in a pkill pattern filled the group quota | lifted: pkill + quota -> CJ §3 |
| rmfit-campaign | 2006-2028 | 2026-10-12 — the width floor belongs on the TOTAL width; a level add only has th | lifted: total width floor with RMFIT_SEED_OMITTED -> LT §4; PO §1 |
| rmfit-campaign | 2029-2050 | 2026-10-12 (b) — sharded convolution follows AZURE2 dad58dc; local ±mΓ table in  | lifted: shard-vs-single after conv change, guard -> PO §5; window table, seed reset -> LT §4 |
| rmfit-campaign | 2051-2060 | 2026-10-13 -- shared-package rule: every rmfit edit passes `tests/test_keys.py`  | lifted: shared-package rule -> SK §8; CJ §4 |
| rmfit-campaign | 2061-2083 | 2026-10-15 — windowed `<targetInt>` convolution in sharded fits; fixed-angle con | lifted: windowed sharded convolution, RMFIT_CONV_PPS=5 -> PO §5; init is single-session -> SK §7 |
| rmfit-campaign | 2084-2099 | 2026-10-16 — recall test of the level add; split rule ignored distant poles; the | lifted: recall test; split rule for distant poles; both seeds -> LT §4-§6 |
| rmfit-campaign | 2100-2113 | 2026-10-20 — centred parameter priors (`param_priors`), first use an ANC prior | lifted: param_priors -> PO §3 |
| rmfit-campaign | 2114-2142 | 2026-10-24 — model-scaled (`relvector`) covariance: refreeze between polish pass | lifted: never refreeze inside a least-squares run; oscillation superseded by reference model -> PO §4 (version note) |
| rmfit-campaign | 2143-2154 | 2026-10-25 — FIX installed: relvector covariance model fixed at a reference (sup | lifted: cov_reference.npz, chi2_self, --cov-reference -> PO §4 |
| rmfit-campaign | 2155-2177 | 2026-10-26 — relvector covariance with a FREE norm: neither reference convention | lifted: relvector + free norm biased -> PO §4; spilled to azure2-uncertainty |
| rmfit-campaign | 2178-2190 | 2026-10-26 — Jacobian-cost shard balancing (rmfit, installed 00:50, DeBoer appro | lifted: Jacobian-cost balancing -> PO §5 |
| rmfit-campaign | 2191-2217 | 2026-10-26 — t0 method for free norms on relvector-covariance rows; guard instal | lifted: t0 method and init guard -> PO §4, SK §2; spilled to azure2-uncertainty |
| rmfit-campaign | 2218-2231 | 2026-10-25 — physical seeds must be TINY, not 1 eV: the Brune transform can fail | lifted: tiny physical seeds -> PO §2; spilled to azure2-fit |
| rmfit-campaign | 2232-2249 | 2026-10-26 -- a frozen build needs AZURE2_ROOT, not just PYTHONPATH (13C+a sessi | lifted: AZURE2_ROOT + PYTHONPATH + assert -> CJ §1, SK intro |
| rmfit-campaign | 2250-2261 | 2026-10-30 -- level labels in models with several compound nuclei carry "@<nucle | lifted: @nucleus labels (version note) -> SK §8 |
| rmfit-campaign | 2262-2280 | 2026-10-30 -- Jacobian shift columns: held ones skipped, the rest forward only ( | lifted: held shift columns skipped (version note) -> PO §5; SK §7 |

## Developer notes (CLAUDE.md files)

| Old skill | Old lines | Old heading / first words | New location |
|---|---|---|---|
| azure2-eval | 3559-3572 | CI (GitHub Actions) rules | tests/CLAUDE.md#Portability (CI runs Linux, macOS and Windows/MinGW) |
| azure2-eval | 3573-3583 | A regression reference is not a correctness check | tests/CLAUDE.md#Which kind of test (paragraph after the table) |
| azure2-gui-review | 153-203 | Changing a dialog or a tab: check it headless before showing it | gui/CLAUDE.md#Changing a dialog or a tab: check it headless before showing it |
| azure2-eval | 639-654 (traps 1-3) | building an m-scheme truth for E1-M1 EC interference | tests/CLAUDE.md#Writing one (trap 3 also in azure2-check) |
