# The old project file through its history

The forms an `.azr` file has had, read from the history of every ref of the repository (`origin/qt5`, `origin/main`,
`origin/master`, `origin/dev` and their branches; read-only `git log`, `git show <rev>:<path>`), and one fixture per
form.  The engine's readers (`include/SegLine.h`, `ExtrapLine.h`, `NucLine.h`, `src/TargetEffect.cpp`, `src/Config.cpp`,
`src/CNuc.cpp`) and the setup program's writers (`gui/src/AZURESetup.cpp`, `SegmentsTab.cpp`, `TargetIntTab.cpp`,
`LevelsTab.cpp`, `PairsTab.cpp`, `FittingTab.cpp`, `MCMCTab.cpp`) decide each form.  The repository's history starts
on 2020-07-03 (`6ecaee22`, a mature code base); the first project file in it is `test/13N.azr` of 2026-04-22.

`make_fixtures.py` writes the constructed forms (the 13N model of `tests/convert/fixtures/tests/13N`, two data files of
`tests/13N/data`) as each era's writer wrote them; `revisions/` holds project files of the repository as they were at
the revision in their name; `glued_*` are the two forms of the 2026-09 writer described below (the `<levels>` one is `glued_levels_2026_09.old`, not
`.azr`: `io_level_row_text` reads every `<levels>` row of the tree's `.azr` files as a valid row).  The ctest
`convert_history` (`convert/tests/history_test.cpp`) converts each and checks the meaning against the same project in
a later era's form; `convert_cli_locale` (`tests/convert/locale_check.sh`) and `gui_casual` (g) convert three of them
under a decimal-comma locale.

## Timeline

| Date, commit | What changed in the file |
|---|---|
| before 2020 (no source in the repository) | `<levels>` rows of 30 columns (no `ecMultMask`); external capture in `<externalCapture>` lines `isActive exitPairKey minJ maxJ multMask`.  The 2020 engine refused such rows; the 2020 setup program imported them (`LevelsTab::readNuclearFile`: a missing 31st column read as 0; `PairsTab::parseOldECSection`: the mask set on the exit pair, `isActive`, `minJ`, `maxJ` not used). |
| 2020-07-03 `6ecaee22` (the state of `origin/qt5`, `origin/master`) | `<config>`: 11 positional lines (A-matrix, output and checks directories as full paths, eight check words); `<levels>`: 31 columns; `<segmentsData>`: `isActive entr exit minE maxE minA maxA isDiff [J l] norm varyNorm normError file`, the file the rest of the line (blanks allowed, unquoted); `<segmentsTest>`: `isActive entr exit minE maxE eStep minA maxA aStep isDiff [J l / order]`; `<targetInt>`: `isActive "segments" points conv sigma integration density "eq" nParams p.. isQ nQ q..`; `<lastRun>`: the setup program's runtime choices (read by no engine).  Numbers in 15-character columns, left-aligned (the config writer set the alignment of the shared stream). |
| 2025-01-26 `17281967` | `<targetInt>`: + `isConvCoefficients nConv c0..` (sigma(E_lab) = c0 + c1 E + c2 E^2 + ...). |
| 2025-01-30 `1c01a007` | `<targetInt>`: + a quoted equation of sigma between the flag and the count (`isConv "eq" nConv c..`). |
| 2025-08-01 `dcfdb344` | `<segmentsData>`: + `energyShift shiftError varyShift` before the file; `<parameterSettings>` section. |
| 2025-08-14 `1225e4b1` | `<mcmc>` section. |
| 2025-08-30 `5697d610` | "advanced" (composite) segments: after the file `1 op n e1 x1 e2 x2 ...` (op 0 sum, 1 ratio; two numbers per component) or ` 0`; test lines the same after their type. |
| 2025-12-22 `de248934` | components with an angle: `e x angle` per component (-999: the segment's own angles). |
| 2026-01-22 `dd6e7a26` | `<targetInt>`: + straggling flag and coefficient. |
| 2026-03-09 `2c406c7a` | `<segmentsData>`: + the secondary-decay (UPOS) field after the composite block: ` 0` or ` 1 L Jf delta`. |
| 2026-03-28 `14753ded` | `<targetInt>`: + resonance-width multiplier and points per width. |
| 2026-05-21 `9b5084d8` | components with a scaling: `-1 n e x angle scaling ...` (the -1 marker only when a component has a scaling; other lines of the file keep three numbers). |
| 2026-07-09 `ee092f9e`, 07-12 `5832bee5` | `<levels>`: + binding energy (32nd column), + `gammaIsRWA` (33rd). |
| 2026-09-17 `922c0537` .. 09-28 `330dc1eb` | isDiff 8 (P(theta) dsigma/dOmega), then removed. |
| 2026-09-25 `d2d30a9a` .. 10-07 `3e7bd976` | round-trip numbers in the fixed columns; a value as wide as its column ran into its neighbour with no blank (`2887d991` fixed `<levels>` and `<segmentsData>` on 09-29, `3e7bd976` the test and target lines on 10-07). |
| 2026-10-12 `a38494e8` .. 10-13 `a60577f6` | dev's `mixture` block on segment lines, then reverted. |
| 2026-10-18 `8f80abde` .. 11-11 | dev's general spin observables written as type 9 with an operator (`9 b:y`); type 10 since, type 9 being Legendre coefficients (`13bed879`, 10-27). |
| 2026-08 .. 11 (`dev`, `main`) | keyword blocks on segment and target lines (sqrtshift, shiftgroup, normgroup, shiftfn, sample, covariance, transmission, ...), `<potential>`, `<thm>`, `<abundances>`, `<isospin>`, `<scatteringLengths>`, `<channelRadii>`; then format 2 (`format=2`) and format 3. |

## Forms and what the converter does

| Form | Fixture | Converter |
|---|---|---|
| before 2020: 30-column `<levels>` and `<externalCapture>` | `pre2020_external_capture` | converts: the mask folded into the 31st column as the 2020 import did (notes; `minJ`, `maxJ` not used, as the import); records and rows equal to `gui_2020`'s, and its calculation's output files byte for byte |
| before 2020: an inactive `<externalCapture>` line; two masks for one pair | (in the test) | refused: the import applied an inactive line, the engine that wrote it is not known to have |
| before 2020: 30-column rows alone | `pre2020_30_columns` | converts: ecMultMask 0 (note) |
| 2020 | `gui_2020` | converts (unchanged); also with CRLF, tabs, blank lines, empty sections, full paths with blanks or backslashes, no `<lastRun>` |
| 2020: `<lastRun>` not closed at the end of the file | (in the test) | converts: closed (note); no engine read that section |
| 2020: unquoted data file name with blanks to the end of the line | (in the test) | converts: read as that name, as the 2020 and 2025-08 engines read the rest of the line (note) |
| 2025-01-26: target line without an equation token | `target_2025_01_no_equation`, `target_2025_01_coefficients` | converts: the polynomial written as the equation (note); records equal to the same line written in the 2025-01-30 form |
| 2025-01-30 .. 2026-03 target lines | (every era fixture) | converts (unchanged) |
| 2025-08: energy shift | `energy_shift_2025_08` | converts (unchanged) |
| 2025-08-30 .. 2026-05: composite lines, two, three or four numbers per component, UPOS, scaling | `components_2025_08`, `components_angle_2025_12`, `secondary_decay_2026_03`, `components_scaling_2026_05` | converts (unchanged): the layout decided as the setup program read it (three numbers when 3n follow, else two; four after -1) |
| 2026-09: a number run into its neighbour | `glued_data_2026_09`, `glued_levels_2026_09` | refused (two meanings): data, test and target lines with the split proposed as a reading; a `<levels>` row naming the token.  Before: a glued `<levels>` row was converted with its columns shifted |
| 2026-09-17: isDiff 8 | `revisions/polarization_product_2026_09_17` | refused: removed by `330dc1eb`, no format-3 quantity for the product (the reason now says so) |
| 2026-10-12: `mixture` | `revisions/mix_data_2026_10_12` | refused: a one-day construct of dev with no format-3 record (the reason now says so) |
| 2026-10-18: `9 <operator>` | `revisions/spin_observables_reaction_2026_10_18` | converts as type 10 (note); records equal to the same lines written as type 10 |
| a target line of inactive lines only | `revisions/12C_pp_ay_2026_08_04` | converts: the nucleus of the lines it lists (note); was a question |
| 2026-04-22 `test/13N.azr` | `revisions/13N_2026_04_22` | converts (unchanged) |
