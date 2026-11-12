# beam_spread_ay — a beam spread folded into analyzing-power data

Named `target_effect_ranges` until 2026-11-12 (when it tested the fold windows, since deleted).

The 7Li_p_ay analyzing-power model with a 30 keV (lab) Gaussian beam spread on both Deineko data sets, which name one
beam record (`beam_targetInt1`, as the old `<targetInt>` line "1,2" did).  Each A_y point is the ratio of two folds,
<sigma A_y> / <sigma>, on the shared mesh of the compound nucleus.  `tests/target_effect_shared` runs this project
against a copy with one beam record per data set.

Until 2026-11-12 the project restricted the spread to the lab-energy window 1.95-2.55 MeV, with a 120 keV blend at
the window's edges and an automatic tolerance of 0.002 (format 3: `fold_windows   1.95-2.55 MeV blend 0.12 MeV
auto_tolerance 0.002`; the old line's trailing `"1.95-2.55" 0.12 0.002`).  Those three options were deleted with the
shared mesh (decision D3, `docs/upgrade/energy_grid.md` section 9): every point of a folded data set is folded, and
the field was removed from this file by hand.  The project now pins the fold everywhere.

Re-pinned 2026-11-12: chi2 3838.66 -> 3760.95 (segments 2313.79 -> 2245.66, 1524.87 -> 1515.29; the engine on the
shared mesh gave 3837.13 with the window).  The moved numbers are the points the windows left unfolded or blended: outside the window the largest change of a model A_y is 0.046
(1.5 data sigma, 15 % relative, E_lab 1.876 MeV), inside it 0.004 (0.13 sigma, the blend and the tolerance); 31 of the
72 points keep their bits.  Independent of the new engine: the engine before the deletion, run on this file and on the
old file with the windows widened to 0-1000 MeV (no blend, no tolerance), gives both output files bit for bit.  The
quadrature: with every interval refined 2x and 4x (`AZURE_REFINE_QUADRATURE`) the A_y move by at most 3.1e-4 (0.01
data sigma) and chi2 by 2.7e-4 relative, inside the default tolerance of `tests/run_tests.sh` (1e-3), so the project
has no tolerance file of its own (the 5e-3 it had absorbed a point flipping across the automatic tolerance).

Earlier pins (with the window): 2026-09-26 3861.43 -> 3861.51 (shift-function fix b6cc41b); 2026-09-28 3861.51 ->
3861.54 (dS/dE step) and 3861.54 -> 3838.66 (each listed segment converts the shared sigma lab -> c.m. once).
