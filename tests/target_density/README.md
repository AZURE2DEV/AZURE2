# target_density — a fitted target density on a `<targetInt>` line

A target-integration line may carry `densityfit vary sigma [max]` and `densitygroup <label>`:
the areal density N (atoms/cm²) of a thin target becomes the fit parameter `density_<n>`
(n = the line's effect number, its place among the active `<targetInt>` lines), with the prior
((N − N0)/sigma)² about the line's value when sigma > 0 (atoms/cm², or `2.3%` of N0). Lines in one
`densitygroup` (one target seen by several detectors) share the parameter of the group's first line,
which carries the `densityfit` block; a member's own density is replaced by the group's (noted).

N moves the energy-loss window only, ΔE with ∫_{E−ΔE}^{E} dE′/ε = N, and the straggling with depth;
the yield is the integral over the window divided by the same N, so the result is the depth average
of the cross section, with no factor of N. The sub-point grids are laid out once, for the largest density
the fit is expected to reach (`max`, by default N0 + 3 sigma but at least 1.1 N0, or 1.5 N0 without a
prior), and the depth integral is cut at E0 − ΔE(N): the yield moves smoothly with N and no grid is
rebuilt. Beyond `max` the engine computes at `max` and warns once; MIGRAD has it as the parameter's limit.

Models (`make_models.py`, from `tests/13N`): 12C(p,p) at 84.3° (and 114.5°) lab from 0.3 to 0.9 MeV
across the 0.457 MeV resonance, every level fixed, a constant stopping cross section of 1e-20 MeV cm²/atom
(10 keV lab at 1e18).

| model | what it is |
|---|---|
| `true.azr` | both segments at N = 1.1e18 (a plain line: the data are made from it) |
| `fit.azr` | 84.3°, `densityfit 1 0` from 1e18 (no prior, grids to 1.5e18) |
| `fit_prior.azr` | the same with `densityfit 1 2%` |
| `fit_group.azr` | both segments in `densitygroup t` (the second line's own 2e18 is replaced) |
| `bad.azr` | `densityfit` on a Gaussian-only line |

`check.sh path/to/AZURE2`:

1. MIGRAD (`--use-gradient`), `--use-lm` and `--use-gsl-lm` recover N = 1.1e18 from 1e18 to 1e-4 (chi2 122
   → 1.5e-4: the truth's grid spans 1.1e18 exactly and the fit's 1.5e18, cut, and the two integrations differ
   by ~1e-5);
2. `AZURE_GRAD_CHECK`: d chi2/dN by finite differences agrees with the Jacobian path;
3. a 2 % prior pulls N between 1e18 and 1.1e18, `Total-Density-Chi-Squared` = ((N − 1e18)/2e16)²;
4. a density group has one parameter (`density_1`, no `density_2`), recovered from both segments;
5. `densityfit` without target integration is refused.

The independent physics check (the engine's yield against the average of the bare cross section over the
window computed from a 0.25 keV grid, at 0.9, 1.0 and 1.1 N0 and with an energy shift; a fitted density
against a plain line at the same N; the group; penalties; the Jacobian column) is
`tests/pyazr/target_density_test.py`. The GUI side (the dialog's Fitted Target Density box, its refusals,
the read/write of both blocks, a group member) is `tests/gui/target_int_tab_test.cpp`, and the Fitting tab's
Target Densities rows `tests/gui/density_fitting_test.cpp`.
