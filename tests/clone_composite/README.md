# clone_composite

The working copies a CLI fit (MIGRAD, LM, GSL) and MCMC evaluate chi2 on are made by `EData::Clone`. Up to
2026-11-05 `Clone` repointed only the points `EDataIterator` visits, and the iterator skips the segments that hold
components (SUM, RATIO, total capture, natural sums with components, a yield's hidden total term). In every working
copy, those segments' points:

- kept the ORIGINAL data as their parent, so they read the original's target effects. A fitted areal density
  (`densityfit`) or a fitted Doppler temperature never moved for them inside a fit. Electron screening was not
  affected: `EData::ApplyScreening` pushes U_e into the copy's own points.
- kept the original's list of mapped component points, next to the copy's own. Every copy then computed, and wrote,
  the original's component points too: twice the component work, and unsynchronised writes between threads.

Model (`make_models.py`): `tests/target_density` (12C(p,p) at 84.3°, thin target, fitted density) with segment 3 made
a SUM of its own term and an identical elastic component, the only active segment.

| check | what |
|---|---|
| (a) | MIGRAD (`--use-gradient`) recovers density_1 = 1.1e18 from 1e18 to 1e-4 with chi2 < 1e-2 (data made at 1.1e18). Before the fix the fit stopped at 1.018e18 with chi2 82. |
| (b) | The gradient check (`AZURE_GRAD_CHECK`) agrees for the density. It passes on the old build as well: the check does not evaluate on the working copies. (a) is the regression test. |

    ./tests/clone_composite/check.sh path/to/AZURE2      # ~5 s
