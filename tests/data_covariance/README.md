# data_covariance — correlated data uncertainties (`covariance <file>`)

A data segment line may end with `covariance <file>`. All segments naming one file share
its covariance matrix C, and their points enter chi2 as r^T M^-1 r, with r the usual
normalized residuals (fit - n d)/(n sigma), M = D^-1 C D^-1 and D = diag(sigma). The file
lists its data files and their data-row counts, then terms that add up to C: `block`,
`diagonal`, `vector`, `relvector` (model-scaled, frozen at the first evaluation) and
`matrix`. The format is in `docs/source/reference/data_formats.rst` and the code in
`include/DataCovariance.h`. AZURE2 never forms C. It whitens with the block-diagonal part
B, factored by connected components, plus a rank-k correction for the vector and
relvector terms. Introduced 2026-10-21 for the MANA 10B(a,n) and 10B(a,p/d/a) data
(10B+a/10-19-26_MANA_Claude): N = 1845 and 14049 rows, about a hundred run blocks, and
per-detector efficiencies as relvector terms.

Model: `tests/angle_mapping` (p+12C and p+13C, thirteen segments on one energy grid, data
taken from there). `make_models.py` writes the four models and three covariance files (all
committed). Every model splits the two-angle file `p12c_120_150.dat` into two segments by
angle window and cuts segment 2 at 2.0 MeV. The covariance covers eight segments: a free
norm, A_y, a shift group with a free shift, and a free sqrt(E) shift. Its header also lists
`p12c_165.dat`, which no covariance segment fits.

| model | covariance |
|---|---|
| `plain.azr` | none |
| `diag.azr` | `cov_diag.txt`: diag(e_i^2), so the chi2 must equal plain |
| `dense.azr` | `cov_dense.txt`: 0.5 diag(e^2), one block per beam energy across all files (0.4 e_i e_j off the diagonal), a 2 % absolute vector per file, one smooth global vector |
| `rel.azr` | `cov_rel.txt`: as dense, with the per-file 2 % terms as `relvector` |
| `fit.azr` | `cov_dense.txt`, with every level energy and width fixed (for the fitter comparison) |

`check.sh`:

- (a) diag = plain, segment by segment (481 of 549 rows fitted).
- (b, c) the engine's chi2 in `covariance.out` equals `verify.py`, an independent numpy
  v^T C^-1 v built from the output files, its own parser of the covariance file and, for
  rel, the model in `covarianceModel.out`. The agreement is 2e-11, which is the precision of
  the output files. The per-segment shares sum r_i (M^-1 r)_i match `chiSquared.out`.
- (d) `AZURE_GRAD_CHECK`: the adjoint gradient and the whitened residual-Jacobian gradient
  agree with central differences for all 21 free parameters, dense and rel, to 3e-7 of the
  largest. The parameters are level energies and widths, norms, a shift-group shift, a
  sqrt(E) shift, a shift function and abundances.
- (e) MIGRAD (`--use-gradient`), `--use-lm` and `--use-gsl-lm` fit `fit.azr` (dense with every
  level fixed: norms, shifts and abundances free) to the same chi2, within 1e-4 relative
  (3e-8 in practice). With the levels free this synthetic model has several minima, and the
  fitters reach different ones depending on the platform: the first CI run on Linux found
  3036 with MIGRAD against 3078 with LM.
- (f) these broken inputs stop the run with an ERROR:
  - a wrong row count in the header
  - a data file missing from the header
  - overlapping angle windows on one file
  - a matrix that is not positive definite
  - an asymmetric block
  - a covariance block on a test segment

About 10 s with 4 threads. `tests/pyazr/data_covariance_test.py` covers the same models
through pyazr:

- the keyword is read and written
- sum w^2 = chi2
- the gradient matches central differences and 2 J^T w
- `refreeze_covariance()` works
