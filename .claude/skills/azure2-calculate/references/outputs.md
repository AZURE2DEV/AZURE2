# Output files

All outputs are centre-of-mass. Full column definitions: `docs/source/reference/output_files.rst`.

| File | Content |
|---|---|
| `AZUREOut_aa=<in>_R=<out>.out` | 9 columns: E_c.m., E_x, θ_c.m., fit σ, fit S, data σ, data σ error, data S, data S error (`TOTAL_CAPTURE` replaces `R=<out>` for summed capture; for A_y segments columns 4 and 6 hold A_y) |
| `AZUREOut_*.extrap` | mode 3: E_c.m., E_x, θ_c.m., σ, S |
| `AZUREOut_*.acoeff` | Legendre coefficients of angular-distribution test segments |
| `*.band` | covariance band (modes 2 and 3 with the band question answered yes) |
| `chiSquared.out` | χ²/N and norm per segment; last line total χ² |
| `scatteringLengths.out` | s-wave scattering lengths and their constraints |
| `param.par`, `param.sav` | initial and best-fit parameters (formal; reload as the external parameter file) |
| `parameters.out` | physical (observed) parameters |
| `param.fit` | autosave during a Minuit run (below) |
| `normalizations.out` | report of fitted norms (one ten-field row per varied norm) |
| `param.errors`, `covariance_matrix.out` | MINOS (mode 4) |
| `reactionrates.out`, `cross_section_grid.dat` | mode 5 |
| `samples.mcmc` | mode 6 |
| `intEC.dat`, `intEC.extrap` | external-capture integral caches (delete freely) |

## The normalization multiplies the data, not the fit

In `AZUREOut_*.out`, column 6 = d·n and column 7 = e·n, while column 4 = c is the bare calculation. The segment's
χ² = Σ((c − d·n)/(e·n))² = Σ((c/n − d)/e)². Lowering n raises χ² when c is too large; "scale the model down by n" is
the wrong formula.

The χ²-optimal norm of one segment at fixed R-matrix parameters, from one output file: with d = col6/n,
e = col7/n, c = col4,

    u* = Σ(c d / e²) / Σ(c² / e²),     n* = 1/u*.

A large χ²/N that no n* removes is a shape problem, not a normalization one: compare the fit column's energy or
angular shape with the data's. A data set's block appears among the active data sets of that entrance/exit pair in file order, separated by
blank lines.

## Data-set numbers

`chiSquared.out` and the `AZUREOut_*` blocks number the data sets by their position among the `dataset` records,
inactive ones included: the `key` that `AZURE2 --describe` prints and `Evaluation.data[i].segment_key`. Quote those,
not a count by eye; inside pyazr use the name (`d.name`).

## `param.fit`: an autosave of an evaluated point, not of the best point

Every 100th objective evaluation rewrites `param.fit` and the `AZUREOut_*`, `normalizations.out` and
`parameters.out` files at the point being evaluated at that call, which may be a line-search or gradient probe. Its
third column is a step size (10 % of the value), not an error. Builds before 2026-09-30 (expires: 2027-02-01) wrote
the unfitted input energies and widths into `param.fit` (only norms and shifts were current). Restarting a killed fit
from it: azure2-fit.

## `normalizations.out`

A report only: nothing reads it back; fitted norms travel in `param.sav`.
builds before 2026-10-27 (expires: 2027-02-01) wrote a malformed six-field row for every fixed norm: read the
ten-field rows only.

## Writing the files from pyazr

`s.write_outputs(x)` writes `AZUREOut_*` and `chiSquared.out` into the existing output directory, after a χ²
evaluation at `x` (a bare forward pass would leave zeros in `chiSquared.out`). It also writes the Park, abundance and
scattering-length penalty totals. It does not write `param.par` (initialisation) or `parameters.out` (end of a CLI
run); use `s.transform(x)` (one value per row of `s.parameters`) for those numbers.

