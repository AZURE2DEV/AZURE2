# A covariance matrix for every data set in a fit

Every data set in an R-matrix fit, own measurement or literature set, carries a covariance matrix C (units of the
data file squared) in AZURE2's `covariance <file>` input format, built before the fit. This file gives the inventory
step, what goes into C, the recipes for own and literature data, the checks and the known traps. The engine side
(the `covariance` keyword, the terms `diagonal`, `block`, `vector`, `relvector`, the χ² it enters) is in
`docs/source/user_guide/chi_squared.rst` and `docs/source/reference/data_formats.rst`. Campaign handling of
`relvector` terms (one fixed reference model) is in `rmfit-campaign`.

## 1. Inventory first

List every data file the fit uses and the covariance file that covers it, as a table in the fit directory's README.
A data analysis can end with every neutron and silicon set covered and the γ-ray sets carrying only
stat ⊕ 3 % errors and normgroups; the table makes such gaps visible.

```bash
# covariance files and the data files they cover
for f in *cov*.txt; do echo "$f: $(grep '^file' $f | awk '{printf "%s(%s) ", $2, $3}')"; done
# active segments of a model WITHOUT a covariance keyword
awk '/<segmentsData>/{f=1;next}/<\/segmentsData>/{f=0}f && $1=="1" && $0!~/ covariance /' model.azr
```

## 2. What goes into C, and what stays out

C holds the uncertainties that move the points of a set relative to each other, or to other sets, in a known
pattern. Two things stay out:

- **The overall normalization of a set** (charge, target thickness, absolute efficiency, the paper's quoted scale
  error). It enters the `.azr` as a free normalization with a prior equal to the stated scale error, or as a
  `normgroup` shared by the segments that share it (all sets of one run: one normgroup, normErr = charge ⊕ areal
  density ⊕ enrichment, e.g. 3.0 %). Putting it in C as well counts it twice. If it must be in C, use a model-scaled
  `relvector`, never a data-scaled `vector` (Peelle's puzzle: the fit is pulled low).
- **Energy-scale uncertainties** (beam-energy calibration, effective-energy offsets). They are uncertainties in x:
  give the segment an energy shift with a prior (`energyShift ... varyEnergyShift`, `sqrtshift`, `shiftgroup`).

Everything else that is documented or physically certain goes in:

| term | AZURE2 term | source |
|---|---|---|
| statistics, independent points | `diagonal` | counting or quoted point errors |
| statistics of quantities fitted together in one spectrum | `block` per spectrum | the spectrum fit's own covariance (Fisher matrix), scaled to the data units |
| feeding or cascade subtraction within one run (γ-ray level cross sections) | `block` per run | J diag(err_line²) Jᵀ, J = ∂(levels)/∂(lines) |
| per-detector solid angle or efficiency, common to all energies and lines of that detector | `relvector` per detector | calibration (gold runs, sources, pair scatter) |
| efficiency shape against energy (one curve shared by all points) | `relvector` per shape parameter | change of the curve's parameter by 1σ, as a fraction of each point |
| analysis choices (line shapes, background models, K-factors, light gain) | `relvector` per variant | half the difference between variant and production (two-sided variants: half the spread), as a fraction shift / max(value, 3σ_stat) |
| per-point scatter beyond statistics (e.g. detector-pair disagreement) | second `diagonal` | per-point term s × value per (line, angle) |

Rows of C are the data rows of the `_cov.dat` named in the `file <path> <rows>` header, not the points left after
cuts: energy windows and angle splits on segment lines need no new matrix. Every fitted row needs its own diagonal or
block variance (the engine factors diagonal + blocks and adds the rank-k terms).

## 3. Own data

- Build C from the same production tables the data files are written from, in the same row order, and write the
  data file's error column as √C_ii (`_cov.dat`). Keep a diagonal-only version of the data file for fits without the
  `covariance` keyword.
- Spectrum fits with deviance/ndf > 1: multiply each fit's Fisher block by max(dev/ndf, 1) (error inflation); the
  correlations inside a block are unchanged.
- Differential and angle-integrated versions of one measurement are not both fitted: the integrals are linear
  combinations of the same points. Give each its own covariance file and use one of them (or the integrals only where
  the angular data are not used).
- An alternate file set for a systematic that cannot be a covariance term (e.g. every silicon angle shifted by the
  effective-angle offset, 3.7°) gets its own matrix of the same structure.
- After any change to the production tables, rebuild C and rerun the read test; keep the previous outputs in a
  `pre_<tag>/` directory.
- A builder script per detector class (neutron, charged-particle, γ-ray, angle-integrated) with a docstring listing
  its terms keeps the matrices reproducible.

## 4. Literature data

One covariance file per experiment (different experiments uncorrelated), covering all of its data files.

1. **Read the error budget.** The EXFOR entry (`pyazr.nds.fetch_exfor(..., source="local")`, or the raw `.x4`):
   `ERR-ANALYS` (meaning of each error column), `ERR-T` (total), `ERR-S` (statistical), `ERR-SYS`, `ERR-1..n`
   (partial systematics), `MONIT` (normalized to another measurement?), the free text, then the paper. Write down per
   set: point-to-point part, correlated parts, scale part, energy-scale part.
2. **Split.** If the quoted point errors include the scale error (an `ERR-T` total), remove the scale in quadrature
   from the diagonal and put it in the norm prior. If only a total is known, use it as the diagonal and say so in the
   README. Never invent correlated terms: add one only when the paper documents it or the method makes it certain
   (several angles measured with separate detectors: a per-angle efficiency term of the size the paper quotes).
3. **Derived quantities.** A 4π value made from an angular distribution, or an excitation function and the angular
   distributions it was integrated from: use one, not both. Keep the derived one only when its source is not in the
   fit.
4. **Linked normalizations.** A set normalized to another set shares that set's `normgroup`; the two are not
   independent scale checks.
5. **Thick-target yields and activation data** need a `<targetInt>` entry (`azure2-targets`); note the stopping
   powers and the target composition next to the segment line.
6. **Deliverables per experiment:** `<tag>_cov.dat` (lab frame, 4 columns, error = √C_ii), `<tag>_cov_rel.txt`,
   ready-to-paste segment lines (norm with its prior, normgroup label, energy-shift prior, energy and angle windows,
   the covariance file name), a correlation plot, and a README table (rows, diagonal source, correlated terms, norm
   prior, shift prior, assumptions).

## 5. Checks before the fit

- **Positive definite** (Cholesky succeeds); report the smallest and largest eigenvalue (`write_covariance.py
  --check` does both).
- **√diag C equals the `_cov.dat` error column** to 1e-6 relative.
- **Correlation plot** per file, inspected: a detector term shows as a block over that detector's rows.
- **AZURE2 read test** (mode 1, any model with the needed exit pairs; run it on a compute node, not a login node):
  the log says `Covariance <file>: N of N rows fitted` for every file (fewer only where the test model lacks an exit
  pair, which must be stated), and no warnings.
- `WARNING: k point(s) of <file> have an error column more than 1 % from sqrt(C_ii)` means the data file and the
  covariance file come from different builds of the tables. Refresh both together, never one.
- `scripts/write_covariance.py --check COV --base DIR` compares √C_ii (with the relvectors evaluated on the data)
  with the error column, tests positive definiteness of B and C and the row counts; confirm separately that the data
  are unchanged against a backup.

## 6. Traps

- **An error column can hide a non-statistical per-detector problem.** Detectors that disagree by 1.4× on an
  isotropic line point to wrong per-detector efficiency anchors. Build C after such checks, not before.
- **A frame or Jacobian error upstream propagates into every derived table and its matrix** (Legendre fits
  multiplied by dΩ_c.m./dΩ_lab instead of dividing). Check a derived angle-integrated set against an independent one
  before building its matrix.
- **A `relvector` fraction (variant − nominal)/nominal computed from the data is meaningless where the value is
  consistent with zero or negative** (fractions up to 70 % with flipped sign on a weak channel). Use the data fraction
  only where the point is above 3σ; elsewhere the fraction the term gives an ideal (e.g. isotropic) value with the
  same detector weights, or shift / max(value, 3σ_stat).
- **`relvector` terms are scaled by a model.** Outside rmfit the engine freezes that model at the first evaluation
  (written to `covarianceModel.out`; pyazr `s.refreeze_covariance()` refreezes it). In a campaign it is one fixed
  reference model (`rmfit-campaign`).
- **Every data file a covariance file names must be present in the fit's `data/`**, not only those the segments use:
  the engine and rmfit match rows by counting each named file's records (`file <path> <rows>`). A missing one stops
  the run with `FileNotFoundError`. A model builder copies the files listed in the covariance headers.
- **Analysis-variant terms must be model-scaled `relvector`s, not data-scaled `vector`s.** A variation with the same
  sign on most rows (a room-background term of −5 % on nearly every row) is a hidden normalization: a residual along
  it costs only (offset/0.05)², about 100 for a 50 % offset over 968 rows, so the set's scale is almost free inside C;
  and because the term is computed from the data it carries Peelle's bias. Symptom: the set's norm moves against
  fit/data (a fit at 0.64 × the data while another set's norm went to 1.32; after conversion to `relvector` the norm
  scan's minimum moved to 0.85). Diagnose any set whose norm and fit/data disagree with a scan of plain against
  covariance χ² as a function of that set's norm alone.
- **A `relvector` covariance plus a free normalization on the same rows is biased** under either freeze of the
  scaling model (the engine's first-evaluation freeze or a fixed reference model), for CLI fits as well as rmfit.
  The `relvector` terms together span the data vector, so scaling the data towards 0 or ∞ costs a bounded amount
  (~1/rel² per term), cheaper than the shape misfit when χ²/N is large; norms then run to their bounds (one ran from
  0.63 up to 3.08 under one convention and down to its 0.2 floor under the other). Before trusting such a fit, check
  the norm against its prior and scan it with everything else fixed. The remedy is the t0 method (Ball et al.,
  JHEP 05 (2010) 075, §6): fix the norm at its nominal value (`varyNorm 0`); add one `relvector` term per former
  norm parameter, (prior) × 1 over the rows that norm scaled; build those terms from a fixed prediction t0 of a
  previous fit, held constant during the fit so χ² stays quadratic; refresh t0 and refit (one iteration usually
  suffices). A covariance built from measured values is biased low (Peelle); one built from the fitted values is
  biased the other way. Neutron standards use the same idea (Chiba–Smith).
