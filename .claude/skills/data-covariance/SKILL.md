---
name: data-covariance
description: Build a covariance matrix for EVERY data set that goes into an AZURE2 R-matrix fit -- our own measurements and literature/EXFOR sets alike -- in AZURE2's `covariance <file>` input format, before the fit is built. Covers the inventory step (which segments still lack a matrix), what goes into C and what stays out (overall normalizations, energy scales), recipes for own data (fit covariances, analysis variants, detector terms, error inflation, angle-integrated sets) and for literature sets (EXFOR ERR-* fields, one file per experiment, derived quantities, linked normalizations), the checks (positive definite, sqrt(diag) vs the error column, AZURE2 read test) and the traps seen so far. Use whenever data are being prepared for an R-matrix fit, a new data set is added to a fit, an EXFOR set is converted for AZURE2, or the task mentions covariance, correlated uncertainties, Peelle's puzzle, relvector, normgroup priors or "error budget for the fit".
---

# A covariance matrix for every data set in a fit

DeBoer, 2026-10-25 (10B+a): "Let's make a new skill entry where we try to create a covariance matrix for all data
sets we are going to put into our R-matrix fits."

**Standing rule: no data set enters an R-matrix fit without a covariance file**, ours or from the literature. A
segment without one is either given one before the fit is built, or the readme says why not (e.g. a set used
only for a plot, or an alternative view of points already in the fit). This file says how. The engine side
(the `covariance` keyword, its terms, `relvector`, `covariance.out`) is in the `azure2-eval` skill,
`.azr` anatomy, "`covariance <file>`"; read that bullet first. Campaign handling of `relvector` terms (one
fixed reference model, `cov_reference.npz`) is in `rmfit-campaign`, log entries of 2026-10-24/25.

## 1. Inventory first: which segments have a matrix?

Before building a fit, list every data file the fit will use and the covariance file that covers it. Do it as a
table in the fit directory's readme. It is easy to finish a data analysis believing everything is covered and
find a gap only when asked: on 10B+a (2026-10-25) every neutron and silicon set had a matrix but the June 2025
HPGe gamma-ray level files had only stat (+) 3 % errors and normgroups.

```bash
# covariance files and the data files they cover
for f in *cov*.txt; do echo "$f: $(grep '^file' $f | awk '{printf "%s(%s) ", $2, $3}')"; done
# active segments of a model WITHOUT a covariance keyword
awk '/<segmentsData>/{f=1;next}/<\/segmentsData>/{f=0}f && $1=="1" && $0!~/ covariance /' model.azr
```

## 2. What goes into C, and what stays out

C holds the uncertainties that make the points of a set move relative to each other or to other sets in a known
pattern. Two things stay OUT of C:

- **The overall normalization of a set** (charge, target thickness, absolute efficiency, the paper's quoted
  scale error). It goes into the `.azr` as a free normalization with a prior equal to the stated scale error,
  or as a `normgroup` shared by the segments that share it (e.g. all sets of one run: our files use normgroup
  `ALL`, normErr 3.0 for charge (+) areal density (+) enrichment). Putting it in C as well counts it twice.
  If it must be in C, use a model-scaled `relvector`, never a data-scaled `vector` (that is Peelle's puzzle:
  the fit is pulled low).
- **Energy-scale uncertainties** (beam-energy calibration, effective-energy offsets). They are x-uncertainties:
  give the segment an energy shift with a prior (`energyShift ... varyEnergyShift`, `sqrtshift`, `shiftgroup`)
  instead.

Everything else that is documented or physically certain goes in:

| term | AZURE2 term | where it comes from |
|---|---|---|
| statistics, independent points | `diagonal` | counting / quoted point errors |
| statistics of quantities fitted together in one spectrum | `block` per spectrum | the spectrum fit's own covariance (Fisher, `--covout`), scaled to the data units |
| feeding / cascade subtraction within one run (gamma-ray level cross sections) | `block` per run | J diag(err_line²) Jᵀ, J = d(levels)/d(lines) |
| per-detector solid angle or efficiency, common to all energies and lines of that detector | `relvector` per detector | calibration (gold runs, sources, pair scatter) |
| efficiency shape vs energy (one curve shared by all points) | `relvector`, one per shape parameter | change of the curve's parameter by 1 sigma, as a fraction of each point |
| analysis choices (line shapes, background models, K-factors, light gain, ...) | `vector`, one per variant | HALF the difference between the variant and production (two-sided variants: half the spread) |
| uncorrelated per-point scatter beyond statistics (e.g. detector-pair disagreement) | second `diagonal` | robust per-point term s x value per (line, angle) |

Rows of the matrix are the DATA ROWS of the `_cov.dat` named in the `file <path> <rows>` header, not the points
left after cuts: energy windows and angle splits on the segment lines need no new matrix. Every row that is
fitted needs its own diagonal or block variance (the engine factors diagonal + blocks and adds rank-k terms).

## 3. Our own data: recipe

Worked examples (read the docstrings; each lists its terms):
- `10B+a/10-19-26_MANA_Claude/azure/make_neutron.py` + `cov_tools.py` (neutrons: Fisher blocks per run, 5 %
  per-detector relvectors, nine analysis vectors, error inflation), `make_charged.py` (silicon + HPGe: blocks
  for lines fitted in one spectrum, per-detector relvectors, solid-angle and Eu relvectors on every row, a
  pair-audit diagonal; HPGe level blocks from the feeding Jacobian), `make_angint.py` (angle-integrated).
- `11B+a/10-22-26_MANA_Claude/step9_azure/make_azure_11b.py`; `10B+a/10-23-26_June5U_Claude/step10_azure/`.

Rules learned on those:
- Build C from the same production tables the data files are written from, in the same row order, and write
  the data file's error column as sqrt(C_ii) (`_cov.dat`). Keep an AZURE2-alone version of the data file
  (diagonal errors only) for fits without the covariance keyword.
- Statistical errors of spectrum fits with dev/ndf > 1: multiply each fit's Fisher block by max(dev/ndf, 1)
  (DeBoer 2026-10-23/24, "error inflation"); the correlations inside a block are unchanged.
- Differential and angle-integrated versions of the same measurement are NOT both put in one fit: the
  integrals are linear combinations of the same points. Give each its own covariance file and use one or the
  other (or the integrals only where the angular data are not used).
- An alternate file set for a systematic choice that cannot be a covariance term (e.g. every silicon angle
  + 3.7 deg effective-angle offset, `_dth37`) gets its own matrix of the same structure.
- After any change to the production tables, rebuild C and rerun the read test; keep the previous outputs
  (`pre_<tag>/`).

## 4. Literature data: recipe

One covariance file per EXPERIMENT (different experiments uncorrelated), covering all of its data files.

1. **Read the error budget.** The EXFOR entry (`pyazr.nds.fetch_exfor(..., source="local")` or the raw `.x4`):
   `ERR-ANALYS` (the meaning of each error column), `ERR-T` (total), `ERR-S` (statistical), `ERR-SYS`,
   `ERR-1..n` (partial systematics), `MONIT` (normalized to another measurement?), and the free text; then the
   paper. Write down, per set: point-to-point part, correlated parts, scale part, energy-scale part.
2. **Split.** If the quoted point errors include the scale error (an `ERR-T` total), remove the scale in
   quadrature from the diagonal and put it in the norm prior. If only a total is known, use it as the diagonal
   and say so in the README. Never invent correlated terms; add one only when the paper documents it or the
   method makes it certain (several angles measured with separate detectors -> a per-angle efficiency term of
   the size the paper quotes).
3. **Derived quantities.** A 4π value made from an angular distribution, or an excitation function and the
   angular distributions it was integrated from: use one of them, not both. Keep the derived one only when its
   source is not in the fit (Van der Zwan 1973's n1/n2+3 4π points, 10B+a 2026-10-25).
4. **Linked normalizations.** A set normalized to another (Liu 2019 was normalized to Van der Zwan at 90 deg)
   shares that set's `normgroup`; they are not independent scale checks.
5. **Thick-target yields and activation data** need a `<targetInt>` entry; note the stopping powers and the
   target composition next to the segment line.
6. **Deliverables** per experiment: `<tag>_cov.dat` (lab, 4 columns, error = sqrt(C_ii)), `<tag>_cov_rel.txt`,
   ready-to-paste segment lines (norm with its prior, normgroup label, energy-shift prior, energy/angle windows,
   the covariance file name), a correlation plot, and a README table (rows, diagonal source, correlated terms,
   norm prior, shift prior, assumptions). Worked example: `10B+a/10-25-26_Rmatrix_Claude/literature/cov/`.

## 5. Checks before the fit

- **Positive definite** (Cholesky succeeds); report the smallest eigenvalue / largest.
- **sqrt(diag C) equals the `_cov.dat` error column** to 1e-6 relative.
- **Correlation plot** per file: look at it; a detector term should show as a block over that detector's rows.
- **AZURE2 read test** (mode 1, any model with the needed exit pairs): the log must say
  `Covariance <file>: N of N rows fitted` for every file (fewer only where the test model lacks an exit pair,
  which must be stated), and no warnings. Run it on the cluster, not the login node.
- **"WARNING: k point(s) of <file> have an error column more than 1 % from sqrt(C_ii)"** means the data file
  and the covariance file come from different builds. On 10B+a (2026-10-25) a test copy had current data files
  next to a 3-day-old neutron covariance: refresh BOTH together, never one.

## 6. Traps seen so far

- A data set's error column may hide a per-detector problem that is not statistical: the June 2025 HPGe
  detectors disagreed by 1.4x for an isotropic line until their per-detector efficiency anchors were fixed
  (10B+a/10-23-26_June5U_Claude readme, 2026-10-25). Build C after such checks, not before.
- A frame or Jacobian error upstream (lab <-> c.m.) propagates into every derived table and its matrix: the
  June 2025 Legendre fits multiplied by dOmega_cm/dOmega_lab instead of dividing (found 2026-10-25). Check a
  derived angle-integrated set against an independent one before building its matrix.
- `relvector` terms are scaled by a model; in a campaign that model is one fixed reference
  (`cov_reference.npz`, rmfit 2026-10-25). Outside rmfit the engine freezes it at the first evaluation.
