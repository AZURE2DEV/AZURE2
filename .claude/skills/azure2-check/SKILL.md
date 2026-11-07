---
name: azure2-check
description: Decide whether an AZURE2 result can be trusted before it is reported or reused - reproduce the total chi2, per-dataset chi2 and pulls, residual structure, chi2/nu, Wigner and Weisskopf limits (theta^2, W.u.), saturated background poles, channel-radius sensitivity, level significance, grid convergence, thread independence, frames and units, nan and threshold failures, data-definition errors. Use after every fit. Fitting is azure2-fit; uncertainties and bands are azure2-uncertainty; missing levels are azure2-resonance-search.
---

# Checking an AZURE2 result

A fit returns a parameter vector p̂ (level energies E_λ in MeV, reduced-width amplitudes γ_λc in MeV^½, normalizations
n_k, energy shifts) and a total χ². This skill lists the checks that decide whether p̂ and the observables computed
from it can be trusted: that the number is reproducible, that the misfit is statistical and not structural, that the
parameters are physical, and that the answer does not depend on choices the data do not make (channel radius a_c,
background poles, grids, thread count). The theory of the χ² is in `docs/source/user_guide/chi_squared.rst`.

## When to use, and when a sibling owns the task

- Use after every fit, after every structural edit (level or segment added or removed, parameter file remapped), and
  before a number or a figure is reported.
- Fitting, starting values, editing levels: `azure2-fit`. Uncertainties, bands, model averaging: `azure2-uncertainty`.
  A structure in the residuals that suggests a missing level: `azure2-resonance-search`. Figures: `azure2-plot`.
  Target and resolution kernels: `azure2-targets`. A general audit of an uncertainty budget: the user-level
  `uncertainty-audit` skill; a referee-style read of the physics: `physics-referee`.

## The checklist

Run the items in order. Each says what it detects.

1. **The number reproduces.** `cat output/chiSquared.out`: a finite total χ² and the expected number of points N per
   segment. From pyazr, `np.sum(m.calculate_chi2_rwa(m.params_rwa))` (data term only) or `m.objective(x)` (data plus
   penalties) must reproduce the model's recorded total; record that total whenever the `.azr` changes. A jump
   concentrated in a few data sets means one of: a stale `output/intEC.extrap`; a `param.sav` whose segment-keyed
   norm or shift names no longer match an edited `<segmentsData>` (renumbering: `azure2-fit`); or ill-conditioning
   in a background level that a `save_fit` round trip's own numerical noise exposes. After a fit, confirm that
   `param.sav` was rewritten and read `parameters.out`.
2. **Per-data-set χ² and pulls.** `m.dataset_chi2(x)` (or `pyazr/examples/per_dataset_chi2.py`) and the pulls
   (data − fit)/σ from `m.residuals(x)`. A total that is acceptable can hide one set at χ²/N = 10. Report the per-set
   table with every result.
3. **Residual structure.** Pulls that keep one sign over a run of adjacent energies or angles are a model defect,
   not noise. Rank such runs and treat them with `azure2-resonance-search` before adding any structure by hand.
4. **χ²/ν.** With ν = N − (number of free parameters, including normalizations and shifts). A χ²/ν well above 1 with
   no structure means underestimated errors; scale the parameter covariance by it (Birge ratio,
   `--scale-covariance`) rather than ignoring it. A χ²/ν far below 1 means overestimated or double-counted errors
   (a scale error in both the covariance and the norm prior: `azure2-uncertainty`).
5. **Physical limits.** Dimensionless reduced widths θ² ≤ 1 for particle channels and γ-ray strengths of at most a
   few Weisskopf units for γ channels (section "Wigner and Weisskopf limits"). θ² > 1 is unphysical.
6. **Background poles.** A pole whose observed widths sit against their ceiling (level-shift factor J → 0) is
   saturated: its widths are unconstrained and the minimum is flat along that direction. Test with caps (section
   "Saturated background poles").
7. **Channel radius.** Refit at two or three radii a_c (one directory per radius, warm-started from the converged
   neighbour). A result that moves by more than its statistical error carries a radius systematic; report it.
8. **Level significance.** Each level's Δχ² (section "Is a level needed?").
9. **Grid convergence.** An extrapolated curve, a rate, or a folded (target-integrated, convolved) point is a
   quadrature on a grid. Halve the grid spacing (or double the sub-points) and recompute; the quantity must move by
   less than the tolerance you quote. Kernel sub-point grids: `azure2-targets`; extrapolation grids: `azure2-calculate`.
10. **Thread independence.** Results do not depend on `OMP_NUM_THREADS`. Recompute the total at 1 and 8 threads;
    they agree to the cross-compiler level (about 1e-5 relative) or better. A larger difference is a bug: report it.
11. **Frames and units.** Section "Frames and units" below; check every number that crosses a file or code boundary.
12. **Data definitions** before structure. Section "Data and their definitions".

## Tolerance for "the numbers do not match"

pyazr against the CLI, and one session against another, differ by up to about 1e-3 relative on an ill-conditioned
model; the largest differences are in channels coupled to broad background poles (measured 6e-5 to 1e-3 per data set
on a 12C+α model). Below 1e-3 treat the difference as noise. Around 1e-2, note it. Investigate only a larger gap, or
one that grows into a failed fit, and then investigate which parameter moved, not the noise floor.

## Failures that look like results

- **An ANC must be physical for the channel radius.** An ANC of 84 fm^-1/2 on a p-wave state bound by 100 keV puts
  130 % of the state outside a = 5 fm; AZURE2 warns `Denominator less than zero while transforming` and the cross
  section stops scaling as ANC². The input is at fault: lower the ANC or enlarge the radius.
- **A `nan` point or `Total-Chi-Squared: nan` is a failed calculation.** The log holds
  `WARNING: the calculation failed at segment K (...), E_lab ... : <error>`, printed once per run. Test hook for such
  paths: `AZURE_TEST_FAIL_E=<E_lab in MeV>` (with `AZURE_TEST_FAIL_EXIT=<pair key>`).
  Builds before 2026-11-01 (expires: 2027-02-01): a failed point was a silent 0 in the CLI objective and pre-pass, a
  failed or NaN term of a SUM or natural-target composite was left out of the sum, and a failed RATIO denominator
  gave 0, so χ² looked normal; on such a build a composite segment with exact zeros or a low SUM is a candidate.
- **One point near a channel threshold can make χ² explode** (10³⁰ to 10⁷²). The Coulomb functions fail when a
  point's c.m. energy, after the segment's energy shift and including every sub-point of a convolution or target
  integration, lies within a few keV of another charged-particle channel's threshold
  (`docs/source/reference/troubleshooting.rst`, "Threshold Effects"). It may appear only for some parameter sets or
  only on some builds. Diagnose: sort `chiSquared.out` by segment χ²; in that segment's `AZUREOut_*.out` find the row
  whose fit column (column 4) is absurd; convert its E_c.m. to E_x and compare with every particle threshold. Fix:
  write a new data file without the point (keep the original; never edit a file a running job reads), repoint the
  segment, verify with a mode-1 run. List the other points within 5 keV of a threshold and watch them in the next fit.
- **Run a plain calculate before submitting a fit** after a structural edit combined with a reparametrization (a
  remapped `param.sav`, a fresh `save_fit` snapshot): CLI mode 1 with the new parameter file, or
  `m.objective(m.params_rwa)`, in a throwaway output directory. Compare with the old total corrected for exactly what
  changed (the removed segment's χ² and N from the old `chiSquared.out`). It costs seconds and catches a corrupted
  start before a multi-day job runs on it.
- **A start far above the baseline is wrong even when the cause is known.** If the calculate check is orders of
  magnitude above the known total because of a deliberate change (a corrected energy shift with norms still fitted
  to the old value), re-seed the affected norms (the n* formula, `docs/source/user_guide/chi_squared.rst`) before
  fitting. Fits started 100× above baseline from one mis-scaled segment converged to a different level running away
  each time; the re-seeded start converged cleanly.

## Wigner and Weisskopf limits

`m.dimensionless_widths(x)` returns θ² for every particle channel and W.u. for every γ channel at any parameter
vector (`docs/source/user_guide/pyazr.rst`, "Dimensionless widths"; `pyazr/examples/dimensionless_widths.py <azr>
[--params output/param.sav]`):

```python
t = m.dimensionless_widths(best)          # WidthTable of ChannelWidth rows
print(t.photons.nonzero.table())          # filters: .particles .photons .free .nonzero
[c for c in t.particles if c.theta2 and c.theta2 > 1]
```

Each row carries `theta2` (= Γ/Γ_W), `theta2_formal` (= γ²/γ²_W), `wigner_width` (Γ_W, eV), `wigner_gamma2` (γ²_W,
MeV) and, for photons, `e_gamma`, `weisskopf`, `wu`. Closed (sub-threshold) channels have `is_open = False` and no
Γ_W: AZURE2 gives an ANC there. `m.wigner_widths(x)` returns `{free_index: Γ_W}`; `weisskopf_width(rad, L, E_γ, A)`
and `m.mass_number` are separate. The Wigner limit exists in two forms; never mix them, and never divide by
`wigner_limit` twice. Details, the Weisskopf estimates and the 1.5 factor between conventions:
`references/widths.md`.

## Saturated background poles

An unbounded background pole drifts to large reduced widths while χ² improves by nothing measurable: the extra
amplitude only drives the level-shift factor J = 1/(1 + Σ_c γ_c² dS_c/dE) toward 0 and the observed widths
Γ_c = 2 P_c γ_c² J against their ceiling. Consequences: the reported widths are meaningless; small numerical
differences are amplified by 1/J (a 1 % change in dS/dE moves the observed widths 1 % and can make a stored `.azr`
fail to transform, `Denominator less than zero while transforming`); fitting slows and becomes unstable. In Park mode
the same flat direction shows as J → 0 (`azure2-exchange`).

Test: one refit per cap k, |γ_c| ≤ k γ_W on the pole's channels, k = 100, 30, 10, 3, from the fit's `param.sav` with
the capped amplitudes clipped to 0.9 of the cap. If the loosest cap costs nothing, prefer the capped solution or
report the pole as unconstrained; the k at which χ² starts to rise measures how much background the data require.
Syntax and a worked example: `references/background-pole-saturation.md`.

## Is a level needed? Is a new one warranted?

- **Existing level.** Rank by frozen Δχ²: recompute χ² with each level switched off at fixed parameters
  (`m.without_level(x, jpi=..., energy=...)`). For the candidates, refit with the level removed
  (`AzrModel.remove_level` or `deactivate_level`, then refit). The frozen Δχ² is an upper bound; the refit Δχ² is the
  value to report. Report both, and which data sets' χ² move: a level supported by one data set is a different claim
  from one supported by all.
- **New level or channel.** Add it (`add_level`, or free one more channel of an existing level), refit the subset it
  can affect, and weigh Δχ² against the parameters spent and the physical limits (θ² ≤ 1; a few W.u.). Candidate
  sweeps run in parallel sessions. The gated search procedure is `azure2-resonance-search`.

## Channel-radius sensitivity

`m.set_channel_radius(pair, r)` on a live session, or `AzrModel.set_channel_radius(pair, r)` and `write()` to
persist it. The radius is the matching surface: AZURE2 rebuilds the compound nucleus and data, recomputes
penetrabilities, shift functions, boundary conditions, Wigner limits and every external-capture integral (the API
does not reuse `output/intEC*`). A reduced width means something different at another radius, so the model is no
longer fitted: refit before reading anything. For a scan, one directory per radius (own `output/`, `data`
symlinked), walk outward, warm-start each radius from its converged neighbour; starting all from one point biases
the ends. `set_channel_radius` raises if the rebuilt model does not carry the requested radius, and overrides persist
through every rebuild (other radii, nuclear-potential changes, `data_mode()`/`extrap_mode()`).
Builds before 2026-11-02 (expires: 2027-02-01): a second `set_channel_radius` undid the first, and so did a
nuclear-potential change or a mode switch; write every radius into the `.azr` and open a fresh session instead.
Sampling a_c inside a likelihood: `azure2-uncertainty`.

## Data and their definitions

Most of the excess χ² of a large fit comes from data definitions and from figures that do not show the fit, not
from the level scheme. Check these before adding structure. Summing rules for totals and shared data files are in
`azure2-setup`; figure-against-output checks in `azure2-plot`.

- **Never fit a derived point together with the data it was derived from.** EXFOR sometimes carries a paper's own
  R-matrix or polynomial extrapolation as a separate subentry, most often S(0). The tell: a lone point or a few
  points labelled "S-factor point" (or "SFC") outside the energy range of that paper's measurements (e.g. one point at
  E_c.m. = 25 keV from a paper whose data start above 1 MeV). Fitting it is circular. Check EXFOR provenance before
  activating any single- or few-point segment.
- **A large, angle-dependent normalization is a symptom.** When a set needs n = 1.7 to 1.9 varying with angle, look
  first for a definition error in another set sharing its channels; such norms fell to 1.0 to 1.2 once the other set
  was corrected.
- **Find outliers by comparing data sets, not only by pulls.** A point 6σ off the fit may be model trouble until a
  second measurement of the same quantity shows that the point is off. Remove it only on the evaluator's decision,
  into a new `_clean` file, and record the row and the reason.
- **Rebuild combined data files with a script, and diff new deliveries row by row.** The script must first
  reproduce the existing combined file exactly from the old inputs (row order, dropped duplicates, added
  systematics). A median ratio of 1.000 between two deliveries hid point-to-point changes of up to 40 %.
- **Judge a data set point by point, not by its χ²/N.** A large χ²/N is often carried by a few points (an outlier,
  one energy off at every angle, the edge of a sharp resonance where a small energy offset gives a large residual).
  Compute r = (model − data)/σ per point from `output/AZUREOut_*.out` (columns 4, 6, 7: fit, data, error; Σ r² over a
  segment equals its `chiSquared.out` entry, check that first) and report the fraction within 1σ and 2σ, the median
  |r|, and how many points carry most of the χ² and where ("15 of 75 points carry 89 % of the χ²: the edge of the
  7/2⁻ resonance, ..."). Flag possible outliers; exclusion is the evaluator's decision.
- **Print χ² per segment and the largest single-point χ² per segment.** A segment whose maximum is most of its total
  is a data question, not a physics one. One point a factor 6.6 below both energy neighbours (pull −61, χ² 3,746)
  held a converged fit: dropping it removed its 3,746 and the refit fell a further 3,743 in data sets that do not
  contain it. Predict the frozen effect of a data change before refitting, so that everything beyond it is
  attributable to the refit; never edit a data file a running job reads (new name, new directory).
- **Secondary-γ detectors viewing a J = 1/2 state must agree in absolute scale.** Its γ emission is isotropic, so
  angle-to-angle differences (10–15 %, χ²/ν = 3.1 between angles) measure detector-scale errors, not physics.
- **Compare data handed over between analyses by value at common energies**, never by file name: a renamed or
  re-binned file with the same name is not the same data.
- **Legendre coefficients are a consistency tool, not the fit target.** Fitting each angular distribution with
  Σ_L a_L P_L(cos θ_c.m.) per energy allows angle-free comparisons between sets: per-detector ratios against a
  reference set; an angle offset scales with ∂ln σ/∂θ; acceptance smearing shows as a_L/a_0 falling like the
  attenuation factors Q_L of a cone. They also bound J: a_6 ≠ 0 requires J ≥ 7/2, a_8 ≠ 0 requires J ≥ 9/2. Fit the
  R-matrix to the angular distributions themselves; the coefficients carry correlated errors and truncation choices.

## Frames and units

The conventions every number crossing a file boundary must respect (the reference for `CLAUDE.md`):

| Quantity | Frame and unit | Where stated |
|---|---|---|
| Data energies, angles, cross sections | lab frame, forward kinematics (light particle is the projectile); E_lab in MeV, θ_lab in degrees, σ in b or b/sr | `docs/source/reference/data_formats.rst` |
| Exceptions | analyzing-power data angles are c.m.; capture data may be lab-angle and are converted with the photon aberration | `docs/source/reference/data_formats.rst` |
| All output files | c.m. frame; σ in b or b/sr, S in MeV b | `docs/source/reference/output_files.rst` |
| E_x column of `AZUREOut_*` | E_c.m. + S + E* of the entrance pair (includes the pair's excitation) | `docs/source/reference/output_files.rst` |
| Particle masses | atomic masses in u (AME), not mass numbers; checked against separation energies at load | `docs/source/user_guide/particle_pairs.rst` |
| `<levels>` width field | open channel: partial width Γ in eV (not a reduced-width amplitude); closed channel: ANC in fm^-1/2 | `docs/source/user_guide/levels_channels.rst` |
| Fit parameters | E_λ in MeV; γ_λc in MeV^½ in the run's basis (Brune, standard, Park) | `docs/source/user_guide/chi_squared.rst` |
| Wigner limit | γ²_W = 3ħ²/(2μa_c²) in MeV (Teichmann-Wigner); Γ_W in eV | `references/widths.md` |
| Rates | N_A⟨σv⟩ in cm³ mol⁻¹ s⁻¹ against T_9 | `azure2-calculate` |

Name the frame and unit in every column header, variable and plot axis you write. Details of target and straggling
quantities (some are c.m.): `azure2-targets`.

## Builds, constants and logs

- **Physical constants are CODATA 2018** (ħc = 197.3269804 MeV fm, u = 931.49410242 MeV, α = 1/137.035999084,
  m_e = 0.51099895 MeV, π to double precision). Builds before 2026-11-05 (expires: 2027-02-01) used π = 3.141592650
  (1.1e-9 low; a 1/√(2π) normalization inherits half of it, so a comparison with a closed form cannot be tighter than
  about 6e-10) and CODATA 2006 or rounded constants. Results move by 1e-9 to 1e-5 relative between the two, more near
  a kinematic edge: a rerun of an archived fit does not reproduce χ² to the last digit, and that is not a model change.
- **Torn progress bars in multi-threaded logs are cosmetic.** Builds before 2026-11-06 (expires: 2027-02-01) let
  every OpenMP thread redraw the external-capture bars, leaving fragments such as `Segment #\r\tSegment #1111 [****`;
  results and `intEC.dat` are identical at 1 and 8 threads. A script that diffs two run logs drops lines that contain
  `\r` or `Segment #` (as `tests/light_excitation` does).

## Pointers

- `docs/source/user_guide/chi_squared.rst` (the objective, penalties, what `chiSquared.out` reports),
  `docs/source/reference/troubleshooting.rst`, `docs/source/user_guide/pyazr.rst`.
- `references/widths.md`, `references/background-pole-saturation.md`.
- Examples: `pyazr/examples/per_dataset_chi2.py`, `dimensionless_widths.py`, `channel_radius_scan.py`.
