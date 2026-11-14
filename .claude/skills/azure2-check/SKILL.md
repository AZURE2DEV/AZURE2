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

## Scripts

| Script | Does |
|---|---|
| `scripts/chi2_report.py` | per data set χ², the largest point, pulls within 1σ/2σ, same-sign runs; χ²/ν, terms |
| `scripts/widths_check.py` | θ² (both forms), W.u. of photons, the norm J of every level; flags |
| `scripts/radius_scan.py` | the channel-radius scan: a project copy and a refit per radius, walked outward |
| `scripts/selftest.py` | runs the three on `tests/upgrade/T1` |

```bash
S=.claude/skills/azure2-check/scripts
python3 $S/chi2_report.py model.azr --params output/param.sav --chisquared-out output/chiSquared.out
python3 $S/widths_check.py model.azr --params output/param.sav --background-above 20
python3 $S/radius_scan.py model.azr --params output/param.sav --pair 1 --radii 4.5 5.0 5.5 6.0 --table radius.txt
```

## The checklist

Run the items in order. Each says what it detects.

1. **The project is what was meant.** `AZURE2 --validate model.azr` (or `pyazr.validate(path)`): every diagnostic as
   `file:line: record.field: severity code: message`, exit 2 with an error; `AZURE2 --describe model.azr` (or
   `--record NAME`): what the engine made of each data set (observable, reaction tree with pair keys, rows kept inside
   the windows, effects, parameters) and the options with their source; `AZURE2 --print-predictions model.azr`: the
   expression each data set evaluates. A warning code (`W_...`) is a decision to make, not noise
   (`docs/source/reference/format3.rst`, "Diagnostics").
2. **The number reproduces.** `chi2_report.py --chisquared-out output/chiSquared.out` compares pyazr's data χ² with
   the program's (agreement 1e-6 on T1). `ev.chi2` is the data χ² only; the objective adds every prior
   (`ev.objective`, `ev.terms`). Record the total whenever the project changes. A jump concentrated in a few data
   sets: a stale capture cache, a `param.sav` whose data-set parameter names no longer match the records (renamed or
   removed data sets), or an ill-conditioned background level.
3. **Per-data-set χ² and pulls.** The table of `chi2_report.py`; report it with every result. A good total can hide
   one set at χ²/N = 10. A data set whose largest point carries most of its χ² is a data question.
4. **Residual structure.** A long run of pulls of one sign (column `run`) over adjacent energies is a model defect,
   not noise: azure2-resonance-search before adding structure by hand.
5. **χ²/ν**, ν = N − free parameters (norms and shifts included). Well above 1 without structure: underestimated
   errors (scale the covariance by it, `--scale-covariance`); far below 1: overestimated or double-counted errors.
6. **Physical limits.** θ² ≤ 1 for particle channels, a few W.u. for γ channels, J not far below 0.5 for levels
   below the background poles (section "Wigner and Weisskopf limits").
7. **Background poles.** Saturated poles (section below).
8. **Channel radius.** The scan with refits (section "Channel radius").
9. **Level significance.** Each level's Δχ² (section "Is a level needed?").
10. **Grid convergence.** An extrapolated curve, a rate or a folded point is a quadrature: halve the step (or raise
    the data set's `quadrature points`) and recompute; the quantity must move less than the tolerance you quote.
11. **Thread independence.** The data χ² at `OMP_NUM_THREADS=1` and 8 agree to rounding (a few units in the last
    place); a larger difference is a defect: report it.
12. **Frames and units** (table below) and **data definitions** (section below) before any new structure.

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
- **One point near a channel threshold can make χ² explode** (10³⁰ to 10⁷²). The Coulomb functions fail when a point's
  c.m. energy, after the segment's energy shift and including every sub-point of a convolution or target integration,
  lies within a few keV of another charged-particle channel's threshold (`docs/source/reference/troubleshooting.rst`,
  "Threshold Effects"). It may appear only for some parameter sets or only on some builds. Diagnose: sort
  `chiSquared.out` by segment χ²; in that segment's `AZUREOut_*.out` find the row whose fit column (column 4) is absurd;
  convert its E_c.m. to E_x and compare with every particle threshold. Fix: write a new data file without the point
  (keep the original; never edit a file a running job reads), repoint the data set (`Project.set(name, file=...)`),
  verify with a `calculate` run. List the other points within 5 keV of a threshold and watch them in the next fit.
- **Run a plain calculate before submitting a fit** after a structural edit combined with a reparametrization (a
  remapped `param.sav`, a fresh `save_fit` snapshot): `AZURE2 calculate` with the new parameter file, or
  `s.objective()(s.parameters.read("new.sav"))`, in a throwaway output directory. Compare with the old total corrected
  for exactly what changed (the removed data set's χ² and N from the old `chiSquared.out`). It costs seconds and catches
  a corrupted start before a multi-day job runs on it.
- **A start far above the baseline is wrong even when the cause is known.** If the calculate check is orders of
  magnitude above the known total because of a deliberate change (a corrected energy shift with norms still fitted
  to the old value), re-seed the affected norms (the n* formula, `docs/source/user_guide/chi_squared.rst`) before
  fitting. Fits started 100× above baseline from one mis-scaled segment converged to a different level running away
  each time; the re-seeded start converged cleanly.

## Wigner and Weisskopf limits

`scripts/widths_check.py` prints, at any parameter vector, θ² = Γ/Γ_W (Γ_W = 2P_l(E_λ)γ²_W, eV) and
θ²_formal = γ²/γ²_W (γ²_W = 3ħ²/(2μa_c²), MeV) of every particle channel, the Weisskopf estimate and W.u. of every
photon channel, and the norm J of every level; it flags θ² > 1, more than 100 W.u., and J < 0.3 below
`--background-above`. In Python:

```python
t = s.dimensionless_widths(x)            # a WidthTable of ChannelWidth rows
print(t.photons.nonzero.table())         # filters: .particles .photons .free .nonzero
[c for c in t.particles if c.theta2 and c.theta2 > 1]
```

Closed channels (`is_open` False) carry an ANC and no Γ_W. The two forms of the limit, the 1.5 between conventions,
the Weisskopf estimates and the level norm J (saturated levels): `references/widths.md`.

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

- **Existing level.** Rank by frozen Δχ²: the χ² with the level's free widths zeroed at fixed parameters
  (`P.put(x, [p for p in P.of_level(label).of_kind("width") if not p.fixed], 0.0)`; azure2-calculate's
  `decompose.py` does the same for the model). The frozen Δχ² is an upper bound; the refit Δχ² (the level removed with
  `Project.remove_level`, or deactivated with `Project.set_level(level, active=False)`, then refit) is the value to
  report, with the data sets whose χ² moved.
- **New level or channel.** Add it, refit, and weigh Δχ² against the parameters spent and the limits above. The gated
  search, significance tests and interference signs: azure2-resonance-search.

## Channel radius

A reduced width means something else at another radius, so a χ² at a new radius means nothing before a refit.
`scripts/radius_scan.py` writes one copy of the project per radius (`Project.set_pair(pair, radius=a)`, its own
output directory, so that no capture cache of another radius is reused), refits each with `least_squares` on the
objective's rows and Jacobian (or the program's fitter, `--fitter engine-lm`), and walks outward from the project's
radius so that every fit starts from its converged neighbour. Report: the objective against a_c, the radii where it is
within 1 (or the Δχ² of your interval) of the minimum, and how the quantities you quote (S(0), a width, a rate) move
over that range: that spread is the radius systematic. A minimum at the end of the scan wants a wider scan; a scan
flat to a fraction of a unit says the data do not choose a_c (state the radius used). Sampling a_c: azure2-uncertainty.

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

## Constants and logs

- **Physical constants are CODATA 2018** (ħc = 197.3269804 MeV fm, u = 931.49410242 MeV, α = 1/137.035999084,
  m_e = 0.51099895 MeV). A rerun of a fit archived with older constants does not reproduce χ² to the last digit
  (1e-9 to 1e-5 relative); that is not a model change.
- **Results of older builds**: what each build range computed differently, and what to recompute:
  `references/older-results.md`.

## Pointers

- `docs/source/user_guide/chi_squared.rst` (the objective, penalties, what `chiSquared.out` reports),
  `docs/source/reference/troubleshooting.rst`, `docs/source/user_guide/pyazr.rst`.
- `references/widths.md`, `references/background-pole-saturation.md`.
- Examples: `pyazr/examples/per_dataset_chi2.py`, `dimensionless_widths.py`, `channel_radius_scan.py`.
