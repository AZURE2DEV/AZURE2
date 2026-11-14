# Normalizations and energy shifts

Every data set carries a normalization n (dimensionless, multiplying the data and their errors) and may carry an
energy shift ΔE (lab frame, added to the data energies: an `offset`, a `sqrt` term b√E, or a `function` of a shared
`energy_scale` record). Each is a nuisance written `<value> fixed|fit [prior <centre> <width>[%]]` on the data set
(`norm`, `shift`) or on a shared record (`normalization`, `energy_scale`) that several data sets name; the syntax is
in `docs/source/reference/format3.rst`. Each enters the objective with its prior. This page says where the values
live, how to audit them, and how to fit an energy calibration. Behaviour reference:
`docs/source/user_guide/chi_squared.rst`, `docs/source/user_guide/segments.rst`.

## Where the values live

| where | what it is |
|---|---|
| the record's value (`norm 1.02 fit prior 1 5%`) | the start of the fit, and the value a calculation uses |
| the record's `prior <centre> <width>` | the reference of the penalty ((n − centre)/width)², whatever the start |
| a parameter file (`param.sav`), by name (`s04_thick.norm`, `s01_int_thin.shift.offset`, `G.offset`) | the fitted value; a run given the file uses it |

A shared record (`normalization`, `energy_scale`) is one parameter, penalised once.

## Consequences

- **The fitted project carries the fit.** `Session.save_fit(path, x)` (default `norms="fitted"`) writes the fitted
  norms and shifts into their records and keeps every prior's centre, so the saved project alone reproduces the
  objective (`verify` checks it). With `norms="nominal"` the records keep the loaded values and the fit is
  project + `.sav`: then a calculation without the `.sav` is not the fit.
- **The prior's centre is the data's own scale** (1.0 for an absolutely normalized set, or the published factor),
  not a previous fit's value. Audit it before trusting any norm penalty:

```python
ev = s.evaluate(x)
for d in ev.data:
    print(f"{d.name:24s} norm {d.norm:.4f}  prior chi2 {d.norm_chi2:.3f}  shift prior chi2 {d.shift_chi2:.3f}")
```

  together with the records' `norm` fields (`pyazr.Project.open(path).describe(name)`).
- **A fixed norm or shift carries no penalty**, and a fixed shift listed in a parameter file is applied from the
  file: when a fixed shift is changed in the project, change it in (or delete it from) every parameter file the run
  reads.
- **Freeing a previously fixed norm or shift starts it from the record's value.** Check for stale values first: six
  data sets of one model carried a forgotten fixed shift of 1.0 MeV from an old conversion; freeing them started
  MIGRAD 50σ from the prior centre (σ = 20 keV).
- **A norm far from its prior's centre** (several σ, or angle-dependent across one experiment) is a symptom, not a
  result (`azure2-check`).

## Fitting a data set's energy calibration

When an excitation function fits at low energy and drifts at high energy, the data's energy calibration is the first
suspect, and a constant offset cannot describe it. In order:

1. **Scan the frozen curve first** (minutes, no refit). Get the model in the data's own observable on a fine grid: an
   extrapolation record for an ordinary data set; for a ratio or composite data set, a copy of the record pointed
   at a pseudo-data file (rows in 0.5 keV steps), checked against the real one at the data energies. Then move the data
   energies, E′ = E + ΔE(E), re-optimizing the norm analytically for every trial
   (1/n = Σ(M y/s²) / Σ(M²/s²), M model, y data, s error), and fit the forms none / a / a + b√E / a + cE /
   a + b√E + cE by Nelder–Mead from a grid of starts. One 201-point set gave χ² = 1953 / 1733 / 1102 / 1123 / 1075:
   ΔE = 0 at E_lab = 1.2 MeV, −15 keV at 3 MeV, with the misfit concentrated at the narrow levels.
2. **Choose the form by the apparatus, not by χ².** An analyzing magnet gives E = kB²: an offset in the field reading
   gives ΔE ∝ √E, an error in k gives ΔE ∝ E (the relativistic term ∝ E², about 8 keV at 4 MeV). Over a factor-3
   energy range √E and linear forms are indistinguishable (1102 against 1123). `shift offset` and `shift sqrt` fit a and
   b; an `energy_scale` record with `function "a0*(x-1.5)" MeV` and `param a0 0 fit prior 0 0.01` fits any form
   (`docs/source/reference/format3.rst`).
   Write the function centred inside the data (`a0*(sqrt(x)-sqrt(E0))`, E0 within the data range): a is then the
   shift at E0 and its penalty means something. The `sqrt` term cannot be centred: its a is the shift at E = 0, so the
   a-penalty fights the fit (one case: scan a = +94 keV, b = −0.059 MeV^½, shift at the data only 22 → 4 keV; a
   16 ± 10 keV a-penalty would cost about 60, more than the whole gain).
3. **Free the shifts, but seed them at the scan optimum.** MIGRAD started from zero shift stops in a local minimum next
   to zero: a and b freed from 0 (σ = 20 keV, 0.02 MeV^½) ended at a = −0.12 keV, b = −10⁻⁴, segment χ² 1,956 → 1,885;
   the same model seeded at the scan values (a = +25.5 keV, b = −0.02335 MeV^½) evaluated to 1,158 before any fitting.
   The frozen one-parameter scans show isolated spikes next to zero (data points crossing narrow structure), enough to
   stop a local optimizer whose initial step is 0.01 × the penalty width. Seed in the records (`shift offset 25.5 keV
   fit prior 0 20 keV`) and verify with `s.evaluate(x)` that the seeded start reproduces the scan's segment χ² before
   submitting. "Give it 20 keV and let it fit" is safe only when the scan says the optimum is a keV or two away.
4. **Checks already done** need not be repeated: the shift is applied correctly on a ratio data set (constant −2 keV:
   engine 1,733.7 against frozen 1,733); `tests/energy_shift_sqrt` recovers a known (a, b) to 5 × 10⁻⁷. Keep penalty
   widths honest: 20 keV and 0.02 MeV^½ cost about 1.5 at the optimum above; a width tight enough to bias the result
   only hides the calibration problem.
5. **One calibration per experiment, not per segment.** A data set split into segments (one per angle, one per figure)
   has one accelerator calibration: one `energy_scale` record named by every member's `shift`. Per-segment
   energy-dependent terms fit angle-by-angle structure instead: in one model b freed per segment gained 91 in χ² for 10
   parameters, but b changed sign between angles of one experiment (−0.016 at 54.7° against +0.006 at 90°), which no
   calibration does. Conversely, grouping digitized sets can cost a lot (+1,913 in a frozen group scan: different
   figures, different digitization offsets); group by what shares a calibration, and test it. When per-segment constants
   are real but the energy-dependent calibration is one, give each member its own `shift offset` beside the shared
   `energy_scale`. The start-up log prints, per function, the condition number of its free terms over the group's
   energies, and warns about a term that only duplicates a (√x on a single-energy angular distribution).
