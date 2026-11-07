# Normalizations and energy shifts

Every data segment carries a normalization n (dimensionless, multiplying the data and their errors) and an energy
shift ΔE (MeV, lab frame, added to the data energies; a constant, a `sqrtshift` a + b√E, or a `shiftfn` function).
Each exists in two places that mean different things, and each enters the objective with a penalty. This page states
the two homes, the penalties, how to audit them, and how to fit an energy calibration. Behaviour reference:
`docs/source/user_guide/chi_squared.rst`, `docs/source/user_guide/segments.rst`, `docs/source/user_guide/fitting.rst`.

## The two homes

| where | field | what it is |
|---|---|---|
| `<segmentsData>` line, tokens 9 and 12 (1-based; `dataNorm`, `energyShift`) | the `.azr` | the nominal value: the starting value and, for the norm only, the reference of the penalty χ²_n = ((n − n₀)/(n₀ · normErr/100))² (`AZURECalc.cpp`, `GetNominalNorm`; `normErr` is a percentage) |
| `param.sav`, `segment_<key>_norm`, `segment_<key>_energy_shift` | written by every fit | the fitted value; a run given this file uses it, while the norm penalty still refers to the `.azr` nominal |

Energy-shift penalties are measured from zero: a free constant shift, `sqrtshift` b and `shiftfn` term a_j add
(ΔE/σ)², (b/σ)², (a_j/σ)²; the value on the segment line is only where the fit starts. pyazr's
`azure2.SHIFT_PENALTY_FROM_ZERO` says which rule an engine uses; rmfit follows it. Norm penalties are measured from the
nominal norm. builds before 2026-10-25 (expires: 2027-02-01) measured the shift penalty from the nominal shift, so an
`.azr` whose shifts were seeded from an earlier fit gives a larger χ² on current builds (one model: shift penalty
57.6 → 317.3 at the same parameters, 45 inherited nominal shifts).

A `shiftgroup` shares one shift among its members; the engine penalises it once, on the first member.
builds before 2026-10-27 (expires: 2027-02-01): pyazr `penalties()` charged each follower a second shift penalty, so
`objective()` differed from the CLI by a constant (gradients and the minimum unaffected).

## Consequences

- **A fit is `.azr` + `.sav`, never the `.azr` alone.** Nothing writes fitted norms back into `<segmentsData>`: not a
  CLI fit (it never rewrites the `.azr`), not `save_fit` (it rewrites `<levels>` and writes the `.sav`, and leaves
  segment lines as found). A calculation with a blank parameter file uses the nominal norms and shifts, not the fit
  (one snapshot: 11,112 against the fit's 9,970; 64 norms differed by up to 21 %). The correct check of a snapshot is
  mode 1 with the `.sav`, which must reproduce `chiSquared.out` including `Total-Norm-Chi-Squared` to about 10⁻⁴
  relative (the residual is the 7-digit `.sav` rounding passed through adaptive target-integration grids).
- **Do not bake fitted norms into `<segmentsData>`.** For a penalized segment (`normErr != 0`) that moves the penalty
  reference onto the fitted value: the blank-file run prints `Total-Norm-Chi-Squared 0`, and a refit from that file
  minimizes a different objective. Baking is harmless only for unpenalized norms.
- **Audit the nominals.** They should be the data's own scale: 1.0 for an absolutely normalized set, or the published
  scale factor. An `.azr` that has been through Qt GUI saves or old snapshots can carry a previous fit's norms as
  nominals (one master file had 0 of its 27 penalized segments at 1.0, values 0.74–1.71), so its 5–30 % systematic
  errors penalize deviations from an arbitrary earlier fit. List `dataNorm` against `normErr` for every penalized
  segment before trusting any `Total-Norm-Chi-Squared`, and say so in the readme when they are not 1.
- **Plotting.** The `.azr` holds the nominals; plots come from `output/`. Data scaled "by
  the fit norm" use the `.sav` value (`normalizations.out` prints it), not the `.azr` field.
- **A norm that is not varied comes from the `.azr` field; its `.sav` value is ignored.** Fixed level parameters do
  take the `.sav` value. To freeze norms at a fit's values (a local refit of a few level parameters with everything
  else held), write the fitted norm into each segment's `dataNorm` as well as setting `varyNorm` to 0; flipping the
  flag alone gave 204,302 instead of the fit's 126,366 in one model. Fixed norms carry no penalty, so
  `Total-Norm-Chi-Squared` reads 0 even when the data χ² reproduces exactly.
- **An energy shift listed in the parameter file is applied whether or not it is varied.** Every `param.sav` written by
  a fit lists `segment_<k>_energy_shift` (and `_shift_sqrt`) for every segment. A fixed +2.75 keV centring shift set
  in the `.azr` was replaced by the 0 of the refit's `param.sav` (that segment's χ² 39,315 instead of 7,905). When a
  fixed shift is set or changed in the `.azr`, change it in, or delete it from, every parameter file the run reads.
- **Freeing a previously fixed norm or shift starts it from the value in the file.** Check for legacy values first:
  six segments of one model carried a forgotten fixed `energyShift` of 1.0 MeV from an old conversion; freeing them
  started MIGRAD 50σ from the prior centre (σ = 20 keV), and the old `param.sav` carried the same value, so it had to
  be reset in both places. Norms fitted around the stale shift (one 20× off its nominal) needed re-seeding.

```python
# nominal vs fitted, penalized segments only
import re, io
S = {}
for l in open("output/param.sav"):
    p = l.split()
    if len(p) >= 2:
        try: S[p[0]] = float(p[1])
        except ValueError: pass
c = io.open("fit.azr", encoding="latin-1").read()
rows = [l.split() for l in re.search(r'<segmentsData>(.*?)</segmentsData>', c, re.DOTALL).group(1).split('\n') if l.split()]
for k, t in enumerate(rows, 1):
    if t[0] == '1' and float(t[10]) > 0:
        nom, fit, err = float(t[8]), S["segment_%d_norm" % k], float(t[10])
        print(k, t[-3], "nominal %.4f fitted %.4f err %g%% penalty %.2f" % (nom, fit, err, ((fit-nom)/(nom*err/100))**2))
```

Whenever you hand over, plot from, or build a variant on a `.azr`, say which home the numbers came from.

## Fitting a data set's energy calibration

When an excitation function fits at low energy and drifts at high energy, the data's energy calibration is the first
suspect, and a constant `energyShift` cannot describe it. In order:

1. **Scan the frozen curve first** (minutes, no refit). Get the model in the data's own observable on a fine grid: a
   `<segmentsTest>` line for an ordinary segment; for a UPOS, ratio or composite segment, a copy of the segment line
   pointed at a pseudo-data file (`E 0 1000 100` rows in 0.5 keV steps) run in mode 1, checked against the real
   segment's `AZUREOut` column 4 at the data energies (agreement 10⁻⁵). Then move the data energies,
   E′ = E + ΔE(E), re-optimizing the norm analytically for every trial
   (1/n = Σ(M y/s²) / Σ(M²/s²), M model, y data, s error), and fit the forms none / a / a + b√E / a + cE /
   a + b√E + cE by Nelder–Mead from a grid of starts. One 201-point set gave χ² = 1953 / 1733 / 1102 / 1123 / 1075:
   ΔE = 0 at E_lab = 1.2 MeV, −15 keV at 3 MeV, with the misfit concentrated at the narrow levels.
2. **Choose the form by the apparatus, not by χ².** An analyzing magnet gives E = kB²: an offset in the field reading
   gives ΔE ∝ √E, an error in k gives ΔE ∝ E (the relativistic term ∝ E², about 8 keV at 4 MeV). Over a factor-3
   energy range √E and linear forms are indistinguishable (1102 against 1123). The `sqrtshift` block fits a and b; a
   `shiftfn` block fits any form (`docs/source/reference/data_formats.rst`), e.g. `shiftfn "a0*(x-1.5)" 1 0 0.01 1`.
   Write the function centred inside the data (`a0*(sqrt(x)-sqrt(E0))`, E0 within the data range): a is then the
   shift at E0 and its penalty means something. `sqrtshift` cannot be centred: its a is the shift at E = 0, so the
   a-penalty fights the fit (one case: scan a = +94 keV, b = −0.059 MeV^½, shift at the data only 22 → 4 keV; a
   16 ± 10 keV a-penalty would cost about 60, more than the whole gain).
3. **Free the shifts, but seed them at the scan optimum.** MIGRAD started from zero shift stops in a local minimum next
   to zero: a and b freed from 0 (σ = 20 keV, 0.02 MeV^½) ended at a = −0.12 keV, b = −10⁻⁴, segment χ² 1,956 → 1,885;
   the same model seeded at the scan values (a = +25.5 keV, b = −0.02335 MeV^½) evaluated to 1,158 before any fitting.
   The frozen one-parameter scans show isolated spikes next to zero (data points crossing narrow structure), enough to
   stop a local optimizer whose initial step is 0.01 × the penalty width. Seed through the `.sav`
   (`segment_K_energy_shift`, `segment_K_energy_shift_sqrt`; a varied shift takes the `.sav` value) and verify by mode
   1 that the seeded start reproduces the scan's segment χ² before submitting. "Give it 20 keV and let it fit" is safe
   only when the scan says the optimum is a keV or two away.
4. **Checks already done** need not be repeated: the shift is applied correctly on a UPOS segment (constant −2 keV:
   engine 1,733.7 against frozen 1,733); `tests/energy_shift_sqrt` recovers a known (a, b) to 5 × 10⁻⁷. Keep penalty
   widths honest: 20 keV and 0.02 MeV^½ cost about 1.5 at the optimum above; a width tight enough to bias the result
   only hides the calibration problem.
5. **One calibration per experiment, not per segment.** A data set split into segments (one per angle, one per figure)
   has one accelerator calibration: share it with `shiftgroup <label>` on every member and put any `shiftfn` on the
   group's first segment. Per-segment energy-dependent terms fit angle-by-angle structure instead: in one model b
   freed per segment gained 91 in χ² for 10 parameters, but b changed sign between angles of one experiment (−0.016
   at 54.7° against +0.006 at 90°), which no calibration does. Conversely, grouping digitized sets can cost a lot
   (+1,913 in a frozen group scan: different figures, different digitization offsets); group by what shares a
   calibration, and test it. When per-segment constants are real but the energy-dependent calibration is one, use
   `shiftfngroup <label>`: members apply the first one's `shiftfn` and keep their own constants. The start-up log
   prints, per function, the condition number of its free terms over the group's energies, and warns about a term
   that only duplicates a (√x on a single-energy angular distribution).
