# The engine's fit as a function

One fit for the program and for a pyazr session: a table of fitters, one function per fitter with one signature, a
`FitResult` record, the error analysis as functions on it, progress as data and a stop request owned by the caller.
The program's run becomes a sequence of named stages. Decision 8 of `pyazr3.md` section 5 (`fit_with_engine`);
findings 2, 3 and 20 of `web_ui.md` section 12.

## 1. The run as it was (`src/AZUREMain.cpp`, 940 lines, `AZUREMain::operator()()` lines 70-941)

| Lines | Stage | Reads | Writes |
|---|---|---|---|
| 71-82 | fill the compound nucleus; its check | `.azr` levels, pairs | `CNuc`; checks file |
| 84-128 | fill the data (data run) or make the extrapolation points; or check the rate's pairs | `.azr` data records, data files; `rateParams` | `WorkingCopy`; data check |
| 130-138 | `CNuc::Initialize` (channels, pathways, capture-A_y terms) | `Config` | `CNuc` |
| 140-147 | the parameter vector (`AZUREParams` from `CNuc`, `WorkingCopy`); `param.par`/`param.sav` with `USE_PREVIOUS_PARAMETERS` | parameter file | `AZUREParams` |
| 148-157 | refine adaptive extrapolation grids at these parameters | `<extrapolations>` | a new `WorkingCopy` |
| 158-174 | match the prior covariance (`priorCovariance=`) onto the vector | the file | `PriorCovariance` on the data |
| 175-178 | write `param.par` (fresh run) | | `param.par` |
| 180-237 | data runs: EC amplitude cache, `WorkingCopy::Initialize`, `ParameterLimitsManager` (`<parameterSettings>` limits and nuisance priors), the radiative-width warning, the `Results` buffer, `AZURECalc`, `ReplaceDiagonalPriors`, the fixed mask | `Config::paramMask`, `<parameterSettings>` | `AZURECalc` |
| 239-363 | lambdas: band covariance capture, parameters at a limit, call limit text, `fitStatus.out`, MIGRAD's report | `maxCalls` | |
| 365-749 | the fit: NLopt (375-498), GSL-LM (500-522), LM (523-546), MIGRAD (547-746) with the refreeze loop (577-598), HESSE (599-623), MINOS and `covariance_matrix.out` (624-741) | `PERFORM_FIT`, the fitter bits, `refreezeCovariance`, `chiVariance`, `CALCULATE_COVARIANCE_BAND`, `PERFORM_ERROR_ANALYSIS` | `param.sav`, `fitStatus.out`, `covariance_matrix.out`, parameter errors |
| 750-782 | the final evaluation; χ²/N per segment and the total on standard output | | `Results` |
| 784-881 | the band: χ² scaling, `covariance.dat`, `BuildBandData` | `SCALE_COVARIANCE_BY_CHI2` | `covariance.dat` |
| 882-900 | output files, scattering lengths, `<modelData>` | | `AZUREOut_*`, `chiSquared.out`, ... |
| 901-922 | the reaction rate (rate runs) | `rateParams` | rate files |
| 924-938 | transform out, `parameters.out`, EC cache cleanup | | `parameters.out` |

## 2. The fitters as they were

| Fitter | Selected by | Start, limits | Call limit, tolerance | Prints | `fitStatus.out` | Fallback, after |
|---|---|---|---|---|---|---|
| MIGRAD (numerical gradient) | default | `MnUserParameters` (limits, fixed) | `--max-calls` or 10^9; Minuit's default tolerance 0.1 (EDM goal 2e-4 Up) | `Iteration: N Chi-Squared: X` every 10th evaluation of the run; "MIGRAD converged ..." or "WARNING: MIGRAD did NOT converge ..." | MIGRAD, status, valid, edm, edm_goal, calls, call_limit | refreeze loop when a data covariance has relvector terms; HESSE with `--covariance-band` (none with no free parameter: Minuit2 crashes); MINOS with `--minos` |
| MIGRAD, analytic gradient | `--use-gradient` (also set for THM data) | as MIGRAD; `EnableG2` seeds Gauss-Newton G2 | as MIGRAD | as MIGRAD, plus "Using analytic gradient ..." and the G2 line | as MIGRAD | as MIGRAD |
| Levenberg-Marquardt | `--use-lm` | free values, limits by projection, diagonal and nuisance priors as rows | 200 iterations; gradient 1e-8, relative step 1e-9 | `LM iteration: N  Chi-Squared: X` per accepted step; curvature warnings | none | no analytic Jacobian: MIGRAD (numerical, no HESSE, no MINOS); covariance (JᵀJ + priors)⁻¹ |
| GSL trust region (`gsl-lm`) | `--use-gsl-lm` | as LM | 200 iterations; xtol = gtol = ftol = 1e-9 | `GSL-LM iteration: N  Chi-Squared: X`; "GSL trust-region solver ... after N iterations." | none | as LM |
| NLopt | `--use-nlopt` (builds with `USE_NLOPT`) | free values; ±1e10 where unbounded | `--max-calls` or none; ftol_rel = xtol_rel = 1e-8 | success / failure / call-limit lines; "Error analysis is not available ..." with `--minos` | NLopt, status, valid, calls | none |

Precedence when several switches are given: NLopt, GSL-LM, LM, MIGRAD with gradient, MIGRAD. Intermediate files
(`param.fit`, `AZUREOut_*`, `parameters.out`) every 100 evaluations (MIGRAD, NLopt) or every 10 iterations (LM,
GSL-LM). `Config::stopFlag` (set by nobody since the Qt GUI left) made `AZURECalc::operator()` return 0 during a fit;
the data fill (`WorkingCopy.cpp`) and the external-capture set-up (`ExternalCaptureSetup.cpp`) skip their remaining
points when it is set.

## 3. Design

**The table** (`src/Fitters.cpp`, `include/Fitters.h`): one row per fitter, in the program's precedence order:
`name` (what `fit_with_engine(method=)` and `pyazr.kinds()["fitters"]` say), `option` (the program's switch, "" for
the default), `configBit` (the `Config::paramMask` bit the switch sets), `gradient`, `residualRows`, `bounds`,
`covariance` (whether it yields one), `minos` (whether MINOS can follow), `progressLabel`, `doc`, and `run`, its
function. The program's choice (`FitterOfConfig`) is the first row whose bit is set, else the row with no switch.
Consistency check 46: every row's switch is parsed by the program and named in `reference/command_line.rst`, and the
fitter table of the user guide lists every row.

**One function per fitter**, `FitResult run(AZURECalc &objective, AZUREParams &params, const FitControl &control)`.
`params` holds the start, the limits and the fixed flags, and on return the end point (with the errors the fitter
sets). `FitControl` (`include/FitControl.h`): `maxCalls` (0: none), `progress` (the sink), `checkpoint` (the
program's intermediate files; empty in a session), `stop` (an `std::atomic<bool>` the caller owns), `log` (the
stream the fitter's lines go to), `refreezePasses`, `callLimitSource` (" (--max-calls)" in the program).
`FitResult`: fitter, full vector `x`, `objective` (the fitter's value at `x`), `status` (`converged`, `call_limit`,
`invalid`, `stopped`, `failed_<reason>`), `converged`, `edm`/`edmGoal` (MIGRAD; NaN otherwise), `calls`, the free
indices, the `covariance` over them (row-major; empty when none), `minos` (per full index), and for the program the
`fitStatus.out` fields, MIGRAD's `FunctionMinimum` (for HESSE and MINOS) and the least-squares band block.

**The error analysis** (`src/FitErrors.cpp`): `RunHesse(objective, result, log)` refines MIGRAD's covariance (none
with no free parameter) and says whether it is valid; `RunMinos(objective, result, up, log)` fills `minos`; the
band's capture stays with the program's band stage.

**Progress as data**: `FitProgress` {kind (evaluation or iteration), label, evaluation (the run's evaluation number),
iteration, calls, objective, best, measure (NaN where none), `x` (the record's point, valid during the call),
errorAnalysis}. `AZURECalc::operator()` hands every evaluation to its `EvaluationObserver`; a fitter hands its
iterations to `control.progress`. The program's sink (`ProgramProgress`, in the fit stage) prints today's lines byte
for byte and calls the checkpoint every 100 evaluations; it is the observer for the whole run, so the evaluations of
HESSE, the final evaluation and the band print as before.

**Stop**: `control.stop`, checked where `Config::stopFlag` was checked (`AZURECalc::operator()` during a fit; now
also `Gradient`/`G2`, the LM loops, the GSL residual, NLopt via `nlopt_force_stop`). A stopped evaluation returns 0
at once (as `stopFlag` did), so MIGRAD ends within its running gradient sweep; the fit returns the best point seen
(tracked by the fit from the evaluations), status `stopped`, no covariance. `Config::stopFlag` keeps its readers in
the data fill and the external-capture set-up (set-up loops of a session's build, not of a fit; one mechanism there
needs the request passed through `WorkingCopy::Fill` and `CNuc::Initialize`: listed, not done); its reader in
`AZURECalc` is replaced by the fit's request.

**The stages** (`AZUREMain::operator()()` calls them in order): `SetUpRun` (`src/RunSetUp.cpp`: lines 71-178 and the
data-run set-up 180-237), `RunFitStage` (`src/ProgramFit.cpp`: the fit, HESSE, MINOS, `fitStatus.out`,
`covariance_matrix.out`, `param.sav`), `FinalEvaluation`, `WriteRunOutputs` (`src/RunOutputs.cpp`: band, output
files, scattering lengths, `<modelData>`), `RunReactionRate`, `FinishRun`. No file above 500 lines.

**State** (`ARCHITECTURE.md` §9): nothing process-wide is added. A fit changes, in the objects it is given: the
working copy's evaluation counter and fit/error-analysis flags (restored), the evaluation buffer (the end point's
model values), and, on the MIGRAD rows with a relvector data covariance, the frozen covariance model (re-frozen at the
fitted point, as in the program). A session's parameter values are not changed by a fit.

**The session** (`AZUREAPI::FitWithEngine`, `api/src/AZUREAPIFit.cpp`): builds the program's `AZUREParams` at the
session's values (x0 for the free ones), applies `<parameterSettings>` as the program does, makes an `AZURECalc` on the
session's compound, data and buffer, runs the row's function, then HESSE (MIGRAD rows) and MINOS when asked, and one
evaluation at the end point. The binding releases the GIL; `pyazr/session_fit.py` gives
`fit_with_engine(method, minos, x0, max_calls, progress, stop)`.
