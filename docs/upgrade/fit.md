# The engine's fit as a function

One fit for the program and for a pyazr session: a table of fitters, one function per fitter with one signature, a
`FitResult` record, the error analysis as functions on it, progress as data and a stop request owned by the caller.
The program's run becomes a sequence of named stages. Decision 8 of `pyazr3.md` section 5 (`fit_with_engine`);
findings 2, 3 and 20 of `web_ui.md` section 12.

## 1. The run as it was (`src/AZUREMain.cpp`, 940 lines, `AZUREMain::operator()()` lines 70-941)

- 71-82: fill the compound nucleus (`.azr` levels and pairs) and write its check.
- 84-128: fill the data (data run: data records and files) or make the extrapolation points; for a rate run check the
  rate's pairs (`rateParams`). Writes the `WorkingCopy` and the data check.
- 130-138: `CNuc::Initialize` (channels, pathways, capture-A_y terms).
- 140-147: the parameter vector (`AZUREParams` from `CNuc` and `WorkingCopy`), `param.par`/`param.sav` with
  `USE_PREVIOUS_PARAMETERS`.
- 148-157: refine adaptive extrapolation grids at these parameters (a new `WorkingCopy`).
- 158-174: match the prior covariance (`priorCovariance=`) onto the vector.
- 175-178: write `param.par` (fresh run).
- 180-237 (data runs): EC amplitude cache, `WorkingCopy::Initialize`, `ParameterLimitsManager` (`<parameterSettings>`:
  limits, nuisance priors), the radiative-width warning, the `Results` buffer, `AZURECalc`, `ReplaceDiagonalPriors`,
  the fixed mask.
- 239-363: lambdas: the band covariance capture, parameters at a limit, the call-limit text, `fitStatus.out`, MIGRAD's
  report.
- 365-749: the fit: NLopt (375-498), GSL-LM (500-522), LM (523-546), MIGRAD (547-746) with the refreeze loop (577-598),
  HESSE (599-623), MINOS and `covariance_matrix.out` (624-741). Reads `PERFORM_FIT`, the fitter bits,
  `refreezeCovariance`, `chiVariance`, `CALCULATE_COVARIANCE_BAND`, `PERFORM_ERROR_ANALYSIS`; writes `param.sav`,
  `fitStatus.out`, `covariance_matrix.out`, `param.errors`.
- 750-782: the final evaluation; χ²/N per segment and the total on standard output.
- 784-881: the band: χ² scaling (`SCALE_COVARIANCE_BY_CHI2`), `covariance.dat`, `BuildBandData`.
- 882-900: output files, scattering lengths, `<modelData>`.
- 901-922: the reaction rate (rate runs).
- 924-938: transform out, `parameters.out`, EC cache cleanup.

## 2. The fitters as they were

- MIGRAD (numerical gradient), the default: `MnUserParameters` (limits, fixed flags); call limit `--max-calls` or
  10^9; Minuit's tolerance 0.1 (EDM goal 2e-4 Up). Prints `Iteration: N Chi-Squared: X` every 10th evaluation of the
  run and "MIGRAD converged ..." or "WARNING: MIGRAD did NOT converge ..."; `fitStatus.out` (MIGRAD, status, valid,
  edm, edm_goal, calls, call_limit). After it: the refreeze loop when a data covariance has relvector terms; HESSE
  with `--covariance-band` (none with no free parameter: Minuit2 crashes in MnHesse); MINOS with menu 4.
- MIGRAD with the analytic gradient (`--use-gradient`, also set for THM data): as MIGRAD, plus "Using analytic
  gradient ..." and the G2 seed line (`EnableG2`).
- Levenberg-Marquardt (`--use-lm`): free values, limits by projection, diagonal and nuisance priors as rows; 200
  iterations, gradient 1e-8, relative step 1e-9; `LM iteration: N  Chi-Squared: X` per accepted step, curvature
  warnings; no `fitStatus.out`. Without an analytic Jacobian: MIGRAD (numerical; no HESSE, no MINOS). Covariance
  (JᵀJ + priors)⁻¹.
- GSL trust region (`--use-gsl-lm`): as LM; 200 iterations, xtol = gtol = ftol = 1e-9; `GSL-LM iteration: ...` and
  "GSL trust-region solver ... after N iterations.".
- NLopt (`--use-nlopt`, builds with `USE_NLOPT`): free values, ±1e10 where unbounded; `--max-calls` or none;
  ftol_rel = xtol_rel = 1e-8; success, failure or call-limit lines; "Error analysis is not available ..." with MINOS;
  `fitStatus.out` (NLopt, status, valid, calls).

Precedence when several switches are given: NLopt, GSL-LM, LM, MIGRAD with gradient, MIGRAD. Intermediate files
(`param.fit`, `AZUREOut_*`, `parameters.out`) every 100 evaluations (MIGRAD, NLopt) or every 10 iterations (LM,
GSL-LM). `Config::stopFlag` (set by nobody since the Qt GUI left) made `AZURECalc::operator()` return 0 during a fit;
the data fill (`WorkingCopy.cpp`) and the external-capture set-up (`ExternalCaptureSetup.cpp`) skip their remaining
points when it is set.

## 3. Design

**The table** (`src/Fitters.cpp`, `include/Fitters.h`): one row per fitter, in the program's precedence order:
`name` (what `fit_with_engine(method=)` and `pyazr.kinds()["fitters"]` say), `option` (the program's switch, "" for
the default), `gradient`, `residualRows`, `bounds`,
`covariance` (whether it yields one), `minos` (whether MINOS can follow), `progressLabel`, `doc`, and `run`, its
function. The program's choice (`FitterOfConfig`) is the row `Config::fitter` names, the earliest row of the switches
given, else the row with no switch (section 5; until 2026-11-17 a `configBit` column and a bit of `Config::paramMask`
per row).
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

## 4. As built

**Files.** `src/AZUREMain.cpp` 940 → 77 lines: `AZUREMain::operator()()` calls `SetUpRun`, `SetUpObjective`
(`src/RunSetUp.cpp`, 181), `RunFitStage` (`src/ProgramFit.cpp`, 285: `ProgramProgress`, `fitStatus.out`,
`covariance_matrix.out`, HESSE for the band, MINOS, `param.sav`), `FinalEvaluation` and `WriteRunOutputs`
(`src/RunOutputs.cpp`, 182), the reaction rate and the transformation out; `ProgramRun` (`include/ProgramRun.h`)
carries what a stage hands to the next. The fitters: `src/Fitters.cpp` (the table, 39), `src/FitMigrad.cpp` (188),
`src/FitLeastSquares.cpp` (439: LM, the free parameters with their priors, the covariance, the row driver with the
MIGRAD fallback), `src/FitGslTrustRegion.cpp` (242), `src/FitNLopt.cpp` (130, compiled with `USE_NLOPT` only),
`src/FitErrors.cpp` (49), `src/FitEvaluations.cpp` (71); headers `FitControl.h`, `Fitters.h`, `FitEvaluations.h`,
`FitErrors.h`, `LeastSquares.h`. `src/AZURECalc.cpp` 1948 → 1418 lines (the least-squares drivers moved out).
The session: `api/src/AZUREAPIFit.cpp` (103), `api/src/pybind_fit.cpp` (140), `pyazr/session_fit.py` (167).

**The table** (precedence order): `nlopt` (`--use-nlopt`), `gsl-lm` (`--use-gsl-lm`), `lm` (`--use-lm`),
`migrad-gradient` (`--use-gradient`), `migrad` (default). Columns as in section 3; `progressLabel` names the
iteration lines of LM and GSL-LM, which the program's printer writes from the record.

**Progress.** `AZURECalc::Observe(const EvaluationObserver *)`: `operator()` hands every evaluation (number, objective,
vector, MINOS or not) to the observer. `ProgramProgress` is the observer for the whole data run and prints the
program's lines (every 10th evaluation outside MINOS, flush on every one; LM and GSL-LM per iteration) and calls the
checkpoint (`WriteIterationOutput`) every 100 evaluations (`HandOverEvaluation`). While a fitter runs,
`FitEvaluations` replaces it, counts, keeps the best point and forwards to the same sink, so the lines are those of
before. `FitProgress::stage` is `fit`, `hesse` or `minos` (the session's error analysis has its own count and best).

**Stop.** `EvaluationObserver::stop`: `operator()` returns 0 at once during a fit (before any work), `Gradient` returns
0 and `G2` 1 (a finite zero step), so MIGRAD ends within its sweep; LM checks before every iteration and every damped
trial; GSL's residual functions fail (`GSL_EBADFUNC`) and the driver returns with its last accepted point; NLopt is
told `nlopt_force_stop`. A stopped fit returns the best point seen with status `stopped` and no covariance.
`Config::stopFlag` lost its reader in `AZURECalc`; its readers in the data fill (`WorkingCopy.cpp`) and the
external-capture set-up (`ExternalCaptureSetup.cpp`) are unchanged (nobody sets the flag; making them read a request
needs it passed through `WorkingCopy::Fill` and `CNuc::Initialize`, not done).

**What a session's fit changes.** The evaluation counter of the session's working copy, the fit and error-analysis
flags (restored), the session's buffer (the end point's model values, as after `evaluate`), and on a MIGRAD row with
a relvector data covariance the frozen model (re-frozen at the fitted point, as in the program). The session's
parameter values are not changed. Evaluating the start after a stopped fit gives the numbers of before, bit for bit
(`tests/pyazr/fit_engine_control_test.py`).

**Not changed on purpose.** No tolerance, call limit, fallback or message. The program's at-limit report after a
MIGRAD fit without a valid HESSE names the parameters' values before the fit (as before; a defect, listed for the
owner). LM and GSL-LM still write no `fitStatus.out`.

**The session against the program** (`tests/pyazr/fit_with_engine_test.py`, one thread). The end point (param.sav's
8 digits), the printed total (12 digits), the band block of the covariance, covariance_matrix.out and MINOS agree
exactly; the errors of param.sav agree to its 8 digits (the largest relative difference, 3e-8, is the rounding of
the print; the test compares the parameters that are not level parameters to 1e-6 relative). The
errors are taken as the program writes them (`MnUserParameters::Error`): after LM the square roots of the diagonal of
the returned covariance differ from them by up to 6e-4 relative for some of those other parameters on T2 (the areal
density; norms about 1e-5); why was not found. Two LM fits from the same start in one session also end 4e-15 apart in the
objective (the pooled working copies keep state between evaluations: applied energy shifts and sub-point grids).
Not found here; for the owner.

## 5. The chosen fitter is a row's name, not a bit (design, 2026-11-17)

The closing plan the PI approved: "a fitter without a `Config.h` bit".

**What was.** Each row of `Fitters()` owned a bit of `Config::paramMask` (`FitterRow::configBit`:
`USE_NLOPT_MINIMIZER` 1<<17, `USE_ANALYTIC_GRADIENT` 1<<19, `USE_LM_MINIMIZER` 1<<20, `USE_GSL_LM_MINIMIZER` 1<<23; the
default row 0). A new fitter needed the next free enumerator of `Config::ParameterFlags` in `include/Config.h`, which
almost every file of the engine, the binding and the tests includes: the recipes trial (`recipes2_report.md`, R3-3,
C-7) rebuilt nearly everything for one table row. The bit carried nothing but "which row".

**The rule of precedence, as it was** (`FitterOfConfig`, `src/Fitters.cpp`): the parser ORs the bit of every fitter
switch it meets (`ApplyFitterSwitch`); the fit runs the *first row of the table* whose bit is set, else the row without
a switch (`migrad`). So with several switches the row earliest in the table wins (nlopt, gsl-lm, lm, migrad-gradient),
whatever their order on the command line or in the options file, and a repeated switch changes nothing. The only other
writer: a project with THM data and the gradient bit unset gets the bit and the line "THM segments present: enabling
the analytic gradient ..." (`src/AZURE2.cpp`, after the project is read).

**The holder.** One member `std::string Config::fitter`: the `name` of the chosen row (`FitterRow::name`), empty when
no switch chose one; `Config::Reset` empties it. Empty means the row without a switch (`migrad`), as before. The
row's key, not its index, because (a) the rows differ by build (`nlopt` is a row only with `USE_NLOPT`), so an index
names a different fitter in another build, while the name is the same everywhere; (b) the name is already the key
of every other reader: `fit_with_engine(method=)`, `pyazr.kinds()["fitters"]`, `fitStatus.out`, the web UI's
`engine_<name>`; (c) a name that is no row is caught (`FitterOfConfig` throws `std::logic_error`), an index past the
end or of a reordered table is not. The cost (a string comparison at the parse and once per fit) is nothing.

**How the parser sets it.** `ApplyFitterSwitch(configure, option)` finds the row of the switch and calls
`ChooseFitter(configure, row.name)`, which sets `fitter` to the row unless the row already chosen comes earlier in the
table. After any sequence of switches `fitter` is therefore the earliest row of those given: the old rule, now applied
when the switch is read instead of when the fit starts. `FitterOfConfig` returns the row named `fitter`, or the row
without a switch. The THM rule becomes `ChooseFitter(configure, "migrad-gradient")` under the same precedence (a
switch of an earlier row still wins, as the bit did). Its line is printed when the chosen row is not migrad-gradient.
This is the old condition, "the gradient bit is unset", except in one case: on a THM project with both a switch of an
earlier row and `--use-gradient` (e.g. `--use-lm --use-gradient`), the line is now printed where it was not before.
The fitter that runs is the same (lm), and so are all the numbers. A single member cannot tell "`--use-gradient` was
given too" apart from "it was not". No test or project of the repository gives two fitter switches.

**Readers of the old bits** (grep of `src/`, `api/`, `include/`, `tests/`, `pyazr/`, `io/`, `experiment/`, before the
change; declarations in `include/Config.h` not counted):

| Bit | Uses | Where |
|---|---|---|
| `USE_NLOPT_MINIMIZER` | 1 | the row (`src/Fitters.cpp`) |
| `USE_GSL_LM_MINIMIZER` | 1 | the row |
| `USE_LM_MINIMIZER` | 1 | the row |
| `USE_ANALYTIC_GRADIENT` | 3 + 2 | the row; the THM rule's test and set (`src/AZURE2.cpp`); check 46 (2) |
| `configBit` / `paramMask` through it | 2 | `FitterOfConfig` (read), `ApplyFitterSwitch` (set) |

No reader in `api/`, `pyazr/` or the tests other than check 46. A session never sets a bit: `fit_with_engine` names
the row. All of them now go through `ChooseFitter` / `FitterOfConfig`.

**The retired bits.** Positions 17, 19, 20 and 23 of `ParameterFlags` stay unused, and a comment in `Config.h` says so.
Nothing writes `paramMask` as a number: no `<<` of it, no file, no `--describe` line (`--describe` and pyazr's
`describe()` print the options of the calculation one by one from `OptionStorages()`, none of them a fitter bit, and the
fitter table from the rows, which never exposed `configBit`), no options-file line (that file holds flags), no
`<lastRun>` writer in the engine. So no saved or printed form depends on these positions. They are still left empty,
so that a mask from an old debugger session or a log keeps its meaning.

**Check 46 and the recipe.** The row regex loses the `configBit` column. The clause "every minimizer bit of
`src/AZURE2.cpp` is a row's" becomes: `include/Config.h` declares no fitter bit (no enumerator `USE_*_MINIMIZER` or
`USE_ANALYTIC_GRADIENT`), and `Config::fitter` is assigned nowhere but `src/Fitters.cpp` (`ChooseFitter`) and
`Config::Reset`. So a fitter is chosen only through the table, and a bit for one cannot come back. In the recipe
"A fitter", step 2 drops the `configBit` column and `Config.h`: a new fitter is its function file, its row, its CMake
line, its docs and its tests.

**Test.** `tests/upgrade/fit_choice_test.cpp` (ctest `upgrade_fit_choice`). With no switch: the row without a switch.
For every non-empty set of switch rows, in every order, each switch given once and given twice: the row chosen is the
earliest of the set in `Fitters()`. The test computes that answer itself from the table's order, not through
`ChooseFitter`. Also: the THM choice under an earlier switch and alone, `Reset`, and an unknown switch, which returns
false and leaves the choice alone. The program's eight fits of `verify_branch.sh` (five fitters) end where they ended.
