# The objective's terms on arrays: design

The fit minimises χ²(θ) = Σ_t χ²_t(θ), the ordered sum of the thirteen term objects of `src/Objective.cpp` (T1-T11 of
`objective.md`; T8 is five objects). Each term today reads the engine's objects (`WorkingCopy`, `Segment`,
`TargetEffect`, `CompoundSystem`, `LevelState`) and does its arithmetic in the same function. This note says, term by
term, which part of that arithmetic needs no engine object and can be a function of plain arrays in `inference/`, with
its value, its gradient contribution, its least-squares rows and its curvature for Minuit's G2; which terms gain
nothing from such a move and stay; whether a `ParameterVector` class should stand beside `AZUREParams`; and how the
`chiSquared.out` writer takes every total from the terms. Physics and every printed number stay identical: code is
moved with its order of floating-point operations.

Status: design only, 2026-11-16, branch `feature/upgrade-7-objective-arrays` at `2b3352a4`. Nothing is built. Line
numbers are those of that commit. "Read" marks what was verified in the code; "inferred" what was not run.

## 1. The cut

Every term becomes three pieces:

1. **The gather** (engine): a walk over the engine's objects that fills a plain struct with numbers in physical units
   and the full-vector index of each parameter. The walk keeps every rule of eligibility (which parameter carries a
   prior, which level is in the R matrix, which channel is a particle channel).
2. **The arrays** (`inference/`): a struct of `std::vector<double>` and `int`, no engine type.
3. **The arithmetic** (`inference/`): functions of the struct that return the value, add the gradient, add the
   curvature, append residual rows. Tested in `inference/tests/` against values by hand and central differences.

The term class stays in `src/Objective.cpp` (recipe "A term of the objective"): its `Add` is the gather and one call.

Five rules keep every output identical.

- **R1, order.** A function adds into the caller's running sum element by element when the term does so today
  (`Segment::AddShiftPenalty` adds each shift prior straight into the objective's running χ²), and forms a subtotal
  from 0 when the term does so today (`WorkingCopy::ScreeningPenalty` sums from 0, then `chi2 +=` the subtotal).
- **R2, the expression.** A diagonal prior's derivative is written in four forms today, each with its own rounding:
  `GaussianPriorGradient` 2(x − x₀)/σ² (adjoint gradients, `AZURECalcGradient.cpp` 207-229, 545-549;
  `AZUREAPIGradient.cpp` 225-228), the row (x − x₀)/σ with Jacobian 1/σ (`AZUREAPIObjective.cpp` 75-89), the normal
  equations (x − x₀)·σ⁻² and σ⁻² (`FitLeastSquares.cpp` 55-70, 97, 306-308), and the GSL row √(σ⁻²)·(x − x₀)
  (`FitGslTrustRegion.cpp` 64-65, 89). Each path keeps its form; only the list it walks changes.
- **R3, the value read.** The value of a T8 prior is formed on the value the working copy holds, clamped into the
  parameter's layout range by the hand-over (`HandOver.cpp` 52-61, 82-91; `ApparatusParameters.cpp` 125-134); its
  gradient on the raw component p[i] (`objective.md` §6, "Where the fit clamps"). The arrays carry both where both
  are read.
- **R4, the same set.** Where two walks of today select the same priors by two written rules, the slice that makes
  them one walk first shows on every project of the tree that the two sets are equal (section 6, item 1).
- **R5, no new state.** The arrays are locals of the call or the caller's; no `mutable`, no `static`, no cache. A list
  that could be built once at set-up (T5) is built per call until a measurement shows the cost (recipe "Add a cache").

A term **gains nothing** from a move when the arithmetic left after the gather is one call of a function that is
already shared (`GaussianPriorChi2`, `PriorCovariance::Chi2`), or when its arrays would be a copy of the model
buffer `Results::model` made at every evaluation, a second data path (`inference.md` §2).

## 2. Term by term

Units: the unit of a parameter's kind is the column `unit` of `ParamKindTable()`. "Paths" are: V the value
(`Objective::Chi2`), A the adjoint gradient (`AZURECalc::Gradient`, `AZUREAPI::CalculateChi2GradRWA`), J the
Jacobian paths (`GradientFromJacobian`, `G2`), R the residual rows (`AppendModelRows`, `ObjectiveRowsRWA`, the LM and
GSL drivers), S the sampler (`AZURECalcMCMCPriors.cpp`).

### T1, the data χ² of a segment

- *Today* (read): `SegmentTerms::Add`, `Objective.cpp` 27-80; the point loop 57-58 calls `azure2::PointResidual`
  (`ModelPass.cpp` 366-369), which reads the point's model value from `Results::model`, the Jacobian
  dΩ_lab/dΩ_c.m. from the segment's state and the datum; the residual is `DataResidual` (`ModelPass.h` 76-89):
  m = f/J, (m − y n)²/(dy n)².
- *Arithmetic without an engine object*: `DataResidual` (4 lines), already a plain inline function.
- *Arrays*: they would be (f_i, J_i, y_i, dy_i) per point, a copy of `Results::model` at every evaluation.
- *Gradient*: A, the adjoint seed 2 r/(err² J) (`AZURECalcGradient.cpp` 186); J and R, the data rows.
- **Verdict: stays.** The point loop is the walk; the arithmetic is one shared inline function. Moving `Residual` and
  `DataResidual` into `inference/` would move 14 lines and change no test.

### T1a, a profiled norm (a THM segment with a free norm)

- *Today* (read): `Segment::ProfileNormChiSquared`, `Segment.cpp` 243-269, called by `SegmentTerms` (55), the writer's
  pass (`DataOutputFiles.cpp` 111), `AZUREAPIChi2.cpp` 110, `AZUREAPIGradient.cpp` 175, `ShiftedFold.cpp` 249,
  `ThmCotangents.cpp` 72. The loop gathers m_i = f_i/J_i, d_i, e_i; the arithmetic is S_mm, S_md, S_dd,
  n* = S_mm/S_md, χ² = max(0, S_dd − S_md²/S_mm), and the degenerate branch (S_md ≤ 0 or S_mm ≤ 0: n = 1,
  χ² = Σ (m_i − d_i)²/e_i²).
- *Proposed* (`inference/include/azure2/inference/LinearProfile.h`):

      // m: the model in the frame and unit of the data (b, b/sr, ...); d, e: the data and their errors, same unit;
      // points with e = 0 carry no weight.  Returns chi2 at the optimum; norm: n* (dimensionless; 1 if degenerate).
      double ProfiledNormChi2(const std::vector<double> &m, const std::vector<double> &d,
                              const std::vector<double> &e, double &norm);

  The degenerate branch must form each point's term as `DataResidual(f, J, y, dy, 1.0).Chi2()` does:
  ((f/J) − y·1)²/(dy·1)². The gather passes m = f/J computed as today; the products by 1.0 are exact.
- *Gradient*: none of its own (the norm is not a parameter; the envelope theorem, `ShiftedFold.h` 216-221).
- **Verdict: moves**, with T1b (one header, one test).

### T1b, a THM experiment (one norm and a polynomial background profiled over its segments)

- *Today* (read): `WorkingCopy::ProfileThmGroup`, `ThmNorm.cpp` 671-693: `GatherThmGroup` (650-) fills the arrays
  m, d, e, E_c.m.; `SolveThmProfile`, `ThmProfile::Residual`, `Background`, `Reported`, `InvertSPD` and
  `ThmProfileDerivative` (`ThmExperiment.cpp` 467-667, about 200 lines) are already functions of those arrays;
  `ThmProfile` (`include/ThmProfile.h`, 47 lines) is plain data. No unit test of `SolveThmProfile` exists under
  `tests/` (grep for the name finds none).
- *Proposed*: `ThmProfile`, `SolveThmProfile`, `InvertSPD`, `ThmProfileDerivative` move unchanged into
  `LinearProfile.h` / `LinearProfile.cpp`; `Results.h` includes the header (it holds `std::vector<ThmProfile>`);
  the engine names the type through `include/InferenceTypes.h`. Units: m, d, e in the unit of the data, E_c.m. in
  MeV, s = 1/n dimensionless, a_k in the unit of the data times MeV^−k.
- *Gradient*: `ThmProfileDerivative` (the rows' Jacobian through the profile), already on arrays.
- **Verdict: moves.** The arithmetic is already on arrays; the gain is a unit test of a closed form that has none
  (the line s m + a₀ + a₁E by hand, a singular normal matrix, s ≤ 0, the derivative against central differences).

### T2, the data covariance

- *Today* (read): `WorkingCopy::CovarianceChi2Correction`, `DataCovarianceGroups.cpp` 152-182: per group the pulls r,
  `Group::EnsureFactored(r)` (freezes the factor at the first evaluation when the group has relvector terms), the
  whitening z = L⁻¹ r, q = z·z, the correction q − Σ r², the seeds L⁻ᵀz into `Results`, the per-segment partition
  Σ r_a (M⁻¹r)_a into `Results::objective.segmentChi2`.
- *Arithmetic without an engine object*: the whitening, already in its own class (`DataCovariance.cpp`, 291 lines;
  `include/DataCovariance.h` includes `experiment/CovarianceFile.h`, a mutex and an atomic for the freezing).
- **Verdict: stays.** The arrays would be the pulls of the group's points (a copy of the model buffer); the factor is
  frozen lazily inside an evaluation, which `inference/` should not own. A move belongs with a decision on when the
  factor is frozen (section 6, item 9).

### T3, the normalization prior of a segment

- *Today* (read): `Objective.cpp` 60-69: for a segment outside a THM experiment and not profiled, with
  `IsVaryNorm`, not replaced by the prior covariance and σ_n = n₀ σ%/100 ≠ 0, `chi2 += GaussianPriorChi2(n, n₀, σ_n)`
  on n = `State(segment).norm`. The same eligibility is written a second time in `ForEachDiagonalPrior`
  (`Objective.h` 210-216, by `pmap.NormIndex(s) >= 0`, which excludes profiled norms: `AZUREGrad.cpp` 305-309), a
  third time in the writer (`DataOutputFiles.cpp` 321-325) and a fourth in the adjoint gradient
  (`AZURECalcGradient.cpp` 216-230, which also skips norm-group followers and THM segments). A follower's σ is set to
  0 (`SegmentParameters.cpp` 68-76), so it carries no prior in any of them (read).
- *Arithmetic*: `GaussianPriorChi2` and its derivatives, already in `inference/`.
- *Proposed*: the eligibility once, in one engine gather (section 2.1 below); n dimensionless.
- *Gradient*: A analytic `GaussianPriorGradient(n, n₀, σ_n)`; J `GaussianPriorGradient` and `GaussianPriorCurvature`;
  R one row each (forms of R2); S a normalised Gaussian log-prior.
- **Verdict: the arithmetic stays where it is (it is one shared function); the four eligibility rules become one
  gather.**

### T4, the energy-shift priors of a segment

- *Today* (read): `Segment::AddShiftPenalty`, `SegmentParameters.cpp` 137-148, adds into the caller's running sum, in
  the order: the constant shift ΔE (MeV), the √E coefficient b (MeV^½), each shift-function parameter a_j (the unit
  the function's expression gives it), each when varied, σ ≠ 0 and not replaced. Callers: `SegmentTerms` (72),
  `SegmentLocalChi2` (`ShiftedFold.cpp` 237), the writer (`DataOutputFiles.cpp` 361, 508, 533).
  `ShiftPenaltyParts` (150-160) gives the three parts for pyazr's breakdown. `ForEachShiftPenalty` (`Objective.h`
  125-143) writes the eligibility a second time.
- *Gradient*: A by a central difference of `SegmentLocalChi2`, which contains the prior (`AZURECalcGradient.cpp`
  313-367), or of the whole χ²; J and R analytic from `ForEachShiftPenalty`.
- **Verdict: as T3** (one gather, the arithmetic already shared).

### 2.1 The one gather of the diagonal priors (T3, T4, T8, T7's per-nucleus priors)

Proposed (`inference/include/azure2/inference/DiagonalPriors.h`):

    // One Gaussian prior ((x - x0)/sigma)^2 on one parameter, every number in the unit of the parameter's kind.
    struct DiagonalPrior {
      int index;      // full-vector index (ParamIndexMap), -1 where the caller has no map
      double value;   // the value the objective's chi2 reads: the working copy's, clamped where the hand-over clamps
      double centre;  // x0
      double sigma;   // > 0
      int group;      // the subtotal it belongs to in an output: the segment (T3, T4), the apparatus row (T8), else 0
      int part;       // T3 norm 0; T4 constant 1, sqrt(E) 2, function 3; T8 and T7 0
    };
    using DiagonalPriors = std::vector<DiagonalPrior>;
    void AddPriorsChi2(const DiagonalPriors &, double &chi2);               // chi2 += r_k^2, k in list order (R1)
    double PriorsChi2(const DiagonalPriors &, int group = kAny, int part = kAny);  // from 0, in list order (R1)
    void AddPriorsGradient(const DiagonalPriors &, const std::vector<double> &p,
                           std::vector<double> &grad);                     // grad[i] += 2 (p[i] - x0)/sigma^2 (R3)
    void AddPriorsCurvature(const DiagonalPriors &, std::vector<double> &g2); // g2[i] += 2/sigma^2

All four call the functions of `GaussianPrior.h`, so the expressions are those of today.

Engine gathers (in `src/Objective.cpp` beside the terms, or a file `src/DiagonalPriors.cpp` of their own):

    // The norm and shift priors of one segment, in AddShiftPenalty's order after the norm (T3, T4).
    void AppendSegmentPriors(const DataSet &, int segment, const SegmentState &, const ParamIndexMap *,
                             DiagonalPriors &out);
    // The priors of one T8 kind, in the order of its parameters k = 1..N (WorkingCopy::<Kind>Penalty's order).
    void AppendTargetPriors(const WorkingCopy &, TargetPriorKind, const ParamIndexMap *, DiagonalPriors &out);
    // Every diagonal prior in ForEachDiagonalPrior's order (norms, shifts, screening, density, Doppler, thickness,
    // apparatus, abundances): the API's rows depend on that order (AZUREAPIObjective.cpp 83-85).
    DiagonalPriors DiagonalPriorsOf(const WorkingCopy &, const ParamIndexMap &);

`ForEachDiagonalPrior` and its six helpers (`Objective.h` 119-232, 114 lines) become a loop over
`DiagonalPriorsOf`; the 9 call sites (`FitLeastSquares.cpp` 47, `AZURECalcGradient.cpp` 207, 210, 213, 545, 575,
`AZURECalcMCMCPriors.cpp` 62, `AZUREAPIGradient.cpp` 225, `AZUREAPIObjective.cpp` 83) keep their lambdas (R2).

### T5, the nuisance priors (`<parameterSettings>` useAsNuisance; level parameters and the THM background)

- *Today* (read): `NuisanceTerm`, `Objective.cpp` 93-104, calls the context's `nuisance` function:
  `AZURECalc::CalculateNuisanceChiSquared` (`AZURECalcNuisance.cpp` 91-104) or `AZUREAPI::NuisanceChi2`
  (`AZUREAPIChi2.cpp` 218-226): for every full index f, `NuisancePriorOf(f)` asks the limits manager and
  `std::find`s f in the prior covariance's list, then adds `GaussianPriorChi2(p[f], x₀, σ)`. The rule is written twice
  (`AZURECalcNuisance.cpp` 70-76, `AZUREAPIChi2.cpp` 207-216), the gradient twice (`AZURECalcNuisance.cpp` 78-89,
  `AZUREAPIChi2.cpp` 228-233), the rows in `AZUREAPIObjective.cpp` 86-89 and the LM penalties in
  `FitLeastSquares.cpp` 59-70.
- *Proposed*: one engine function `DiagonalPriors NuisancePriorsOf(const ParameterLimitsManager &, const
  PriorCovariance *, const vector_r &p)` (value = p[f], index = f); the value is `PriorsChi2` (sum from 0 in f order,
  as today), the gradient `AddPriorsGradient`.
- *Gradient*: A analytic; J none (G2 is switched off when a nuisance prior is present: `AZURECalcGradient.cpp` 589,
  `HasNuisanceParameters`); R rows; S normalised Gaussian log-priors.
- **Verdict: the arithmetic stays shared; the two copies of the rule become one gather.** Cost per evaluation today:
  n_full calls of the limits manager and n_full `std::find`s over the covariance's list (not measured).

### T6, Park's wall J > 0

- *Today* (read): three walks of the R-matrix levels with three forms of the arithmetic.
  - Value: `ParkTerm::Wall`, `Objective.cpp` 127-140: Σ over levels in the R matrix with J < 0 of
    `pow(J/s, 2.0)`, s = `kParkNormScale` = 10⁻³, J = `state.ParkNorm(j, la)`.
  - Gradient: `AddParkPenaltyGradient`, `AZUREGrad.cpp` 82-103: for J < 0 (skipped when J ≥ 0 or NaN since the
    repair of 2026-11-16, section 6 item 6; before it a NaN J passed),
    factor = 2J/s², dχ²/dγ_c = factor · (−2 γ_c S'_c), dχ²/dE = factor · (−Σ_c γ_c² S''_c), over the channels with
    radiation type 'P' of a pair of type 0.
  - Rows: `AppendModelRows`, `Objective.cpp` 424-451: one row per level, J/s for J < 0 (else 0, also for a NaN J),
    dr/dγ_c = −2 γ_c S'_c/s, dr/dE = −Σ γ_c² S''_c / s.
- *Arrays*: per level J (dimensionless), and per particle channel γ_c (MeV^½), S'_c = dS_c/dE (MeV⁻¹), S''_c (MeV⁻²).
- *Proposed* (`ParkWall.h`):

      struct ParkLevel {
        double J;                                  // 1 - sum_c gamma_c^2 S'_c
        std::vector<int> channels;                 // the particle channels (1-based, the J group's numbering)
        std::vector<double> gamma, dS, d2S;        // MeV^1/2, MeV^-1, MeV^-2
      };
      double ParkWallChi2(const std::vector<ParkLevel> &, double scale);  // sum over J < 0 of pow(J/scale, 2.0)
      void ParkWallGradient(const ParkLevel &, double scale, std::vector<double> &dGamma, double &dE);  // factor form
      double ParkWallRow(const ParkLevel &, double scale, std::vector<double> &dGamma, double &dE);     // row form

  The engine keeps one gather `ParkLevelsOf(compound, state, bool secondDerivatives)` (the value path has no S'':
  `AZURECalc.cpp` 71) and maps (j, la, c) to `GradAccum::AddGamma`/`AddE`.
- *Curvature*: none outside the rows (J: J ᵀJ of the rows).
- **Verdict: moves.** Three forms side by side in one header, tested against central differences of J(γ, E) with a
  linear S'(E) given in the test; three walks become one. The three skip a NaN J alike (section 6, item 6).

### T7, the abundances

- *Today* (read): value `AbundanceTerm::Penalty`, `Objective.cpp` 191-206 (per sample: each nucleus's prior with
  σ > 0 and not replaced, then the sum-to-one term ((Σw) − c)²/s², s from `AbundanceSamples::SumSigma` and
  `Config::abundanceSumSigma`); the sampler's part, sums only (179-181); gradient `AddAbundancePenaltyGradient`,
  `AZUREGrad.cpp` 110-125 (per ordinal k: the sum part first, then the prior); rows `AppendModelRows`, `Objective.cpp`
  455-470 (one row per sample with s > 0, d/dw_k = 1/s). Σw is `AbundanceSamples::Sum` (`AbundanceSamples.cpp`
  105-109), in nucleus order.
- *Arrays*: w_k (dimensionless, ordinal k = (s − 1) N_nuclei + i), per ordinal the prior (centre, σ, active), per
  sample the sum's centre c and width s.
- *Proposed* (`AbundancePriors.h`):

      struct AbundanceBlock {
        int numSamples = 0, numNuclei = 0;
        std::vector<double> w, centre, sigma;   // per ordinal; sigma used only where active
        std::vector<char> active;               // sigma > 0 and not replaced by the prior covariance
        std::vector<double> sumCentre, sumSigma;  // per sample; sumSigma <= 0: no sum term
      };
      double AbundanceChi2(const AbundanceBlock &, bool withPriors);         // T7; false: the sampler's part
      void AbundanceGradient(const AbundanceBlock &, std::vector<double> &dw);  // d chi2/d w_k, the order of today
      void AbundanceSumRows(const AbundanceBlock &, std::vector<double> &r,
                            std::vector<std::vector<std::pair<int, double>>> &dr);  // (ordinal, 1/s)

  The per-nucleus priors are also entries of `DiagonalPriorsOf` (J, R, S); one gather fills both from the same
  `AbundanceSamples`, so the rule "σ > 0 and not replaced" is written once.
- **Verdict: moves.** Three walks and three forms become one gather and one header.

### T8, the target and apparatus priors (five term objects)

- *Today* (read): `TargetPriorTerm::Add`, `Objective.cpp` 155-163, adds the subtotal of
  `WorkingCopy::ScreeningPenalty` (`HandOver.cpp` 161-166; U_e in eV), `DensityPenalty` (118-123; N in atoms/cm²),
  `DopplerPenalty` (66-71; T in K), `ThicknessPenalty` (34-41; n in atoms/barn), `ApparatusPenalty`
  (`ApparatusParameters.cpp` 143-151; the kernel parameter's engine unit). Each sums from 0, in k order, the value of
  `TargetEffect::ScreeningPenalty` (`TargetEffect.h` 86-89), `DensityPenalty` (120-123), `DopplerPenalty` (321-324),
  `Segment::ThicknessPenalty` (`Segment.h` 476-479), `ApparatusParam::Penalty` (`ApparatusParameters.cpp` 27-30): one
  `GaussianPriorChi2` each, on the value the working copy holds (clamped, R3). The gradient's eligibility is
  `ForEachScreeningPrior` ... `ForEachApparatusPrior` (`Objective.h` 147-204), a second rule (section 6, item 1).
- *Gradient*: A analytic for screening, Doppler, thickness (`AZURECalcGradient.cpp` 207-215), central differences of
  the whole χ² for density and apparatus (277-375); J, R, S from `ForEachDiagonalPrior`.
- *Proposed*: `AppendTargetPriors` (2.1) and `PriorsChi2`; the term adds `PriorsChi2(list)` to its running sum as
  today (R1); the writer's per-row apparatus total is `PriorsChi2(list, row)`, bit for bit `ApparatusPenalty(row)`.
  Nine functions of 4-9 lines go (the five sums and the four per-object penalties; `ApparatusParam::Penalty` too).
- **Verdict: the arithmetic stays shared; two eligibility rules become one gather; the writer's five re-formed
  totals (section 4) read it.**

### T9, the model data

- *Today* (read): `ModelDataTerm`, `Objective.cpp` 211-225, `ModelData::Chi2` and its rows and gradient
  (`src/ModelData.cpp`): per row (q_k(θ) − d_k)/e_k, q_k from `ScatteringLength::ConstraintValue` or from
  `CompoundSystem::ObservedOfFitValues` and `physics/ResonanceQuantities.h`; the gradient by the level-matrix adjoint
  (scattering lengths) or central differences with `CompoundSystem::SetFitLevelParameter`.
- **Verdict: stays.** Its content is the engine walk and physics functions already in `physics/`; the arithmetic left
  is Σ r_k².

### T10, the prior covariance

- *Today* (read): `PriorCovarianceTerm`, `Objective.cpp` 243-255, calls `PriorCovariance::Chi2(p)`, already
  `inference/`'s; gradient `AddGradient`, rows `AppendRows`.
- **Verdict: nothing to move.**

### T11, the `gamma_decay` priors and the wall on a negative remainder

- *Today* (read): `GammaDecayTerm`, `Objective.cpp` 229-240; `GammaDecaySet::PenaltyRows` (`GammaDecayMap.cpp`
  672-712) forms rows (r, sparse dr/dq) from the record's numbers q (the level state's fill) and the record's priors
  (`experiment::GammaDecaySpec`); `Penalty` = Σ r² (714-717); `AddPenaltyGradient` = 2 r dr/dq (719-722).
- The rows are already plain (`PenaltyRow`); the recipe "A derived parametrisation of a level" keeps the record's
  priors in `GammaDecayMap.cpp`, and its input type is an `experiment/` type, which `inference/` may not include.
- **Verdict: stays.**

### Summary

| term | arithmetic | eligibility walks today → after | verdict |
|---|---|---|---|
| T1 points | `DataResidual` (shared) | 1 → 1 | stays |
| T1a, T1b | profiles, ~230 lines | 1 → 1 | moves (`LinearProfile.h`) |
| T2 | whitening (own class) | 1 → 1 | stays |
| T3, T4 | `GaussianPrior.h` | 4 (T3), 2 (T4) → 1 | gather only |
| T5 | `GaussianPrior.h` | 2 → 1 | gather only |
| T6 | three forms | 3 → 1 | moves (`ParkWall.h`) |
| T7 | three forms | 3 → 1 | moves (`AbundancePriors.h`) |
| T8 | `GaussianPrior.h` | 2 → 1 | gather only |
| T9, T10, T11 | engine, library, record | — | stay |

## 3. `ParameterVector`

**What `AZUREParams` holds** (read, `include/AZUREParams.h`, 81 lines; `src/AZUREParams.cpp`, 267): Minuit2's
`MnUserParameters` (per parameter its name, value, error used as Minuit's step and afterwards as the fit's error,
lower and upper limits, fixed flag), one `ParameterRecord` per parameter, a counter per kind and the basis tag of the
last parameter file. Its own functions: `Add`, `Record`, `IndexOfName` (a linear search), the parameter-file reader
and writer, `ReconcileBasis`.

**What depends on Minuit2** (grep over `src/`, `include/`, `api/`): `GetMinuitParams()` is called at 117 places in 32
files: 43 `.Params()` (the values), 28 `.Parameter(i)` (fixed flag, limits, value), 8 `GetName`, 6 `Value`, 2
`SetValue`, 2 `Fix`, 2 `Error`, 1 `SetError`, 1 `Release`, 26 passing the object on. `MnUserParameters` is named in
54 lines of 30 files; 18 files include a Minuit2 header. Outside `src/Fit*`, 31 lines set limits or fix and release
(`ParameterLimitsManager`, `ApplyParameterSettings`, the `AddParameters` of abundances, apparatus and `gamma_decay`,
`SyncShiftGroupParams`, `SyncNormGroupParams`, `ConvertAmplitudeBasis`).

**What a plain class would hold**: names, values, steps, lower and upper limits (±∞ for none), fixed flags, records,
in the full vector's order; the packed map is `ParamIndexMap`'s.

**Who would own it**: `AZUREParams`, as the only holder. Minuit's `MnUserParameters` would be made from it at the
Minuit boundary (`FitMigrad.cpp`, `FitErrors.cpp`, `ProgramFit.cpp`) and the end point written back. A class beside
`AZUREParams` holding the same values is a second holder kept equal by hand (`inference.md` §4) and is not proposed.

**What it buys the terms**: nothing. No term reads `AZUREParams` (read: `Objective.cpp`, `Objective.h` and
`AZUREGrad.h` name no Minuit type); every term reads `const vector_r &p` and the working copy. The fixed flags the
gradient and the band need are read today by building a whole `AZUREParams` through `FillMnParams` (12 such sites,
`AZURECalc::FixedMask`, `HasNuisanceParameters`, `CovarianceBand.cpp` 198); the cost there is the fill, which a plain
holder does not remove.

**Cost**: the 30 files that name `MnUserParameters`, about 300-500 changed lines (inferred from the counts above), the
limits manager's semantics re-written against the plain class (which of `SetLimits`, `SetLowerLimit`,
`SetUpperLimit` was called decides Minuit's internal transformation), and the gate's eight fits as the check that
MIGRAD starts from the same internal state.

**Recommendation: not worth it in this job, and not as a class beside `AZUREParams`.** If a fitter without Minuit2 or
the web UI needs the vector without Minuit2, the replacement (one holder, Minuit's object made at its boundary) is a
job of its own after the slices of section 5.

## 4. The writer, `WorkingCopy::WriteChiSquaredFile` (`DataOutputFiles.cpp` 304-418)

| field of `chiSquared.out` | source today (read) | after |
|---|---|---|
| segment Chi-Squared, Total-Chi-Squared | `objective.segmentChi2`, stored by T1, T1a, T1b, T2 | unchanged |
| per segment Norm | `NormUsed` (stored n* for a profiled norm) | unchanged |
| per segment Norm-Chi-Squared, Total-Norm | **re-formed**, 321-328, own eligibility | `PriorsChi2(list, s, norm)` |
| Total-Park | stored (`objective.park`) | unchanged |
| Total-Shift | **re-formed**, 354-363, `AddShiftPenalty` loop | `PriorsChi2(list of shift parts)` |
| Total-Abundance | stored (`objective.abundance`) | unchanged |
| Total-Screening | **re-formed**, 368-373, `ScreeningPenalty()` | `PriorsChi2(AppendTargetPriors)` |
| Total-Density | **re-formed**, 375-380 | as screening |
| Total-<row>-Chi-Squared (apparatus) | **re-formed**, 382, `WriteApparatusChiSquared` | `PriorsChi2(list, row)` |
| Total-Doppler | **re-formed**, 384-389 | as screening |
| Total-Thickness | **re-formed**, 391-398 | as screening |
| Total-Nuisance | stored (`objective.nuisance`) | unchanged |
| Total-Model-Data | stored | unchanged |
| Total-Gamma-Decay | stored | unchanged |
| Total-Prior-Covariance | stored | unchanged |

Seven totals are re-formed. The per-segment norm prior and Total-Shift are the only ones whose arithmetic is written
in the writer itself; the five T8 totals call the same function as the term.

**The design.** The re-formed totals are formed by the term's own gather and function, called by the writer: the
writer holds no eligibility rule and no formula. They are not stored in `Results::objective`, because a stored total
can be stale: pyazr's `write_outputs` (`AZUREAPI::WriteOutputFiles`, `AZUREAPIReports.cpp` 273-287) evaluates only
the model-only terms before writing (`ModelChi2`: T6, T7, T9, T11), not T3, T4 and T8, while the parameters may have
moved since the last `Objective::Chi2` of the session; the writer's call of today reads the working copy at the
present parameters, and so does the design. The totals that need the compound, the level state or the limits
manager stay stored by their terms, as today.

**Could a printed number change?** No, by construction, for these reasons (read):

- the gather visits the elements in the order the writer's loop does (segments in order; for each, the constant
  shift, the √E coefficient, the function's parameters; T8 in k order) and `PriorsChi2` sums from 0 as the writer
  does; adding a 0 for a segment without a norm prior, as the writer does, leaves the sum's bits unchanged;
- the eligibility of the norm prior is the writer's (σ_n ≠ 0, varied, not profiled, not replaced), which equals the
  term's (section 2, T3: a segment of a THM experiment with a varied norm is profiled; a follower has σ = 0);
- the value read is the working copy's, as today (R3);
- the print conditions (`anyShiftVaried`, `anyPrior` per kind: σ > 0, the thickness's "not a follower") stay in the
  writer: they decide whether a field is printed, not its number.

The stale nuisance total through pyazr's `write_outputs` (`inference.md` §6) is kept: correcting it changes a printed
number and needs the PI's word. The function penalty of `shifts.out` (`DataOutputFiles.cpp` 530-540) is not a field of
`chiSquared.out`; it was repaired on its own, with the PI's word (section 6, item 7).

## 5. Slices

Each slice is one commit with no output changed. Gate of every code slice (as `inference.md` §5): build with 0
warnings; every project in the three parametrizations against the reference outputs (887, 851, 887 files identical,
`work/director/verify_branch.sh`); the eight fits identical to the base build; `ctest -LE slow`; `tests/run_tests.sh`;
user-space instructions at one thread of a T2 calculation and a T2 gradient fit within 1 % of the base build
(`work/director/infer_speed.sh`); `repo_consistency` (checks 8b, 11, 26, 38). Line counts are estimates.

1. **T8 on one list, and the writer's five T8 totals.** `inference/`: `DiagonalPriors.h` (about 60 lines), test
   `inference_diagonal_priors` (about 90): values by hand; `AddPriorsChi2` into a running sum equal bit for bit to the
   loop of `GaussianPriorChi2` it replaces; `PriorsChi2` by group; the gradient and the curvature against central
   differences; an empty list adds nothing. Engine: `AppendTargetPriors` (about 45 lines), `TargetPriorTerm::Add` and
   the writer's five totals on it; the five `WorkingCopy::*Penalty`, the four per-object penalties and
   `ApparatusParam::Penalty` go (about 70 lines). Checked besides the gate: `pyazr_screening`, `pyazr_target_density`,
   `pyazr_doppler`, `pyazr_thickness`, `pyazr_output_penalties`, `pyazr_objective_terms`, `sampler_objective`; scripts
   `electron_screening`, `screening_6li`, `target_density`, `doppler_fgm`, `transmission_yield`, `stopping_yield`,
   `rectangular_resolution`, `target_effect_shared`, `prior_covariance`.
2. **One enumeration of the diagonal priors (T3, T4, T8, T7's per-nucleus priors).** `AppendSegmentPriors`,
   `DiagonalPriorsOf`; `ForEachDiagonalPrior` and its six helpers a loop over the list (114 → about 40 lines);
   `Segment::AddShiftPenalty` and `ShiftPenaltyParts` go (24 lines) and their five callers use the list; the
   writer's norm and shift totals on it (about −20). Before the commit, a scratch program (not committed) compares on
   every project of `tests/` the sets of the value path and of the old enumeration (R4) and the commit message gives
   the count. Checked besides the gate: `pyazr_objective`, `pyazr_objective_object`, `pyazr_output_penalties`;
   scripts `shift_fn`, `shift_fn_groups`, `shift_groups`, `energy_shift_sqrt`, `shift_penalty_zero`, `prior_centre`,
   `norm_groups`, `normalizations_out`, `thm_fit`, `18O_p_a_thm`, `prior_covariance`; `AZURE_GRAD_CHECK` and
   `AZURE_SHIFT_GRAD_CHECK` on T2.
3. **T5 on the list.** `NuisancePriorsOf`; the two copies of the rule, of the value and of the gradient in
   `AZURECalcNuisance.cpp` and `AZUREAPIChi2.cpp` call it (about −30 lines). No new library code. Checked besides the
   gate: `pyazr_objective_terms`, `sampler_objective`; scripts `nuisance_prior`, `nuisance_prior_centre`,
   `nuisance_kinds`, `nuisance_level_index`, `parameter_settings_nan`.
4. **T7 on arrays.** `AbundancePriors.h` (about 70 lines), test `inference_abundance_priors` (about 80): two samples of
   two nuclei by hand, one prior replaced, one sample without a sum term; the gradient against central differences;
   the rows' squares = the sums' part of the value; the sampler's part = value minus the priors. Engine: one gather,
   three call sites (about −30 lines). Checked: scripts `natural_target`, `natural_neon`; `AZURE_GRAD_CHECK` on
   `natural_target`; `pyazr_objective_terms` (T2 has a non-zero abundance prior).
5. **T6 on arrays.** `ParkWall.h` (about 60 lines), test `inference_park_wall` (about 80): one level of two particle
   channels, J by hand, J ≥ 0 gives no value, no gradient and a zero row; value, gradient and row against central
   differences in γ and E with S'(E) linear; the factor form and the row form give rows whose Gauss-Newton gradient
   equals the gradient to 10⁻¹². Engine: `ParkLevelsOf`, three walks become one (about −40 lines in `AZUREGrad.cpp`,
   which is listed in `file_sizes.tsv` at 1108 and must not grow). Checked: scripts `park_formalism`, `thm_park`,
   `mcmc_park_start`; ctest `park_chi2_fallback`, `pyazr_park_gradient`, `pyazr_park_dsde`, `pyazr_thm_park`.
6. **T1a and T1b on arrays.** `LinearProfile.h`, `.cpp` (about 230 lines moved from `ThmExperiment.cpp` and
   `ThmProfile.h`, which goes), `ProfiledNormChi2` (the arithmetic of `ProfileNormChiSquared`, about 20 lines), test
   `inference_linear_profile` (about 150): the line s m + a₀ + a₁ E by hand with 0, 1, 2, 3 background terms; χ² equal
   to the minimum found by a scan; a singular normal matrix; s ≤ 0; the profiled norm's closed form against the
   scan; `ThmProfileDerivative` against central differences. Checked: scripts `thm_experiment`, `thm_fit`,
   `thm_coherent_background`, `thm_options`, `18O_p_a_thm`, `mcmc_classification` (its THM chain).
7. **Documents.** `inference.md` (a section "Done" for this work), `objective.md` §4, `CLAUDE.md` (the layout row of
   `inference/` and the recipe "A term of the objective": the gather in the engine, the arithmetic in `inference/`),
   the Developer Guide's class tables, `output_files.rst` if a field's description names a function that went.

About 330 lines of library code and 470 of tests are added; about 280 engine lines go.

## 6. Risks, and what reading did not settle

1. **Two eligibility rules for one prior** (read). The value path and the enumeration differ in writing:
   `TargetEffect::ScreeningPenalty` also asks `hasScreening_`, `DensityPenalty` `hasDensityFit_`, `DopplerPenalty`
   `dopplerFit_`, which `ForEachScreeningPrior` etc. do not; `ThicknessPenalty` and `ApparatusParam::Penalty` test
   `!(σ > 0)`, the enumerations `σ <= 0`, which differ for σ = NaN. Not settled: whether a screening, density or
   Doppler parameter k can exist on an effect whose flag is false, and whether a NaN σ can reach either (the format-3
   grammar may refuse it). Slice 2's comparison (R4) settles it on the tree's projects; a difference stops the slice.
2. **The value read** (read): T8's value is on the clamped value, its gradient and rows on p (R3). Equal for T3 (n is
   set from p, `HandOver.cpp` 172-187). For T4 not settled: whether `State(segment).energyShift` equals p bit for bit
   for every segment (shift-group followers).
3. **Order of additions** (read): any change of the element order or of the target of an addition moves the last
   bits of χ², which the comparison of 887 files detects only where a printed digit moves; the eight fits detect it in
   MIGRAD's path. The tests of slices 1 and 4 compare bit for bit with the loop they replace.
4. **Allocation per evaluation** (inferred): the gathers allocate small vectors at every evaluation (T8: one entry
   per fitted target parameter; T3 and T4: up to four per segment). Against 1.24 × 10⁹ instructions of one T2
   calculation (`inference.md` §7) this is expected below 10⁻⁴ relative; the gate measures it.
5. **The API's own data χ²** (read): `AZUREAPI::EvaluateFilledChi2` (`AZUREAPIChi2.cpp` 68-125) forms T1, T1a, T1b and
   T2 a second time, with the sum of `Pull()²` where `SegmentTerms` sums `Chi2()` ((d/e)² against d²/e²: different
   rounding). pyazr's `calculate_chi2` and its `objective()` (`Objective::Chi2`) can therefore differ in the last bits
   of the data part. Not changed by any slice; unifying it changes pyazr numbers.
6. **Park and a NaN J** (read): the value and the rows skip a NaN J, the gradient did not (`AZUREGrad.cpp` 90 tested
   `>= 0.0`). **Repaired** 2026-11-16 (commit "fix(grad): the Park wall's gradient skips a NaN J as its value and
   rows do", the PI's word): the gradient tests `!(J < 0)`; ctest `park_chi2_fallback`, variant E.
7. **`shifts.out`'s function penalty** (read): `pen − penOwn` (`DataOutputFiles.cpp` 532-535) subtracts the constant
   shift's prior without asking whether the prior covariance replaced it, and does not subtract the √E coefficient's
   prior, although the comment says "minus those of its constant and sqrt(E) terms". For a segment that owns a shift
   function and has a varied √E shift with σ ≠ 0, or whose constant prior is replaced, the printed penalty is not the
   function's. **Repaired** 2026-11-16 (commit "fix(output): shifts.out prints the shift function's own priors", the
   PI's word): the line prints `PriorsChi2(list, kShiftFn, s)`, the function's own priors by the objective's rule;
   `tests/shift_fn/check.sh` case 8 (a replaced offset prior: 5 before, 9 after, 9 by hand). No `shifts.out` of the
   projects of `tests/` changed (24 with a function line, three bases); the √E case cannot occur (a function and a
   √E term on one data set are refused).
8. **Stale totals through pyazr** (read): `objective.nuisance` is set only by an `Objective::Chi2` with a nuisance
   function; `objective.priorCovariance` by `CalculateChi2RWA` / `Physical` (`AZUREAPIChi2.cpp` 142, 202), not by
   `write_outputs`. Unchanged by the design.
9. **The covariance factor** (read): `CovarianceChi2Correction` freezes a group's factor inside the first evaluation
   (`Group::EnsureFactored`), on a group reached through the set-up data, guarded by a mutex. Whether this is state in
   one of the three owned places of `CLAUDE.md` was not settled; it is why T2 stays.
10. **Coverage** (not settled): which projects of `tests/` exercise each term with a non-zero value in each path was
    not counted; slices 4 and 5 depend on `natural_target` and `park_formalism` reaching J < 0 and a sum term.

## 7. Where the documents and the code disagree (read)

- `include/Objective.h` 28-37 lists the order as T1-T4, T2, T5, T6, "T8 ... (four terms)", T7, T9 "measured scattering
  lengths", T10. The constructor (`Objective.cpp` 259-273) has five T8 objects (apparatus added), T9 covers all model
  data, and T11 is last.
- `inference.md` §6 says a session's `write_outputs` "refreshes the model-only terms and the prior covariance"; the code
  refreshes `ModelChi2` only (T6, T7, T9, T11); the prior covariance total is set by the last `CalculateChi2RWA`.
- The comment of `AZUREAPI::WriteOutputFiles` (`AZUREAPIReports.cpp` 275-279) says pyazr adds the Park, abundance and
  scattering-length terms in Python; `objective.md` §4 step 4 says pyazr computes no prior since pyazr 2.11.0.
- `CLAUDE.md` check 11 says no term is added outside `src/Objective.cpp`; `AZUREAPI::EvaluateFilledChi2` re-forms the
  data terms (item 5 of section 6).
- In `include/Objective.h` 234-250 the comment of `ReplaceDiagonalPriors` stands above that of `AppendModelRows`, and
  the latter names Park's wall and the abundance sums but not the `gamma_decay` rows the function also appends
  (`Objective.cpp` 452-454).
- `api/src/AZUREAPIObjective.cpp` includes `InferenceTypes.h` twice (lines 13-14).
