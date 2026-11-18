// AZURECalc, the evaluation: the chi-squared of a parameter vector as Minuit2 asks for it (operator()), and the
// same value on a given level state and working copy (Chi2On, Chi2OnWith, Chi2Value), which the gradient's finite
// differences, the least-squares drivers and the checks call.  One unit because these are the one forward path:
// hand the vector over (levels, norms, energy shifts), run the point pass, ask the objective (Objective.h).
// The other jobs of the class are files of their own: AZURECalcGradient.cpp (the gradient and G2),
// AZURECalcJacobian.cpp (the residual rows and their Jacobian), AZURECalcNuisance.cpp (the limits manager's
// priors), AZURECalcPool.cpp (the idle states of the caller's FitterStates), AZURECalcOutput.cpp (the
// intermediate files), AZURECalcDiagnostics.cpp (AZURE_DUMP_NONFINITE).
#include "AZURECalc.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "ModelPass.h"
#include "Objective.h"
#include "Profile.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#include "Config.h"
#include "CompoundSystem.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "ParameterLimitsManager.h"

namespace {

// Value returned to Minuit2 for a non-finite or absurdly large chi-squared.
constexpr double kChi2Ceiling = 1.0e100;

}  // namespace

double AZURECalc::operator()(const vector_r &p) const {
  AZUREProfile::ScopedReport profileReport(configure().profile, "evaluation");
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kEvaluation);
  if (callCounter_) callCounter_->fetch_add(1, std::memory_order_relaxed);  // CountCallsInto
  bool isFit = data()->IsFit();
  // A stopped fit (FitControl::stop): every later evaluation returns 0 at once, so the minimizer ends within the
  // sweep it is running; the fit keeps the best point it saw (docs/upgrade/fit.md section 3).
  if (isFit && StopRequested()) return 0.;
  // One atomic read-and-increment: this evaluation's number, distinct across the threads of a parallel gradient.
  const int thisIteration = data()->Iterate();

  WorkingCopy *localData = NULL;
  WorkingData *working = nullptr;
  std::unique_ptr<LevelState> pooled;
  if (isFit) {
    // An evaluation of its own (made when no idle one is left), so evaluations running at once never share one.
    pooled = GetPooledState();
    working = GetPooledData();
    localData = working->data.get();
  } else {
    localData = data();
  }
  const CompoundSystem *localCompound = compound();
  // The level state of this evaluation: the pooled one, or the caller's for an evaluation on the data set itself.
  LevelState &ls = pooled ? *pooled : state_;
  // The buffer of this evaluation: the working copy's, or the caller's for an evaluation on the data set itself.
  azure2::Results &results = working ? working->results : results_;

  // Fill Compound Nucleus From Minuit Parameters
  const auto te0 = std::chrono::steady_clock::now();
  localCompound->FillCompoundFromParams(p, ls);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p, results, localCompound, &ls, &inputs_, &configure());
  // A chi-squared alone: no Park second derivatives (CompoundSystem::CalcShiftFunctions).
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), ls, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
  const auto te1 = std::chrono::steady_clock::now();

  // Calculate all points in one pass over the compound energies: with a team of its own, or on this thread alone
  // when the caller is already inside a parallel region (Minuit2's numerical gradient).
  const bool pre = PrecalculatePoints(ls, localData, results, configure());
  const azure2::support::Environment &environment = azure2::support::ProcessEnvironment();
  if (environment.Set(azure2::support::Switch::kGradTiming)) {
    const int c = ++configure().diagnostics->evaluationTimings;
    if (c <= 20 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      std::cerr << "[eval-timing] evaluation " << c << ": parameters + energy shifts " << ms(te0, te1)
                << " ms, points " << ms(te1, std::chrono::steady_clock::now()) << " ms" << std::endl;
    }
  }

  // The model value of every point, then the objective (Objective.h: data, priors and constraints, in one order).
  azure2::EvaluateModel(*localData, localCompound, ls, configure(), results,
                        pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
  ObjectiveContext context{localCompound, ls, inputs_, localData, results, p, configure()};
  if (limitsManager_) context.nuisance = [this](const vector_r &q) { return CalculateNuisanceChiSquared(q); };
  double chiSquared = Objective::Engine().Chi2(context);

  if (environment.Set(azure2::support::Switch::kDumpNonFinite)) DumpNonFinite(p, chiSquared, localCompound, localData, results);
  if (environment.Set(azure2::support::Switch::kTraceFcn)) {
    const std::string &trace = environment.Text(azure2::support::Switch::kTraceFcn);
    // Diagnostic (docs/upgrade/fit_threads.md): each objective value with its parameter vector, bit for bit
    // (hexadecimal floating point), so two runs can be compared value by value.  With "fresh", a fit's value is
    // also computed on a fresh clone of the master objects: a pooled working copy must give the same bits.
    std::ostringstream os;
    os << std::hexfloat << "[fcn] " << thisIteration << " " << chiSquared << " |";
    for (double v : p) os << " " << v;
    if (isFit && trace == "fresh") {
      std::vector<double> segPool, segFresh;
      Chi2On(ls, localData, results, p, &segPool);
      LevelState fresh(state_);
      std::unique_ptr<WorkingCopy> fd(data()->Clone());
      azure2::Results freshResults;
      os << " || fresh " << Chi2On(fresh, fd.get(), freshResults, p, &segFresh);
      for (size_t s = 1; s < segPool.size() && s < segFresh.size(); s++)
        if (segPool[s] != segFresh[s]) os << " segment " << s << " " << segPool[s] << " " << segFresh[s];
    }
    static std::mutex traceMutex;
    std::lock_guard<std::mutex> lk(traceMutex);
    std::cerr << os.str() << std::endl;
  }

  // The evaluation as data, to the one who runs this objective (the program prints it and writes its intermediate
  // files from it; a fit tracks its best point).
  if (observer_ && observer_->evaluated) observer_->evaluated(thisIteration, chiSquared, p, localData->IsErrorAnalysis());
  if (isFit) {
    // New multithreading with object pools
    ReturnPooledState(std::move(pooled));
    ReturnPooledData(working);
  }

  // A non-finite (or absurdly large) chi-squared must reach Minuit2 as a large
  // FINITE number.  MnLineSearch fits a parabola through the function values;
  // +inf (which this used to return for NaN) turns that parabola into NaN, the
  // line search returns garbage and MIGRAD ends with no status message -- the
  // 11B+a pol_iso --use-gradient fit stopped that way after 30 evaluations
  // (2026-10-11).  A finite ceiling makes the step simply look very bad, so the
  // line search backtracks.  The ceiling stays far above any real chi-squared
  // and far below overflow in the parabola arithmetic.
  if (!std::isfinite(chiSquared) || chiSquared > kChi2Ceiling) chiSquared = kChi2Ceiling;

  return chiSquared;
}

bool AZURECalc::PrecalculatePoints(const LevelState &ls, WorkingCopy *ld, azure2::Results &results, const Config &config,
                                   bool thmOnly) const {
  // thmOnly: the THM segments alone (the THM-only finite differences of Chi2On).
  std::vector<const Segment *> segments;
  for (int i = 1; i <= ld->Data().NumSegments(); i++) {
    const Segment *seg = ld->Data().GetSegment(i);
    if (!seg || (thmOnly && !seg->IsTHM())) continue;
    segments.push_back(seg);
  }
  return ld->CalculatePointsInParallel(segments, compound(), ls, config, results) >= 0;
}

/*!
 * The chi-squared of p on the working copies lc and ld (mirrors the chi-squared
 * math of operator()).  With thmOnly, only THM (HOES) segments are computed and
 * summed: Gradient() finite-differences the THM contribution that way, at a
 * cost that scales with the THM points alone.  The nuisance, Park and abundance
 * penalties are left out there: Gradient() adds them analytically.
 */
double AZURECalc::Chi2On(LevelState &ls, WorkingCopy *ld, azure2::Results &results, const vector_r &p, std::vector<double> *segChis,
                         bool thmOnly) const {
  return Chi2OnWith(configure(), ls, ld, results, p, segChis, thmOnly);
}

double AZURECalc::Chi2OnWith(const Config &config, LevelState &ls, WorkingCopy *ld, azure2::Results &results,
                             const vector_r &p, std::vector<double> *segChis, bool thmOnly) const {
  const CompoundSystem *lc = compound();
  lc->FillCompoundFromParams(p, ls);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, results, lc, &ls, &inputs_, &config);
  // A chi-squared alone: no Park second derivatives (CompoundSystem::CalcShiftFunctions).
  if (config.paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(config, ls, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
  const bool pre = PrecalculatePoints(ls, ld, results, config, thmOnly);

  for (int i = 1; i <= ld->Data().NumSegments(); i++)
    if (!thmOnly || ld->Data().GetSegment(i)->IsTHM())
      azure2::EvaluateModelSegment(*ld, *ld->Data().GetSegment(i), lc, ls, config, results,
                                   pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
  // The objective of operator() (Objective.h), on these working copies: the finite-difference columns of
  // Gradient(), the AZURE_GRAD_CHECK reference and the LM cost differentiate or minimize the same function.
  // thmOnly: the THM segments' data part alone, which adds to an analytic gradient that has every other term.
  ObjectiveContext context{lc, ls, inputs_, ld, results, p, config};
  context.segmentChi2 = segChis;
  context.thmOnly = thmOnly;
  if (limitsManager_) context.nuisance = [this](const vector_r &q) { return CalculateNuisanceChiSquared(q); };
  return Objective::Engine().Chi2(context);
}

double AZURECalc::Chi2Value(const vector_r &p, bool thmOnly) const {
  // Side-effect-free chi-squared (mirrors the chi-squared math of operator()).
  // Pooled working copies keep their applied energy shifts and mapping between
  // calls; a fresh Clone() re-shifts every shifted segment and rebuilds the
  // mapping each time.
  const CompoundSystem *lc = compound();
  std::unique_ptr<LevelState> pooled = GetPooledState();
  LevelState &ls = *pooled;
  WorkingData *working = GetPooledData();
  std::vector<double> segPooled;
  const bool poolCheck = azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kPoolCheck);
  const double chi = Chi2On(ls, working->data.get(), working->results, p, poolCheck ? &segPooled : nullptr, thmOnly);
  ReturnPooledState(std::move(pooled));
  ReturnPooledData(working);
  if (poolCheck) {
    if (++configure().diagnostics->poolChecks <= 5) {
      LevelState fresh(state_);
      WorkingCopy *fd = data()->Clone();
      std::vector<double> segFresh;
      // The reference: the per-point (serial) path on a fresh copy, as AZURE_SERIAL_CHI2 runs it.
      Config serial(configure());
      serial.environment.serialPointPass = true;
      azure2::Results freshResults;
      const double chiFresh = Chi2OnWith(serial, fresh, fd, freshResults, p, &segFresh, thmOnly);
      delete fd;
      std::cerr << "[pool-check] pooled+parallel " << chi << " fresh+serial " << chiFresh << std::endl;
      std::vector<std::pair<double, int>> d;
      for (size_t i = 1; i < segFresh.size() && i < segPooled.size(); i++)
        d.push_back({std::fabs(segPooled[i] - segFresh[i]), (int)i});
      std::sort(d.rbegin(), d.rend());
      for (int k = 0; k < 6 && k < (int)d.size() && d[k].first > 1e-6; k++)
        std::cerr << "[pool-check]   segment " << d[k].second << ": pooled " << segPooled[d[k].second] << " fresh "
                  << segFresh[d[k].second] << std::endl;
    }
  }
  return chi;
}
