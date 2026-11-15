// The MIGRAD fitters (Minuit2), with Minuit's numerical gradient or the objective's analytic one, and the rerun with a
// re-frozen relvector data covariance (Fitters.h; docs/upgrade/fit.md).
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserCovariance.h>

#include "AZURECalc.h"
#include "CompoundSystem.h"
#include "AZUREParams.h"
#include "FitEvaluations.h"
#include "MigradSeed.h"
#include "Fitters.h"
#include "WorkingCopy.h"

namespace {

// Minuit2 scales the call limit by 1.3 in an int, so "none" is 10^9 calls.
constexpr unsigned int kNoCallLimit = 1000000000u;

std::string CallLimitText(const FitControl &control) {
  return control.maxCalls > 0 ? std::to_string(std::min(control.maxCalls, kNoCallLimit)) + control.callLimitSource
                              : std::string("none");
}

// The covariance of the free parameters from Minuit's packed matrix, row-major over result.freeIndices.
void TakeCovariance(const ROOT::Minuit2::FunctionMinimum &m, const ROOT::Minuit2::MnUserParameters &start,
                    FitResult &result) {
  result.freeIndices.clear();
  for (unsigned int i = 0; i < start.Params().size(); i++)
    if (!start.Parameter(i).IsFixed()) result.freeIndices.push_back((int)i);
  const ROOT::Minuit2::MnUserCovariance &cov = m.UserCovariance();
  const unsigned int n = (unsigned int)result.freeIndices.size();
  result.covariance.clear();
  if (cov.Nrow() != n || n == 0) return;
  result.covariance.assign((size_t)n * n, 0.0);
  for (unsigned int a = 0; a < n; a++)
    for (unsigned int b = 0; b < n; b++) result.covariance[(size_t)a * n + b] = cov(a, b);
}

// Say how MIGRAD ended; a fit that stopped at the call limit or with an invalid minimum is flagged.
void ReportMigrad(const ROOT::Minuit2::FunctionMinimum &m, long calls, const FitControl &control, FitResult &result) {
  std::ostream &log = *control.log;
  const std::string callLimitText = CallLimitText(control);
  const double goal = 0.002 * 0.1 * m.Up();  // MnMigrad's default tolerance 0.1 (VariableMetricBuilder)
  std::string status;
  std::ostringstream why;
  why << std::scientific << std::setprecision(3);
  if (m.IsValid()) {
    status = "converged";
    log << std::endl << "MIGRAD converged: valid minimum, EDM " << std::scientific << std::setprecision(3) << m.Edm()
        << " (goal " << goal << ")" << std::defaultfloat << std::setprecision(6) << ", " << calls
        << " function calls, call limit " << callLimitText << "." << std::endl;
  } else {
    if (m.HasReachedCallLimit()) {
      status = "call_limit";
      why << "stopped at the call limit (" << calls << " function calls, limit " << callLimitText << ")";
    } else {
      status = "invalid";
      why << "invalid minimum";
      if (m.IsAboveMaxEdm()) why << ", EDM above its goal";
      if (!m.HasValidCovariance()) why << ", no valid covariance";
      if (m.HasMadePosDefCovar()) why << ", covariance forced positive definite";
      if (m.HesseFailed()) why << ", HESSE failed";
      why << "; " << calls << " function calls";
    }
    why << "; EDM " << m.Edm() << " (goal " << goal << ")";
    // At the call limit the fit ends at the lowest objective MIGRAD evaluated (FitMigradRows); otherwise at
    // Minuit2's end state.
    log << std::endl << "WARNING: MIGRAD did NOT converge: " << why.str() << ".  The parameters and outputs written are "
        << (status == "call_limit" ? "the lowest point it evaluated" : "its end point")
        << ", not a minimum; rerun the fit from param.sav" << (status == "call_limit" ? " or raise --max-calls" : "")
        << "." << std::endl;
  }
  result.status = status;
  result.converged = m.IsValid();
  result.edm = m.Edm();
  result.edmGoal = goal;
  result.calls = calls;
}

FitResult FitMigradRows(AZURECalc &objective, AZUREParams &params, const FitControl &control, bool gradient) {
  FitResult result = MigradPass(objective, params.GetMinuitParams(), gradient, control);
  // relvector data-covariance terms use the model frozen at the FIRST evaluation of the run,
  // i.e. at the starting parameters, so MIGRAD minimized a chi^2 whose covariance belongs to
  // the start.  Re-freeze it at the fitted parameters: the outputs (chiSquared.out, ...) then
  // match a calculation from the param.sav written below; with --refreeze-covariance[=N]
  // MIGRAD is rerun with the re-frozen covariance until the objective stops changing.
  // (12C+a_onefile 10-29-26: fits that started 15k-90k above their minimum reported a chi^2
  // 50-600 away from a mode-1 calculation of their own param.sav.)
  if (result.status != "stopped" && objective.data()->Data().HasRelativeCovariance()) {
    std::ostream &log = *control.log;
    for (int pass = 0;; pass++) {
      const double frozen = result.objective;
      objective.data()->RefreezeCovariance();
      const double refrozen = objective(result.x);
      log << "Data covariance: relvector model re-frozen at the fitted parameters; objective " << frozen
          << " (model frozen at the start of this pass) -> " << refrozen << "." << std::endl;
      if (std::fabs(refrozen - frozen) < 0.5) break;
      if (pass >= control.refreezePasses) {
        log << "WARNING: the fitted parameters minimize the chi^2 with the start's covariance "
               "model, not with their own"
            << (control.refreezePasses > 0 ? " (pass limit reached)" : "")
            << "; rerun from param.sav or use --refreeze-covariance." << std::endl;
        break;
      }
      log << "Refitting with the re-frozen covariance (pass " << pass + 1 << " of at most " << control.refreezePasses
          << ")..." << std::endl;
      ROOT::Minuit2::MnUserParameters restart = result.minimum->UserParameters();
      for (unsigned int i = 0; i < result.x.size(); i++) restart.SetValue(i, result.x[i]);
      result = MigradPass(objective, restart, gradient, control);
      if (result.status == "stopped") break;
    }
  }
  AdoptMigradEndPoint(params, result);
  return result;
}

}  // namespace

namespace {

/// One MIGRAD from `start`.  `seed`: the analytic gradient with its bounded Gauss-Newton seed, watched (MigradSeed.h);
/// nullptr with `analyticGradient`: Minuit2's numerical seed.  `*stalled`: the seed's first step found no improvement
/// and the pass ended on the seed's own stop request; the result is then not a fit.
FitResult MigradOnce(AZURECalc &objective, const ROOT::Minuit2::MnUserParameters &start, bool analyticGradient,
                     const FitControl &control, MigradSeed *seed, bool *stalled) {
  FitResult result;
  result.fitter = analyticGradient ? "migrad-gradient" : "migrad";
  const unsigned int migradCalls = control.maxCalls > 0 ? std::min(control.maxCalls, kNoCallLimit) : kNoCallLimit;
  // MIGRAD's function calls as reported (fitStatus.out, the log): counted by the objective itself, from 0 for each
  // MIGRAD.  Minuit2's NFcn() is not used: its counter is a plain int that its OpenMP numerical gradient increments
  // from several threads (lost increments, a count a few low, varying from run to run); Minuit2 still stops at
  // --max-calls by that counter, so under threads the limit is approximate by as much.
  std::atomic<long> calls{0};
  objective.CountCallsInto(&calls);
  std::unique_ptr<ROOT::Minuit2::FunctionMinimum> minimum;
  bool haveBest = false;  // the lowest objective MIGRAD evaluated, and its point
  double bestObjective = 0.0;
  std::vector<double> bestPoint;
  {
    FitEvaluations evaluations(objective, control, result.fitter.c_str());
    // With the analytic gradient MIGRAD is given the objective as an FCNGradientBase (AZURECalc::Gradient);
    // otherwise as a plain FCNBase, so that Minuit uses its own numerical gradient.
    if (analyticGradient) {
      // The seed's watch: an observer in front of `evaluations`' own, with the seed's stop request (which follows
      // the caller's).  Put back before `evaluations` ends.
      const EvaluationObserver *inner = objective.Observer();
      EvaluationObserver watching;
      const long freeParameters = start.VariableParameters();
      watching.evaluated = [&](long evaluation, double value, const std::vector<double> &p, bool errorAnalysis) {
        seed->Forward(control.stop);
#ifdef _OPENMP
        const bool inParallel = omp_in_parallel() != 0;  // inside Minuit2's numerical gradient
#else
        const bool inParallel = false;
#endif
        if (seed->Evaluated(value, inParallel))
          *control.log << std::endl
                       << "MIGRAD is computing a numerical Hessian (HESSE): about "
                       << freeParameters * (freeParameters - 1) / 2 + 4 * freeParameters << " function evaluations for "
                       << freeParameters << " free parameters, with no gradient calls; chi^2 does not change while it runs."
                       << std::endl;
        if (inner && inner->evaluated) inner->evaluated(evaluation, value, p, errorAnalysis);
      };
      if (seed) {
        watching.stop = seed->Stop();
        objective.Observe(&watching);
        objective.SeedWith(seed);
      } else {
        objective.DisableG2();
      }
      if (seed && objective.EnableG2(start.Params())) {
        *control.log << "Seeding MIGRAD with the analytic Gauss-Newton G2." << std::endl;
        const MigradSeed::Bounded bounded = seed->Seed();
        if (bounded.raised > 0) {
          *control.log << "MIGRAD seed: the first step of " << bounded.raised << " of " << bounded.free
                       << " free parameters was bounded (weakly constrained or near a limit; largest: parameter "
                       << bounded.worst + 1 << ", ";
          if (std::isfinite(bounded.worstSteps))
            *control.log << bounded.worstSteps << " initial step sizes";
          else
            *control.log << "no curvature";
          *control.log << ")." << std::endl;
        }
      } else if (seed) {
        *control.log << "Analytic G2 unavailable for this model; MIGRAD seeds numerically." << std::endl;
      }
      ROOT::Minuit2::MnMigrad migrad(objective, start);
      minimum.reset(new ROOT::Minuit2::FunctionMinimum(migrad(migradCalls)));
      if (seed) {
        objective.SeedWith(nullptr);
        objective.Observe(inner);
        if (stalled) *stalled = seed->Stalled() && !evaluations.Stopped();
      }
    } else {
      const ROOT::Minuit2::FCNBase &fb = objective;
      ROOT::Minuit2::MnMigrad migrad(fb, start);
      minimum.reset(new ROOT::Minuit2::FunctionMinimum(migrad(migradCalls)));
    }
    if (stalled && *stalled) {
      objective.CountCallsInto(nullptr);
      result.calls = calls.load();
      return result;
    }
    if (evaluations.Stopped()) {
      // The evaluations after the request returned 0: Minuit's minimum is not a point of the objective.  The fit
      // ends at the best point it saw.
      result.status = "stopped";
      result.calls = calls.load();
      result.x = evaluations.HasBest() ? evaluations.BestPoint() : start.Params();
      result.objective = evaluations.HasBest() ? evaluations.BestObjective() : objective.Chi2Value(result.x);
      *control.log << std::endl
                   << "MIGRAD stopped on request after " << result.calls << " function calls; the best point seen "
                   << "is kept." << std::endl;
    }
    haveBest = evaluations.HasBest();
    if (haveBest) {
      bestObjective = evaluations.BestObjective();
      bestPoint = evaluations.BestPoint();
    }
  }
  objective.CountCallsInto(nullptr);
  if (result.status != "stopped") {
    ReportMigrad(*minimum, calls.load(), control, result);
    result.x = minimum->UserParameters().Params();
    result.objective = minimum->Fval();
    TakeCovariance(*minimum, start, result);
    // At the call limit Minuit2 returns its last accepted state, which can be the start itself (the limit reached
    // while it seeds the gradient).  The fit's end point is then the lowest objective it evaluated.
    if (result.status == "call_limit" && haveBest && bestObjective < result.objective &&
        bestPoint.size() == result.x.size()) {
      std::ostringstream notice;
      notice << "MIGRAD stopped at its call limit; the end point is the lowest objective it evaluated, "
             << std::setprecision(10) << bestObjective << ", not Minuit2's last state (" << result.objective
             << "), whose covariance is the one reported.";
      *control.log << notice.str() << std::endl;
      result.notices.push_back(notice.str());
      result.x = bestPoint;
      result.objective = bestObjective;
    }
  }
  result.minimum = std::shared_ptr<ROOT::Minuit2::FunctionMinimum>(std::move(minimum));
  return result;
}

}  // namespace

FitResult MigradPass(AZURECalc &objective, const ROOT::Minuit2::MnUserParameters &start, bool analyticGradient,
                     const FitControl &control) {
  if (!analyticGradient) return MigradOnce(objective, start, false, control, nullptr, nullptr);
  *control.log << "Using analytic gradient for Minuit2." << std::endl;
  MigradSeed seed(start, objective.Up());
  bool stalled = false;
  FitResult result = MigradOnce(objective, start, true, control, &seed, &stalled);
  if (!stalled) return result;
  // Minuit2 answers a failed first line search with a numerical HESSE of n(n-1)/2 evaluations.  Its own numerical
  // seed (second derivatives by differences, the limit transformations included) costs a few evaluations per
  // parameter: tried once before that.
  *control.log << "MIGRAD's first step from the analytic seed found no improvement (" << result.calls
               << " function calls); restarting with Minuit2's numerical seed." << std::endl;
  return MigradOnce(objective, start, true, control, nullptr, nullptr);
}

void AdoptMigradEndPoint(AZUREParams &params, const FitResult &result) {
  if (result.status != "stopped" && result.minimum) params.GetMinuitParams() = result.minimum->UserParameters();
  for (unsigned int i = 0; i < result.x.size(); i++) params.GetMinuitParams().SetValue(i, result.x[i]);
}

FitResult FitMigrad(AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  return FitMigradRows(objective, params, control, false);
}

FitResult FitMigradGradient(AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  return FitMigradRows(objective, params, control, true);
}
