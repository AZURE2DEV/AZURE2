// The NLopt fitter (builds with USE_NLOPT): a derivative-free algorithm of NLopt on the free parameters
// (Fitters.h; docs/upgrade/fit.md).
#ifdef USE_NLOPT
#include <algorithm>
#include <string>
#include <vector>

#include <nlopt.h>

#include "AZURECalc.h"
#include "WorkingCopy.h"
#include "CNuc.h"
#include "AZUREParams.h"
#include "Config.h"
#include "FitEvaluations.h"
#include "Fitters.h"

namespace {

// What the objective wrapper needs: the objective, the full vector (fixed parameters included) and where the free
// ones go, and the stop request.
struct NLoptData {
  const AZURECalc *theFunc;
  std::vector<double> *full_params;
  std::vector<unsigned> *variable_indices;
  const std::atomic<bool> *stop;
  nlopt_opt opt;
};

// NLopt hands over the free parameters only; the full vector is rebuilt with the fixed ones.
double nlopt_objective(unsigned n, const double *x, double *grad, void *f_data) {
  NLoptData *data = static_cast<NLoptData *>(f_data);
  if (data->stop && data->stop->load(std::memory_order_relaxed)) nlopt_force_stop(data->opt);
  for (unsigned j = 0; j < n; j++) {
    unsigned i = (*data->variable_indices)[j];
    (*data->full_params)[i] = x[j];
  }
  const double chi_squared = (*data->theFunc)(*data->full_params);
  (void)grad;  // derivative-free algorithms only: NLopt never asks for it
  return chi_squared;
}

}  // namespace

FitResult FitNLopt(AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  std::ostream &log = *control.log;
  FitResult result;
  result.fitter = "nlopt";
  log << "Using NLopt minimizer..." << std::endl;
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const unsigned n_total = mp.Params().size();
  std::vector<double> full_params(n_total);
  std::vector<unsigned> variable_indices;
  for (unsigned i = 0; i < n_total; i++) {
    full_params[i] = mp.Value(i);
    if (!mp.Parameter(i).IsFixed()) variable_indices.push_back(i);
  }
  const unsigned n_variable = variable_indices.size();
  std::vector<double> x(n_variable), lb(n_variable), ub(n_variable);
  for (unsigned j = 0; j < n_variable; j++) {
    const unsigned i = variable_indices[j];
    x[j] = mp.Value(i);
    lb[j] = mp.Parameter(i).HasLowerLimit() ? mp.Parameter(i).LowerLimit() : -1e10;  // large for unbounded
    ub[j] = mp.Parameter(i).HasUpperLimit() ? mp.Parameter(i).UpperLimit() : 1e10;
  }
  log << "Optimizing " << n_variable << " parameters (out of " << n_total << " total)..." << std::endl;

  const nlopt_algorithm algorithms[] = {NLOPT_LN_SBPLX, NLOPT_LN_COBYLA, NLOPT_LN_BOBYQA,
                                        NLOPT_LN_NEWUOA, NLOPT_LN_PRAXIS, NLOPT_LN_NELDERMEAD};
  const char *algorithm_names[] = {"SBPLX", "COBYLA", "BOBYQA", "NEWUOA", "PRAXIS", "Nelder-Mead"};
  const int alg_index = objective.configure().nloptAlgorithm;
  log << "Using NLopt " << algorithm_names[alg_index] << " algorithm..." << std::endl;
  nlopt_opt opt = nlopt_create(algorithms[alg_index], n_variable);

  NLoptData nlopt_data{&objective, &full_params, &variable_indices, control.stop, opt};
  nlopt_set_min_objective(opt, nlopt_objective, &nlopt_data);
  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  const unsigned int noCallLimit = 1000000000u;
  nlopt_set_maxeval(opt, control.maxCalls > 0 ? static_cast<int>(std::min(control.maxCalls, noCallLimit)) : 0);
  nlopt_set_ftol_rel(opt, 1e-8);
  nlopt_set_xtol_rel(opt, 1e-8);

  const std::string callLimitText =
      control.maxCalls > 0 ? std::to_string(std::min(control.maxCalls, noCallLimit)) + control.callLimitSource
                           : std::string("none");
  double minf;
  nlopt_result status;
  bool stopped;
  {
    FitEvaluations evaluations(objective, control, "nlopt");
    status = nlopt_optimize(opt, x.data(), &minf);
    stopped = evaluations.Stopped();
    if (stopped && evaluations.HasBest()) {
      const std::vector<double> best = evaluations.BestPoint();
      for (unsigned j = 0; j < n_variable; j++) x[j] = best[variable_indices[j]];
      minf = evaluations.BestObjective();
    }
  }
  const long nloptCalls = nlopt_get_numevals(opt);
  result.calls = nloptCalls;
  result.statusMinimizer = "NLopt";
  if (stopped) {
    result.status = "stopped";
    log << std::endl << "NLopt stopped on request after " << nloptCalls << " function calls." << std::endl;
  } else if (status == NLOPT_MAXEVAL_REACHED) {
    log << std::endl
        << "WARNING: NLopt did NOT converge: stopped at the call limit (" << nloptCalls << " function calls, limit "
        << callLimitText << "); final chi-squared " << minf << ".  Rerun the fit from param.sav or raise --max-calls."
        << std::endl;
    result.status = "call_limit";
  } else if (status < 0) {
    log << std::endl;
    log << "NLopt optimization failed: " << nlopt_result_to_string(status) << std::endl;
    result.status = std::string("failed_") + nlopt_result_to_string(status);
  } else {
    result.status = std::string("converged_") + nlopt_result_to_string(status);
    result.converged = true;
    result.statusValid = true;
    log << std::endl;
    log << "NLopt optimization succeeded: " << nlopt_result_to_string(status) << std::endl;
    log << "Final chi-squared: " << minf << std::endl;
  }
  for (unsigned j = 0; j < n_variable; j++) mp.SetValue(variable_indices[j], x[j]);
  nlopt_destroy(opt);
  result.x = mp.Params();
  result.objective = minf;
  return result;
}
#endif
