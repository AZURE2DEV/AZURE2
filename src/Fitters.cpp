// The table of the engine's fitters (Fitters.h; docs/upgrade/fit.md section 3), the driver that runs a row and says
// which settings it ignored, and the program's switches and usage lines made from the rows.
#include "Fitters.h"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "Config.h"

const std::vector<FitterRow> &Fitters() {
  // In the program's precedence: with several switches the first row given runs.  Columns: name, option, configBit,
  // gradient, residualRows, bounds, covariance, minos, callLimit, fallback, progressLabel, doc, run.
  static const std::vector<FitterRow> rows = {
#ifdef USE_NLOPT
      {"nlopt", "--use-nlopt", Config::USE_NLOPT_MINIMIZER, false, false, true, false, false, true, "", "",
       "NLopt's derivative-free algorithm (nlopt_algorithm=, default SBPLX) on the objective; no covariance.",
       &FitNLopt},
#endif
      {"gsl-lm", "--use-gsl-lm", Config::USE_GSL_LM_MINIMIZER, false, true, true, true, false, false, "migrad",
       "GSL-LM iteration",
       "GSL's geodesic-accelerated trust-region least squares on the residual rows; covariance (J^T J + priors)^-1.",
       &FitGslTrustRegion},
      {"lm", "--use-lm", Config::USE_LM_MINIMIZER, false, true, true, true, false, false, "migrad", "LM iteration",
       "Levenberg-Marquardt on the residual rows and their analytic Jacobian; covariance (J^T J + priors)^-1.",
       &FitLevenbergMarquardt},
      {"migrad-gradient", "--use-gradient", Config::USE_ANALYTIC_GRADIENT, true, false, true, true, true, true, "", "",
       "Minuit2's MIGRAD with the analytic gradient and Gauss-Newton G2; HESSE covariance, MINOS.",
       &FitMigradGradient},
      {"migrad", "", 0u, false, false, true, true, true, true, "", "",
       "Minuit2's MIGRAD with its numerical gradient (the default); HESSE covariance, MINOS.", &FitMigrad},
  };
  return rows;
}

const std::vector<AbsentFitterRow> &AbsentFitters() {
  static const std::vector<AbsentFitterRow> rows = {
#ifndef USE_NLOPT
      {"nlopt", "--use-nlopt", "NLopt installed and the CMake option -DUSE_NLOPT=ON"},
#endif
  };
  return rows;
}

const FitterRow *FitterNamed(const std::string &name) {
  for (const FitterRow &row : Fitters())
    if (name == row.name) return &row;
  return nullptr;
}

const FitterRow &FitterOfConfig(const Config &configure) {
  for (const FitterRow &row : Fitters())
    if (row.configBit != 0u && (configure.paramMask & row.configBit)) return row;
  return *FitterNamed("migrad");
}

const FitterRow &FitterThatRan(const FitterRow &row, const FitResult &result) {
  if (result.fellBackToMigrad && *row.fallback)
    if (const FitterRow *ran = FitterNamed(row.fallback)) return *ran;
  return row;
}

FitResult RunFitter(const FitterRow &row, AZURECalc &objective, AZUREParams &params, const FitControl &control) {
  FitResult result = row.run(objective, params, control);
  const FitterRow &ran = FitterThatRan(row, result);
  if (control.maxCalls > 0 && !ran.callLimit) {
    const std::string notice = std::string(ran.name) + " does not honour a call limit: the limit " +
                               std::to_string(control.maxCalls) + control.callLimitSource +
                               " was ignored; it ended by its own criterion after " + std::to_string(result.calls) +
                               (ran.residualRows ? " residual passes." : " calls.");
    *control.log << "Note: " << notice << std::endl;
    result.notices.push_back(notice);
  }
  return result;
}

std::string FitOutcome(const FitResult &result) {
  if (result.converged) return "converged";
  if (result.status == "call_limit" || result.status == "iteration_limit" || result.status == "stopped" ||
      result.status == "invalid")
    return result.status;
  return "failed";
}

std::string FitterSwitchesHelp() {
  std::ostringstream out;
  std::string honour;
  for (const FitterRow &row : Fitters()) {
    if (row.callLimit) honour += (honour.empty() ? "" : ", ") + std::string(row.name);
    if (!*row.option) continue;
    out << std::setw(25) << std::left << ("\t" + std::string(row.option) + ":") << std::setw(0) << "Fit with "
        << row.name << ": " << row.doc;
    if (*row.fallback) out << " Runs " << row.fallback << " for a model without an analytic Jacobian.";
    out << std::endl;
  }
  out << std::setw(25) << std::left << "\t--max-calls=N:" << std::setw(0) << "Stop a fit after N objective calls ("
      << honour << "; the other fitters ignore it and say so; default: no limit, the fit ends at convergence)."
      << std::endl;
  return out.str();
}

bool ApplyFitterSwitch(Config &configure, const std::string &option) {
  for (const FitterRow &row : Fitters())
    if (*row.option && option == row.option) {
      configure.paramMask |= row.configBit;
      return true;
    }
  for (const AbsentFitterRow &row : AbsentFitters())
    if (option == row.option) {
      std::string allowed;
      for (const FitterRow &r : Fitters())
        if (*r.option) allowed += (allowed.empty() ? "" : ", ") + std::string(r.option);
      throw std::invalid_argument(option + " selects the fitter " + row.name +
                                  ", which this build of AZURE2 does not have.  The fitter switches of this build: " +
                                  allowed + ", or none for migrad.  To use " + row.name + ", build AZURE2 with " +
                                  row.build + ".");
    }
  return false;
}
