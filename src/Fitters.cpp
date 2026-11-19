// The table of the engine's fitters (Fitters.h; docs/upgrade/fit.md section 3), the driver that runs a row and says
// which settings it ignored, and the program's switches and usage lines made from the rows.
#include "Fitters.h"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "Config.h"

const std::vector<FitterRow> &Fitters() {
  // In the program's precedence: with several switches the first row given runs.  Columns: name, option, gradient,
  // residualRows, bounds, covariance, minos, callLimit, fallback, progressLabel, doc, run.
  static const std::vector<FitterRow> rows = {
#ifdef USE_NLOPT
      {"nlopt", "--use-nlopt", false, false, true, false, false, true, "", "",
       "NLopt's derivative-free algorithm (nlopt_algorithm=, default SBPLX) on the objective; no covariance.",
       &FitNLopt},
#endif
      {"lm", "--use-lm", false, true, true, true, false, false, "migrad", "LM iteration",
       "Levenberg-Marquardt on the residual rows and their analytic Jacobian; covariance (J^T J + priors)^-1.",
       &FitLevenbergMarquardt},
      {"migrad-gradient", "--use-gradient", true, false, true, true, true, true, "", "",
       "Minuit2's MIGRAD with the analytic gradient and Gauss-Newton G2; HESSE covariance, MINOS.",
       &FitMigradGradient},
      {"migrad", "", false, false, true, true, true, true, "", "",
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

const std::vector<RemovedFitterRow> &RemovedFitters() {
  static const std::vector<RemovedFitterRow> rows = {
      // The owner's answer of 2026-11-18 (docs/upgrade/OWNER_ANSWERS.md, section 10 item 8: "Remove both").
      {"gsl-lm", "--use-gsl-lm", "2026-11-18",
       "GSL's trust-region least squares minimised the same residual rows as lm and reached the same minima.  Use "
       "--use-lm (method \"lm\": Levenberg-Marquardt on the same rows and analytic Jacobian, the same covariance "
       "(J^T J + priors)^-1) or --use-gradient (method \"migrad-gradient\": MIGRAD with the analytic gradient, HESSE "
       "and MINOS)"},
  };
  return rows;
}

const RemovedFitterRow *RemovedFitterNamed(const std::string &name) {
  for (const RemovedFitterRow &row : RemovedFitters())
    if (name == row.name || name == row.option) return &row;
  return nullptr;
}

std::string RemovedFitterMessage(const RemovedFitterRow &row, const std::string &asked) {
  return asked + " asks for the fitter " + row.name + ", which was removed from AZURE2 on " + row.removed + ".  " +
         row.instead + ".";
}

const FitterRow *FitterNamed(const std::string &name) {
  for (const FitterRow &row : Fitters())
    if (name == row.name) return &row;
  return nullptr;
}

namespace {

/// The place of the row named `name` in Fitters(); throws for a name that is no row.
std::size_t FitterPlace(const std::string &name) {
  const std::vector<FitterRow> &rows = Fitters();
  for (std::size_t k = 0; k < rows.size(); k++)
    if (name == rows[k].name) return k;
  throw std::logic_error("'" + name + "' is no fitter of this build (Fitters(), src/Fitters.cpp)");
}

}  // namespace

const FitterRow &FitterOfConfig(const Config &configure) {
  if (!configure.fitter.empty()) return Fitters()[FitterPlace(configure.fitter)];
  for (const FitterRow &row : Fitters())
    if (!*row.option) return row;
  throw std::logic_error("Fitters() has no row without a switch (src/Fitters.cpp)");
}

void ChooseFitter(Config &configure, const std::string &name) {
  const std::size_t place = FitterPlace(name);
  if (configure.fitter.empty() || place < FitterPlace(configure.fitter)) configure.fitter = name;
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
      ChooseFitter(configure, row.name);
      return true;
    }
  for (const RemovedFitterRow &row : RemovedFitters())
    if (option == row.option) throw std::invalid_argument(RemovedFitterMessage(row, option));
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
