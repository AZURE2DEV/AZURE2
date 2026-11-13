// The table of the engine's fitters (Fitters.h; docs/upgrade/fit.md section 3).
#include "Fitters.h"

#include "Config.h"

const std::vector<FitterRow> &Fitters() {
  // In the program's precedence: with several switches the first row given runs.
  static const std::vector<FitterRow> rows = {
#ifdef USE_NLOPT
      {"nlopt", "--use-nlopt", Config::USE_NLOPT_MINIMIZER, false, false, true, false, false, "",
       "NLopt's derivative-free algorithm (nlopt_algorithm=, default SBPLX) on the objective; no covariance.",
       &FitNLopt},
#endif
      {"gsl-lm", "--use-gsl-lm", Config::USE_GSL_LM_MINIMIZER, false, true, true, true, false, "GSL-LM iteration",
       "GSL's geodesic-accelerated trust-region least squares on the residual rows; covariance (J^T J + priors)^-1.",
       &FitGslTrustRegion},
      {"lm", "--use-lm", Config::USE_LM_MINIMIZER, false, true, true, true, false, "LM iteration",
       "Levenberg-Marquardt on the residual rows and their analytic Jacobian; covariance (J^T J + priors)^-1.",
       &FitLevenbergMarquardt},
      {"migrad-gradient", "--use-gradient", Config::USE_ANALYTIC_GRADIENT, true, false, true, true, true, "",
       "Minuit2's MIGRAD with the analytic gradient and Gauss-Newton G2; HESSE covariance, MINOS.",
       &FitMigradGradient},
      {"migrad", "", 0u, false, false, true, true, true, "",
       "Minuit2's MIGRAD with its numerical gradient (the default); HESSE covariance, MINOS.", &FitMigrad},
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
