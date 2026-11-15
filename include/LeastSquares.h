#ifndef LEASTSQUARES_H
#define LEASTSQUARES_H

// What the two least-squares fitters (LM: src/FitLeastSquares.cpp, GSL-LM: src/FitGslTrustRegion.cpp) share: the free
// parameters with their limits and priors, the covariance (J^T J + priors)^-1 at the end, and the row driver with its
// MIGRAD fallback.

#include <ostream>
#include <string>
#include <vector>

#include "Constants.h"
#include "FitControl.h"

class AZURECalc;
class AZUREParams;
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

/// What one least-squares run hands back besides the end point it leaves in its AZUREParams.
struct LeastSquaresOutcome {
  double objective = -1.0;          ///< the objective at the end; negative: the analytic Jacobian is unsupported
  std::string status = "converged"; ///< converged, iteration_limit, stopped, failed_<reason>
  long calls = 0;                   ///< residual passes (forward or with the Jacobian)
  std::vector<int> freeIndices;     ///< full indices of the covariance's rows
  std::vector<double> covariance;   ///< (J^T J + priors)^-1, row-major; empty when singular
  BandCovariance band;              ///< its band block (BandPackedColumns)
};

/*!
 * One residual Jacobian at `full` fixes the packed column order (`p2f`, in `r`/`J`); fills the free start values
 * `x`, the projection limits `lo`/`hi` and the Gaussian priors `pen_nom`/`pen_inv2` (0 where none: the diagonal
 * priors, then the nuisance priors of levels and of the THM coherent background).  Returns the number of free
 * parameters, or -1 when the analytic Jacobian is unsupported.
 */
int PrepareFreeParams(const AZURECalc &self, const vector_r &full, const ROOT::Minuit2::MnUserParameters &mp,
                      vector_r &r, vector_r &J, std::vector<int> &p2f, vector_r &x, std::vector<double> &lo,
                      std::vector<double> &hi, std::vector<double> &pen_nom, std::vector<double> &pen_inv2);

/*!
 * At the end point `full`: H = J^T J (+ priors), inverted (with a ridge when singular); the errors go into `mp`, the
 * poorly constrained directions are named on `log`, the covariance and its band block into `out`.
 */
void FinalizeLeastSquaresCovariance(const AZURECalc &self, const vector_r &full, const std::vector<int> &p2f,
                                    const std::vector<double> &pen_inv2, ROOT::Minuit2::MnUserParameters &mp,
                                    LeastSquaresOutcome &out, std::ostream &log);

/// A least-squares solver: from `params` to its end point (left in `params`), at most `maxIter` iterations.
using LeastSquaresSolver = LeastSquaresOutcome (*)(AZURECalc &self, AZUREParams &params, int maxIter,
                                                   const char *fitter, const FitControl &control);
LeastSquaresOutcome RunLevenbergMarquardt(AZURECalc &self, AZUREParams &params, int maxIter, const char *fitter,
                                          const FitControl &control);
LeastSquaresOutcome RunGSLNonlinear(AZURECalc &self, AZUREParams &params, int maxIter, const char *fitter,
                                    const FitControl &control);

/// The row function of a least-squares fitter: the solver, or MIGRAD (numerical) where it has no analytic Jacobian.
FitResult FitByLeastSquares(AZURECalc &objective, AZUREParams &params, const FitControl &control, const char *fitter,
                            LeastSquaresSolver solver);

#endif
