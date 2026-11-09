// Numerical derivatives (Derivatives.h).  central_derivative follows GSL's gsl_deriv_central (deriv/deriv.c,
// Copyright (C) 2004, 2007 Brian Gough, GPL-3) operation by operation, so that it reproduces its bits.
#include "azure2/physics/Derivatives.h"

#include <cfloat>
#include <cmath>

namespace azure2 {
namespace physics {

namespace {

void CentralDeriv(const std::function<double(double)> &f, double x, double h, double *result, double *abserr_round,
                  double *abserr_trunc) {
  // The 5-point rule from the 3-point rule at h and at h/2 (Richardson); the error estimates as GSL's.
  const double fm1 = f(x - h);
  const double fp1 = f(x + h);
  const double fmh = f(x - h / 2);
  const double fph = f(x + h / 2);
  const double r3 = 0.5 * (fp1 - fm1);
  const double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;
  const double e3 = (std::fabs(fp1) + std::fabs(fm1)) * DBL_EPSILON;
  const double e5 = 2.0 * (std::fabs(fph) + std::fabs(fmh)) * DBL_EPSILON + e3;
  // Rounding in x + h and x - h.
  const double dy = std::fmax(std::fabs(r3 / h), std::fabs(r5 / h)) * (std::fabs(x) / h) * DBL_EPSILON;
  *result = r5 / h;
  *abserr_trunc = std::fabs((r5 - r3) / h);
  *abserr_round = std::fabs(e5 / h) + dy;
}

}  // namespace

double central_derivative(const std::function<double(double)> &f, double x, double h, double *abserr) {
  double r_0, round, trunc, error;
  CentralDeriv(f, x, h, &r_0, &round, &trunc);
  error = round + trunc;
  if (round < trunc && (round > 0 && trunc > 0)) {
    double r_opt, round_opt, trunc_opt, error_opt;
    // The step that minimizes the total error, when rounding is not dominant.
    const double h_opt = h * std::pow(round / (2.0 * trunc), 1.0 / 3.0);
    CentralDeriv(f, x, h_opt, &r_opt, &round_opt, &trunc_opt);
    error_opt = round_opt + trunc_opt;
    // Kept only if its error is smaller and it agrees within the old error.
    if (error_opt < error && std::fabs(r_opt - r_0) < 4.0 * error) {
      r_0 = r_opt;
      error = error_opt;
    }
  }
  if (abserr) *abserr = error;
  return r_0;
}

double second_derivative(const std::function<double(double)> &f, double x, double h, double g, bool one_sided,
                         Side side) {
  if (one_sided) {
    const bool below = (side == Side::Below);
    const double s0 = f(x);
    const double s1 = below ? f(x - g) : f(x + g);
    const double s2 = below ? f(x - 2.0 * g) : f(x + 2.0 * g);
    const double s3 = below ? f(x - 3.0 * g) : f(x + 3.0 * g);
    return (2.0 * s0 - 5.0 * s1 + 4.0 * s2 - s3) / (g * g);
  }
  const double s0 = f(x);
  const double d1 = (f(x + h) - 2.0 * s0 + f(x - h)) / (h * h);
  const double d2 = (f(x + 2.0 * h) - 2.0 * s0 + f(x - 2.0 * h)) / (4.0 * h * h);
  return (4.0 * d1 - d2) / 3.0;
}

}  // namespace physics
}  // namespace azure2
