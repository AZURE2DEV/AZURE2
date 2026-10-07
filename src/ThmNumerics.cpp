#include "ThmNumerics.h"

#include <cmath>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

std::vector<complex> ThmCoulombPhases(double eta, int lmax) {
  std::vector<complex> out(lmax + 1, complex(1.0, 0.0));
  if (eta == 0.0) return out;
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e(1.0, eta, &lnr, &arg);
  double sigma = arg.val;
  for (int l = 0; l <= lmax; l++) {
    if (l > 0) sigma += std::atan(eta / l);
    out[l] = std::polar(1.0, sigma);
  }
  return out;
}

std::string ThmNumberText(double x) {
  std::ostringstream s;
  s.precision(6);
  s << x;
  return s.str();
}

double ThmBoundStateTail(double r, double kappa, double etaB, bool yukawa) {
  if (yukawa || etaB == 0.0) return std::exp(-kappa * r) / r;
  gsl_sf_result U;
  int status = gsl_sf_hyperg_U_e(1.0 + etaB, 2.0, 2.0 * kappa * r, &U);
  return status == GSL_SUCCESS || status == GSL_EUNDRFLW ? 2.0 * kappa * std::exp(-kappa * r) * U.val : 0.0;
}
