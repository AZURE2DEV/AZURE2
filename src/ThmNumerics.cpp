#include "ThmNumerics.h"

#include <cmath>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include "azure2/physics/coulomb/Whittaker.h"

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
  // W_{-eta_b, 1/2}(2 kappa r) / r = 2 kappa e^{-kappa r} U(1 + eta_b, 2, 2 kappa r) (physics/coulomb, no GSL); outside
  // the double range the IEEE value of the scaled form (0 on underflow), 0 for invalid arguments.
  namespace co = azure2::physics::coulomb;
  const co::WhittakerWave w = co::whittaker(0, etaB, kappa * r);
  if (w.status == co::Status::ok) return w.W / r;
  const co::WhittakerScaled ws = co::whittaker_scaled(0, etaB, kappa * r);
  return ws.status == co::Status::ok ? ws.mantissa * std::pow(10.0, (double)ws.log10_scale) / r : 0.0;
}
