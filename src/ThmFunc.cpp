#include "ThmFunc.h"
#include "Constants.h"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>

//! forward finite-difference step for dj_l/drho (as in mrmpy's model.py)
static const double kThmFDStep = 1.0e-6;

double ThmSphericalBessel(int l, double x) {
  return gsl_sf_bessel_jl(l, x);
}

double ThmRho(double mu, double E, double B, double radius) {
  // p = sqrt(2 mu (E + B)) in MeV/c; rho = p r / (hbar c), dimensionless.
  return std::sqrt(2.0 * mu * (E + B)) * radius / hbarc;
}

double ThmFormFactor(int l, double b, double mu, double E, double B,
                     double radius) {
  double rho = ThmRho(mu, E, B, radius);
  double jl = gsl_sf_bessel_jl(l, rho);
  // dj_l/drho by forward finite difference (matches mrmpy), then rho*dj_l/drho.
  double rho_djl = rho * (gsl_sf_bessel_jl(l, rho + kThmFDStep) - jl) / kThmFDStep;
  return (b - 1.0) * jl - rho_djl;
}
