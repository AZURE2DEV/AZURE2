// IntegratedFermiFunc near the beta end point.  The integrand multiplied exp(pi eta), which overflows, by
// |Gamma(gamma0 + i eta)|^2, which underflows, once the electron momentum is below ~2.6e-4 m_e c; the quadrature
// samples such momenta whenever the end-point kinetic energy T0 is below ~0.2 keV, and inf * 0 made the integrated
// Fermi function NaN there.  A beta-delayed spectrum point, or a resolution sub-point, within 0.2 keV of the decay's
// upper end then came out NaN: the 12C(a,g) MCMC's 16N spectrum with 5 sigma resolution windows (2026-10-29).
//
// Checks, for beta-minus to Z = 8 with R = 6 fm (the 16N -> 16O + e- pair of tests/beta_delayed_alpha):
//   finite and positive from T0 = 1e-7 keV to 1 MeV, and increasing with T0;
//   the values where the old form was finite are unchanged (reference values from the old form, 7 digits);
//   f ~ T0^3 at the end point: f / T0^3 changes by less than 2 % between T0 = 1e-5 and 1e-1 keV;
//   no end point (T0 <= 0) gives 0.
#include "IntegratedFermiFunc.h"
#include <gsl/gsl_errno.h>
#include <cmath>
#include <cstdio>
static int failures = 0;
static void check(bool ok, const char *what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}
static double W0of(double keV) { return 1. + keV * 1e-3 / 0.510998903; }
int main() {
  gsl_set_error_handler_off();
  IntegratedFermiFunc f(-1);
  const double Z = 8., R = 6.;
  bool finite = true, increasing = true;
  double last = 0.;
  for (double lg = -7.; lg <= 3.0001; lg += 0.25) {
    const double v = f(W0of(std::pow(10., lg)), Z, R);
    if (!(std::isfinite(v) && v > 0.)) {
      finite = false;
      std::printf("        T0 = %.3g keV: f = %g\n", std::pow(10., lg), v);
    }
    if (!(v > last)) increasing = false;
    last = v;
  }
  check(finite, "finite and positive for T0 = 1e-7 keV to 1 MeV (NaN below ~0.2 keV before)");
  check(increasing, "increasing with the end-point energy");
  const struct { double keV, value; } ref[] = {{0.3, 2.542379e-11}, {1.0, 9.423127e-10}, {10.0, 9.791124e-07},
                                               {100.0, 1.440678e-03}, {1000.0, 5.646789e+00}};
  bool same = true;
  for (const auto &r : ref) {
    const double v = f(W0of(r.keV), Z, R);
    if (std::fabs(v / r.value - 1.) > 1e-6) {
      same = false;
      std::printf("        T0 = %g keV: %.7e, was %.7e\n", r.keV, v, r.value);
    }
  }
  check(same, "unchanged where the old form was finite (0.3 keV - 1 MeV, to 1e-6)");
  const double a = f(W0of(1e-5), Z, R) / 1e-15, b = f(W0of(1e-1), Z, R) / 1e-3;
  check(std::fabs(a / b - 1.) < 0.02, "f ~ T0^3 at the end point (f/T0^3 within 2 % from 1e-5 to 0.1 keV)");
  check(f(1.0, Z, R) == 0. && f(0.999, Z, R) == 0., "no end point: 0");
  std::printf(failures ? "FAILED: %d\n" : "PASS\n", failures);
  return failures ? 1 : 0;
}
