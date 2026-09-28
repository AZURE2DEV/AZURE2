#include "ShftFunc.h"
#include "Constants.h"
#include <algorithm>
#include <cmath>
#include <math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>

double ShftFunc::thisShftFunc(double x, void *p) {
  ShftFunc *shift = (ShftFunc *)p;
  int l = shift->params_.lValue;

  return shift->operator()(l, x);
}

/*!
 * Zero-energy limit of the shift function, the value both the positive- and the
 * negative-energy branches tend to at threshold.  At E = 0 the decaying
 * solution of the Coulomb equation is u(r) = sqrt(r) K_{2l+1}(sqrt(8 c r)),
 * c = eta k = Z1 Z2 alpha mu c^2/(hbar c) (fm^-1), so with x = sqrt(8 c a)
 *
 *   S_l(0) = a u'(a)/u(a) = -l - (x/2) K_{2l}(x) / K_{2l+1}(x);
 *
 * without Coulomb (c = 0) u = r^-l and S_l(0) = -l.
 */
double ShftFunc::ZeroEnergyLimit(int l) const {
  const WhitFunc &w = *params_.whitFunc;
  const double c = w.z1() * w.z2() * fstruc * w.redmass() * uconv / hbarc;
  if (!(c > 0.0)) return -(double)l;
  const double x = std::sqrt(8.0 * c * radius());
  gsl_sf_result kl, kl1;
  gsl_error_handler_t *oldHandler = gsl_set_error_handler_off();
  int s1 = gsl_sf_bessel_Kn_scaled_e(2 * l, x, &kl);
  int s2 = gsl_sf_bessel_Kn_scaled_e(2 * l + 1, x, &kl1);
  gsl_set_error_handler(oldHandler);
  if (s1 != GSL_SUCCESS || s2 != GSL_SUCCESS || kl1.val == 0.0) return -(double)l;
  return -(double)l - 0.5 * x * kl.val / kl1.val;
}

/*!
 * Below threshold, S = z W'(z)/W(z) at z = 2 kappa a for the Whittaker function
 * W_{-eta, l+1/2}, and z W'_{k,m} = (z/2 - k) W_{k,m} - W_{k+1,m} (DLMF 13.15.23)
 * turns that into a ratio of two confluent hypergeometric functions,
 *
 *   S = z/2 + eta - U(A-1, 2l+2, z) / U(A, 2l+2, z),   A = l + 1 + eta,
 *
 * evaluated with their powers of ten carried separately (gsl_sf_hyperg_U_e10_e).
 * Returns false if GSL cannot evaluate it.
 */
bool ShftFunc::WhittakerShift(int l, double binding, double &s) const {
  const WhitFunc &w = *params_.whitFunc;
  const double mu = w.redmass() * uconv;
  const double kappa = std::sqrt(2.0 * mu * binding) / hbarc;  // fm^-1
  const double eta = w.z1() * w.z2() * fstruc * mu / hbarc / kappa;
  const double z = 2.0 * kappa * radius();
  const double A = l + 1.0 + eta;
  const double b = 2.0 * l + 2.0;

  gsl_sf_result_e10 uA, uA1;
  gsl_error_handler_t *oldHandler = gsl_set_error_handler_off();
  int s1 = gsl_sf_hyperg_U_e10_e(A, b, z, &uA);
  int s2 = gsl_sf_hyperg_U_e10_e(A - 1.0, b, z, &uA1);
  gsl_set_error_handler(oldHandler);
  if (s1 != GSL_SUCCESS || s2 != GSL_SUCCESS || uA.val == 0.0 || !std::isfinite(uA.val) ||
      !std::isfinite(uA1.val))
    return false;
  s = 0.5 * z + eta - uA1.val / uA.val * std::pow(10.0, (double)(uA1.e10 - uA.e10));
  return std::isfinite(s);
}

/*!
 * The negative-energy shift function, continuous down to and at threshold.
 *
 * This used to be a numerical derivative of the plain-double Whittaker
 * function, which underflows to zero close to threshold, where eta grows like
 * 1/sqrt(binding): for 6Li+d it returned S = 0 for every binding below ~1e-5
 * MeV and NaN at round-off distance from threshold.  A THM sub-point landing
 * there (the Gaussian-folding grid of a THM segment crosses E = 0) put a NaN
 * into the level matrix and ~1e20 into the folded cross section.
 *
 * Now: the analytic log-derivative (WhittakerShift) where GSL's U is trusted
 * (bindingMin below).  The two U's there nearly cancel against eta, and beyond
 * it GSL's values lose the digits that survive the cancellation, so closer to
 * threshold S is interpolated linearly between S(-bindingMin) and the exact
 * zero-energy limit S(0): S is smooth through threshold, and the neglected
 * curvature is below 1e-8 absolute on the systems tried (p+7Li, 6Li+d,
 * 12C+12C; the slope within 1e-4 relative of dS/dE(0)).  At exactly zero
 * binding S(0) itself.
 */
double ShftFunc::operator()(int l, double energy) {
  const double binding = fabs(energy - totalSepE());
  if (!(binding > 0.0)) return ZeroEnergyLimit(l);

  // Where GSL's U is trusted: eta <= 300, and either eta <= 50 or z = 2 kappa a
  // >= 0.1.  Light pairs reach eta = 300 only at z ~ 1e-3, where the values
  // carry errors of ~1e-7 (p+7Li l=1 against mpmath: 4.4e-7 at B = 1e-6 MeV,
  // eta = 443; 5e-9 at eta = 81-140, z ~ 0.01; scattered by up to ~1e-9 with
  // the last bits of B at eta = 44).  An S good to 1e-7 is harmless, but the
  // slope of the interpolation below came out 0.60 instead of 0.373 MeV^-1
  // there -- and dS/dE at a level at threshold is that slope (now within 1e-5).
  // Heavy pairs keep eta = 300 (12C+12C at eta = 254, z = 0.35: 7.5e-11).
  const double kEtaMax = 300.0, kEtaSafe = 50.0, kZSafe = 0.1;
  const WhitFunc &w = *params_.whitFunc;
  const double mu = w.redmass() * uconv;
  const double c = w.z1() * w.z2() * fstruc * mu / hbarc;  // eta * kappa, fm^-1
  auto bindingAtKappa = [mu](double kappa) { return std::pow(hbarc * kappa, 2.0) / (2.0 * mu); };
  const double bindingMin =
      (c > 0.0) ? std::max(bindingAtKappa(c / kEtaMax),
                           std::min(bindingAtKappa(c / kEtaSafe), bindingAtKappa(kZSafe / (2.0 * radius()))))
                : 0.0;

  double s;
  if (binding >= bindingMin) {
    if (WhittakerShift(l, binding, s)) return s;
    return ZeroEnergyLimit(l);
  }
  const double s0 = ZeroEnergyLimit(l);
  if (!WhittakerShift(l, bindingMin, s)) return s0;
  return s0 + (s - s0) * (binding / bindingMin);
}

double ShftFunc::EnergyDerivative(int l, double energy, double step) {
  double result;
  double error;

  params_.lValue = l;

  gsl_function F;
  F.function = &thisShftFunc;
  F.params = this;

  gsl_deriv_central(&F, energy, step, &result, &error);

  return result;
}
