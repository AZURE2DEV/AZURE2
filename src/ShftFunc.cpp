#include "ShftFunc.h"
#include "Config.h"
#include "Constants.h"
#include "CoulombTable.h"
#include "GSLException.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
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

double ShftFunc::thisShftFuncDeriv(double x, void *p) {
  ShftFunc *shift = (ShftFunc *)p;
  int l = shift->params_.lValue;

  return shift->EnergyDerivative(l, x);
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
  int s1, s2;
  {
    GslQuiet quiet;  // statuses checked here; thread-safe (GSLException.h)
    s1 = gsl_sf_bessel_Kn_scaled_e(2 * l, x, &kl);
    s2 = gsl_sf_bessel_Kn_scaled_e(2 * l + 1, x, &kl1);
  }
  if (s1 != GSL_SUCCESS || s2 != GSL_SUCCESS || kl1.val == 0.0) return -(double)l;
  return -(double)l - 0.5 * x * kl.val / kl1.val;
}

double ShftFunc::theWhitFunc(double x, void *p) {
  Params *params = (Params *)p;
  return params->whitFunc->operator()(params->lValue, x, params->bindingEnergy);
}

/*!
 * Below threshold, S = z W'(z)/W(z) at z = 2 kappa a for the Whittaker function
 * W_{-eta, l+1/2} = exp(-z/2) z^(l+1) U(a, b, z), a = l + 1 + eta, b = 2l + 2,
 * and U' = -a U(a+1, b+1, z):
 *
 *   S = l + 1 - z/2 - a z U(a+1, b+1, z) / U(a, b, z)
 *
 * (CoulombTable::ExactBound), with the powers of ten of the two U's carried
 * separately (gsl_sf_hyperg_U_e10_e).  No large terms cancel, unlike the
 * U(a-1, b, z)/U(a, b, z) form used on thm before the merge of October 2026
 * (z/2 + eta - ratio, both ~eta), which was 10-100 times less accurate near
 * threshold.  Returns false if GSL cannot evaluate it.
 */
bool ShftFunc::WhittakerShift(int l, double binding, double &s) const {
  const double etaB = sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / binding);
  const double rhoB = sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
  s = CoulombTable::ExactBound(l, etaB, rhoB);
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
 * Now, in order: S(0) at exactly zero binding; the closed form of a neutral
 * channel (NeutralBoundShift); the Coulomb table's bound sheet when a table is
 * in use and covers (eta_b, rho_b); the exact U ratio (WhittakerShift) where
 * GSL's U is trusted (bindingMin below), and closer to threshold a linear
 * interpolation between S(-bindingMin) and S(0), so S is smooth through
 * threshold.
 */
double ShftFunc::operator()(int l, double energy) {
  const double binding = fabs(energy - totalSepE());
  if (!(binding > 0.0)) return ZeroEnergyLimit(l);
  params_.bindingEnergy = binding;
  params_.lValue = l;

  // Neutral channel: the Whittaker function is rho k_l(rho) up to a constant,
  // whose logarithmic derivative has the closed form below.
  if (neutral_ && UseNeutralClosedForm()) {
    const double rho = sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
    if (rho > 0.0) return NeutralBoundShift(l, rho);
  }
  // Charged channel: the table's Whittaker branch, at eta_b and rho_b = kappa a.
  if (!neutral_ && table_) {
    const double etaB = sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / binding);
    const double rhoB = sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
    double s;
    if (table_->LookupBound(l, etaB, rhoB, s)) return s;
  }
  // AZURE_BOUND_SHIFT_NUMERICAL=1: the route used before 2026-10-19, a numerical
  // radial derivative of W, whose noise makes a 1 eV central difference dS/dE
  // wrong by 3e-4 - 4e-3 (p+12C); for comparisons only.
  if (!UseURatio()) return NumericalShift();

  // Where GSL's U is trusted (the exact U ratio, WhittakerShift): eta <= 1000,
  // and either eta <= 50 or z = 2 kappa a >= 0.05.  Against mpmath (October
  // 2026; p+p, p+7Li, 6Li+d, p+12C, p+18O, a+12C, a+16O, 12C+12C, l = 0-3,
  // binding 1e-10 - 3 MeV) the ratio is good to ~1e-9 there, and outside it it
  // can be off by 1e-6 (p+p, eta = 111, z = 1.2e-3) or fail; closer to
  // threshold S is interpolated linearly between S(-bindingMin) and the exact
  // zero-energy limit S(0), so it is smooth through threshold.  Largest error
  // of the whole function on that set: 8.7e-9 absolute (it was 2.2e-8 with
  // the U(A-1)/U(A) form and eta <= 300).
  const double kEtaMax = 1000.0, kEtaSafe = 50.0, kZSafe = 0.05;
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

/*!
 * The numerical radial derivative of the Whittaker function (the route before
 * 2026-10-19), at params_.bindingEnergy and params_.lValue.
 */
double ShftFunc::NumericalShift() {
  double result;
  double error;

  gsl_function F;
  F.function = &theWhitFunc;
  F.params = &params_;

  gsl_deriv_central(&F, radius(), 1e-4, &result, &error);

  // Prevent division by zero
  if (theWhitFunc(radius(), &params_) == 0.) return 0.0;

  return radius() * result / theWhitFunc(radius(), &params_);
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

double ShftFunc::EnergySecondDerivative(int l, double energy) {
  // Richardson-extrapolated central second difference (O(h^4), 10 keV step),
  // kept below the threshold where the Whittaker form holds: a backward
  // stencil when the level is within 4h of it.
  const double h = 1.0e-2;
  double binding = totalSepE() - energy;
  if (binding < 4.0 * h) {
    const double f = 1.0e-3;
    double s0 = operator()(l, energy);
    double s1 = operator()(l, energy - f);
    double s2 = operator()(l, energy - 2.0 * f);
    double s3 = operator()(l, energy - 3.0 * f);
    return (2.0 * s0 - 5.0 * s1 + 4.0 * s2 - s3) / (f * f);
  }
  double s0 = operator()(l, energy);
  double d1 = (operator()(l, energy + h) - 2.0 * s0 + operator()(l, energy - h)) / (h * h);
  double d2 = (operator()(l, energy + 2.0 * h) - 2.0 * s0 + operator()(l, energy - 2.0 * h)) / (4.0 * h * h);
  return (4.0 * d1 - d2) / 3.0;
}

const CoulombTable *ShftFunc::SharedTable(bool neutral) {
  if (neutral || !g_config || !g_config->useCoulombTable) return nullptr;
  return CoulombTable::Shared(g_config->coulombTablePath);
}

bool ShftFunc::UseURatio() {
  static const bool use = (std::getenv("AZURE_BOUND_SHIFT_NUMERICAL") == nullptr);
  return use;
}

bool ShftFunc::UseNeutralClosedForm() {
  static const bool use = (std::getenv("AZURE_NO_NEUTRAL_CLOSED_FORM") == nullptr);
  return use;
}

double ShftFunc::NeutralBoundShift(int l, double rho) {
  double s = -rho;
  for (int ll = 1; ll <= l; ll++) s = -rho * rho / (ll - s) - ll;
  return s;
}
