#include "ShftFunc.h"
#include "Config.h"
#include "Constants.h"
#include "CoulombTable.h"
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_deriv.h>

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

double ShftFunc::theWhitFunc(double x, void *p) {
  Params *params = (Params *)p;
  int l = (params->lValue);
  double bindingenergy = (params->bindingEnergy);
  WhitFunc *whitFunc = (params->whitFunc.get());

  return whitFunc->operator()(l, x, bindingenergy);
}

double ShftFunc::operator()(int l, double energy) {
  params_.bindingEnergy = fabs(energy - totalSepE());
  params_.lValue = l;

  // Neutral channel: the Whittaker function is rho k_l(rho) up to a constant,
  // whose logarithmic derivative has the closed form below.  The general route
  // takes a numerical derivative of GSL's hyperg_U (6-10 calls per value).
  if (neutral_ && UseNeutralClosedForm()) {
    const double rho = kin_.rel ? kin_.RhoBound(radius(), params_.bindingEnergy)
                               : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * params_.bindingEnergy);
    if (rho > 0.0) return NeutralBoundShift(l, rho);
  }
  // Charged channel: the table's Whittaker branch, at eta_b and rho_b = kappa a,
  // else the exact S = z W'/W from U' = -a U(a+1, b+1, z) (CoulombTable::ExactBound).
  // The numerical radial derivative of W below is the fallback, and with
  // AZURE_BOUND_SHIFT_NUMERICAL=1 the route used before 2026-10-19: its noise makes
  // the 1 eV central difference dS/dE wrong by 3e-4 - 4e-3 (p+12C).
  if (!neutral_ && params_.bindingEnergy > 0.0 && (table_ || UseURatio())) {
    const double etaB = kin_.rel ? kin_.EtaBound(params_.bindingEnergy)
                                 : sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / params_.bindingEnergy);
    const double rhoB = kin_.rel ? kin_.RhoBound(radius(), params_.bindingEnergy)
                                 : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * params_.bindingEnergy);
    double s;
    if (table_ && table_->LookupBound(l, etaB, rhoB, s)) return s;
    if (UseURatio()) {
      s = CoulombTable::ExactBound(l, etaB, rhoB);
      if (std::isfinite(s)) return s;
    }
  }

  double result;
  double error;

  gsl_function F;
  F.function = &theWhitFunc;
  F.params = &params_;

  gsl_deriv_central(&F, radius(), 1e-4, &result, &error);

  // std::cout << "DEBUG: " << radius()*result/theWhitFunc(radius(),&params_) << " " << radius() << " " << result << "  " << theWhitFunc(radius(),&params_) << std::endl;

  // Prevent division by zero
  if (theWhitFunc(radius(), &params_) == 0.) return 0.0;

  return radius() * result / theWhitFunc(radius(), &params_);
}

double ShftFunc::EnergyDerivative(int l, double energy) {
  double result;
  double error;

  params_.lValue = l;

  gsl_function F;
  F.function = &thisShftFunc;
  F.params = this;

  gsl_deriv_central(&F, energy, 1e-6, &result, &error);

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
