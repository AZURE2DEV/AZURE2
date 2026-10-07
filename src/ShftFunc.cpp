#include "ShftFunc.h"
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
