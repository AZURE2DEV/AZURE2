#include "ChannelFunc.h"
#include "Constants.h"
#include "CoulFunc.h"
#include "PPair.h"
#include "ShftFunc.h"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_deriv.h>

ChannelFunc::ChannelFunc(PPair *pair, bool useGSLFunctions)
    : pair_(pair),
      useGSL_(useGSLFunctions),
      threshold_(pair->GetSepE() + pair->GetExE()),
      radius_(pair->GetChRad()) {}

ChannelFunc::~ChannelFunc() = default;

CoulFunc &ChannelFunc::coul() {
  if (!coul_) coul_.reset(new CoulFunc(pair_, useGSL_));
  return *coul_;
}

ShftFunc &ChannelFunc::shft() {
  if (!shft_) shft_.reset(new ShftFunc(pair_));
  return *shft_;
}

bool ChannelFunc::IsCoulombThreshold(PPair *pair, double e) {
  const double z1z2 = pair->GetZ(1) * pair->GetZ(2);
  if (!(e > 0.0) || z1z2 == 0.0) return false;
  const double eta = std::sqrt(uconv / 2.) * fstruc * z1z2 * std::sqrt(pair->GetRedMass() / e);
  return eta > kEtaClosed;
}

bool ChannelFunc::IsClosed(double e) const {
  return !(e > 0.0) || IsCoulombThreshold(pair_, e);
}

double ChannelFunc::Shift(int l, double e) {
  if (e < 0.0) return shft()(l, e + threshold_);
  if (e == 0.0) return shft().ZeroEnergyLimit(l);
  if (IsCoulombThreshold(pair_, e))
    return 2.0 * shft().ZeroEnergyLimit(l) - shft()(l, threshold_ - e);
  return coul().PEShift(l, radius_, e);
}

double ChannelFunc::Penetrability(int l, double e) {
  if (IsClosed(e)) return 0.0;
  return coul().Penetrability(l, radius_, e);
}

double ChannelFunc::shiftAdaptor(double e, void *p) {
  ChannelFunc *self = static_cast<ChannelFunc *>(p);
  return self->Shift(self->lDeriv_, e);
}

double ChannelFunc::DerivativeStep(double e) {
  return std::max(kStep, std::min(kMaxStep, 0.25 * std::fabs(e)));
}

double ChannelFunc::ShiftDerivative(int l, double e) {
  // gsl_deriv_central probes e +- h and e +- h/2.
  const double h = DerivativeStep(e);
  // Whole stencil in one branch: the routine AZURE2 always used there.
  if (e + h < 0.0) return shft().EnergyDerivative(l, e + threshold_, h);
  if (e - h > 0.0 && !IsCoulombThreshold(pair_, e - h))
    return coul().PEShift_dE(l, radius_, e, h);
  // Coulomb threshold: S(e) = 2 S(0) - S(-e), so dS/dE(e) = dS/dE(-e).
  if (e > 0.0 && IsCoulombThreshold(pair_, e)) return ShiftDerivative(l, -e);
  // A stencil across threshold (or across the eta = 100 line): difference the
  // continuous S.
  lDeriv_ = l;
  gsl_function F;
  F.function = &shiftAdaptor;
  F.params = this;
  double result = 0.0, error = 0.0;
  gsl_deriv_central(&F, e, kStep, &result, &error);
  return result;
}
