#include "ChannelFunc.h"
#include "Constants.h"
#include "CoulFunc.h"
#include "PPair.h"
#include "ShftFunc.h"
#include <algorithm>
#include <cmath>
#include "azure2/physics/Derivatives.h"

ChannelFunc::ChannelFunc(const PPair *pair, bool useGSLFunctions)
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

bool ChannelFunc::IsCoulombThreshold(const PPair *pair, double e) {
  const double z1z2 = pair->GetZ(1) * pair->GetZ(2);
  if (!(e > 0.0) || z1z2 == 0.0) return false;
  const double eta = std::sqrt(uconv / 2.) * fstruc * z1z2 * std::sqrt(pair->GetRedMass() / e);
  return eta > kEtaClosed;
}

bool ChannelFunc::IsClosed(double e) const {
  return !(e > 0.0) || IsCoulombThreshold(pair_, e);
}

double ChannelFunc::Shift(int l, int two_j, double e) {
  if (e < 0.0) return shft()(l, two_j, e + threshold_);
  if (e == 0.0) return shft().ZeroEnergyLimit(l, two_j);
  if (IsCoulombThreshold(pair_, e))
    return 2.0 * shft().ZeroEnergyLimit(l, two_j) - shft()(l, two_j, threshold_ - e);
  return coul().PEShift(l, two_j, radius_, e);
}

double ChannelFunc::Penetrability(int l, int two_j, double e) {
  if (IsClosed(e)) return 0.0;
  return coul().Penetrability(l, two_j, radius_, e);
}

double ChannelFunc::DerivativeStep(double e) {
  return std::max(kStep, std::min(kMaxStep, 0.25 * std::fabs(e)));
}

double ChannelFunc::ShiftDerivative(int l, int two_j, double e) {
  // central_derivative (GSL's gsl_deriv_central) probes e +- h and e +- h/2.
  const double h = DerivativeStep(e);
  // Whole stencil in one branch: the routine AZURE2 always used there.
  if (e + h < 0.0) return shft().EnergyDerivative(l, two_j, e + threshold_, h);
  if (e - h > 0.0 && !IsCoulombThreshold(pair_, e - h))
    return coul().PEShift_dE(l, two_j, radius_, e, h);
  // Coulomb threshold: S(e) = 2 S(0) - S(-e), so dS/dE(e) = dS/dE(-e).
  if (e > 0.0 && IsCoulombThreshold(pair_, e)) return ShiftDerivative(l, two_j, -e);
  // A stencil across threshold (or across the eta = 100 line): difference the
  // continuous S (central_derivative: GSL's gsl_deriv_central, ported to physics/).
  return azure2::physics::central_derivative([this, l, two_j](double x) { return Shift(l, two_j, x); }, e, kStep);
}
