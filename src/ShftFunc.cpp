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
#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Derivatives.h"
#include "azure2/physics/HybridPotential.h"
#include <gsl/gsl_errno.h>

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
double ShftFunc::ZeroEnergyLimit(int l, int two_j) const {
  const WhitFunc &w = *params_.whitFunc;
  if (potential_)  // the zero-energy decaying solution of the radial equation with the potential
    return azure2::physics::HybridBoundShift(*potential_, w.z1(), w.z2(), w.redmass() * uconv, l, pair_->SpinOrbitFactor(l, two_j), radius(), 0.0);
  const double c = w.z1() * w.z2() * fstruc * w.redmass() * uconv / hbarc;
  return azure2::physics::zero_energy_shift(c, radius(), l);  // K_n without GSL (coulomb.md section 1.4)
}

double ShftFunc::theWhitFunc(double x, void *p) {
  Params *params = (Params *)p;
  return params->whitFunc->operator()(params->lValue, params->twoJ, x, params->bindingEnergy);
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
 * in use and covers (eta_b, rho_b); the exact ratio rho W'/W (physics/ coulomb_bound_shift) inside
 * eta_b <= 1000 (bindingMin below), and closer to threshold a linear
 * interpolation between S(-bindingMin) and S(0), so S is smooth through
 * threshold.
 */
double ShftFunc::operator()(int l, int two_j, double energy) {
  const double binding = fabs(energy - totalSepE());
  if (!(binding > 0.0)) return ZeroEnergyLimit(l, two_j);
  params_.bindingEnergy = binding;
  params_.lValue = l;
  params_.twoJ = two_j;
  if (potential_) {
    const WhitFunc &w = *params_.whitFunc;
    return azure2::physics::HybridBoundShift(*potential_, w.z1(), w.z2(), w.redmass() * uconv, l, pair_->SpinOrbitFactor(l, two_j), radius(), binding);
  }

  // Neutral channel: the Whittaker function is rho k_l(rho) up to a constant,
  // whose logarithmic derivative has the closed form below.
  if (neutral_ && UseNeutralClosedForm()) {
    const double rho = kin_.rel ? kin_.RhoBound(radius(), binding)
                               : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
    if (rho > 0.0) return NeutralBoundShift(l, rho);
  }
  // Charged channel: the table's Whittaker branch, at eta_b and rho_b = kappa a (relativistic under
  // "relativistic=1", Kinematics.h), unless EnergyDerivative's stencil leaves the table (forceExact_).
  if (!neutral_ && table_ && !forceExact_) {
    const double etaB = kin_.rel ? kin_.EtaBound(binding) : sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / binding);
    const double rhoB = kin_.rel ? kin_.RhoBound(radius(), binding)
                                 : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
    double s;
    if (table_->LookupBound(l, etaB, rhoB, s, pair_->FunctionSettings().tableUse)) return s;
  }
  // AZURE_BOUND_SHIFT_NUMERICAL=1: the route used before 2026-10-19, a numerical
  // radial derivative of W, whose noise makes a 1 eV central difference dS/dE
  // wrong by 3e-4 - 4e-3 (p+12C); for comparisons only.
  if (!UseURatio()) return NumericalShift();

  // The exact ratio within eta_b <= 1000 and the linear interpolation to S(0) closer to threshold: physics/
  // coulomb_bound_shift (coulomb.md section 1.4), with this channel's eta_b and rho_b.
  const WhitFunc &w = *params_.whitFunc;
  azure2::physics::BoundShiftArgs args;
  args.mu_MeV = w.redmass() * uconv;
  args.c_fm = w.z1() * w.z2() * fstruc * args.mu_MeV / hbarc;  // eta * kappa, fm^-1
  args.a_fm = radius();
  args.eta_rho = [this](double B, double &etaB, double &rhoB) {
    etaB = kin_.rel ? kin_.EtaBound(B) : sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / B);
    rhoB = kin_.rel ? kin_.RhoBound(radius(), B) : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * B);
  };
  return azure2::physics::coulomb_bound_shift(l, binding, args);
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

bool ShftFunc::BoundTableOk(int l, double energy) const {
  if (!table_ || neutral_) return false;
  const double binding = fabs(energy - totalSepE());  // as operator() sets params_.bindingEnergy
  if (!(binding > 0.0)) return false;
  const double etaB = kin_.rel ? kin_.EtaBound(binding) : sqrt(uconv / 2.) * fstruc * zz_ * sqrt(redmass_ / binding);
  const double rhoB = kin_.rel ? kin_.RhoBound(radius(), binding) : sqrt(2. * uconv) / hbarc * radius() * sqrt(redmass_ * binding);
  double s;
  return table_->LookupBound(l, etaB, rhoB, s, pair_->FunctionSettings().tableUse);
}

double ShftFunc::EnergyDerivative(int l, int two_j, double energy, double step) {
  params_.lValue = l;
  params_.twoJ = two_j;

  // central_derivative (GSL's gsl_deriv_central, ported to physics/) samples E +- h and E +- h/2: one route for all four (the table only when it covers
  // every one of them), so the table's edge cannot put a jump into the derivative (dev, audit 2026-10-25,
  // item 32).  The step is the caller's (ChannelFunc::DerivativeStep for a closed channel; 1 eV by default).
  const double h = step;
  forceExact_ = table_ && !(BoundTableOk(l, energy + h) && BoundTableOk(l, energy - h) &&
                            BoundTableOk(l, energy + 0.5 * h) && BoundTableOk(l, energy - 0.5 * h));
  const double result = azure2::physics::central_derivative([this, l, two_j](double e) { return operator()(l, two_j, e); }, energy, h);
  forceExact_ = false;

  return result;
}

double ShftFunc::EnergySecondDerivative(int l, int two_j, double energy) {
  // Richardson-extrapolated central second difference (O(h^4), 10 keV step),
  // kept below the threshold where the Whittaker form holds: a backward
  // stencil when the level is within 4h of it.
  const double h = 1.0e-2;
  double binding = totalSepE() - energy;
  // The whole stencil by the exact route: the table's Hermite interpolant is C1 across its cells, and a 10 keV
  // stencil crosses them, which put a kink of the interpolation error into the second difference (audit
  // 2026-10-25, 03 F8; 2026-10-28).  Five exact bound-state evaluations per level are cheap.
  struct ExactScope {
    bool &flag; bool saved;
    ExactScope(bool &f) : flag(f), saved(f) { flag = true; }
    ~ExactScope() { flag = saved; }
  } exact(forceExact_);
  return azure2::physics::second_derivative([this, l, two_j](double e) { return operator()(l, two_j, e); }, energy, h, 1.0e-3,
                                            binding < 4.0 * h, azure2::physics::Side::Below);
}

const CoulombTable *ShftFunc::SharedTable(const PPair *pair, bool neutral) {
  const PairFunctionSettings &settings = pair->FunctionSettings();
  if (neutral || !settings.useCoulombTable) return nullptr;
  return CoulombTable::Shared(settings.coulombTablePath, settings.coulombTableDefault);
}

bool ShftFunc::UseURatio() const { return !pair_->FunctionSettings().boundShiftNumerical; }

bool ShftFunc::UseNeutralClosedForm() const { return !pair_->FunctionSettings().noNeutralClosedForm; }

double ShftFunc::NeutralBoundShift(int l, double rho) {
  return azure2::physics::neutral_bound_shift(l, rho);
}
