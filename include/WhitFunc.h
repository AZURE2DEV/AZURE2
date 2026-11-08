#ifndef WHITFUNC_H
#define WHITFUNC_H

#include "PPair.h"
#include <cmath>
#include "azure2/physics/coulomb/Whittaker.h"
#include "Constants.h"
#include "GSLException.h"
#include "Kinematics.h"

#include <iostream>

extern double gsl_whit_function(int, double, double, double, int, int);

/// A function class to calculate Whittaker functions for negative energy channels.

/*!
 * The function class WhitFunc evaluates the Whittaker function W_{-eta_b, l+1/2}(2 kappa r) of a negative-energy
 * channel at the engine's eta_b and kappa r (non-relativistic, or relativistic under "relativistic=1"); the function
 * itself is physics/coulomb's whittaker / whittaker_scaled (no GSL; docs/upgrade/coulomb.md section 1.4).
 */

class WhitFunc {
 public:
  /*!
   * The WhitFunc object is created with reference to a PPair object.
   */
  WhitFunc(PPair *pPair) {
    z1_ = pPair->GetZ(1);
    z2_ = pPair->GetZ(2);
    redmass_ = (double)pPair->GetRedMass();
    kin_ = Kinematics::Channel(pPair);
  };
  /*!
   * Returns the atomic number of the first particle in the pair.
   */
  int z1() const {
    return z1_;
  };
  /*!
   * Returns the atomic number of the second particle in the pair.
   */
  int z2() const {
    return z2_;
  };
  /*!
   * Returns the reduced mass of the particle pair.
   */
  double redmass() const {
    return redmass_;
  };
  /*!
   * The parenthesis operator is defined to make the class instance callable as a function.  The orbital
   * angular momentum, binding energy, and radius are the dependent variables.
   * The function returns the value of the Whittaker function.
   */
  double operator()(int l, double radius, double energy) const {
    double etaB, rho;
    if (!Arguments(radius, energy, etaB, rho)) return 0;  // NaN arguments: 0, as before (FIXME upstream)
    namespace co = azure2::physics::coulomb;
    const co::WhittakerWave w = co::whittaker(l, etaB, rho);
    if (w.status == co::Status::ok) return w.W;
    // Outside the double range: the IEEE value of mantissa * 10^scale (0 on underflow, inf on overflow).
    const co::WhittakerScaled ws = co::whittaker_scaled(l, etaB, rho);
    if (ws.status != co::Status::ok) return 0;
    return ws.mantissa * std::pow(10.0, (double)ws.log10_scale);
  };
  /*!
   * Whittaker function in scaled form: returns the mantissa and sets
   * log10Scale so that W = mantissa * 10^log10Scale.
   *
   * Just below a channel threshold the confluent hypergeometric argument a
   * grows like sqrt(1/E) (a ~ 255 only 200 eV below threshold), and U(a,b,z)
   * behaves like 1/Gamma(a) -- around 1e-501 here.  A plain double flushes
   * that to exactly zero, which is what makes the sub-threshold overlap
   * integral evaluate 0/0.  The ratios these values appear in are perfectly
   * well conditioned, so the exponent has to be carried rather than discarded:
   * whittaker_scaled sums in logarithms and returns the power of ten separately.
   */
  double Scaled(int l, double radius, double energy, double &log10Scale) const {
    log10Scale = 0.0;
    double etaB, rho;
    if (!Arguments(radius, energy, etaB, rho) || !(rho > 0.0)) return 0.0;
    namespace co = azure2::physics::coulomb;
    const co::WhittakerScaled ws = co::whittaker_scaled(l, etaB, rho);
    if (ws.status != co::Status::ok) return 0.0;
    log10Scale = (double)ws.log10_scale;
    return ws.mantissa;  // in [1, 10): W > 0 for every valid argument
  };

 private:
  /// eta_b and rho = kappa r at binding energy `energy` (MeV) and `radius` (fm), as the engine always formed them
  /// (z = 2 kappa r); false if either is NaN.
  bool Arguments(double radius, double energy, double &etaB, double &rho) const {
    etaB = kin_.rel ? kin_.EtaBound(energy) : sqrt(uconv / 2.) * fstruc * z1() * z2() * sqrt(redmass() / energy);
    const double z = kin_.rel ? 2.0 * kin_.RhoBound(radius, energy) : 2.0 * sqrt(2.0 * uconv) / hbarc * radius * sqrt(redmass() * energy);
    rho = z / 2.0;
    return etaB == etaB && rho == rho;
  }
  int z1_;
  int z2_;
  double redmass_;
  Kinematics::Channel kin_;
};

#endif
