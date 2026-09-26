#ifndef SHIFTFUNC_H
#define SHIFTFUNC_H

#include "PPair.h"
#include "WhitFunc.h"
#include <memory>

class ShftFunc;

/// A function class for negative energy shift functions.

/*!
 * A shift function for negative energy channels is calculated as
 * \f$ S=\rho \frac{W_c'(k\rho)}{W_c(k\rho)} \f$, where the prime indicates
 * the derivative with respect to \f$ \rho \f$.  The derivative is taken
 * analytically through the Whittaker recurrence (see operator()), and at
 * threshold the zero-energy limit is used (ZeroEnergyLimit()).
 */

class ShftFunc {
 public:
  /*!
   * The ShftFunc object is created with reference to a particle pair.
   */
  ShftFunc(PPair *pPair) {
    totalSepE_ = pPair->GetSepE() + pPair->GetExE();
    radius_ = (double)pPair->GetChRad();
    params_.whitFunc = std::make_unique<WhitFunc>(pPair);
  };
  /*!
   *
   */
  ~ShftFunc() = default;
  /*!
   * The parenthesis operator is defined to make the class instance callable as a function.  The orbital
   * angular momentum and energy in the compound system are the dependent variables.
   * The function returns the value of the shift function.
   */
  double operator()(int l, double energy);
  /*!
   * Returns the energy derivative of the shift function at the specified orbital
   * angular momentum and energy in the compound system.
   */
  double EnergyDerivative(int l, double energy);
  /*!
   * The threshold (E = 0) value S_l(0), common limit of the positive- and
   * negative-energy shift functions: -l - (x/2) K_{2l}(x)/K_{2l+1}(x) with
   * x = sqrt(8 eta k a), or -l without Coulomb.
   */
  double ZeroEnergyLimit(int l) const;
  /*!
   * S = z W'/W below threshold from the Whittaker recurrence, at a binding
   * energy (MeV, > 0).  False if GSL cannot evaluate it.
   */
  bool WhittakerShift(int l, double binding, double &s) const;

 private:
  /// Separation plus excitation energy of the pair, the threshold this is measured from.
  double totalSepE() const { return totalSepE_; };
  /// Channel radius, fm.
  double radius() const { return radius_; };
  /// GSL adaptor for the shift function, for differentiation.
  static double thisShftFunc(double, void *);
  typedef struct Params {
    int lValue;
    double bindingEnergy;
    std::unique_ptr<WhitFunc> whitFunc;
  } Params;
  Params params_;
  double totalSepE_;
  double radius_;
};

#endif
