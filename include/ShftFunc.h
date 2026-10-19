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
 * the derivative with respect to \f$ \rho \f$. The AZURE function class ShftFunc
 * uses the GSL package to calculates the numerical derivative.
 */

class ShftFunc {
 public:
  /*!
   * The ShftFunc object is created with reference to a particle pair.
   */
  ShftFunc(PPair *pPair) {
    totalSepE_ = pPair->GetSepE() + pPair->GetExE();
    radius_ = (double)pPair->GetChRad();
    neutral_ = (pPair->GetZ(1) * pPair->GetZ(2) == 0);
    redmass_ = pPair->GetRedMass();
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
   * Returns the second energy derivative of the shift function (a central
   * difference of EnergyDerivative with a 1 keV step).
   */
  double EnergySecondDerivative(int l, double energy);
  /*!
   * Exact shift function of a closed neutral channel (eta = 0) at rho = kappa a:
   * S_0 = -rho, S_l = -rho^2 / (l - S_{l-1}) - l, i.e. rho u'/u for the decaying
   * solution u = rho k_l(rho).  Every term is negative, so the recursion is stable.
   */
  static double NeutralBoundShift(int l, double rho);
  /// Is the closed form in use (off with AZURE_NO_NEUTRAL_CLOSED_FORM, for comparisons)?
  static bool UseNeutralClosedForm();

 private:
  /// Separation plus excitation energy of the pair, the threshold this is measured from.
  double totalSepE() const { return totalSepE_; };
  /// Channel radius, fm.
  double radius() const { return radius_; };
  /// GSL adaptor for the shift function, for differentiation.
  static double thisShftFunc(double, void *);
  static double thisShftFuncDeriv(double, void *);
  /// GSL adaptor for the Whittaker function.
  static double theWhitFunc(double, void *);
  typedef struct Params {
    int lValue;
    double bindingEnergy;
    std::unique_ptr<WhitFunc> whitFunc;
  } Params;
  Params params_;
  double totalSepE_;
  double radius_;
  bool neutral_;
  double redmass_;
};

#endif
