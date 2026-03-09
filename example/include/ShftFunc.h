#ifndef SHIFTFUNC_H
#define SHIFTFUNC_H

#include "PPair.h"
#include "WhitFunc.h"

class ShftFunc;

///A function class for negative energy shift functions.

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
  ShftFunc(PPair* pPair) {
    totalSepE_=pPair->GetSepE()+pPair->GetExE();
    radius_=(double) pPair->GetChRad();
    params_.whitFunc= new WhitFunc(pPair);
  };
  /*!
   *
   */
  ~ShftFunc() {
    delete params_.whitFunc;
  };
  /*!
   * The parenthesis operator is defined to make the class instance callable as a function.  The orbital
   * angular momentum and energy in the compound system are the dependent variables.
   * The function returns the value of the shift function.
   */
  double operator()(int l,double energy);
  /*!
   * Returns the energy derivative of the shift function at the specified orbital 
   * angular momentum and energy in the compound system.
   */
  double EnergyDerivative(int l,double energy);
 private:
  double totalSepE() const {return totalSepE_;};
  double radius() const {return radius_;};
  static double thisShftFunc(double,void*);
  static double theWhitFunc(double,void*);
  typedef struct Params {
    int lValue;
    double bindingEnergy;
    WhitFunc *whitFunc;    
  } Params;
  Params params_;
  double totalSepE_;
  double radius_;
};

#endif
