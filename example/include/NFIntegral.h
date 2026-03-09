#ifndef NFINTEGRAL_H
#define NFINTEGRAL_H

#include "PPair.h"
#include "WhitFunc.h"
#include <math.h>

///A function class to calculate the channel integrals in the denominator of the \f$ N_f^{1/2} \f$ term.

/*!
 * The NFIntegral class returns the channel integral given by 
 * \f$ \int_a^\infty \left[ \frac{W_c(kr)}{W_c{ka_c}} \right]^2 \f$.
 */

class NFIntegral {
 public:
  /*! 
   * The NFIntegral object is created with reference to a PPair object. A WhitFunc object is also created.  
   */
  NFIntegral(PPair* pPair) {
    params_.whitFunc = new WhitFunc(pPair);
    chanrad_ = pPair->GetChRad();
    total_sep_e_ = pPair->GetSepE()+pPair->GetExE();
  };
  /*!
   * The WhitFunc object is destroyed with the NFIntegral object.
   */
  ~NFIntegral() {
    delete params_.whitFunc;
  };
  /*!
   * The parenthesis operator is defined so the instance can be callable as
   * a function.  The final channel orbital angular momentum and final state 
   * energy in the compound system are passed as dependent variables.
   */
  double operator()(int lFinal,double levelEnergy);
  /*!
   * Returns the channel radius of the particle pair.
   */
  double chanRad() const {return chanrad_;};
  /*! 
   * Returns the total seperation energy of the particle pair.
   */
  double totalSepE() const {return total_sep_e_;};
 private:
  static double Integrand(double,void*);
  typedef struct Params {
    WhitFunc *whitFunc;
    int lfValue;
    double bindingEnergy;
    double whitChRadSquaredValue;
  } Params;
  Params params_;
  double chanrad_;
  double total_sep_e_;
};

#endif
