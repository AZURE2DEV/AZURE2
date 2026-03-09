#ifndef WHITFUNC_H
#define WHITFUNC_H

#include "PPair.h"
#include <gsl/gsl_sf_hyperg.h>
#include "Constants.h"

extern double gsl_whit_function(int,double,double,double,int,int);

/// A function class to calculate Whittaker functions for negative energy channels.

/*!
 * The function class WhitFunc uses the GSL package to calculate Whittaker functions for negative energy channels from the GSL confluent hypergeometric functions.  
 */

class WhitFunc {
  public:
  /*!
   * The WhitFunc object is created with reference to a PPair object.
   */
  WhitFunc(PPair *pPair) {
    z1_=pPair->GetZ(1);
    z2_=pPair->GetZ(2);
    redmass_=(double)pPair->GetRedMass();
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
    return  z2_;
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
    const double k=-sqrt(uconv/2.)*fstruc*z1()*z2()*sqrt(redmass()/energy);
    const double m=l+0.5;
    const double z=2.0*sqrt(2.0*uconv)/hbarc*radius*sqrt(redmass()*energy);

    const double a=m-k+0.5;
    const double b=1.0+2.0*m;
    
    return exp(-z/2.0)*pow(z,m+0.50)*gsl_sf_hyperg_U(a,b,z);
  };
  private:
    int z1_;
    int z2_;
    double redmass_;  
};

#endif
