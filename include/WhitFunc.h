#ifndef WHITFUNC_H
#define WHITFUNC_H

#include "PPair.h"
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "Constants.h"

#include <iostream>

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
    
    // FIXME: sometimes a and z are NaN so we make a check before
    if( a != a || b != b || z != z ) return 0;

    return exp(-z/2.0)*pow(z,m+0.50)*gsl_sf_hyperg_U(a,b,z);
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
   * well conditioned, so the exponent has to be carried rather than discarded.
   * gsl_sf_hyperg_U_e10_e returns U with its power of ten kept separate.
   */
  double Scaled(int l, double radius, double energy, double& log10Scale) const {
    log10Scale = 0.0;
    const double k=-sqrt(uconv/2.)*fstruc*z1()*z2()*sqrt(redmass()/energy);
    const double m=l+0.5;
    const double z=2.0*sqrt(2.0*uconv)/hbarc*radius*sqrt(redmass()*energy);
    const double a=m-k+0.5;
    const double b=1.0+2.0*m;
    if( a != a || b != b || z != z || !(z>0.0) ) return 0.0;

    gsl_sf_result_e10 u;
    gsl_error_handler_t* oldHandler = gsl_set_error_handler_off();
    int status = gsl_sf_hyperg_U_e10_e(a,b,z,&u);
    gsl_set_error_handler(oldHandler);
    if(status != GSL_SUCCESS || u.val == 0.0) return 0.0;

    const double log10W = (-z/2.0)/M_LN10 + (m+0.5)*log10(z)
                        + log10(fabs(u.val)) + (double)u.e10;
    double intPart;
    double fracPart = modf(log10W, &intPart);
    if(fracPart < 0.0) { fracPart += 1.0; intPart -= 1.0; }  // mantissa in [1,10)
    log10Scale = intPart;
    return ((u.val < 0.0) ? -1.0 : 1.0)*pow(10.0, fracPart);
  };
  private:
    int z1_;
    int z2_;
    double redmass_;  
};

#endif
