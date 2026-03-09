#ifndef ANGCOEFF_H
#define ANGCOEFF_H

///A container class for angular coupling coefficient functions.

/*!
 * The AngCoeff class serves as a container class for the angular momentum coupling 
 * coefficients.  
 */

class AngCoeff {
 public:
/*!
 * Returns the Clebsh-Gordan coefficient for the given angular momentum quantum numbers.
 */
  static double ClebGord(double,double,double,double,double,double);
/*!
 * Returns the Racah coefficient for the given angular momentum quantum numbers.
 */
  static double Racah(double,double,double,double,double,double);
};

#endif
