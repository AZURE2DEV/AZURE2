#ifndef ANGCOEFF_H
#define ANGCOEFF_H

#include "Constants.h"

/// A container class for angular coupling coefficient functions.

/*!
 * The AngCoeff class serves as a container class for the angular momentum coupling
 * coefficients.
 */

class AngCoeff {
 public:
  /*!
   * Returns the Clebsh-Gordan coefficient for the given angular momentum quantum numbers.
   */
  static double ClebGord(double, double, double, double, double, double);
  /*!
   * Returns the Racah coefficient for the given angular momentum quantum numbers.
   */
  static double Racah(double, double, double, double, double, double);
  /*!
   * Spherical harmonic \f$Y_l^m(\theta,\phi)\f$ with the Condon-Shortley
   * phase. AZURE2 otherwise only needs \f$m=0\f$, for which points carry
   * Legendre polynomials; polarization observables need \f$m \neq 0\f$.
   */
  static complex SphericalHarmonic(int, int, double, double phi = 0.0);
  /*!
   * The associated Legendre function \f$P_l^1(x)\f$ in the convention of
   * Seyler and Weller's Eq. (12) -- i.e. *without* the Condon-Shortley phase,
   * so that \f$P_l^1(\cos\theta) = \sin\theta\, dP_l/d(\cos\theta) \ge 0\f$
   * near \f$\theta = 0\f$. GSL's \c gsl_sf_legendre_Plm carries the phase, so
   * this is its negative.
   */
  static double LegendreP1(int, double);
};

#endif
