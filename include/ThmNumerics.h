#ifndef THM_NUMERICS_H
#define THM_NUMERICS_H

#include "Constants.h"

#include <string>
#include <vector>

/*!
 * Small numerical pieces shared by the THM distortion factor R(E)
 * (ThmDistortion) and the DW vertex (ThmDwVertex).  Each one is the single
 * copy of what both units computed in the same way.
 */

/// e^{i sigma_l}, sigma_l = arg Gamma(l + 1 + i eta), l = 0..lmax (eta = 0: 1).
std::vector<complex> ThmCoulombPhases(double eta, int lmax);

/// A number as the THM messages print it (6 significant digits).
std::string ThmNumberText(double x);

/*!
 * The s-x bound state phi(r) without its normalization: exp(-kappa r)/r for a
 * Yukawa tail (or eta_B = 0), else the Whittaker tail
 * 2 kappa e^{-kappa r} U(1 + eta_B, 2, 2 kappa r) = W(2 kappa r)/r up to a
 * constant; 0 if GSL fails (other than by underflow).  Call it inside a
 * GslQuiet scope (GSLException.h).
 */
double ThmBoundStateTail(double r, double kappa, double etaB, bool yukawa);

/*!
 * The cubic Lagrange weights through the nodes s = -1, 0, 1, 2 at s:
 * f(s) = c[0] f(-1) + c[1] f(0) + c[2] f(1) + c[3] f(2).
 */
inline void ThmCubicLagrange(double s, double c[4]) {
  c[0] = -s * (s - 1.0) * (s - 2.0) / 6.0;
  c[1] = (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0;
  c[2] = -(s + 1.0) * s * (s - 2.0) / 2.0;
  c[3] = (s + 1.0) * s * (s - 1.0) / 6.0;
}

#endif
