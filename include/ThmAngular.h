#ifndef THM_ANGULAR_H
#define THM_ANGULAR_H

#include <string>
#include <vector>

#include "Constants.h"

/*!
 * Fixed-angle (differential) HOES observable of a THM experiment
 * (`<thm> experiment[...] theta=thmin-thmax`; docs/source/theory/
 * thm_implementation.rst, "Fixed-angle observable").
 *
 * With the quantization axis along p_xA the THM entrance carries only
 * m_l = 0, exactly like the plane (or Coulomb) wave of an ordinary unpolarized
 * beam along z.  The spin-summed angular distribution of the exit pair is then
 * the Blatt-Biedenharn form that AZURE2 uses for differential cross sections
 * (CNuc::CalcAngularDists, GenMatrixFunc), with the T-matrix element replaced
 * by the HOES partial amplitude x of each J^pi group, entrance (s, l) and exit
 * (s', l'):
 *
 *   sum_{spins} |F(theta)|^2 = (1/pi) sum_L b_L P_L(cos theta),
 *   b_L = sum_{ij} (-1)^{s'-s}/4 Zb(l_i J_i l_j J_j; s L) Zb(l'_i J_i l'_j J_j; s' L) Re(x_i x_j^*),
 *   Zb(l1 J1 l2 J2; s L) = sqrt((2l1+1)(2l2+1)(2J1+1)(2J2+1)) (l1 0 l2 0|L 0) W(l1 J1 l2 J2; s L),
 *
 * the sum over pairs with the same s and the same s'.  It is the direct sum
 * F_{nu nu'} = sum sqrt(2l+1) <s nu l 0|J nu><s' nu' l' m'|J nu> x Y_l'^m'(theta, 0)
 * in AZURE2's phase convention, and 4 pi times the angular average is
 * 4 b_0 = sum (2J+1) |x|^2, the angle-integrated HOES cross section.
 */

/// Quantum numbers of one HOES partial amplitude.
struct ThmPartialWave {
  double J = 0.0;   ///< J of the J^pi group
  double s = 0.0;   ///< entrance channel spin
  int l = 0;        ///< entrance orbital momentum
  double sp = 0.0;  ///< exit channel spin
  int lp = 0;       ///< exit orbital momentum
  bool operator==(const ThmPartialWave &o) const {
    return J == o.J && s == o.s && l == o.l && sp == o.sp && lp == o.lp;
  }
};

/*!
 * The Legendre coefficients b_L (L = 0 .. size-1) of the spin-summed angular
 * distribution of the amplitudes x (one per wave), as above.  The geometric
 * coefficients depend only on the waves; a per-thread cache keeps them while
 * the list of waves does not change (every point of a model has the same).
 */
void ThmLegendreCoefficients(const std::vector<ThmPartialWave> &waves, const std::vector<complex> &x,
                             std::vector<double> &b);

/*!
 * The angular window of a THM experiment: theta (degrees) is the c.m. angle
 * of particle 1 of the exit pair relative to particle 2, measured from the
 * relative momentum of particle 1 of the entrance pair (the transferred
 * particle x, when it is particle 1) relative to particle 2 -- theta_cm =
 * arccos(k_xA . k_bB) of the Catania analyses (Tribble et al., RPP 77 (2014)
 * 106901, section 4.3).  meanP[L] is the average of P_L(cos theta) over the
 * window, uniform in cos theta (solid angle); for thetaMin == thetaMax it is
 * P_L at that angle.
 */
struct ThmAngleWindow {
  static const int kMaxL = 80;  ///< highest Legendre order the window carries
  std::string experiment;
  double thetaMin = 0.0, thetaMax = 180.0;  ///< degrees
  std::vector<double> meanP;                ///< L = 0 .. kMaxL
  /// <dsigma/dOmega> over the window from the b_L: (1/pi) sum_L b_L meanP[L].
  double Mean(const std::vector<double> &b) const;
};

/// Fills the window for [thetaMin, thetaMax] (degrees, 0 <= min <= max <= 180):
/// Gauss-Legendre in cos theta, 48 nodes, exact for every L <= kMaxL.
void BuildThmAngleWindow(double thetaMin, double thetaMax, ThmAngleWindow &out);

#endif
