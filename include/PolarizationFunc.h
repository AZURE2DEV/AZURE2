#ifndef POLARIZATIONFUNC_H
#define POLARIZATIONFUNC_H

#include <vector>

#include "Constants.h"

class CNuc;
class EPoint;

/*!
 * Polarization observables from the channel-spin amplitude matrix.
 *
 * The scheme follows R. G. Seyler, Nucl. Phys. A124 (1969) 253, Eq. (4), which
 * writes the transition matrix in the channel-spin representation directly in
 * terms of the Lane-Thomas collision matrix -- the formalism AZURE2 already
 * uses:
 *
 *   M_{s'v'sv}(theta) = sqrt(pi)/k [ -C(theta) delta_ss' delta_vv'
 *       + i sum_{J,l,l'} sqrt(2l+1) (s l v 0|J v) (s' l' v' v-v'|J v)
 *              exp{i(w_l + w_l')} (delta_ss' delta_ll' - U^J_{s'l'sl})
 *              Y_{l'}^{v-v'}(theta, 0) ]
 *
 * The bracketed factor is exactly what AZURE2 already forms as its "tmatrix"
 * element (see AMatrixFunc/RMatrixFunc): with
 * U = uphase (1 + umatrix) and uphase carrying exp{i(w_l+w_l')} together with
 * the hard-sphere phases, the code's tmatrix is
 * delta exp(2 i w_l) - U, i.e. the bracket with its exponential already
 * applied. So no phase has to be reconstructed here.
 *
 * Once M exists every polarization observable is a density-matrix trace,
 * rho_out = M rho_in M^dagger, with no further angular-momentum algebra.
 *
 * Scope: particle channels with a spin-1/2 projectile, which covers the
 * vector analyzing power A_y. Capture channels have no amplitude matrix of this
 * form and take the other route -- the Legendre coefficients of Seyler and
 * Weller, PRC 20 (1979) 453, built by CNuc::CalcCaptureAnalyzingPower and
 * evaluated by GenMatrixFunc::CalculateCaptureAnalyzingPower. Tensor
 * observables are not implemented on either side.
 */
namespace Polarization {

// The spherical harmonics this needs live in AngCoeff, beside the
// Clebsch-Gordan and Racah coefficients.

/*!
 * One element of the channel-spin amplitude matrix, indexed by the entrance
 * and exit channel spins and their projections.
 */
struct Amplitude {
  double s, v;    //!< entrance channel spin and projection
  double sp, vp;  //!< exit channel spin and projection
  complex value;
};

/*!
 * The amplitude matrix for one point, for the reaction connecting particle
 * pairs \p aa (entrance) and \p ir (exit).
 *
 * \p tmatrix supplies AZURE2's T-matrix element for a given (J-group, entrance
 * channel, exit channel); it returns false when that combination carries no
 * pathway, which is how the caller signals a channel pair the R-matrix does not
 * connect.
 */
class AmplitudeMatrix {
 public:
  AmplitudeMatrix(CNuc *compound, EPoint *point, int aa, int ir);

  //! Add the contribution of one pathway: T element for (jNum, ch, chp).
  void AddPathway(int jNum, int chNum, int chpNum, complex tMatrixElement);

  //! Coulomb amplitude term, added once per (s, v) diagonal element.
  void AddCoulomb(complex coulombAmplitude);

  /*!
   * Spin-averaged, spin-summed |M|^2 -- the unpolarized differential cross
   * section. This is the correctness gate: it has to reproduce what
   * GenMatrixFunc::CalculateCrossSection produces by the Blatt-Biedenharn
   * route, and a disagreement means the coupling order or the choice of l
   * versus l' in the spherical harmonic is wrong.
   */
  double UnpolarizedCrossSection() const;

  /*!
   * Vector analyzing power A_y for a spin-1/2 projectile, in the Madison
   * convention with y along k_in x k_out.
   *
   * Returns 0 when the entrance channel spin cannot support it (no spin-1/2
   * projectile), which is the correct value rather than an error.
   */
  double AnalyzingPowerAy() const;

  /*!
   * Reverse mode for A_y. \c AnalyzingPowerBar returns, per amplitude slot,
   * 2 dA_y/dM* -- the factor of two matching the convention AMatrixFunc uses
   * for its cotangents, where the parameter gradient is finally taken as
   * Re(conj(bar) dz/dp) with no further factor.
   *
   * \c PathwayAdjoint then walks exactly the loop \c AddPathway walks and
   * contracts those cotangents with the same coefficients, returning
   * dA_y/dT* for that one (J-group, entrance channel, exit channel). Since
   * M is linear in T, that contraction is the whole derivative.
   */
  std::vector<complex> AnalyzingPowerBar() const;
  complex PathwayAdjoint(int jNum, int chNum, int chpNum,
                         const std::vector<complex> &bar) const;

  /*!
   * Vector polarization P_y of the OUTGOING particle produced with an
   * unpolarized beam and target, in the same Madison convention (y along
   * k_in x k_out).
   *
   * Where A_y is Tr(M sigma_y M+)/Tr(M M+), with the Pauli matrix acting on
   * the entrance spin, this is Tr(sigma_y M M+)/Tr(M M+): the same amplitude
   * matrix, the same trace, taken on the exit index instead. The exit channel
   * spin is decomposed into (ejectile, residual) exactly as AnalyzingPowerAy
   * decomposes the entrance one, so a residual nucleus of any spin is handled.
   *
   * Returns 0 when the ejectile is not spin-1/2, which is the correct value.
   *
   * Note the sign: carrying the Pauli matrix through on the exit index gives
   * -2 Im(...) where the entrance-index trace gives +2 Im(...). For elastic
   * scattering P_y and A_y must come out equal, which is the test.
   */
  double OutgoingPolarizationPy() const;

  /*!
   * Reverse mode for P_y, the exit-index counterpart of AnalyzingPowerBar:
   * per amplitude slot, 2 dP_y/dM*, in the same cotangent convention.
   * Feed the result to PathwayAdjoint, which is observable-agnostic.
   */
  std::vector<complex> OutgoingPolarizationBar() const;

  /*!
   * Numerator of P_y on its own: N = -2 Im(sum u' conj(d')) over the exit-spin
   * decomposition, with no denominator.
   *
   * This is worth having separately because the denominator of P_y is, by
   * unitarity of the Clebsch-Gordan transformation, exactly the spin sum that
   * UnpolarizedCrossSection divides by nEntrance. So the published observable
   *
   *     P dsigma/dOmega  =  P_y * sigma  =  (N/D) * (C D / nEntrance)
   *                      =  C N / nEntrance
   *
   * with C the kinematic factor the caller supplies. The product carries no
   * ratio at all: it is linear in the numerator, so its adjoint needs no
   * quotient rule and cannot blow up where the cross section is small -- which
   * is precisely where A_y's adjoint is worst conditioned.
   */
  double OutgoingPolarizationNumerator() const;
  std::vector<complex> OutgoingPolarizationNumeratorBar() const;

  /*!
   * Unit check of OutgoingPolarizationNumeratorBar: perturbs every amplitude
   * slot in turn and compares the analytic cotangent against a central
   * difference of OutgoingPolarizationNumerator. The convention is
   * bar = 2 dN/dM*, and for a real N of complex M that is exactly
   * dN/dRe(M) + i dN/dIm(M), so both parts are directly checkable.
   * Returns the worst relative disagreement over all slots.
   */
  double SelfCheckNumeratorBar(double h = 1.e-6) const;

  //! Largest |M| with v != v' -- the spin-flip strength, which is what makes
  //! a vector analyzing power non-zero. Diagnostic.
  double MaxSpinFlip() const;
  void DumpSpinHalf() const;

  //! Number of amplitudes accumulated; zero means no pathway contributed.
  std::size_t size() const { return amplitudes_.size(); };

 private:
  complex &At(double s, double v, double sp, double vp);
  int IndexOf(double s, double v, double sp, double vp) const;
  complex Get(double s, double v, double sp, double vp) const;

  CNuc *compound_;
  EPoint *point_;
  int aa_, ir_;
  double theta_;
  std::vector<double> entranceSpins_, exitSpins_;
  std::vector<Amplitude> amplitudes_;
};

}  // namespace Polarization

#endif  // POLARIZATIONFUNC_H
