#ifndef POLARIZATIONFUNC_H
#define POLARIZATIONFUNC_H

#include <utility>
#include <vector>

#include "Constants.h"

class CompoundSystem;
class Site;
namespace azure2 {
namespace physics {
struct SpinObservable;
}
}  // namespace azure2

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
 * Weller, PRC 20 (1979) 453, built at set-up into the pathway table and
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
/*!
 * Y_l^mu(theta, 0) for l = 0 .. lMax and mu = -muMax .. muMax, at (l, mu) index l (2 muMax + 1) + mu + muMax: the
 * angular factors of the amplitude matrix, one table per angle (a reader's row, RowSet.h; the matrix's own).
 */
void HarmonicTable(double theta, int lMax, int muMax, vector_c &table);
/// The largest |v - v'| of the amplitudes of a reaction between pairs of these spins.
int MaxProjection(double j1, double j2, double jp1, double jp2);

class AmplitudeMatrix {
 public:
  /// The matrix of the reaction from pair `aa` to pair `ir` at the c.m. angle `angle` (deg) of a reader.
  AmplitudeMatrix(const CompoundSystem *compound, double angle, int aa, int ir);
  /// The matrix of the reaction from pair `aa` to pair `ir` without an angle: for AddPathwayCoefficients alone.
  AmplitudeMatrix(const CompoundSystem *compound, int aa, int ir);
  /// An empty matrix of a reaction between pairs of these spins (j1 the projectile's), filled with SetAmplitude.
  AmplitudeMatrix(double j1, double j2, double jp1, double jp2);

  //! Add the contribution of one pathway: T element for (jNum, ch, chp).
  void AddPathway(int jNum, int chNum, int chpNum, complex tMatrixElement);

  /*!
   * The same loop as AddPathway without the spherical harmonic: adds i sqrt(2l+1) (s l v 0|J v) (s' l' v' mu|J v) T
   * x identical to coefficients[a nLp + l'] (a the amplitude's index, l' the exit's orbital momentum, l' < nLp), so
   * that M_a = sum_l' coefficients[a nLp + l'] Y_l'^mu (a site's record of an analysing power, SiteObservable.h).
   */
  void AddPathwayCoefficients(int jNum, int chNum, int chpNum, complex tMatrixElement, complex *coefficients,
                              int nLp) const;
  //! The amplitudes in their order (entrance spin and projection, then exit spin and projection).
  const Amplitude &AmplitudeAt(std::size_t a) const { return amplitudes_[a]; }
  void SetAmplitude(std::size_t a, complex value) { amplitudes_[a].value = value; }
  //! The angular factor of amplitude coefficients: Y_l'^mu from the harmonic table `table` of half-width muMax.
  static complex HarmonicOf(const complex *table, int muMax, int lp, int mu) {
    return table[(std::size_t)lp * (2 * muMax + 1) + (mu + muMax)];
  }

  /*!
   * Coulomb amplitude term, added once per (s, v) diagonal element. The
   * amplitude is the point's own, taken per channel spin: for two identical
   * particles with spin the exchange term enters as (-1)^s f_C(pi - theta)
   * (see Unit::GetCoulombAmplitude(double)).
   */
  /// (A reader's amplitudes at its slot: the point's own, f_d, f_x of Site::CoulombAt.)
  /// With a reader's Coulomb amplitudes f_d, f_x (Site::CoulombAt): f_C(s) = f_d + (-1)^s f_x.
  void AddCoulomb(complex direct, complex exchange);

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

  /*!
   * A general spin observable (azure2/physics/SpinObservable.h):  Re(c Tr(F A F^dagger B)) / Tr(F F^dagger),  with F
   * the amplitude in the product basis of the particles' spin projections (both channel spins decomposed with the
   * Lane-Thomas coefficients A_y uses), A the entrance and B the exit operator of \p obs.  \c SpinObservableBar is
   * its reverse mode, 2 dO/dM* per amplitude slot, for \c PathwayAdjoint.
   */
  double SpinObservableValue(const azure2::physics::SpinObservable &obs) const;
  std::vector<complex> SpinObservableBar(const azure2::physics::SpinObservable &obs) const;
  complex PathwayAdjoint(int jNum, int chNum, int chpNum,
                         const std::vector<complex> &bar) const;

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

  //! Exchange symmetrization of one nuclear pathway: 1 + (-1)^(l'+s') for
  //! elastic scattering of identical particles, 1 otherwise.
  double IdenticalFactor(int lp, double sp) const;

  const CompoundSystem *compound_ = nullptr;
  //! What At() hands out for a combination the matrix does not hold: written and discarded (zeroed at every such call).
  complex discarded_ = complex(0.0, 0.0);
  int aa_ = 0, ir_ = 0;
  double theta_ = 0.0;
  double j1_ = 0.0, j2_ = 0.0;  //!< the entrance pair's spins (projectile, target)
  double jp1_ = 0.0, jp2_ = 0.0;  //!< the exit pair's spins (light, heavy)
  //! F in the product basis (exit index, entrance index; row-major) and, per amplitude slot, the (element of F,
  //! coefficient) pairs it maps to.
  std::vector<complex> ProductBasis(std::vector<std::vector<std::pair<int, double>>> *slots) const;
  int muMax_ = 0;
  vector_c harmonics_;  //!< HarmonicTable at theta_
  std::vector<double> entranceSpins_, exitSpins_;
  std::vector<Amplitude> amplitudes_;
};

}  // namespace Polarization

#endif  // POLARIZATIONFUNC_H
