#ifndef AZURE2_PHYSICS_EXTERNAL_CAPTURE_H
#define AZURE2_PHYSICS_EXTERNAL_CAPTURE_H

/*!
 * external_capture: the formulas of the external part of radiative capture (docs/upgrade/external_capture.md).
 *
 *   bound_state_normalization      N_f of a final bound state, Eq. (1)
 *   inverse_sqrt_normalization     N_f^(-1/2), the factor of the ANC scale s_fc = gamma~_fc N_f^(-1/2) K_c, Eq. (2)
 *   couples                        the multipole selection rule of a pair (entrance channel, final channel)
 *   external_reduced_width         gamma^ext of a level's photon channel, Eq. (4), with the radial integral injected
 *   capture_mask                   the --ignore-externals mask of channel capture in the formal basis, Eq. (5)
 *
 * Plain arrays in physical units; nothing of the engine.  The per-point unit amplitudes of Eq. (3) stay with the data
 * points (Unit), the radial integral with ECIntegral until the Whittaker functions are in physics/ (step 3b (i)).
 */

#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace azure2 {
namespace physics {

/// One particle channel of a final bound state f.
struct BoundChannel {
  double radius_fm = 0.0;        ///< a_c (fm)
  double reduced_mass_u = 0.0;   ///< mu_c (u)
  double gamma = 0.0;            ///< Brune amplitude gamma~_fc (MeV^1/2); Park: gamma~_fc J_f^1/2
  double nf_integral_fm = 0.0;   ///< I_fc = int_a^inf [W(2 kappa r) / W(2 kappa a)]^2 dr (fm); 0 for a channel without
};

/// 2 mu_c a_c / hbar^2 (MeV^-1 fm^-1) of a channel of radius a_c (fm) and reduced mass mu_c (u): the coefficient of
/// gamma_fc^2 I_fc in N_f, so dN_f / d(gamma_fc) = 2 bound_state_coefficient gamma_fc I_fc.
double bound_state_coefficient(double radius_fm, double reduced_mass_u);

/// N_f = n0 + sum_c (2 mu_c a_c / hbar^2) gamma_fc^2 I_fc (dimensionless); n0 = 1 (Brune) or Park's J_f.
double bound_state_normalization(double n0, const std::vector<BoundChannel> &channels);

/// N_f^(-1/2).  Throws std::domain_error for N_f <= 0.
double inverse_sqrt_normalization(double n0, const std::vector<BoundChannel> &channels);

/// Multipole kind of a photon channel.
enum class Multipole : unsigned char { E, M };

/// Whether the entrance channel (l, s) couples to the final channel (l_f, s_f) of a state of spin J_f by the multipole
/// XL, as in external_reduced_width: orbital rule |l - L| <= l_f <= l + L, |s - l_f| <= J_f <= s + l_f, s = s_f; for
/// M also the spin rule |s - L| <= s_f <= s + L, |l - s_f| <= J_f <= l + s_f, l = l_f.  Spins in units of hbar.
bool couples(Multipole X, int L, int l, double s, int l_f, double s_f, double J_f);

/// A particle channel of the level lambda (entrance side of Eq. (4)).
struct CaptureChannel {
  int partition = 0;   ///< alpha (any label shared with the final channels)
  int l = 0;
  double s = 0.0;      ///< channel spin (hbar)
  double gamma = 0.0;  ///< gamma_lambda c (MeV^1/2)
};

/// A particle channel of the final state f.
struct FinalChannel {
  int partition = 0;
  int l = 0;
  double s = 0.0;
  double scale = 0.0;  ///< s_fc = gamma~_fc N_f^(-1/2) K_c, Eq. (2)
};

/// The radial integral I_{c c_f XL}(E_lambda, E_f) of (entrance index, final index), complex (ECIntegral in the
/// engine, channel-capture form).
using RadialIntegral = std::function<std::complex<double>(std::size_t c, std::size_t c_f)>;

/*!
 * gamma^ext_lambda gamma = - sum_c sum_{c_f} s_{f c_f} gamma_lambda c I_{c c_f XL}, over the pairs of one partition that
 * couple (couples, J_f the final state's spin); zero unless XL is E1, E2 or M1.  MeV^1/2, complex.  The integral is
 * called only for the pairs that couple, in the order c, then c_f.
 */
std::complex<double> external_reduced_width(Multipole X, int L, double J_f, const std::vector<CaptureChannel> &level,
                                            const std::vector<FinalChannel> &final_state,
                                            const RadialIntegral &integral);

/*!
 * The projector Pi = T^-1 P T of Eq. (5) (n x n, row-major): T the Brune-from-formal transform of a block
 * (LevelSet::transform), P = diag(1 - masked).  The masked Brune bilinear gamma~^T P A~ P gamma~' is (Pi gamma)^T A
 * (Pi gamma') in the formal basis.  Throws std::domain_error for a singular T.
 */
std::vector<double> capture_mask(const std::vector<double> &T, std::size_t n, const std::vector<unsigned char> &masked);

}  // namespace physics
}  // namespace azure2

#endif
