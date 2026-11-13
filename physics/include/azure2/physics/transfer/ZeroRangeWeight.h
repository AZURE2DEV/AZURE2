#ifndef AZURE2_PHYSICS_TRANSFER_ZERORANGEWEIGHT_H
#define AZURE2_PHYSICS_TRANSFER_ZERORANGEWEIGHT_H

#include <complex>
#include <vector>

#include "azure2/physics/transfer/DistortedWaves.h"
#include "azure2/physics/transfer/Overlap.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The zero-range weight of the engine's THM distortion factor (docs/upgrade/transfer.md sections 13 and 24): the
 * orbital integral of Eq. (6.7) for an l = 0 overlap u(r) (the s-x bound state, phi = u/r without Y_00),
 *
 *   O_0(x) = sqrt(4 pi)/(k_i k_f c_i c_f) S(x),
 *   S(x)   = sum_l (2l+1) e^{i(sigma_l^i + sigma_l^f)} P_l(x) int_{r_min} dr u_l^f(c_f r) u(r) u_l^i(c_i r)/r,     (24.2)
 *
 * at the directions x = cos(k_f, k_i), Simpson's rule on the overlap's grid from its first node (r_min of a cut tail)
 * to its last.  The engine's THM weight is M = sqrt(4 pi) O_0 = 4 pi S/(k_sF beta k_aA) (c_i = beta, c_f = 1).  The
 * sum over l is (24.2) term by term with P_l(x) (no Clebsch-Gordan table, so no limit on l), each exit wave solved
 * when its term is reached (SolvePartialWave on the exit channel tabulated once per energy), stopped by a bound on
 * the terms left:
 *
 *   |sum_{l' > l} ...| <= max(2, max_r |u_l^f|) A_{>l},   A_l = (2l+1) int |u u_l^i| / r,   A_{>l} = sum_{l' > l} A_l',
 *
 * once it is below sum_tolerance times the smallest |S(x)| (|u_l'^f| is bounded by its largest value and grows
 * slowly with l outside the turning point).  The entrance waves (at c_i r on the overlap's grid) are solved once,
 * for l = 0, 1, ... until A_l is below entrance_cut times its largest value and decreasing.  Every wave is matched
 * beyond its turning point (MatchingBeyondTurningPoint): a sum that stops at an l it does not know in advance needs
 * waves that are right for any l.
 */

/// The numerical settings (the values reproduce the engine's distortion factor before the wiring).
struct ZeroRangeWeightRequest {
  double c_i = 1.0, c_f = 1.0;  ///< the scale factors of Eq. (6.7): entrance waves at c_i r, exit waves at c_f r
  MatchingBeyondTurningPoint matching;  ///< where each wave is matched
  double entrance_cut = 1e-17;   ///< entrance waves until A_l < entrance_cut x the largest A_l
  double sum_tolerance = 1e-13;  ///< the sum stops once the bound on the terms left is below this x min |S|
  int extra_partial_waves = 0;   ///< terms added after the stop (a refinement)
  int l_cap = 400;               ///< refused beyond (the sum does not converge)
  double max_megabytes = 512.0;  ///< refused when the entrance waves would need more
};

/// The entrance side of the sum, once per channel: the waves at c_i r_p on the overlap's grid and the bounds.
struct ZeroRangeEntrance {
  double k_fm = 0.0;
  std::vector<double> sigma;                            ///< sigma_l^i (rad)
  std::vector<std::vector<std::complex<double>>> u;     ///< u_l^i(c_i r_p), l = 0 ... L
  std::vector<double> remaining;                        ///< A_{>=l} = sum_{l' >= l} A_l', l = 0 ... L + 1
};

/// Solves the entrance waves (step c_i h of the overlap's grid h).  Throws std::runtime_error when the sum would
/// not converge by l_cap or the waves would need more than max_megabytes (the message says which), or a wave cannot
/// be normalised.
ZeroRangeEntrance PrepareZeroRangeEntrance(const RadialOverlap &overlap, const WaveChannel &entrance,
                                           const ZeroRangeWeightRequest &request);

/// The sum at the directions asked.
struct ZeroRangeWeight {
  std::vector<std::complex<double>> sum;  ///< S(x) of Eq. (24.2) per direction (fm^-1/2 times the unit of u)
  std::vector<std::complex<double>> O;    ///< O_0(x) = sqrt(4 pi) S/(k_i k_f c_i c_f)
  int l_max = -1;                         ///< the last l summed
  double end_integrand = 0.0;             ///< |u u_0^f u_0^i / r| at the last node (the size of the radial cutoff)
};

/// Eq. (24.2) with the exit channel solved at its energy (step c_f h).  Throws std::runtime_error (with the l) when
/// an exit wave cannot be normalised.
ZeroRangeWeight ZeroRangeWeightSum(const RadialOverlap &overlap, const ZeroRangeEntrance &entrance,
                                   const WaveChannel &exit, const std::vector<double> &cos_theta,
                                   const ZeroRangeWeightRequest &request);

/// The plane-wave limit of the weight for an l = 0 overlap: phi~(q) = 4 pi int_{r_min} r u(r) j_0(qr) dr (the
/// Fourier transform of phi = u/r; fm^3/2 for u in fm^-1/2), q in fm^-1, Simpson's rule on the overlap's grid.
double ZeroRangeMomentumAmplitude(const RadialOverlap &overlap, double q_fm);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
