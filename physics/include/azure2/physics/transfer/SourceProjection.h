#ifndef AZURE2_PHYSICS_TRANSFER_SOURCEPROJECTION_H
#define AZURE2_PHYSICS_TRANSFER_SOURCEPROJECTION_H

#include <complex>
#include <functional>
#include <vector>

#include "azure2/physics/transfer/DistortedWaves.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The projection of the source of the prior-form distorted-wave vertex on the harmonics of the x + A channel
 * (docs/upgrade/transfer.md section 25).  With r = r_xA, u = r_sx, R_s = alpha r + u, R_a = r + beta u, an l = 0
 * overlap phi(u) of s + x and spinless waves chi_i (entrance, a + A) and chi_f (exit, s + F):
 *
 *   S(r) = int d^3u phi(u) chi_f^(-)*(k_f, R_s) chi_i^(+)(k_i, R_a),   S_lm(r) = int dOmega_r Y*_lm(r^) S(r),   (25.1)
 *
 * and, in the frame k_i^ = z^, k_f^ at the angle theta_fi in the xz plane,
 *
 *   S_lm(r) = sum_{Ls La} H^l_{Ls La}(r) <Ls m La 0|l m> sqrt((2La+1)/4pi) Y_{Ls m}(theta_fi, 0),               (25.2)
 *   H = (4pi)^2/(k_f k_i) (-i)^Ls i^La e^{i(sigma^f_Ls + sigma^i_La)} sqrt(4pi/(2l+1)) h,
 *   h^l_{Ls La}(r) = 2pi int u^2 du phi(u) int dcos(theta_u) f_Ls(R_s) f_La(R_a) K^l_{Ls La}(theta_s, theta_a),
 *   K = sum_M <Ls M La -M|l 0> Y_{Ls M}(theta_s, 0) Y_{La -M}(theta_a, 0),   f_L(R) = u_L(kR)/R.
 *
 * h (and dh/dr at fixed u) does not depend on the directions of k_i and k_f: it is computed once per exit energy at
 * the radius a; the directions cost the sum over the triples (l; Ls, La) only.  The surface term of EntranceVertex
 * (ThmVertex.h) for the wave w = r S_lm is (hbar^2/2mu) u(a) [(B - 1) S_lm(a) - a S_lm'(a)] (transfer.md (25.3)).
 *
 * Numerics: u from r_min to r_min + tail_per_kappa/kappa in Gauss-Legendre panels, cos(theta_u) on Gauss-Legendre
 * nodes, the waves on uniform grids read by six-point Lagrange interpolation (value and derivative), partial waves
 * L <= k R_max + margin in each channel.  The sum over the points (u, cos theta_u) is cut into blocks of consecutive
 * points, each summed in order and the blocks added in order, so that a caller that sums the blocks on several threads
 * gets the numbers of one thread.  Nothing is kept between calls: the plan and every result belong to the caller.
 */

/// The overlap of the Trojan horse a = s + x (l = 0): phi(r) = u(r)/r (fm^-3/2 times the unit of u), 0 below r_min;
/// kappa (fm^-1) sets where the quadrature in u ends.
struct SourceOverlap {
  std::function<double(double r_fm)> phi;
  double r_min_fm = 0.0;
  double kappa_fm = 0.0;
};

/// What is projected: the overlap, the entrance channel, the grids of the wave tables, the exit wave numbers to come,
/// the mass ratios and the radius.
struct SourceSetup {
  SourceOverlap overlap;
  WaveChannel entrance;            ///< i (a + A) at its c.m. energy
  double entrance_step_fm = 0.0;   ///< step of the entrance wave table
  double exit_step_fm = 0.0;       ///< step of the exit wave tables
  double exit_k_max_fm = 0.0;      ///< the largest exit wave number of the energies to come (partial waves, panels)
  double alpha = 0.0;              ///< m_A/m_F: R_s = alpha r + u
  double beta = 0.0;               ///< m_s/m_a: R_a = r + beta u
  double radius_fm = 0.0;          ///< a > 0
  std::vector<int> l_values;       ///< the channel's orbital momenta l, ascending and distinct
};

/// The numerical settings (the values reproduce the engine's vertex before the projection moved here).
struct SourceProjectionRequest {
  double tail_per_kappa = 26.0;    ///< u_max = r_min + tail_per_kappa/kappa
  int nodes_per_panel = 16;        ///< Gauss-Legendre nodes per u panel
  double panel_phase_rad = 4.0;    ///< a panel spans at most this of (k_f,max + beta k_i) u ...
  double panel_max_fm = 5.0;       ///< ... and at most this
  int partial_wave_margin = 14;    ///< L <= k R_max + margin (R_max = a + beta u_max, alpha a + u_max)
  double angle_nodes_per_l = 1.2;  ///< cos(theta_u) nodes: angle_nodes_per_l max(L_i,max, L_f,max) + extra
  int angle_nodes_extra = 20;
  int l_cap = 200;                 ///< refused beyond (kappa too small)
  MatchingBeyondTurningPoint matching;  ///< where each wave is matched
  int energies_per_pass = 8;       ///< exit energies summed together (the entrance side once per point and pass)
  int sum_blocks = 64;             ///< blocks of the sum over the points
  double sum_megabytes = 32.0;     ///< at most this for their buffers (fewer blocks if needed)
};

/// One triple (l; Ls, La) with its coupling coefficients.
struct SourceTriple {
  int l_index = 0;                 ///< index of l in l_values
  int l = 0, l_exit = 0, l_entrance = 0;  ///< l, Ls, La
  /// sqrt((2Ls+1)/(2l+1)) (-1)^nu <Ls 0 La nu|l nu> (x 2 for nu > 0), nu <= min(l, La): the kernel K's coefficients
  std::vector<double> kernel;
  std::vector<double> direction;   ///< <Ls m La 0|l m> sqrt((2La+1)/4pi), m <= min(l, Ls)
};

/// Everything that does not depend on the exit energy: built once, owned by the caller.
struct SourceProjectionPlan {
  SourceSetup setup;
  SourceProjectionRequest request;
  double u_max_fm = 0.0;
  int l_entrance_max = 0, l_exit_max = 0;
  std::vector<double> u_nodes, u_weights, phi_at_u;  ///< the u quadrature and phi at its nodes
  std::vector<double> cos_nodes, cos_weights;        ///< the cos(theta_u) quadrature
  std::vector<SourceTriple> triples;
  double entrance_k_fm = 0.0;
  int entrance_nodes = 0, exit_nodes = 0;            ///< nodes of the wave tables
  std::vector<std::vector<std::complex<double>>> entrance_waves;  ///< u_L(k_i R), L <= l_entrance_max
  std::vector<std::complex<double>> entrance_phases;               ///< e^{i sigma_L}
};

/*!
 * The plan: the quadratures, phi at the u nodes, the triples and the entrance waves.  Throws std::invalid_argument for
 * a bad setup, std::length_error ("the partial-wave sums would need L up to ... (kappa too small)") beyond l_cap, and
 * passes on the exception of an entrance wave that cannot be solved (std::runtime_error).
 */
SourceProjectionPlan PrepareSourceProjection(const SourceSetup &setup, const SourceProjectionRequest &request);

/// The exit waves at one energy.
struct SourceExitWaves {
  double k_fm = 0.0;
  std::vector<std::vector<std::complex<double>>> waves;  ///< u_L(k_f R), L <= l_exit_max
  std::vector<std::complex<double>> phases;              ///< e^{i sigma_L}
};
/// Solves them (the channel's energy is the exit c.m. energy).  Throws as SolveStoredWaves.
SourceExitWaves SolveSourceExit(const SourceProjectionPlan &plan, const WaveChannel &exit);

/// The reduced amplitudes H of (25.2) and dH/dr at r = a, per triple of the plan (the unit of S: fm^3 times that of
/// phi; the derivative per fm).
struct SourceAmplitudes {
  std::vector<std::complex<double>> value, derivative;
};

/// The number of blocks of the sum over the points for `energies` exit energies summed together, and the size of one
/// block's buffer (complex numbers: triples x energies x 2).
int SourceSumBlocks(const SourceProjectionPlan &plan, int energies);
std::size_t SourceBlockSize(const SourceProjectionPlan &plan, int energies);
/// Adds block `block` of `blocks` (the points of that block in order) for the exit energies `exits` to `sums`
/// (SourceBlockSize values, zeroed by the caller).
void SumSourceBlock(const SourceProjectionPlan &plan, const std::vector<SourceExitWaves> &exits, int block, int blocks,
                    std::complex<double> *sums);
/// The blocks' buffers (blocks x SourceBlockSize, in block order) added in block order, with the factors of (25.2).
std::vector<SourceAmplitudes> CombineSourceBlocks(const SourceProjectionPlan &plan,
                                                  const std::vector<SourceExitWaves> &exits,
                                                  const std::vector<std::complex<double>> &block_sums, int blocks);

/// All of it on one thread, in passes of request.energies_per_pass exit energies (the numbers of a threaded caller
/// that runs the same passes and blocks).
std::vector<SourceAmplitudes> ProjectSource(const SourceProjectionPlan &plan, const std::vector<WaveChannel> &exits);

/// S_lm(a) and S_lm'(a) (fm^-1 the derivative) for m = 0 ... l at the direction cos(theta_fi) of k_f to k_i, per l
/// index: value[li][m], derivative[li][m] (S_l,-m = (-1)^m S_lm).
struct SourceMultipoles {
  std::vector<std::vector<std::complex<double>>> value, derivative;
};
SourceMultipoles SourceMultipolesAt(const SourceProjectionPlan &plan, const SourceAmplitudes &amplitudes,
                                    double cos_theta_fi);

/// phi~(q) = 4 pi int u^2 phi(u) j_0(qu) du on the plan's u nodes (q in fm^-1): the plane-wave source,
/// S = phi~(q) e^{i p.r} without distortion (p = k_i - alpha k_f, q = k_f - beta k_i).
double SourceMomentumAmplitude(const SourceProjectionPlan &plan, double q_fm);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
