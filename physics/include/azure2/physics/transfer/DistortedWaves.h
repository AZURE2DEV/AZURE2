#ifndef AZURE2_PHYSICS_TRANSFER_DISTORTEDWAVES_H
#define AZURE2_PHYSICS_TRANSFER_DISTORTEDWAVES_H

#include <complex>
#include <functional>
#include <vector>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/transfer/OpticalPotential.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * Distorted waves of a projectile of spin s (0, 1/2 or 1) on a target whose spin is ignored, in a complex optical
 * potential (docs/upgrade/transfer.md section 4).  The radial function u_lj(r) of orbital momentum l and total
 * angular momentum j of the projectile (l + s = j) solves, at the c.m. energy E > 0 (MeV),
 *
 *   u'' + [ k^2 - (2 mu/hbar^2) U_lj(r) - l(l+1)/r^2 ] u = 0,   U_lj = U(r) of Eq. (2.1) with <2 l.s> of (l, j, s), (4.1)
 *
 * k = sqrt(2 mu E)/hbar (fm^-1), mu = mu_u x u (non-relativistic, as FRESCO's default), u(0) = 0, and beyond the
 * matching radius r_m (the last node of the grid) it is
 *
 *   u_lj(r) = (i/2) [ H^-_l(eta, kr) - S_lj H^+_l(eta, kr) ] = F_l + T_lj H^+_l,  T = (S - 1)/(2i),            (4.2)
 *
 * H^+- = G +- iF the Coulomb functions of physics/coulomb (RadialNumerics.h's CoulombForMatching), eta = Z1 Z2 e^2 mu/(hbar^2 k).  S_lj is the
 * nuclear S-matrix element (the Coulomb phase e^{2 i sigma_l} is not in it), as FRESCO prints it.  The distorted
 * wave of Eq. (6.3) of the note is chi^(+)(k, r) = (4 pi/(k r)) sum_lm i^l e^{i sigma_l} u_l(r) Y*_lm(k^) Y_lm(r^)
 * for s = 0.  Numerov's method (RadialNumerics.h) outward from the origin; S from the values at r_m and at the
 * node nearest r_m - 1 fm (two-point matching, no derivative).  A projectile with spin but no spin-orbit term has
 * waves independent of j: one per l is kept (two_j = 2l), and the observables add the 2s + 1 projections.
 */

/// What to solve: the pair, the projectile spin, the energy and the grid.
struct ScatteringRequest {
  ChannelKinematics pair;   ///< z1, z2, mu_u (u); non-relativistic only
  int two_s = 0;            ///< 2s of the projectile: 0, 1 or 2
  double energy_MeV = 0.0;  ///< E (c.m.) > 0
  double step_fm = 0.02;    ///< h
  double matching_radius_fm = 20.0;  ///< r_m, the last node (rounded to a multiple of h)
  int l_max = 20;           ///< highest l solved
  bool keep_waves = true;   ///< false: S only (u is not stored)
};

/// One partial wave.
struct PartialWave {
  int l = 0;
  int two_j = 0;
  std::complex<double> S;            ///< S_lj
  std::vector<std::complex<double>> u;  ///< u_lj(r_i), normalised as Eq. (4.2) (empty unless kept)
};

/// The waves of one channel at one energy.
struct DistortedWaves {
  double k_fm = 0.0;     ///< k (fm^-1)
  double eta = 0.0;      ///< Sommerfeld parameter
  double step_fm = 0.0;  ///< h of the grid r_i = i h
  int two_s = 0;
  bool spin_orbit = false;           ///< whether the potential has a spin-orbit term (waves per j)
  std::vector<double> sigma;         ///< Coulomb phases sigma_l, l = 0 ... l_max (rad)
  std::vector<PartialWave> waves;    ///< by l, then j ascending
};

/// Solves Eq. (4.1) for l = 0 ... l_max (and every j when the potential has a spin-orbit term and s > 0).
/// Throws std::invalid_argument for E <= 0, a bad spin or grid, std::runtime_error when a Coulomb function fails.
DistortedWaves SolveDistortedWaves(const OpticalPotential &potential, const ScatteringRequest &request);

/*!
 * The same for any potential: U(r) = central(r) + spin_orbit_form(r) <2 l.s> (MeV, r > 0 in fm), with
 * `has_spin_orbit` saying whether the second term exists; `point_charge_MeV_fm` is C of a point-Coulomb term C/r in
 * central (0 when central is finite at the origin), which the series start of the regular solution needs.
 * What SolveDistortedWaves calls.
 */
DistortedWaves SolveDistortedWavesInPotential(const std::function<std::complex<double>(double r_fm)> &central_MeV,
                                              const std::function<std::complex<double>(double r_fm)> &spin_orbit_form_MeV,
                                              bool has_spin_orbit, double point_charge_MeV_fm,
                                              const ScatteringRequest &request);

/*!
 * One partial wave at a time (what SolveDistortedWaves is made of, and what a sum over l that stops at an l it does
 * not know in advance needs).  RadialChannel is the potential of one channel at one energy on r_i = i h, tabulated
 * once; SolvePartialWave integrates Eq. (4.1) for one (l, j) outward to the outer matching node and matches at that
 * node and the inner one (S of Eq. (4.2) from the two values, as SolveDistortedWaves), keeping the first
 * `stored_nodes` values.  A plane-wave channel (no potential, no charge) gives u_l = F_l(0, kr) = kr j_l(kr) exactly,
 * S = 1 (Riccati-Bessel functions by Miller's downward recurrence normalised with sum (2l+1) j_l^2 = 1, upward where
 * l < kr).
 */
struct RadialChannel {
  double k_fm = 0.0, eta = 0.0, step_fm = 0.0;
  int two_s = 0;
  bool spin_orbit = false;   ///< waves per j
  bool plane_wave = false;   ///< no potential and no Coulomb: the waves are Riccati-Bessel functions
  double potential_range_fm = 0.0;   ///< beyond it the potential is point Coulomb (nuclear terms below e^-12)
  double k2mu_per_MeV_fm2 = 0.0;     ///< 2 mu / hbar^2
  double point_charge_MeV_fm = 0.0;  ///< C of a point-Coulomb term C/r in `central` (0 when finite at the origin)
  std::complex<double> central_at_origin_MeV;     ///< the finite part of `central` at r -> 0 (C/r removed)
  std::complex<double> spin_orbit_at_origin_MeV;  ///< the spin-orbit form at the first node (series start)
  std::vector<std::complex<double>> central_MeV;     ///< central U(r_i) (MeV), index 0 unused
  std::vector<std::complex<double>> spin_orbit_MeV;  ///< the spin-orbit form at r_i, empty without a spin-orbit term
};

/// The potential's channel on r_i = i h up to the node nearest radius_fm (the outer matching node of every wave
/// solved on it must lie inside).  Throws std::invalid_argument as SolveDistortedWaves.
RadialChannel TabulateChannel(const OpticalPotential &potential, const ScatteringRequest &request, double radius_fm);
/// The plane-wave channel of the pair's mass at the request's energy (the pair's charges are ignored).
RadialChannel PlaneWaveChannel(const ScatteringRequest &request, double radius_fm);

/// Where a wave is matched: two nodes, outer > inner.
struct WaveMatching {
  int outer_node = 0;
  int inner_node = 0;
};
/// SolveDistortedWaves' matching: the last node of the channel and the node nearest 1 fm inside (at least 5 nodes).
WaveMatching MatchingAtLastNode(const RadialChannel &channel);
/*!
 * The matching of a wave whose sum over l is cut by a bound (the zero-range THM weight, ZeroRange.h): beyond the
 * stored range, the potential and the turning point (eta + sqrt(eta^2 + l(l+1)))/k by max(n_q quarter wavelengths,
 * min_fm), the inner node a quarter wavelength inside (at least 2 nodes) and beyond the stored nodes.  The channel must
 * reach the outer node (MatchingRadius gives the radius to tabulate for the highest l).
 */
struct MatchingBeyondTurningPoint {
  double quarter_wavelengths = 2.0;
  double min_fm = 1.0;
};
WaveMatching MatchingFor(const RadialChannel &channel, int l, int stored_nodes, const MatchingBeyondTurningPoint &rule);
/// The radius (fm) a channel must reach for the waves l <= l_max under the rule.
double MatchingRadius(const RadialChannel &channel, int l_max, int stored_nodes, const MatchingBeyondTurningPoint &rule);

/// One partial wave (l, j) normalised as Eq. (4.2), u on r_i, i < stored_nodes.  Throws std::runtime_error when the
/// wave cannot be normalised (S or the scale not finite), std::invalid_argument for a matching outside the channel.
PartialWave SolvePartialWave(const RadialChannel &channel, int l, int two_j, const WaveMatching &matching,
                             int stored_nodes);
/// The Coulomb phases sigma_l, l = 0 ... l_max (rad) of the channel (0 for a plane wave).
std::vector<double> CoulombPhases(const RadialChannel &channel, int l_max);

/// Plane waves on the request's grid: u_l = F_l(0, kr) for l = 0 ... l_max, S = 1, sigma = 0.
DistortedWaves PlaneWaves(const ScatteringRequest &request);

/// A spinless channel to solve: a potential, or a plane wave (no potential, no charge), the pair and the c.m.
/// energy (MeV).
struct WaveChannel {
  bool plane_wave = false;
  OpticalPotential potential;  ///< ignored for a plane wave
  ChannelKinematics pair;
  double energy_MeV = 0.0;
};
/// The channel on r_i = i h (h = step_fm), tabulated as far as the waves l <= l_max stored on `stored_nodes` nodes
/// need under the rule.
RadialChannel ChannelFor(const WaveChannel &channel, double step_fm, int stored_nodes, int l_max,
                         const MatchingBeyondTurningPoint &rule);
/// The waves l = 0 ... l_max of the channel on `stored_nodes` nodes of step h, each matched by the rule, with the
/// Coulomb phases: what the THM distortions read.
DistortedWaves SolveStoredWaves(const WaveChannel &channel, double step_fm, int stored_nodes, int l_max,
                                const MatchingBeyondTurningPoint &rule);

/// The wave of (l, j) (j ignored without spin-orbit); throws std::out_of_range when it was not solved.
const PartialWave &WaveOf(const DistortedWaves &w, int l, int two_j);

/*!
 * The elastic cross section over the Rutherford cross section at the c.m. angle theta (deg, 0 < theta <= 180),
 * averaged over the 2s + 1 projections of the projectile spin (target spin ignored):
 *
 *   f_m'm(theta) = f_C(theta) delta_m'm + sum_lj sqrt(4 pi (2l+1))/(2ik) e^{2i sigma_l} (S_lj - 1)
 *                  <l 0 s m | j m> <l m-m' s m' | j m> Y_l,m-m'(theta, 0),                                  (4.3)
 *   f_C = -eta/(2k sin^2(theta/2)) exp(-i eta ln sin^2(theta/2) + 2 i sigma_0),
 *   sigma/sigma_R = sum_mm' |f_m'm|^2 / ((2s + 1) |f_C|^2)     (charged pair).
 *
 * For a neutral pair (eta = 0) it returns dsigma/dOmega in mb/sr instead (there is no Rutherford cross section).
 */
double ElasticOverRutherford(const DistortedWaves &w, double theta_deg);
/// dsigma/dOmega (mb/sr) of the elastic scattering at theta (deg), Eq. (4.3).
double ElasticCrossSection(const DistortedWaves &w, double theta_deg);
/// The reaction cross section (mb): (pi/k^2) sum_lj (2j+1)/(2s+1) (1 - |S_lj|^2) x 10 mb/fm^2.
double ReactionCrossSection(const DistortedWaves &w);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
