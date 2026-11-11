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
