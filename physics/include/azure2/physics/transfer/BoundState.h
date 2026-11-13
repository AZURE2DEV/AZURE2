#ifndef AZURE2_PHYSICS_TRANSFER_BOUNDSTATE_H
#define AZURE2_PHYSICS_TRANSFER_BOUNDSTATE_H

#include <functional>
#include <vector>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/transfer/OpticalPotential.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The bound state of two clusters in a real potential well (docs/upgrade/transfer.md section 3).  The radial
 * function u(r) (fm^-1/2) of orbital momentum l and total angular momentum j of the bound particle (spin s) solves
 *
 *   u'' + [ (2 mu/hbar^2) (-B - V(r)) - l(l+1)/r^2 ] u = 0,   V(r) = -V0 f(r) + V_fixed(r),                  (3.1)
 *
 * at the binding energy B > 0 (MeV, c.m. of the pair), mu the reduced mass of the pair (non-relativistic,
 * mu = mu_u x u of physics/Constants.h), with n interior nodes (0 < r < infinity; FRESCO's NN = n + 1).  The depth
 * V0 (MeV) of the form f is searched; V_fixed holds the terms that are not scaled: the Coulomb potential and the
 * spin-orbit term times <2 l.s> (OpticalPotential.h).  u is normalised, int_0^r_max u^2 dr = 1, and positive next
 * to the origin; beyond the nuclear range u(r) = C W_{-eta, l+1/2}(2 kappa r) with kappa = sqrt(2 mu B)/hbar (fm^-1),
 * eta = Z1 Z2 e^2 mu/(hbar^2 kappa) and the Whittaker function normalised to e^{-kappa r} (2 kappa r)^{-eta} at large
 * r (physics/coulomb/Whittaker.h): C (fm^-1/2) is the single-particle asymptotic normalisation coefficient (ANC),
 * of sign (-1)^n.  This is FRESCO's convention ("ANC" of its bound-state listing).
 *
 * Numerics: Numerov (RadialNumerics.h) outward from the origin and inward from r_max, where u = W; the depth is
 * bracketed by the Sturm count (nodes of the outward solution) and found where the discrete Wronskian (Casoratian)
 * of the two solutions at the matching radius vanishes, to 1e-13 relative.  The error is that of Numerov's method,
 * O(h^4) for a smooth potential, O(h^2) for one with a step (on a node, where the mean of the two sides is used).
 */

/// What to solve for.
struct BoundStateRequest {
  ChannelKinematics pair;    ///< z1, z2 and mu_u (u); non-relativistic only
  int l = 0;                 ///< orbital angular momentum
  int two_j = 0;             ///< 2j, j = l + s coupled (2l when two_s = 0)
  int two_s = 0;             ///< 2s, s the spin the spin-orbit term couples to l
  double binding_MeV = 0.0;  ///< B > 0
  int nodes = 0;             ///< interior nodes
  double step_fm = 0.02;     ///< h
  double r_max_fm = 30.0;    ///< end of the grid; the inward start
  double matching_radius_fm = 0.0;  ///< where the two solutions meet; 0: the radius of the searched term
};

/// The solution, on r_i = i h, i = 0 ... n - 1, r_{n-1} = r_max.
struct BoundState {
  int l = 0, two_j = 0, two_s = 0;
  double binding_MeV = 0.0;
  double kappa_fm = 0.0;     ///< sqrt(2 mu B)/hbar (fm^-1)
  double eta = 0.0;          ///< Z1 Z2 e^2 mu/(hbar^2 kappa)
  double depth_MeV = 0.0;    ///< V0 found (0 for an analytic tail)
  double step_fm = 0.0;
  std::vector<double> u;     ///< u(r_i) (fm^-1/2)
  std::vector<double> potential_MeV;  ///< V(r_i) of Eq. (3.1) at V0 (the interaction of a post-form transfer)
  std::vector<double> coulomb_MeV;    ///< its Coulomb part V_C(r_i) (SolveBoundState; 0 from SolveBoundStateOfForm)
  double rms_radius_fm = 0.0;         ///< sqrt(int r^2 u^2 dr)
  double anc_fm_half = 0.0;           ///< C (fm^-1/2)
  double r_min_fm = 0.0;              ///< u = 0 below it (a cut analytic tail); radial integrals start there
};

/*!
 * Solves Eq. (3.1) with f the real volume Woods-Saxon form of `well` (its depth is the starting value and is
 * replaced) and V_fixed = V_C + Re(spin-orbit form) <2 l.s> of `well` (FRESCO's ISC = 1: only the central depth is
 * scaled).  The matching radius defaults to the radius of the real volume term.  Throws std::invalid_argument for
 * an imaginary term, B <= 0 or bad grid settings, std::runtime_error when no depth below 1e5 MeV gives the state.
 */
BoundState SolveBoundState(const OpticalPotential &well, const BoundStateRequest &request);

/*!
 * The same for any form f(r) (dimensionless) and V_fixed(r) (MeV): what SolveBoundState calls.  At a node where f
 * jumps the caller returns the mean of its two limits.  request.matching_radius_fm must be positive.
 */
BoundState SolveBoundStateOfForm(const std::function<double(double r_fm)> &form,
                                 const std::function<double(double r_fm)> &fixed_MeV, const BoundStateRequest &request);

/// D0 = int d^3r V(r) phi(r) = sqrt(4 pi) int_0^r_max r V(r) u(r) dr (MeV fm^3/2) of an l = 0 state, with
/// phi = (u/r) Y_00: the zero-range constant (FRESCO's "D0" of its bound-state listing).
double ZeroRangeConstant(const BoundState &state);

/// The tail alone, u(r) = C W_{-eta, l+1/2}(2 kappa r) (Whittaker) or C e^{-kappa r} (Yukawa, eta ignored), for
/// r >= r_min and 0 below, on the grid of `request` (its pair, l, B, step, r_max): the two analytic bound states of
/// the engine's THM distortion code.  Not normalised; depth 0, no potential.
/// r_min is a node when it is within 1e-12 fm of one (the node then belongs to the tail).
enum class TailShape { Whittaker, Yukawa };
BoundState AnalyticTail(TailShape shape, const BoundStateRequest &request, double anc_fm_half, double r_min_fm);
/// The same tail at one radius r (fm), off the grid (u in fm^-1/2; 0 below r_min).
double AnalyticTailAt(TailShape shape, const BoundStateRequest &request, double anc_fm_half, double r_min_fm,
                      double r_fm);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
