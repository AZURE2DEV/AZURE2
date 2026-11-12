#ifndef AZURE2_PHYSICS_TRANSFER_CONTINUUM_H
#define AZURE2_PHYSICS_TRANSFER_CONTINUUM_H

#include <functional>
#include <vector>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/transfer/OpticalPotential.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * Scattering states of x + A in a real potential and FRESCO's continuum bins (docs/upgrade/transfer.md section 14).
 *
 * The scattering state of orbital momentum l and total j (spin s of x) at the c.m. energy E > 0 (MeV) solves
 * Eq. (3.1) at -B = E, regular at the origin, with the real potential V = V_C - V f + V_so <2 l.s> (Eq. (2.1), the
 * imaginary terms must be 0).  It is real and normalised to asymptotic amplitude 1 of the sine:
 *
 *   u_E(r) = cos(delta) F_l(eta, kr) + sin(delta) G_l(eta, kr)  ->  sin(kr - l pi/2 - eta ln 2kr + sigma_l + delta)
 *                                                                                            (r beyond V_N), (14.1)
 *
 * delta the nuclear phase shift (S = e^{2 i delta}), k = sqrt(2 mu E)/hbar (fm^-1).  It is the distorted wave of
 * Eq. (4.2) times e^{-i delta}.  The other two normalisations differ by constant factors:
 *
 *   delta(k - k'):  u_k = sqrt(2/pi) u_E                      (int u_k u_k' dr = delta(k - k'), u_k in fm^1/2),
 *   delta(E - E'):  u_E^(E) = sqrt(2 mu/(pi hbar^2 k)) u_E    (int ... = delta(E - E'), in MeV^-1/2 fm^-1/2),  (14.2)
 *
 * since int_0^inf sin(kr + phi) sin(k'r + phi') dr = (pi/2) delta(k - k') + (terms without a delta function) and
 * delta(k - k') = (dE/dk) delta(E - E') = (hbar^2 k/mu) delta(E - E').
 *
 * A bin (FRESCO's &OVERLAP with BE < 0; the FRESCO namelist manual, "CONTINUUM BINS") is the weighted, normalised
 * average of the scattering states over k_min <= k <= k_max (E_min, E_max: FRESCO's E - |ER|/2, E + |ER|/2 with
 * E_min >= 1e-5 MeV):
 *
 *   u_bin(r) = sqrt(2/(pi N)) int dk w(k) u_E(k)(r),   N = int dk |w(k)|^2,                                 (14.3)
 *
 * so that int_0^inf u_bin^2 dr = 1 (by (14.2)) and two bins of the same potential over disjoint intervals, and a bin
 * and a bound state of the same potential, are orthogonal.  u_bin is in fm^-1/2, as a bound state.  The weight
 * (FRESCO's ISC): w = 1 on the real u_E (ISC = 2: FRESCO weights its complex wave by e^{-i delta}, which is u_E),
 * w = sin(delta) (ISC = 4: T* = e^{-i delta} sin(delta)), and either times k (ISC = 12, 14).  The phase shift is
 * continued across the bin from its value in (-pi/2, pi/2] at k_min, as FRESCO does (a jump by more than pi/2
 * between two nodes adds pi), so that u_E is continuous in k.
 *
 * The k integral: FRESCO's rule (NK equally spaced nodes from k_min to k_max, each with the weight
 * (k_max - k_min)/(NK - 1), in both the integral and N) or Simpson's rule on the same nodes (NK odd).  FRESCO's
 * rule is not a consistent quadrature (its end points carry a full weight), but the same rule in (14.3) and in N
 * makes a bin of the same nodes; both converge to (14.3) as NK grows.
 */

/// What to solve: the pair, the quantum numbers, the potential's grid.
struct ContinuumRequest {
  ChannelKinematics pair;    ///< z1, z2 and mu_u (u); non-relativistic only
  int l = 0;                 ///< orbital angular momentum
  int two_j = 0;             ///< 2j (2l when two_s = 0)
  int two_s = 0;             ///< 2s of x (0, 1 or 2)
  double step_fm = 0.05;     ///< h
  double r_max_fm = 40.0;    ///< the last node: the matching radius of every scattering state, the end of a bin
};

/// One scattering state of Eq. (14.1) on r_i = i h, i = 0 ... n - 1.
struct ScatteringState {
  int l = 0, two_j = 0, two_s = 0;
  double energy_MeV = 0.0;  ///< E (c.m.)
  double k_fm = 0.0;        ///< k (fm^-1)
  double eta = 0.0;
  double phase_shift_rad = 0.0;  ///< delta in (-pi/2, pi/2]
  double step_fm = 0.0;
  std::vector<double> u;         ///< u_E(r_i), amplitude 1 (dimensionless)
  std::vector<double> potential_MeV;  ///< V(r_i) with <2 l.s>
  std::vector<double> coulomb_MeV;    ///< its Coulomb part
};

/// The normalisations of Eq. (14.2).
enum class ContinuumNormalisation { AmplitudeOne, DeltaK, DeltaE };

/// The factor that turns u_E of Eq. (14.1) into the normalisation asked (1, fm^1/2 or MeV^-1/2 fm^-1/2).
double NormalisationFactor(const ScatteringState &state, double mu_u, ContinuumNormalisation n);

/// Solves Eq. (14.1) at E > 0 (MeV).  Throws std::invalid_argument for an imaginary term, E <= 0 or a bad grid.
ScatteringState SolveScatteringState(const OpticalPotential &potential, const ContinuumRequest &request,
                                     double energy_MeV);

/// The same for any real potential V(r) (MeV, for the request's l and j; Coulomb included) with Coulomb part V_C(r)
/// and the point charge C (MeV fm) of a C/r term in V (0 when V is finite at the origin): what SolveScatteringState
/// calls.
struct RealPotential {
  std::function<double(double r_fm)> total_MeV;
  std::function<double(double r_fm)> coulomb_MeV;
  double point_charge_MeV_fm = 0.0;
};
ScatteringState SolveScatteringStateIn(const RealPotential &potential, const ContinuumRequest &request,
                                       double energy_MeV);
RealPotential RealPotentialOf(const OpticalPotential &potential, const ContinuumRequest &request);

/// The bin weight: FRESCO's ISC = 2, 4, 12, 14.
enum class BinWeight { Flat, TMatrix, FlatTimesK, TMatrixTimesK };
/// The rule of the k integral of Eq. (14.3).
enum class BinRule { Fresco, Simpson };

/// A bin of Eq. (14.3) on r_i = i h, with what FRESCO lists for it.
struct ContinuumBin {
  int l = 0, two_j = 0, two_s = 0;
  double energy_min_MeV = 0.0, energy_max_MeV = 0.0;
  double step_fm = 0.0;
  std::vector<double> k_fm;             ///< the nodes of the k integral (fm^-1)
  std::vector<double> phase_shift_rad;  ///< delta at each node, continued across the bin
  std::vector<double> u;                ///< u_bin(r_i) (fm^-1/2)
  std::vector<double> potential_MeV;    ///< V(r_i) with <2 l.s> (as the scattering states')
  std::vector<double> coulomb_MeV;
};

/// The bin of Eq. (14.3) at the centre E (MeV) and width |ER| (MeV) with NK nodes (FRESCO's E = -BE, ER < 0, NK).
/// Throws std::invalid_argument for NK < 3 (Simpson: NK odd), a width <= 0 or the reasons of SolveScatteringState.
ContinuumBin MakeBin(const OpticalPotential &potential, const ContinuumRequest &request, double centre_MeV,
                     double width_MeV, int nodes, BinWeight weight, BinRule rule);
ContinuumBin MakeBinIn(const RealPotential &potential, const ContinuumRequest &request, double centre_MeV,
                       double width_MeV, int nodes, BinWeight weight, BinRule rule);

/// What FRESCO prints for a bin ("Norm = ..., RMS = ... & D0 = ..."), on the bin's grid r_1 ... r_{n-2} by the
/// rectangle rule as FRESCO sums them: the norm (int u^2 dr)^1/2, the rms radius (int r^2 u^2 dr)^1/2 / norm and,
/// for l = 0, D0 = sqrt(4 pi) int r V u dr (MeV fm^3/2).
struct BinListing {
  double norm = 0.0, rms_fm = 0.0, d0_MeV_fm32 = 0.0;
};
BinListing ListBin(const ContinuumBin &bin);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
