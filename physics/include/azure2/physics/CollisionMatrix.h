#ifndef AZURE2_PHYSICS_COLLISION_MATRIX_H
#define AZURE2_PHYSICS_COLLISION_MATRIX_H

// The collision matrix of one J^pi block at one energy, from arrays (docs/upgrade/physics_functions.md):
//
//   the inverse level matrix  (A^-1)_lm = (E_l - E) delta_lm - sum_c gamma_lc gamma_mc (S_c(E) - B_c + i P_c)
//   (standard; fill_level_matrix has the Brune and the Park forms), its factorization (LevelMatrixSolver), the
//   width bilinear forms g_c^T A g_c', and from them per pair of channels
//
//   U_cc' = Omega_c Omega_c' (delta_cc' + 2 i sqrt(P_c P_c') g_c^T A g_c'),   Omega_c = exp(i (omega_c - phi_c)),
//   T_cc' = exp(2 i omega_c) delta_cc' - U_cc'.
//
// Nothing here knows a project, a data point or the engine's classes: a caller hands over the level parameters of
// the block and the channel functions at the energy.  The engine (AMatrixFunc) is one such caller; a unit test, a
// Python binding or another program are others.  Levels and channels count from 0; matrices are row-major.
//
// The expressions and their order are those of the engine's routines they were moved from; the output files of
// every project depend on them to the last digit.

#include "azure2/physics/LevelMatrixSolver.h"

#include <complex>
#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {

/// A reduced-width amplitude below this (MeV^1/2) couples nothing: the level matrix and the bilinear forms skip it.
constexpr double kWidthCutoff = 1.0e-12;

/// The level parameters of one J^pi block and its channel functions at one energy, as fill_level_matrix reads
/// them.  n levels (those in the R-matrix), nc channels; the pointers are borrowed.
struct LevelMatrixTables {
  int n = 0;                                 ///< levels
  int nc = 0;                                ///< channels
  const double *E = nullptr;                 ///< level energies E_l (MeV), n
  const double *gamma = nullptr;             ///< reduced-width amplitudes gamma_lc (MeV^1/2), n x nc
  const double *shift = nullptr;             ///< S_c(E_l), n x nc (read for the particle channels of Brune and Park)
  const double *shift_d1 = nullptr;          ///< dS_c/dE at E_l, n x nc (Brune and Park, between two levels)
  const double *shift_d2 = nullptr;          ///< d2S_c/dE2 at E_l, n x nc (the same)
  const double *park_norm = nullptr;         ///< J_ll of Park (2021) Eq. (28), n (Park)
  const std::complex<double> *lo = nullptr;  ///< L_c(E) - B_c = S_c(E) - B_c + i P_c, nc
  const double *boundary = nullptr;          ///< B_c, nc (read for the particle channels of Brune and Park)
  const char *rad_type = nullptr;            ///< 'P' particle, 'E' or 'M' photon, per channel, nc
  bool brune = false;                        ///< the level parameters are Brune's (2002): the on-resonance ones
  bool park = false;                         ///< the level parameters are Park's (2021); takes precedence over brune
  bool rmc = false;                          ///< R-matrix capture: i gamma^2 added to the photon channels' diagonal
};

/*!
 * The inverse level matrix of the block at energy E (MeV), n x n, into M (resized).  An element that is not finite is
 * left zero.  A level with no width above kWidthCutoff has a zero row and column but for its diagonal, which
 * vanishes at E = E_l; it contributes nothing at any energy (everything that uses A weighs it by its widths), so a
 * unit diagonal there keeps the matrix regular.
 */
void fill_level_matrix(const LevelMatrixTables &tables, double E_MeV, LevelMatrixSolver::matrix_c &M);

/// g^T x over the levels: with x = A g_c' (LevelMatrixSolver::Solve) the bilinear form g_c^T A g_c'.
inline std::complex<double> width_bilinear(const double *g, const std::complex<double> *x, std::size_t n) {
  std::complex<double> sum(0.0, 0.0);
  for (std::size_t i = 0; i < n; i++) sum += g[i] * x[i];
  return sum;
}

/// The resonant amplitude of the pathway c -> c': 2 i sqrt(P_c) sqrt(P_c') g_c^T A g_c'.
inline std::complex<double> pathway_amplitude(double sqrt_pen_in, double sqrt_pen_out, std::complex<double> bilinear) {
  const std::complex<double> sqrtPenEn = sqrt_pen_in, sqrtPenEx = sqrt_pen_out;
  return 2.0 * std::complex<double>(0.0, 1.0) * sqrtPenEn * sqrtPenEx * bilinear;
}

/// The phases of one channel at one energy: exp(i omega_c) (Coulomb) and exp(-i phi_c) (hard sphere).
struct ChannelPhases {
  std::complex<double> coulomb;
  std::complex<double> hard_sphere;
};

/// The T-matrix element of an internal pathway c -> c' (`diagonal`: c = c') from its resonant amplitude:
/// exp(2 i omega_c) - Omega_c^2 (1 + u) on the diagonal, - Omega_c Omega_c' u off it.
inline std::complex<double> t_matrix_element(bool diagonal, const ChannelPhases &in, const ChannelPhases &out,
                                             std::complex<double> umatrix) {
  const std::complex<double> uphase = in.coulomb * in.hard_sphere * out.coulomb * out.hard_sphere;
  if (diagonal) {
    const std::complex<double> tphase = in.coulomb * in.coulomb;
    return tphase - uphase * (1.0 + umatrix);
  }
  return -uphase * umatrix;
}

/// The collision-matrix element U_cc' from the resonant amplitude: Omega_c Omega_c' (delta_cc' + u).
inline std::complex<double> collision_element(bool diagonal, const ChannelPhases &in, const ChannelPhases &out,
                                              std::complex<double> umatrix) {
  const std::complex<double> uphase = in.coulomb * in.hard_sphere * out.coulomb * out.hard_sphere;
  return diagonal ? uphase * (1.0 + umatrix) : uphase * umatrix;
}

/// The T-matrix element of a pathway of R-matrix capture (the photon as a channel): -i u.
inline std::complex<double> rmc_t_matrix_element(std::complex<double> umatrix) {
  return std::complex<double>(0.0, -1.0) * umatrix;
}

/// The T-matrix element of an external-capture pathway: the normalisation of the final state (its reduced width
/// amplitude, the square root of its normalisation factor, the conversion factor of the multipole) times the
/// amplitude of the external integral.
inline std::complex<double> external_capture_element(double ec_gamma, double sqrt_nf, double conversion,
                                                     std::complex<double> amplitude) {
  const double ecNormParam = ec_gamma * sqrt_nf * conversion;
  return ecNormParam * amplitude;
}

/// The factor of channel capture (the external integral fed through the resonant entrance channel):
/// 2 i sqrt(P_c) u^T A v.
inline std::complex<double> channel_capture_factor(double sqrt_pen, std::complex<double> bilinear) {
  return 2.0 * std::complex<double>(0.0, 1.0) * sqrt_pen * bilinear;
}

/// The channel functions of the block's channels at the energy, as collision_matrix reads them (nc of each).
struct ChannelTablesAtE {
  const std::complex<double> *coulomb_phase = nullptr;      ///< exp(i omega_c)
  const std::complex<double> *hard_sphere_phase = nullptr;  ///< exp(-i phi_c)
  const double *sqrt_penetrability = nullptr;               ///< sqrt(P_c)
};

/// The level matrix A and the collision and T-matrices of one block at one energy.
struct CollisionMatrix {
  int n = 0;                            ///< levels
  int nc = 0;                           ///< channels
  bool singular = false;                ///< the inverse level matrix has no inverse: A, U and T are zero
  std::vector<std::complex<double>> A;  ///< n x n
  std::vector<std::complex<double>> U;  ///< nc x nc, U[c * nc + c'] for the pathway c -> c'
  std::vector<std::complex<double>> T;  ///< nc x nc
  std::complex<double> u(int c, int cp) const { return U[(std::size_t)c * nc + cp]; }
  std::complex<double> t(int c, int cp) const { return T[(std::size_t)c * nc + cp]; }
};

/*!
 * The whole block at one energy on the level route: the inverse level matrix, one factorization, one solve per
 * channel, and every pair of channels.  The engine does the same per pathway of a data point (and, for a block
 * with fewer coupled channels than levels, through the channel route of ChannelSolve.h, which agrees to round-off).
 */
CollisionMatrix collision_matrix(const LevelMatrixTables &tables, double E_MeV, const ChannelTablesAtE &channels);

}  // namespace physics
}  // namespace azure2

#endif
