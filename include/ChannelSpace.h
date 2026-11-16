#ifndef CHANNELSPACE_H
#define CHANNELSPACE_H

#include "Constants.h"
#include <cstdint>
#include "LevelState.h"
#include "azure2/physics/ChannelSolve.h"
#include "azure2/physics/CollisionMatrixAdjoint.h"
#include <memory>
#include <unordered_map>
#include <vector>

class JGroup;
class Config;

/// Channel-space solution of the level matrix (Woodbury identity).

/*!
 * The level matrix of a J-group is, in every formalism AZURE2 offers, of the form
 *
 *   \f$ M(E) = K - E\,N - G\,\Lambda(E)\,G^T \f$,
 *
 * with \f$K, N\f$ real symmetric n x n matrices that depend on the level
 * parameters only, \f$G\f$ the n x c reduced widths and \f$\Lambda\f$ diagonal in
 * the channels (\f$\Lambda_c = S_c - B_c + iP_c\f$, plus \f$B_c\f$ for the particle
 * channels of the Brune and Park level matrices, whose boundary condition cancels):
 *
 *   - standard:  \f$K = \mathrm{diag}(E_\lambda)\f$, \f$N = 1\f$;
 *   - Brune:     \f$K_{\lambda\mu} = \tilde E_\lambda\delta_{\lambda\mu} + \sum_c \gamma_{\lambda c}\gamma_{\mu c} s^0_{\lambda\mu c}\f$,
 *                \f$N_{\lambda\mu} = \delta_{\lambda\mu} - \sum_c \gamma_{\lambda c}\gamma_{\mu c} s^1_{\lambda\mu c}\f$,
 *                where the shift-function term of the Brune level matrix is
 *                \f$s^0 + E s^1\f$ (\f$s^0_{\lambda\lambda} = S_{\lambda c}\f$, \f$s^1_{\lambda\lambda} = 0\f$);
 *   - Park:      \f$K_{\lambda\mu} = E_\lambda J_{\lambda\mu} + \sum_c \gamma_{\lambda c}\gamma_{\mu c} S_{\lambda c}\f$, \f$N = J\f$.
 *
 * The generalized eigenproblem \f$K v = e N v\f$ (\f$V^T N V = 1\f$) turns it into a
 * standard (formal) R-matrix with energies \f$e_\lambda\f$ and widths
 * \f$\Gamma = V^T G\f$:
 *
 *   \f$ M^{-1} = V\,[D - \Gamma\Lambda\Gamma^T]^{-1} V^T, \quad D = \mathrm{diag}(e_\lambda - E) \f$,
 *
 * and by the Woodbury identity the bilinear forms the T-matrix needs are
 *
 *   \f$ G^T M^{-1} G = (1 - R\Lambda)^{-1} R, \quad R = \Gamma^T D^{-1} \Gamma \f$,
 *
 * a c x c solve instead of an n x n one.  FormalLevels holds the
 * energy-independent part of this (once per parameter set and J-group); the
 * per-energy part lives in AMatrixFunc.
 */
namespace ChannelSpace {

/// False when the run's AZURE_LEVEL_SPACE is set (Config::environment.levelSpaceOnly): every J-group then takes the
/// level-space route.
bool Enabled(const Config &configure);

/// The limits of the channel route (azure2/physics/ChannelSolve.h says what each one bounds).
using azure2::physics::channel_space::kMaxPoleRatio;
using azure2::physics::channel_space::kDeflateRatio;
using azure2::physics::channel_space::kMaxDeflated;
using azure2::physics::channel_space::kMaxPoleError;
using azure2::physics::channel_space::kEnergyUlps;
using azure2::physics::channel_space::kMinPivotRatio;

/// Energy-independent formal representation of one J-group's level matrix: the arrays the per-energy solve reads
/// (azure2::physics::channel_space::FormalBlock) and the key of the parameters they were built from.
struct FormalLevels : azure2::physics::channel_space::FormalBlock {
  std::vector<std::uint64_t> key;       ///< the state, the levels and their LevelState::Version it was built from
  std::vector<std::uint64_t> scratch_;  ///< key of the current parameters (reused)
};

/*!
 * The formal representation of J-group `jg` kept for this thread, shared by
 * every AMatrixFunc the thread creates (the gradient creates one per point).
 * Valid only after Update().
 */
struct Memo;
std::shared_ptr<FormalLevels> Cached(Memo &memo, const JGroup *jg);

/// The level data LevelMatrixAdjoint reads, per J-group (levels 1..n), kept in the caller's Memo.
struct LevelArrays {
  std::vector<std::uint64_t> key;  ///< the state, the levels and their LevelState::Version
  int n = 0, nch = 0;
  std::vector<double> E;    ///< fitted level energies (n)
  std::vector<double> G;    ///< fitted widths, n x nch
  std::vector<double> S;    ///< shift functions at the level energies, n x nch
  std::vector<double> SG;   ///< S times G
  std::vector<double> inv;  ///< 1/(E_r - E_c) off the diagonal, n x n
};

/*!
 * The LevelArrays of J-group `jg` (all its levels; number `j` in `state`, the
 * levels' fit values) for the current parameters: rebuilt when a level's
 * version changed, so the n^2 divisions of `inv` are paid once per parameter
 * set, not once per level-matrix adjoint.
 */
const LevelArrays &CurrentLevelArrays(Memo &memo, const JGroup *jg, const LevelState &state, int j);

/*!
 * Derivatives of K and N of a J-group's level matrix with respect to its free
 * level parameters (AMatrixFunc::LevelBasis): for each, dK = e_l k^T + k e_l^T
 * and dN = e_l q^T + q e_l^T over all levels (standard, Brune and Park level
 * matrices).
 */
struct LevelDerivatives {
  std::vector<std::uint64_t> key;  ///< state, levels, versions, formalism and free mask
  /// One free parameter: level, channel (-1: the level energy), the vectors k and q
  /// (azure2/physics/CollisionMatrixAdjoint.h, level_parameter_adjoint).
  typedef azure2::physics::LevelParameterDerivative Entry;
  std::vector<Entry> entries;
};

/*!
 * The LevelDerivatives of J-group `jg` for the free parameters (`freeE`,
 * `freeG`, by level and level/channel) and the current parameters.  `mode`
 * selects the level matrix: 0 standard, 1 Brune (shift-function derivatives
 * from `shiftDeriv` of the group, [level-1][ch-1]), 2 Park (the first and
 * second derivatives and the overlap J of `state`,
 * CompoundSystem::CalcShiftFunctions); `particle` marks the channels with the shift
 * terms (by channel, 0-based).  `j` is the group's number in `state`.  Kept
 * per thread, rebuilt when a level's version or the mask changes.
 */
const LevelDerivatives &CurrentLevelDerivatives(Memo &memo, const JGroup *jg, const LevelState &state, int j, int mode,
                                                const std::vector<char> &particle,
                                                const std::vector<char> &freeE,
                                                const std::vector<std::vector<char>> &freeG,
                                                const std::vector<std::vector<double>> *shiftDeriv);

/*!
 * Brings `f` up to date with the parameters of J-group `jg` (number `j` of
 * `state`; active levels `levels`, 1-based original numbers).  Returns true if it was rebuilt.  The
 * standard level matrix is rebuilt on every call (it costs what comparing
 * would); the Brune and Park ones only when a parameter they depend on changed.
 */
bool Update(const JGroup *jg, const LevelState &state, int j, const std::vector<int> &levels, const Config &configure,
            FormalLevels &f);

/*!
 * What the channel-space route keeps between the energies of one thread's tasks: per J group the formal-basis
 * solver, the level arrays and the level-matrix derivatives of the present parameters.  Owned by the caller (the
 * scratch of an energy task, azure2::TaskScratch, or an AMatrixFunc of its own), never process-wide.  Each entry is
 * keyed by the J group's address and the stamps of its levels (ALevel::Stamp, unique in the process), so an entry is
 * used only for the levels it was built from; at most kMaxGroups groups (cleared when full).
 */
struct Memo {
  std::unordered_map<const JGroup *, std::shared_ptr<FormalLevels>> formal;
  std::unordered_map<const JGroup *, LevelArrays> arrays;
  std::unordered_map<const JGroup *, LevelDerivatives> derivatives;
  static constexpr std::size_t kMaxGroups = 1024;
};

/// The LU routines of the channel route, in azure2/physics/ChannelSolve.h (names kept for the callers here).
inline double LUFactor(complex *a, int *piv, int m) { return azure2::physics::channel_space::lu_factor(a, piv, m); }
inline void LUSolve(const complex *lu, const int *piv, int m, complex *b) {
  azure2::physics::channel_space::lu_solve(lu, piv, m, b);
}
inline double HessenbergLU(double *a, char *swapped, int n) {
  return azure2::physics::channel_space::hessenberg_lu(a, swapped, n);
}
inline void HessenbergSolve(const double *lu, const char *swapped, int n, double *b, int nrhs) {
  azure2::physics::channel_space::hessenberg_solve(lu, swapped, n, b, nrhs);
}
inline double LUFactorReal(double *a, int *piv, int m) {
  return azure2::physics::channel_space::lu_factor_real(a, piv, m);
}
inline void LUSolveReal(const double *lu, const int *piv, int m, double *b, int nrhs) {
  azure2::physics::channel_space::lu_solve_real(lu, piv, m, b, nrhs);
}

}  // namespace ChannelSpace

#endif
