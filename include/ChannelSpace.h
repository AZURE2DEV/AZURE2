#ifndef CHANNELSPACE_H
#define CHANNELSPACE_H

#include "Constants.h"
#include <cstdint>
#include <memory>
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

/// False when AZURE_LEVEL_SPACE is set (to anything but "0"): every J-group then takes the level-space route.
bool Enabled();

/// Largest \f$\kappa_\lambda = |e_\lambda - E|^{-1} \sum_c \Gamma_{\lambda c}^2 |\Lambda_c|\f$ the channel route accepts.
extern const double kMaxPoleRatio;
/// Levels with a larger |Gamma_k|^2_Lambda / |e_k - E| are deflated (kept out of R, in the bordered system).
extern const double kDeflateRatio;
/// At most this many deflated levels; more send the J-group to the level route.
extern const int kMaxDeflated;
/// Largest relative error of a level's term in R, from the round-off of \f$e_\lambda - E\f$, the channel route accepts.
extern const double kMaxPoleError;
/// Round-off of \f$e_\lambda - E\f$, in units of epsilon times the largest energy involved.
extern const double kEnergyUlps;
/// Smallest pivot ratio of \f$1 - R\Lambda\f$ the channel route accepts.
extern const double kMinPivotRatio;

/// Energy-independent formal representation of one J-group's level matrix.
struct FormalLevels {
  int mode = -1;         ///< 0 standard, 1 Brune, 2 Park
  int n = 0;             ///< active levels
  int nch = 0;           ///< channels
  bool usable = false;   ///< the channel route may be taken (set-up succeeded, fewer active channels than levels)
  bool identity = true;  ///< V = 1 (standard level matrix)
  /*!
   * No formal representation (N not positive definite: Brune or Park
   * parameters without a real standard equivalent): the channel route then
   * solves the real pencil \f$P = K - E N\f$ at each energy,
   * \f$R = G^T P^{-1} G\f$, which still spares the complex level-matrix fill and
   * factorization.  K and N are kept for it, and `gamma` holds G.
   */
  bool pencil = false;
  /// Pencil route: Qt K Z = H (upper Hessenberg), Qt N Z = T (upper triangular).
  std::vector<long double> H, T;
  std::vector<double> Qt, Z;  ///< n x n row-major
  std::vector<double> gq, gz; ///< Qt G and Z^T G over the active channels, n x m
  bool adjointOK = false;  ///< formal energies separated enough for the eigen-decomposition adjoint
  bool complexE = false;   ///< formal energies complex (RMC photon widths)
  int status = 0;        ///< why not usable: 0 usable, 1 size, 2 set-up (eigen-decomposition)
  std::vector<complex> e;  ///< formal energies (n)
  /// Low part of the formal energies (extended-precision eigenvalue minus e):
  /// e - E is formed as (e - E) + elo, so that a narrow resonance at a data
  /// energy sees its distance to E as accurately as on the level route.
  std::vector<double> elo;
  double emax = 0.0;       ///< largest |e|
  std::vector<double> V;   ///< n x n row-major; column k is eigenvector k
  std::vector<double> gcut;  ///< level-space widths with the 1e-12 cutoff, n x nch row-major
  std::vector<double> gamma;  ///< formal widths \f$\Gamma = V^T G\f$ over the active channels, n x m row-major
  std::vector<int> act;    ///< active channels (0-based), m of them
  std::vector<int> pos;    ///< channel (0-based) -> position in act, or -1
  /*!
   * Per channel constant added to the Lo-matrix element to give the kernel:
   * \f$\Lambda_c = L_c + \lambda_c\f$.  For the Brune and Park particle channels
   * \f$\lambda_c = B_c - \bar S_c\f$: their boundary condition cancels, and the
   * reference shift \f$\bar S_c\f$ (the mean of \f$S_c(E_\lambda)\f$ over the coupled
   * levels, moved into K as \f$-\bar S_c G_c G_c^T\f$) keeps the real part of the
   * kernel small, so \f$1 - R\Lambda\f$ cancels less.
   */
  std::vector<double> lamShift;
  std::vector<std::uint64_t> key;       ///< levels (and their ALevel::Stamp) it was built from
  std::vector<std::uint64_t> scratch_;  ///< key of the current parameters (reused)
};

/*!
 * The formal representation of J-group `jg` kept for this thread, shared by
 * every AMatrixFunc the thread creates (the gradient creates one per point).
 * Valid only after Update().
 */
std::shared_ptr<FormalLevels> Cached(const JGroup *jg);

/// The level data LevelMatrixAdjoint reads, per J-group (levels 1..n), kept per thread.
struct LevelArrays {
  std::vector<std::uint64_t> key;  ///< levels and their ALevel::Stamp
  int n = 0, nch = 0;
  std::vector<double> E;    ///< fitted level energies (n)
  std::vector<double> G;    ///< fitted widths, n x nch
  std::vector<double> S;    ///< shift functions at the level energies, n x nch
  std::vector<double> SG;   ///< S times G
  std::vector<double> inv;  ///< 1/(E_r - E_c) off the diagonal, n x n
};

/*!
 * The LevelArrays of J-group `jg` (all its levels) for the current
 * parameters: rebuilt when a level's parameters changed, so the n^2
 * divisions of `inv` are paid once per parameter set, not once per
 * level-matrix adjoint.
 */
const LevelArrays &CurrentLevelArrays(JGroup *jg);

/*!
 * Brings `f` up to date with the parameters of J-group `jg` (active levels
 * `levels`, 1-based original numbers).  Returns true if it was rebuilt.  The
 * standard level matrix is rebuilt on every call (it costs what comparing
 * would); the Brune and Park ones only when a parameter they depend on changed.
 */
bool Update(JGroup *jg, const std::vector<int> &levels, const Config &configure, FormalLevels &f);

/*!
 * In-place LU factorization with partial pivoting of the m x m row-major
 * complex matrix `a`.  Returns the ratio of the smallest to the largest pivot
 * modulus (0 if a pivot vanishes or is not finite).
 */
double LUFactor(complex *a, int *piv, int m);
/// Solves with the factors of LUFactor, overwriting `b`.
void LUSolve(const complex *lu, const int *piv, int m, complex *b);
/// LU factorization with partial pivoting of an upper Hessenberg n x n row-major
/// matrix, in place; returns the pivot ratio (0 if a pivot vanishes).
double HessenbergLU(double *a, char *swapped, int n);
/// Solves with the factors of HessenbergLU for `nrhs` right-hand sides (b is n x nrhs row-major).
void HessenbergSolve(const double *lu, const char *swapped, int n, double *b, int nrhs);
/// Real LUFactor.
double LUFactorReal(double *a, int *piv, int m);
/// Solves with the factors of LUFactorReal for `nrhs` right-hand sides (b is m x nrhs row-major).
void LUSolveReal(const double *lu, const int *piv, int m, double *b, int nrhs);

}  // namespace ChannelSpace

#endif
