#ifndef AZURE2_PHYSICS_CHANNEL_SOLVE_H
#define AZURE2_PHYSICS_CHANNEL_SOLVE_H

// The channel-space solution of the level matrix of one J^pi block at one energy (the Woodbury identity;
// docs/source/developer/level_matrix.rst, "Channel-Space Solution").  The level matrix is, in every formalism,
//
//   M(E) = K - E N - G Lambda(E) G^T,
//
// K, N real symmetric over the levels (level parameters only), G the reduced widths, Lambda diagonal over the
// channels.  With the formal levels of the block (K v = e N v, V^T N V = 1: energies e, widths Gamma = V^T G; a
// FormalBlock) the bilinear forms the T-matrix needs are
//
//   G^T M^-1 G = (1 - R Lambda)^-1 R,   R = Gamma^T D^-1 Gamma,   D = diag(e - E),
//
// a solve over the channels in place of one over the levels.  The functions here are the part that runs at each
// energy: R, the deflation of a level near a pole, the factors of the bordered system, and the products with
// A = M^-1 that follow from them.  They take the formal levels and the kernel as arrays; which block takes this
// route, the memo of an energy and the set-up of the formal levels are the caller's.

#include <complex>
#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {
namespace channel_space {

/// Largest kappa_l = |e_l - E|^-1 sum_c Gamma_lc^2 |Lambda_c| the channel route accepts.
constexpr double kMaxPoleRatio = 1.0e3;
/// Smallest pivot ratio of 1 - R Lambda the channel route accepts.
constexpr double kMinPivotRatio = 1.0e-10;
/// Largest relative error of a level's term in R, from the round-off of e_l - E, the channel route accepts.
constexpr double kMaxPoleError = 1.0e-11;
/// Levels with a larger |Gamma_k|^2_Lambda / |e_k - E| are deflated (kept out of R, in the bordered system).
constexpr double kDeflateRatio = 1.0e2;
/// At most this many deflated levels; more send the block to the level route.
constexpr int kMaxDeflated = 8;
/// Round-off of e_l - E, in units of epsilon times the largest energy involved.
constexpr double kEnergyUlps = 2.0;

/// Energy-independent formal representation of one block's level matrix.
struct FormalBlock {
  int mode = -1;         ///< 0 standard, 1 Brune, 2 Park
  int n = 0;             ///< active levels
  int nch = 0;           ///< channels
  bool usable = false;   ///< the channel route may be taken (set-up succeeded, fewer active channels than levels)
  bool identity = true;  ///< V = 1 (standard level matrix)
  /*!
   * No formal representation (N not positive definite: Brune or Park
   * parameters without a real standard equivalent): the channel route then
   * solves the real pencil P = K - E N at each energy,
   * R = G^T P^-1 G, which still spares the complex level-matrix fill and
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
  std::vector<std::complex<double>> e;  ///< formal energies (n)
  /// Low part of the formal energies (extended-precision eigenvalue minus e):
  /// e - E is formed as (e - E) + elo, so that a narrow resonance at a data
  /// energy sees its distance to E as accurately as on the level route.
  std::vector<double> elo;
  double emax = 0.0;       ///< largest |e|
  std::vector<double> V;   ///< n x n row-major; column k is eigenvector k
  std::vector<double> gcut;  ///< level-space widths with the 1e-12 cutoff, n x nch row-major
  std::vector<double> gamma;  ///< formal widths Gamma = V^T G over the active channels, n x m row-major
  std::vector<int> act;    ///< active channels (0-based), m of them
  std::vector<int> pos;    ///< channel (0-based) -> position in act, or -1
  /*!
   * Per channel constant added to the Lo-matrix element to give the kernel:
   * Lambda_c = L_c + lambda_c.  For the Brune and Park particle channels
   * lambda_c = B_c - Sbar_c: their boundary condition cancels, and the
   * reference shift Sbar_c (the mean of S_c(E_l) over the coupled
   * levels, moved into K as -Sbar_c G_c G_c^T) keeps the real part of the
   * kernel small, so 1 - R Lambda cancels less.
   */
  std::vector<double> lamShift;
};

/// What the channel route holds of one block at one energy: the kernel, R, the factors, and the caller's memo of
/// the columns and solves it has asked for (the members below `F` belong to the caller; the functions here use
/// w1, w2, w4, w3b, wp, wz and rreal as scratch).
struct ChannelSolve {
  std::vector<std::complex<double>> dinv;   ///< 1/(e_lambda - E), n
  std::vector<std::complex<double>> lam;    ///< Lambda per channel, nch
  std::vector<std::complex<double>> lamAct; ///< Lambda over the active channels, m
  std::vector<std::complex<double>> R;      ///< R over the active channels (without deflated levels), m x m
  std::vector<std::complex<double>> lu;     ///< LU factors of the bordered system (1 - Lambda R without deflated levels)
  std::vector<int> piv;
  std::vector<std::complex<double>> F;      ///< (1 - R Lambda)^-1 R, column b at F[b*m]
  std::vector<char> fcol;      ///< columns of F computed
  std::vector<std::vector<std::complex<double>>> usolve;  ///< formal channel solves A_f Gamma_c (per channel)
  std::vector<char> uvalid;
  std::vector<std::vector<std::complex<double>>> inverse;  ///< A, materialized on request
  bool inverseValid = false;
  std::vector<std::complex<double>> w1, w2, w3, w4;     ///< scratch
  std::vector<std::complex<double>> w3b, wp, wz;        ///< scratch of the bordered solves
  std::vector<int> defl;       ///< deflated formal levels (near a pole), kept in the bordered system
  std::vector<char> isDefl;
  std::vector<double> rreal;   ///< scratch of the real R
  std::vector<double> P;       ///< pencil route: LU factors of H - E T (upper Hessenberg)
  std::vector<char> hswap;     ///< pencil route: row interchanges of that LU
  std::vector<double> X;       ///< pencil route: (H - E T)^-1 Qt G over the active channels, n x m
};

/// The outcome of the per-energy part: the route is taken, or the level route must be (near a formal pole, where
/// R is dominated by one level and the Woodbury form loses digits the level matrix keeps; or 1 - R Lambda
/// ill-conditioned).
enum class Route { kTaken, kPole, kPivot };

/*!
 * The per-energy part of the channel route for one block: the kernel Lambda_c = lo[c] + lamShift[c] (lo[c] the
 * Lo-matrix element L_c(E) - B_c of each of the block's nch channels), D^-1, R = Gamma^T D^-1 Gamma over the
 * coupled channels and the factors of 1 - R Lambda (with the deflated levels, of the bordered system).  Clears
 * the memo of `st`.
 */
Route solve_channel_space(const FormalBlock &f, const std::complex<double> *lo, double E_MeV, ChannelSolve &st);

/// R of the formal route, R = Gamma^T D^-1 Gamma; false near a pole (see Route).  Needs st.lamAct.
bool solve_formal(const FormalBlock &f, ChannelSolve &st, double E_MeV);
/// R of the pencil route, R = G^T (K - E N)^-1 G from a real LU of the pencil; false near a formal pole.
bool solve_pencil(const FormalBlock &f, ChannelSolve &st, double E_MeV);

/// Solves the bordered system for the right-hand side (wp over the deflated levels, wz over the channels), in
/// place: on return wp holds y_p and wz holds z = Lambda Gamma^T y.
void solve_bordered(const FormalBlock &f, ChannelSolve &st, std::complex<double> *wp, std::complex<double> *wz);

/// Computes column b (a position in f.act) of F = G^T A G into st.F and marks it.
void compute_f_column(const FormalBlock &f, ChannelSolve &st, int b);

/// Column b of F = Gamma^T A_f Gamma = G^T A G (m values, by position in f.act), computed once per energy.
inline const std::complex<double> *f_column(const FormalBlock &f, ChannelSolve &st, int b) {
  if (!st.fcol[b]) compute_f_column(f, st, b);
  return &st.F[(std::size_t)b * f.act.size()];
}

/// y = A_f w in the coordinates of the formal levels (formal route).
void formal_apply(const FormalBlock &f, ChannelSolve &st, const std::complex<double> *w, std::complex<double> *y);
/// out = A u over the block's active levels (formal or pencil route).
void apply_a(const FormalBlock &f, ChannelSolve &st, const std::complex<double> *u, std::complex<double> *out);
/// u^T A v over the block's active levels (vectors that are not columns of G: the masked capture pathways).
std::complex<double> channel_bilinear(const FormalBlock &f, ChannelSolve &st, const std::complex<double> *u,
                                      const std::complex<double> *v);

/*!
 * In-place LU factorization with partial pivoting of the m x m row-major
 * complex matrix `a`.  Returns the ratio of the smallest to the largest pivot
 * modulus (0 if a pivot vanishes or is not finite).
 */
double lu_factor(std::complex<double> *a, int *piv, int m);
/// Solves with the factors of lu_factor, overwriting `b`.
void lu_solve(const std::complex<double> *lu, const int *piv, int m, std::complex<double> *b);
/// LU factorization with partial pivoting of an upper Hessenberg n x n row-major
/// matrix, in place; returns the pivot ratio (0 if a pivot vanishes).
double hessenberg_lu(double *a, char *swapped, int n);
/// Solves with the factors of hessenberg_lu for `nrhs` right-hand sides (b is n x nrhs row-major).
void hessenberg_solve(const double *lu, const char *swapped, int n, double *b, int nrhs);
/// Real lu_factor.
double lu_factor_real(double *a, int *piv, int m);
/// Solves with the factors of lu_factor_real for `nrhs` right-hand sides (b is m x nrhs row-major).
void lu_solve_real(const double *lu, const int *piv, int m, double *b, int nrhs);

}  // namespace channel_space
}  // namespace physics
}  // namespace azure2

#endif
