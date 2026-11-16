#ifndef AZURE2_PHYSICS_CHANNEL_SOLVE_ADJOINT_H
#define AZURE2_PHYSICS_CHANNEL_SOLVE_ADJOINT_H

// The reverse mode of the channel route (ChannelSolve.h) for one J^pi block (docs/upgrade/physics_functions.md,
// stage 5), and the solves A g_c it starts from.  In the formal representation of the block's level matrix,
//
//   M_f = D - Gamma Lambda Gamma^T,   D = diag(e - E),   F_ab = Gamma_a^T M_f^-1 Gamma_b,
//
// the cotangent of a bilinear form F_ab goes to the formal energies e and widths Gamma in O(n c) operations
// (formal_bilinear_adjoint).  For the standard level matrix (V = 1) those are the level parameters.  For Brune's and
// Park's they are the eigenvalues and the projected widths of K v = e N v, V^T N V = 1, Gamma = V^T G: a pass sums
// their cotangents over its energies and turns them into cotangents of K, N and G once
// (resolve_formal_cotangents), which the structural adjoint of CollisionMatrixAdjoint.h takes to E_l and gamma_lc.
//
// Conventions as in CollisionMatrixAdjoint.h: zbar = dL/dRe z + i dL/dIm z; contributions are added through the
// caller's sink (energy, width, and widths(level): the row of a level's width cotangents, nc values) in the order
// of the engine's routines these were moved from (AMatrixFunc::FormalSolve, GetRawChannelSolve,
// FormalBilinearAdjoint, ResolveFormalCotangents).  The memo of an energy (which solves and columns exist) is the
// caller's, in ChannelSolve.

#include "azure2/physics/ChannelSolve.h"
#include "azure2/physics/CollisionMatrixAdjoint.h"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {
namespace channel_space {

/// u = A_f Gamma_c over the formal levels, for channel ch0 (from 0) of the block: zero for a channel without widths
/// (not among f.act).  Uses st.w2.
inline void formal_channel_solve(const FormalBlock &f, ChannelSolve &st, int ch0, std::vector<std::complex<double>> &u) {
  const int n = f.n, m = (int)f.act.size();
  u.assign(n, std::complex<double>(0.0, 0.0));
  const int b = f.pos[ch0];
  if (b >= 0) {
    std::vector<std::complex<double>> &w = st.w2;
    w.resize(n);
    for (int k = 0; k < n; k++) w[k] = f.gamma[(std::size_t)k * m + b];
    formal_apply(f, st, w.data(), u.data());
  }
}

/// A g_c over the levels from the formal solve u_c = A_f Gamma_c: V u_c (u_c itself for the standard matrix).
inline void level_solve_from_formal(const FormalBlock &f, const std::vector<std::complex<double>> &uf,
                                    std::vector<std::complex<double>> &solve) {
  const int n = f.n;
  if (f.identity)
    solve = uf;
  else
    for (int r = 0; r < n; r++) {
      const double *vr = &f.V[(std::size_t)r * n];
      std::complex<double> sum(0.0, 0.0);
      for (int k = 0; k < n; k++) sum += vr[k] * uf[k];
      solve[r] = sum;
    }
}

/// A g_c over the levels on the pencil route, for the channel at position b of f.act (b < 0: zero):
/// A g_c = Z X (e_c + z),  (1 - Lambda R) z = Lambda R e_c.  `solve` has n values.  Uses st.wz and st.w4.
inline void pencil_channel_solve(const FormalBlock &f, ChannelSolve &st, int b, std::vector<std::complex<double>> &solve) {
  typedef std::complex<double> complex;
  const int n = f.n, m = (int)f.act.size();
  std::fill(solve.begin(), solve.end(), complex(0.0, 0.0));
  if (b >= 0) {
    std::vector<complex> &z = st.wz;
    z.resize(m);
    for (int a = 0; a < m; a++) z[a] = st.lamAct[a] * st.R[(std::size_t)a * m + b];
    solve_bordered(f, st, nullptr, z.data());
    z[b] += 1.0;
    std::vector<complex> &y = st.w4;
    y.resize(n);
    for (int k = 0; k < n; k++) {
      complex sum(0.0, 0.0);
      for (int a = 0; a < m; a++) sum += st.X[(std::size_t)k * m + a] * z[a];
      y[k] = sum;
    }
    for (int r = 0; r < n; r++) {
      const double *zr = &f.Z[(std::size_t)r * n];
      complex sum(0.0, 0.0);
      for (int k = 0; k < n; k++) sum += zr[k] * y[k];
      solve[r] = sum;
    }
  }
}

/// Column c0 (a channel of the block, from 0) of Y = 1 + Lambda F, over all nch channels: Y_dc = delta_dc +
/// Lambda_d F_dc (F = Gamma^T A_f Gamma; zero outside the coupled channels).  y has nch values.
inline void kernel_column(const FormalBlock &f, ChannelSolve &st, int c0, std::vector<std::complex<double>> &y) {
  std::fill(y.begin(), y.end(), std::complex<double>(0.0, 0.0));
  y[c0] = 1.0;
  const int b = f.pos[c0];
  if (b < 0) return;
  const std::complex<double> *Fb = f_column(f, st, b);
  const int m = (int)f.act.size();
  for (int a = 0; a < m; a++) y[f.act[a]] += st.lamAct[a] * Fb[a];
}

/*!
 * The adjoint of one term beta F_{ch,chp} in the formal representation.  With u_c = M_f^-1 Gamma_c (ua for ch, ub
 * for chp; n values) and the columns ya, yb of Y = 1 + Lambda F (kernel_column; nch values):
 *
 *   ebar_k        += -Re(conj(beta) ua_k ub_k)                       (from dD),
 *   Gammabar_{kd} +=  Re(conj(beta) (ua_k Y_{d,chp} + ub_k Y_{d,ch}))  (the widths as they stand and -Gamma Lambda Gamma^T,
 *                                                                      with Gamma^T u_b = F e_b).
 * O(n c).  The sink's levels are the formal levels k.
 */
template <class Sink>
void formal_bilinear_adjoint(const BilinearTerm &term, const std::complex<double> *ua, const std::complex<double> *ub,
                             const std::complex<double> *ya, const std::complex<double> *yb, int n, int nch,
                             Sink &&sink) {
  typedef std::complex<double> complex;
  const complex cb = std::conj(term.beta);
  for (int k = 0; k < n; k++) {
    const complex pa = cb * ua[k], pb = cb * ub[k];
    sink.energy(k, -std::real(pa * ub[k]));
    if (pa == complex(0.0, 0.0) && pb == complex(0.0, 0.0)) continue;
    double *g = sink.widths(k);
    for (int d = 0; d < nch; d++) g[d] += std::real(pa * yb[d] + pb * ya[d]);
  }
}

/// The cotangents of K and N from resolve_formal_cotangents, as the two matrices the structural adjoint contracts:
/// level_matrix_adjoint(X, E) with the channel kernel set aside contracts X with dK - E dN, so
/// <Kbar, dK> + <Nbar, dN> = LMA(X1, 0) + LMA(X2, 1) with X1 = Kbar + Nbar, X2 = -Nbar.
struct FormalResolve {
  LevelMatrixSolver::matrix_c X1, X2;
};

/*!
 * Formal cotangents (ebar_k, Gammabar_kc) -> cotangents of K and N.  With K V = N V diag(e), V^T N V = 1 and
 * Gamma = V^T G, a variation gives
 *   de_l = v_l^T (dK - e_l dN) v_l,   dV = V Omega,
 *   Omega_kl = v_k^T (dK - e_l dN) v_l / (e_l - e_k) (k != l),   Omega_ll = -v_l^T dN v_l / 2,
 * so with H = V^T Vbar = Gamma Gammabar^T (Vbar = G Gammabar^T)
 *   Kbar = V Phi V^T,   Phi_ll = ebar_l,   Phi_kl = H_kl / (e_l - e_k),
 *   Nbar = V Psi V^T,   Psi_ll = -e_l ebar_l - H_ll / 2,   Psi_kl = -e_l H_kl / (e_l - e_k),
 * taken symmetric (K and N are).  gbar[k][c] over all channels c of the block; the result is out.X1, out.X2.
 * The cotangent of G, V Gammabar, is add_level_widths_from_formal.  Formal energies closer than the caller's limit
 * (FormalBlock::adjointOK) must not come here.
 */
inline void resolve_formal_cotangents(const FormalBlock &f, const std::vector<double> &ebar,
                                      const std::vector<std::vector<double>> &gbar, FormalResolve &out) {
  typedef std::complex<double> complex;
  const int n = f.n, m = (int)f.act.size();
  std::vector<double> e(n);
  for (int k = 0; k < n; k++) e[k] = std::real(f.e[k]);
  // H = Gamma Gammabar^T over the coupled channels.
  std::vector<double> H((std::size_t)n * n, 0.0);
  for (int k = 0; k < n; k++)
    for (int l = 0; l < n; l++) {
      double s = 0.0;
      for (int a = 0; a < m; a++) s += f.gamma[(std::size_t)k * m + a] * gbar[l][f.act[a]];
      H[(std::size_t)k * n + l] = s;
    }
  std::vector<double> Phi((std::size_t)n * n, 0.0), Psi((std::size_t)n * n, 0.0);
  for (int l = 0; l < n; l++) {
    Phi[(std::size_t)l * n + l] = ebar[l];
    Psi[(std::size_t)l * n + l] = -e[l] * ebar[l] - 0.5 * H[(std::size_t)l * n + l];
    for (int k = 0; k < n; k++) {
      if (k == l) continue;
      const double num = H[(std::size_t)k * n + l] - H[(std::size_t)l * n + k];
      const double den = e[l] - e[k];
      const double ratio = (num == 0.0 || den == 0.0) ? 0.0 : num / (2.0 * den);
      Phi[(std::size_t)k * n + l] = ratio;
      Psi[(std::size_t)k * n + l] = -e[l] * ratio - 0.5 * H[(std::size_t)l * n + k];
    }
  }
  // Psi as built is symmetric only up to round-off; use its symmetric part.
  for (int k = 0; k < n; k++)
    for (int l = k + 1; l < n; l++) {
      const double s = 0.5 * (Psi[(std::size_t)k * n + l] + Psi[(std::size_t)l * n + k]);
      Psi[(std::size_t)k * n + l] = Psi[(std::size_t)l * n + k] = s;
    }
  // V X V^T for X = Phi, Psi.
  auto congruence = [&](const std::vector<double> &X) {
    std::vector<double> T((std::size_t)n * n, 0.0), res((std::size_t)n * n, 0.0);
    for (int r = 0; r < n; r++)
      for (int l = 0; l < n; l++) {
        double s = 0.0;
        for (int k = 0; k < n; k++) s += f.V[(std::size_t)r * n + k] * X[(std::size_t)k * n + l];
        T[(std::size_t)r * n + l] = s;
      }
    for (int r = 0; r < n; r++)
      for (int c = 0; c < n; c++) {
        double s = 0.0;
        for (int l = 0; l < n; l++) s += T[(std::size_t)r * n + l] * f.V[(std::size_t)c * n + l];
        res[(std::size_t)r * n + c] = s;
      }
    return res;
  };
  const std::vector<double> Kbar = congruence(Phi), Nbar = congruence(Psi);
  out.X1.assign(n, std::vector<complex>(n));
  out.X2.assign(n, std::vector<complex>(n));
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) {
      out.X1[r][c] = Kbar[(std::size_t)r * n + c] + Nbar[(std::size_t)r * n + c];
      out.X2[r][c] = -Nbar[(std::size_t)r * n + c];
    }
}

/// The cotangent of the level widths from the formal ones: Gbar = V Gammabar (gbar[k][c], all nch channels), added
/// through the sink level by level.
template <class Sink>
void add_level_widths_from_formal(const FormalBlock &f, const std::vector<std::vector<double>> &gbar, Sink &&sink) {
  const int n = f.n, nch = f.nch;
  for (int r = 0; r < n; r++)
    for (int c = 0; c < nch; c++) {
      double s = 0.0;
      for (int k = 0; k < n; k++) s += f.V[(std::size_t)r * n + k] * gbar[k][c];
      sink.width(r, c, s);
    }
}

}  // namespace channel_space
}  // namespace physics
}  // namespace azure2

#endif
