#ifndef AZURE2_PHYSICS_COLLISION_MATRIX_ADJOINT_H
#define AZURE2_PHYSICS_COLLISION_MATRIX_ADJOINT_H

// The reverse mode of CollisionMatrix.h for one J^pi block at one energy (docs/upgrade/physics_functions.md,
// stage 5): from the cotangent of a T-matrix element, of a width bilinear form g_c^T A g_c', of the level matrix A
// or of its inverse M = A^-1 to the cotangents of the level energies E_l and the reduced-width amplitudes gamma_lc.
//
// Convention, everywhere: for a complex intermediate z and the real scalar L (a chi-squared, a model value) the
// cotangent is zbar = dL/dRe(z) + i dL/dIm(z), so that dL = Re(conj(zbar) dz).  Then
//   z = c a (c a complex constant)   gives   abar += conj(c) zbar,
//   a real parameter p               gets    dL/dp += Re(conj(zbar) dz/dp),
//   A = M^-1                         gives   Mbar = -A^H Abar A^H.
//
// Nothing here knows the engine: a caller hands over the arrays of the block and a *sink*, an object with
//   void energy(int level, double value);              // dL/dE_level += value
//   void width(int level, int channel, double value);  // dL/dgamma_{level, channel} += value
// (levels and channels of the block from 0).  The functions add through it term by term, in the order of the
// engine's routines they were moved from (AMatrixFunc::LevelMatrixAdjoint, BilinearAdjoint, LevelAdjoint): the
// gradients of every project depend on that order to the last digit.  LevelCotangents is a sink over two arrays.
// Scratch is the caller's (a struct it keeps between calls); nothing is kept here.

#include "azure2/physics/CollisionMatrix.h"
#include "azure2/physics/LevelPairShift.h"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {

/// One term beta g_ch^T A g_chp of the width bilinear forms of a block: beta is the cotangent of the form.  The
/// channels count from 1 (the engine's pathway table).
struct BilinearTerm {
  int ch;
  int chp;
  std::complex<double> beta;
};

/// A sink over two arrays: e[l] and gamma[l * nc + c].
struct LevelCotangents {
  int n = 0, nc = 0;
  std::vector<double> e, gamma;
  LevelCotangents(int levels, int channels) : n(levels), nc(channels), e(levels, 0.0), gamma((std::size_t)levels * channels, 0.0) {}
  void energy(int l, double v) { e[l] += v; }
  void width(int l, int c, double v) { gamma[(std::size_t)l * nc + c] += v; }
};

/// What level_matrix_adjoint reads of one block: n levels, nc channels; the pointers are borrowed.  Arrays over
/// (level, channel) are row-major, l * nc + c.
struct LevelMatrixAdjointTables {
  int n = 0;
  int nc = 0;
  const double *E = nullptr;            ///< level energies E_l (MeV), n
  const double *gamma = nullptr;        ///< reduced-width amplitudes (MeV^1/2), n x nc
  const double *shift = nullptr;        ///< S_c(E_l), n x nc
  const double *shift_gamma = nullptr;  ///< S_c(E_l) gamma_lc, n x nc (standard and Brune)
  const double *inv_gap = nullptr;      ///< 1 / (E_r - E_c) off the diagonal (inverse_level_gap), n x n (the same)
  const double *shift_d1 = nullptr;     ///< dS_c/dE at E_l as the levels hold it, n x nc (Park; coincident levels)
  const double *shift_d2 = nullptr;     ///< d2S_c/dE2 at E_l, n x nc (the same)
  const double *park_norm = nullptr;    ///< J_ll of Park (2021) Eq. (28), n (Park)
  /// dS_c/dE at E_l of a Brune gradient pass, [l][c]; null: the level energies get no cotangent through S(E_l).
  const std::vector<std::vector<double>> *brune_shift_d1 = nullptr;
  const std::complex<double> *kernel = nullptr;  ///< Lambda_c: L_c(E) - B_c, or what the caller differentiates with
  const double *boundary = nullptr;     ///< B_c, added to the kernel of a shift channel; null: the kernel has it
  const char *shift_channel = nullptr;  ///< channels whose kernel carries the level's shift (Brune, Park: particle), nc
  bool brune = false;
  bool park = false;  ///< takes precedence over brune
};

/// Scratch of level_matrix_adjoint (the caller keeps it between calls).
struct LevelMatrixAdjointScratch {
  std::vector<double> re_sym, im_sym, dE, acc;
  std::vector<char> no_channel;
  std::vector<std::complex<double>> kernel;
};

namespace adjoint_detail {
inline LevelPairShift<double> pair_shift(const LevelMatrixAdjointTables &t, int la, int lb, int ch) {
  auto at = [&](int l) {
    const std::size_t i = (std::size_t)l * t.nc + ch;
    return ShiftAtLevel<double>{t.E[l], t.shift[i], t.shift_d1[i], t.shift_d2[i]};
  };
  return level_pair_shift(at(la), at(lb));
}
}  // namespace adjoint_detail

/*!
 * The structural adjoint of the inverse level matrix M(E) of one block (fill_level_matrix): given its cotangent
 * Mbar (n x n), adds the cotangents of the level energies and the reduced-width amplitudes of
 *
 *   standard:  M_rc = (E_r - E) delta_rc - sum_d g_rd g_cd Lambda_d,
 *   Brune:     the same with Lambda_d + B_d - D_d[r][c] in the shift channels, D = S_r on the diagonal and the chord
 *              of S between the two levels at E off it (LevelPairShift),
 *   Park:      M_rc = (E_r - E) J_rc - sum_d g_rd g_cd (Lambda_d + B_d - S_rd), J of Park (2021) Eqs. (9), (28).
 *
 * Lambda_d is `kernel` as the caller gives it (the Lo-matrix element at E; zero, or a constant, where the caller
 * contracts with dK - E dN only: AMatrixFunc::ResolveFormalCotangents).  Moved from AMatrixFunc::LevelMatrixAdjoint
 * with its order of operations.
 */
template <class Sink>
void level_matrix_adjoint(const LevelMatrixAdjointTables &t, const LevelMatrixSolver::matrix_c &Mbar, double inEnergy,
                          LevelMatrixAdjointScratch &scratch, Sink &&sink) {
  typedef std::complex<double> complex;
  const int n = (int)Mbar.size();
  if (n == 0) return;
  const int nch = t.nc;
  // The (E_r - E) delta term: dM_aa/dE_a = 1.
  for (int a = 0; a < n; a++) sink.energy(a, std::real(Mbar[a][a]));

  if (t.park) {
    // M is real-coefficient in E_l, so only Re(Mbar) enters.  The delta term above gave dM_rr/dE_r = 1; here it is
    // J_rr, so J_rr - 1 is added, plus the derivatives of J and of S_rd.
    for (int r = 0; r < n; r++) {
      const double Er = t.E[r];
      const double *gr = &t.gamma[(std::size_t)r * nch];
      for (int c = 0; c < n; c++) {
        const double reM = std::real(Mbar[r][c]);
        if (reM == 0.0) continue;
        if (r == c) {
          double dJ = 0.0, sumS1 = 0.0;
          for (int d = 0; d < nch; d++) {
            if (!t.shift_channel[d]) continue;
            const double g2 = gr[d] * gr[d];
            if (g2 == 0.0) continue;
            dJ -= g2 * t.shift_d2[(std::size_t)r * nch + d];
            sumS1 += g2 * t.shift_d1[(std::size_t)r * nch + d];
          }
          sink.energy(r, reM * ((t.park_norm[r] - 1.0) + (Er - inEnergy) * dJ + sumS1));
        } else {
          const double *gcl = &t.gamma[(std::size_t)c * nch];
          double Jrc = 0.0, dJdEr = 0.0, dJdEc = 0.0, sumS1 = 0.0;
          for (int d = 0; d < nch; d++) {
            if (!t.shift_channel[d]) continue;
            const double g = gr[d] * gcl[d];
            if (g == 0.0) continue;
            const LevelPairShift<double> chord = adjoint_detail::pair_shift(t, r, c, d);
            const double q = chord.slope() / chord.gap;
            Jrc -= g * q;
            dJdEr -= g * chord.slope_by_Ea(q) / chord.gap;
            dJdEc -= g * chord.slope_by_Eb(q) / chord.gap;
            sumS1 += g * chord.a.dS;
          }
          sink.energy(r, reM * (Jrc + (Er - inEnergy) * dJdEr + sumS1));
          sink.energy(c, reM * (Er - inEnergy) * dJdEc);
        }
      }
    }
    // The widths.  For M_ac and M_ca the kernel is row-indexed, K_d[a][c] = Lambda_d + B_d - S_ad, and J gives
    //   dJ_aa/dg_ad = -2 g_ad S'_ad,   dJ_ac/dg_ad = dJ_ca/dg_ad = -g_cd (S_ad - S_cd) / (E_a - E_c).
    for (int d = 0; d < nch; d++) {
      const bool pd = t.shift_channel[d];
      const complex Ld = t.kernel[d];
      const double Bd = (pd && t.boundary) ? t.boundary[d] : 0.0;
      for (int a = 0; a < n; a++) {
        const double Ea = t.E[a];
        const double ga = t.gamma[(std::size_t)a * nch + d];
        const double Sa = pd ? t.shift[(std::size_t)a * nch + d] : 0.0;
        double grad = 0.0;
        for (int c = 0; c < n; c++) {
          const double gc = t.gamma[(std::size_t)c * nch + d];
          if (c == a) {
            const complex Kaa = pd ? (Ld + Bd - Sa) : Ld;
            grad += -2.0 * ga * std::real(std::conj(Mbar[a][a]) * Kaa);
            if (pd) grad += std::real(Mbar[a][a]) * (Ea - inEnergy) * (-2.0 * ga * t.shift_d1[(std::size_t)a * nch + d]);
          } else {
            if (gc == 0.0) continue;
            const double Sc = pd ? t.shift[(std::size_t)c * nch + d] : 0.0;
            const complex Kac = pd ? (Ld + Bd - Sa) : Ld;
            const complex Kca = pd ? (Ld + Bd - Sc) : Ld;
            grad += -gc * std::real(std::conj(Mbar[a][c]) * Kac + std::conj(Mbar[c][a]) * Kca);
            if (pd) {
              const double Ec = t.E[c];
              const LevelPairShift<double> chord = adjoint_detail::pair_shift(t, a, c, d);
              const double delta = chord.slope() / chord.gap;
              grad += -gc * delta * (std::real(Mbar[a][c]) * (Ea - inEnergy) + std::real(Mbar[c][a]) * (Ec - inEnergy));
            }
          }
        }
        sink.width(a, d, grad);
      }
    }
    return;
  }

  // Standard and Brune.  The loops are n^2 per channel over flat arrays.
  const double *E = t.E, *G = t.gamma, *S = t.shift, *inv = t.inv_gap;
  std::vector<double> &reSym = scratch.re_sym, &imSym = scratch.im_sym;
  // The shift channels of Brune; none for the standard matrix, whatever the caller marked.
  const char *bruneP = t.shift_channel;
  if (!t.brune) {
    scratch.no_channel.assign(nch, 0);
    bruneP = scratch.no_channel.data();
  }
  bool anyBruneP = false;
  for (int ch = 0; ch < nch; ch++)
    if (bruneP[ch]) anyBruneP = true;

  // Brune: the level energies enter through the kernel too, K_ch[r][c] = Lambda + B - D_ch[r][c] with D = S_r on the
  // diagonal and Q_rc = (S_r (E - E_c) - S_c (E - E_r)) / (E_r - E_c) off it, so
  //   dL/dE_nu += sum_rc Re(Mbar_rc) sum_ch g_r g_c dD/dE_nu;   dD/dE couples to E_r and E_c only.
  if (t.brune && t.brune_shift_d1 && anyBruneP) {
    const std::vector<std::vector<double>> &sd = *t.brune_shift_d1;
    std::vector<double> &dE = scratch.dE;
    dE.assign(n, 0.0);
    for (int ch = 0; ch < nch; ch++) {
      if (!bruneP[ch]) continue;
      for (int r = 0; r < n; r++) {
        const double gr = G[(std::size_t)r * nch + ch];
        if (gr == 0.0) continue;
        const double Sr = S[(std::size_t)r * nch + ch];
        const double dSr = sd[r][ch];
        const complex *mrow = &Mbar[r][0];
        // D = S_r;  dD/dE_r = S'_r.
        dE[r] += std::real(mrow[r]) * gr * gr * dSr;
        for (int c = 0; c < n; c++) {
          if (c == r) continue;
          const double g = gr * G[(std::size_t)c * nch + ch];
          if (g == 0.0) continue;
          const double iv = inv[(std::size_t)r * n + c];
          const LevelPairShift<double> chord = level_pair_shift<double>(
              {E[r], Sr, dSr, t.shift_d2[(std::size_t)r * nch + ch]},
              {E[c], S[(std::size_t)c * nch + ch], sd[c][ch], t.shift_d2[(std::size_t)c * nch + ch]});
          const double Q = chord.chord(inEnergy) * iv;
          const double dQdEr = chord.chord_by_Ea(inEnergy, Q) * iv;
          const double dQdEc = chord.chord_by_Eb(inEnergy, Q) * iv;
          const double w = std::real(mrow[c]) * g;
          dE[r] += w * dQdEr;
          dE[c] += w * dQdEc;
        }
      }
    }
    for (int r = 0; r < n; r++) sink.energy(r, dE[r]);
  }

  // The widths:  dL/dg_ad += -sum_c Re(conj(Mbar_ac + Mbar_ca) K_d[a][c]) g_cd.  With sym = Mbar + Mbar^T and D real,
  //   -Re(conj(sym) K) = -Re(sym) Re(Lambda + B) - Im(sym) Im(Lambda) + Re(sym) D,
  // and the D term splits into
  //   S_a sum_{c != a} Re(sym_ac) (E - E_c) / (E_a - E_c) g_c
  //   - (E - E_a) sum_{c != a} Re(sym_ac) / (E_a - E_c) S_c g_c  + Re(sym_aa) S_a g_a,
  // so every channel is a few row-times-matrix products.
  reSym.resize((std::size_t)n * n);
  imSym.resize((std::size_t)n * n);
  for (int a = 0; a < n; a++)
    for (int c = 0; c < n; c++) {
      const complex sym = Mbar[a][c] + Mbar[c][a];
      reSym[(std::size_t)a * n + c] = std::real(sym);
      imSym[(std::size_t)a * n + c] = std::imag(sym);
    }
  const double *SG = t.shift_gamma;
  std::vector<double> &acc = scratch.acc;
  acc.resize((std::size_t)4 * nch);
  double *sRe = &acc[0], *sIm = &acc[nch], *t1 = &acc[2 * (std::size_t)nch], *t2 = &acc[3 * (std::size_t)nch];
  std::vector<complex> &Lconst = scratch.kernel;
  Lconst.resize(nch);
  for (int d = 0; d < nch; d++) {
    complex Ld = t.kernel[d];
    if (t.boundary && bruneP[d]) Ld += t.boundary[d];
    Lconst[d] = Ld;
  }
  for (int a = 0; a < n; a++) {
    std::fill(acc.begin(), acc.end(), 0.0);
    const double *rs = &reSym[(std::size_t)a * n], *is = &imSym[(std::size_t)a * n];
    for (int c = 0; c < n; c++) {
      const double *gc = &G[(std::size_t)c * nch];
      const double wr = rs[c], wi = is[c];
      for (int d = 0; d < nch; d++) {
        sRe[d] += wr * gc[d];
        sIm[d] += wi * gc[d];
      }
      if (anyBruneP && c != a) {
        const bool coincident = levels_coincide(E[a], E[c]);  // 1/gap is 1, D = S_a + S'_m (E - E_a)
        const double iv = inv[(std::size_t)a * n + c];
        const double w1 = coincident ? wr : wr * (inEnergy - E[c]) * iv, w2 = wr * iv;
        const double *sgc = &SG[(std::size_t)c * nch];
        for (int d = 0; d < nch; d++) {
          t1[d] += w1 * gc[d];
          t2[d] += coincident ? -w2 * adjoint_detail::pair_shift(t, a, c, d).slope() * gc[d] : w2 * sgc[d];
        }
      }
    }
    for (int d = 0; d < nch; d++) {
      double g = -(std::real(Lconst[d]) * sRe[d] + std::imag(Lconst[d]) * sIm[d]);
      if (bruneP[d]) {
        const double Sa = S[(std::size_t)a * nch + d];
        g += Sa * t1[d] - (inEnergy - E[a]) * t2[d] + rs[a] * Sa * G[(std::size_t)a * nch + d];
      }
      sink.width(a, d, g);
    }
  }
}

/// The inverse lemma: with A = M^-1 and the cotangent Abar of A, Mbar = -A^H Abar A^H (n x n, two products).
inline LevelMatrixSolver::matrix_c inverse_adjoint(const LevelMatrixSolver::matrix_c &A,
                                                   const LevelMatrixSolver::matrix_c &Ab) {
  typedef std::complex<double> complex;
  typedef std::vector<complex> vector_c;
  const int n = (int)A.size();
  // Conjugate transpose A^H.
  LevelMatrixSolver::matrix_c Ah(n, vector_c(n, complex(0.0, 0.0)));
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) Ah[r][c] = std::conj(A[c][r]);
  LevelMatrixSolver::matrix_c tmp(n, vector_c(n, complex(0.0, 0.0)));  // Abar A^H
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) {
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += Ab[r][k] * Ah[k][c];
      tmp[r][c] = s;
    }
  LevelMatrixSolver::matrix_c Mbar(n, vector_c(n, complex(0.0, 0.0)));
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++) {
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += Ah[r][k] * tmp[k][c];
      Mbar[r][c] = -s;
    }
  return Mbar;
}

/*!
 * The widths as they stand in b = g_ch^T A g_chp: db/dg_{l,ch} = (A g_chp)_l and db/dg_{l,chp} = (A^T g_ch)_l, and A
 * is complex symmetric, so both are the solves u_c = A g_c (u_in for ch, u_out for chp, n values each).
 */
template <class Sink>
void add_explicit_widths(const BilinearTerm &term, const std::complex<double> *u_in, const std::complex<double> *u_out,
                         int n, Sink &&sink) {
  const std::complex<double> cb = std::conj(term.beta);
  for (int l = 0; l < n; l++) {
    sink.width(l, term.ch - 1, std::real(cb * u_out[l]));
    sink.width(l, term.chp - 1, std::real(cb * u_in[l]));
  }
}

/// The cotangent of M from one term: Abar = beta g_ch g_chp^T gives Mbar = -A^H Abar A^H = -beta conj(u_in) conj(u_out)^T.
inline void add_rank_one_cotangent(std::complex<double> beta, const std::complex<double> *u_in,
                                   const std::complex<double> *u_out, int n, LevelMatrixSolver::matrix_c &Mbar) {
  for (int r = 0; r < n; r++) {
    const std::complex<double> a = -beta * std::conj(u_in[r]);
    std::vector<std::complex<double>> &row = Mbar[r];
    for (int c = 0; c < n; c++) row[c] += a * std::conj(u_out[c]);
  }
}

/*!
 * One term on a pass that sums a scalar over many energies: Mbar = -beta conj(u_in) conj(u_out)^T enters the
 * structural adjoint through K - E N linearly in E, so only Re(Mbar) and E Re(Mbar) are summed here (into K and N,
 * n x n row-major, which the caller resolves once per pass with level_matrix_adjoint), and the kernel
 * G Lambda(E) G^T is differentiated now, through the rank-one term:
 *   gammabar_{l d} -= Re(conj(Lambda_d) ((Mbar + Mbar^T) g_d)_l).
 * G: the widths, n x nc; kernel: Lambda_d, nc; ga, gb: scratch of nc values each.
 */
template <class Sink>
void deferred_bilinear_adjoint(const BilinearTerm &term, const std::complex<double> *uEn, const std::complex<double> *uEx,
                               const double *G, const std::complex<double> *lam, double inEnergy, int n, int nch,
                               double *K, double *N, std::vector<std::complex<double>> &ga,
                               std::vector<std::complex<double>> &gb, Sink &&sink) {
  typedef std::complex<double> complex;
  // conj(u)^T g_d for both solves.
  std::fill(ga.begin(), ga.end(), complex(0.0, 0.0));
  std::fill(gb.begin(), gb.end(), complex(0.0, 0.0));
  for (int la = 0; la < n; la++) {
    const complex ca = std::conj(uEn[la]), cbx = std::conj(uEx[la]);
    const double *g = &G[(std::size_t)la * nch];
    for (int d = 0; d < nch; d++) {
      ga[d] += ca * g[d];
      gb[d] += cbx * g[d];
    }
  }
  for (int la = 0; la < n; la++) {
    const complex a = -term.beta * std::conj(uEn[la]), b = -term.beta * std::conj(uEx[la]);
    // (Mbar g_d)_l = a_l (conj(uEx) . g_d);  (Mbar^T g_d)_l = b_l (conj(uEn) . g_d).
    for (int d = 0; d < nch; d++) sink.width(la, d, -std::real(std::conj(lam[d]) * (a * gb[d] + b * ga[d])));
    double *kr = &K[(std::size_t)la * n], *nr = &N[(std::size_t)la * n];
    for (int c = 0; c < n; c++) {
      const double v = std::real(a * std::conj(uEx[c]));
      kr[c] += v;
      nr[c] += inEnergy * v;
    }
  }
}

/// The derivative of K and N of a block's level matrix M = K - E N - G Lambda G^T by one level parameter (the energy
/// of `level`, ch = -1, or its width in channel ch): dK = e_l k^T + k e_l^T, dN = e_l q^T + q e_l^T over the levels.
struct LevelParameterDerivative {
  int level;  ///< 0-based
  int ch;     ///< 0-based channel, or -1 for the level energy
  bool zero;  ///< k = q = 0 (only the kernel and the explicit terms remain)
  std::vector<double> k, q;
};

/// gu[s * nc + d] = g_d . u_s over the levels, for the solves u_s (nu of them, n values each); G is n x nc.
inline void width_solve_products(const double *G, const std::vector<const std::vector<std::complex<double>> *> &u, int n,
                                 int nch, std::vector<std::complex<double>> &gu) {
  typedef std::complex<double> complex;
  const int nu = (int)u.size();
  gu.assign((std::size_t)nu * nch, complex(0.0, 0.0));
  for (int i = 0; i < n; i++) {
    const double *g = &G[(std::size_t)i * nch];
    for (int s = 0; s < nu; s++) {
      const complex ui = (*u[s])[i];
      complex *row = &gu[(std::size_t)s * nch];
      for (int d = 0; d < nch; d++) row[d] += g[d] * ui;
    }
  }
}

/*!
 * The derivatives of sum_t beta_t g_ch^T A g_chp by the level parameters of `entries`, one at a time: with
 * conj(Mbar) = -sum_t conj(beta_t) u_ch u_chp^T (u = A g) a parameter's derivative is
 *   -Re(sum_t conj(beta_t) u_ch^T dM u_chp),   dM = dK - E dN - (dG Lambda G^T + G Lambda dG^T),
 * and for dK = e_l k^T + k e_l^T:  u^T dK v = u_l (k . v) + (k . u) v_l.  A width also enters the form as it
 * stands: conj(beta) u_chp and conj(beta) u_ch.  O(n) per parameter.  slot[c]: the place of channel c (0-based)
 * among the solves u; gu: width_solve_products; kernel: Lambda_d; ku: scratch.  add(level, ch, value) receives each
 * derivative (ch = -1: the level energy).  Moved from AMatrixFunc::LevelAdjoint.
 */
template <class Add>
void level_parameter_adjoint(const std::vector<LevelParameterDerivative> &entries, const std::vector<BilinearTerm> &terms,
                             const std::vector<int> &slot, const std::vector<const std::vector<std::complex<double>> *> &u,
                             const std::vector<std::complex<double>> &gu, const std::complex<double> *lam, double E,
                             int n, int nch, std::vector<std::complex<double>> &ku, Add &&add) {
  typedef std::complex<double> complex;
  const int nu = (int)u.size();
  ku.resize(nu);
  for (const LevelParameterDerivative &e : entries) {
    const int l = e.level;
    if (!e.zero) {
      std::fill(ku.begin(), ku.end(), complex(0.0, 0.0));
      for (int i = 0; i < n; i++) {
        const double kt = e.k[i] - E * e.q[i];
        if (kt == 0.0) continue;
        for (int s = 0; s < nu; s++) ku[s] += kt * (*u[s])[i];
      }
    }
    double v = 0.0;
    for (const BilinearTerm &t : terms) {
      const int sa = slot[t.ch - 1], sb = slot[t.chp - 1];
      const complex ua = (*u[sa])[l], ub = (*u[sb])[l];
      complex dm(0.0, 0.0);
      if (!e.zero) dm += ua * ku[sb] + ku[sa] * ub;
      if (e.ch >= 0) dm -= lam[e.ch] * (ua * gu[(std::size_t)sb * nch + e.ch] + gu[(std::size_t)sa * nch + e.ch] * ub);
      const complex cb = std::conj(t.beta);
      v -= std::real(cb * dm);
      if (e.ch == t.ch - 1) v += std::real(cb * ub);
      if (e.ch == t.chp - 1) v += std::real(cb * ua);
    }
    add(l, e.ch, v);
  }
}

}  // namespace physics
}  // namespace azure2

#endif
