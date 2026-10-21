#include "ChannelSpace.h"
#include "AChannel.h"
#include "ALevel.h"
#include "Config.h"
#include "JGroup.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>

namespace ChannelSpace {

const double kMaxPoleRatio = 1.0e3;
const double kMinPivotRatio = 1.0e-10;
const double kMaxPoleError = 1.0e-11;
const double kEnergyUlps = 2.0;

bool Enabled() {
  static const bool on = [] {
    const char *e = std::getenv("AZURE_LEVEL_SPACE");
    return !(e && *e && std::string(e) != "0");
  }();
  return on;
}

std::shared_ptr<FormalLevels> Cached(const JGroup *jg) {
  // Keyed by address only: Update() checks the levels themselves.  Bounded,
  // for processes that build many compound nuclei (pyazr sessions); entries in
  // use stay alive through their shared pointers.
  thread_local std::unordered_map<const JGroup *, std::shared_ptr<FormalLevels>> cache;
  auto it = cache.find(jg);
  if (it != cache.end()) return it->second;
  if (cache.size() >= 1024) cache.clear();
  std::shared_ptr<FormalLevels> f = std::make_shared<FormalLevels>();
  cache.emplace(jg, f);
  return f;
}

double LUFactor(complex *a, int *piv, int m) {
  double pmin = 0.0, pmax = 0.0;
  for (int k = 0; k < m; k++) {
    int p = k;
    double best = std::norm(a[(size_t)k * m + k]);
    for (int i = k + 1; i < m; i++) {
      const double v = std::norm(a[(size_t)i * m + k]);
      if (v > best) {
        best = v;
        p = i;
      }
    }
    piv[k] = p;
    if (p != k)
      for (int j = 0; j < m; j++) std::swap(a[(size_t)k * m + j], a[(size_t)p * m + j]);
    const complex pivot = a[(size_t)k * m + k];
    const double mod = std::abs(pivot);
    if (!(mod > 0.0) || !std::isfinite(mod)) return 0.0;
    if (k == 0) pmin = pmax = mod;
    pmin = std::min(pmin, mod);
    pmax = std::max(pmax, mod);
    const complex inv = 1.0 / pivot;
    const complex *rowK = &a[(size_t)k * m];
    for (int i = k + 1; i < m; i++) {
      complex *rowI = &a[(size_t)i * m];
      const complex f = (rowI[k] *= inv);
      if (f == complex(0.0, 0.0)) continue;
      for (int j = k + 1; j < m; j++) rowI[j] -= f * rowK[j];
    }
  }
  return m > 0 ? pmin / pmax : 1.0;
}

void LUSolve(const complex *lu, const int *piv, int m, complex *b) {
  for (int k = 0; k < m; k++)
    if (piv[k] != k) std::swap(b[k], b[piv[k]]);
  for (int i = 1; i < m; i++) {
    const complex *row = &lu[(size_t)i * m];
    complex s = b[i];
    for (int j = 0; j < i; j++) s -= row[j] * b[j];
    b[i] = s;
  }
  for (int i = m - 1; i >= 0; i--) {
    const complex *row = &lu[(size_t)i * m];
    complex s = b[i];
    for (int j = i + 1; j < m; j++) s -= row[j] * b[j];
    b[i] = s / row[i];
  }
}

double LUFactorReal(double *a, int *piv, int m) {
  double pmin = 0.0, pmax = 0.0;
  for (int k = 0; k < m; k++) {
    int p = k;
    double best = std::fabs(a[(size_t)k * m + k]);
    for (int i = k + 1; i < m; i++) {
      const double v = std::fabs(a[(size_t)i * m + k]);
      if (v > best) {
        best = v;
        p = i;
      }
    }
    piv[k] = p;
    if (p != k)
      for (int j = 0; j < m; j++) std::swap(a[(size_t)k * m + j], a[(size_t)p * m + j]);
    const double pivot = a[(size_t)k * m + k];
    const double mod = std::fabs(pivot);
    if (!(mod > 0.0) || !std::isfinite(mod)) return 0.0;
    if (k == 0) pmin = pmax = mod;
    pmin = std::min(pmin, mod);
    pmax = std::max(pmax, mod);
    const double inv = 1.0 / pivot;
    const double *rowK = &a[(size_t)k * m];
    for (int i = k + 1; i < m; i++) {
      double *rowI = &a[(size_t)i * m];
      const double f = (rowI[k] *= inv);
      if (f == 0.0) continue;
      for (int j = k + 1; j < m; j++) rowI[j] -= f * rowK[j];
    }
  }
  return m > 0 ? pmin / pmax : 1.0;
}

void LUSolveReal(const double *lu, const int *piv, int m, double *b, int nrhs) {
  // b is m x nrhs row-major.
  for (int k = 0; k < m; k++)
    if (piv[k] != k)
      for (int r = 0; r < nrhs; r++) std::swap(b[(size_t)k * nrhs + r], b[(size_t)piv[k] * nrhs + r]);
  for (int i = 1; i < m; i++) {
    const double *row = &lu[(size_t)i * m];
    double *bi = &b[(size_t)i * nrhs];
    for (int j = 0; j < i; j++) {
      const double l = row[j];
      if (l == 0.0) continue;
      const double *bj = &b[(size_t)j * nrhs];
      for (int r = 0; r < nrhs; r++) bi[r] -= l * bj[r];
    }
  }
  for (int i = m - 1; i >= 0; i--) {
    const double *row = &lu[(size_t)i * m];
    double *bi = &b[(size_t)i * nrhs];
    for (int j = i + 1; j < m; j++) {
      const double u = row[j];
      if (u == 0.0) continue;
      const double *bj = &b[(size_t)j * nrhs];
      for (int r = 0; r < nrhs; r++) bi[r] -= u * bj[r];
    }
    const double inv = 1.0 / row[i];
    for (int r = 0; r < nrhs; r++) bi[r] *= inv;
  }
}

namespace {

/// Smallest Cholesky pivot of N accepted (N has unit diagonal for Brune).
const double kMinOverlapPivot = 1.0e-4;

bool Debug() {
  static const bool on = std::getenv("AZURE_CHANNEL_DEBUG") != nullptr;
  return on;
}

/// Width with the cutoff the level-matrix fill applies (AMatrixFunc::FillMatrices).
inline double Cut(double g) { return std::fabs(g) < 1.0e-12 ? 0.0 : g; }

typedef long double ldouble;

/*!
 * Eigen-decomposition of the symmetric n x n matrix C (row-major, destroyed)
 * by cyclic Jacobi rotations: C = W diag(e) W^T.  Slower than Householder
 * reduction but accurate to the working precision in every eigenvector
 * however close the eigenvalues are.
 */
bool Jacobi(std::vector<ldouble> &C, int n, std::vector<ldouble> &e, std::vector<ldouble> &W) {
  W.assign((size_t)n * n, 0.0L);
  for (int i = 0; i < n; i++) W[(size_t)i * n + i] = 1.0L;
  const ldouble eps = std::numeric_limits<ldouble>::epsilon();
  for (int sweep = 0; sweep < 60; sweep++) {
    ldouble off = 0.0L, diag = 0.0L;
    for (int p = 0; p < n; p++) {
      diag += C[(size_t)p * n + p] * C[(size_t)p * n + p];
      for (int q = p + 1; q < n; q++) off += C[(size_t)p * n + q] * C[(size_t)p * n + q];
    }
    if (!std::isfinite((double)off) || !std::isfinite((double)diag)) return false;
    if (off <= eps * eps * 1.0e-4L * diag || off == 0.0L) {
      e.resize(n);
      for (int i = 0; i < n; i++) e[i] = C[(size_t)i * n + i];
      return true;
    }
    for (int p = 0; p < n; p++)
      for (int q = p + 1; q < n; q++) {
        const ldouble apq = C[(size_t)p * n + q];
        if (apq == 0.0L) continue;
        const ldouble theta = (C[(size_t)q * n + q] - C[(size_t)p * n + p]) / (2.0L * apq);
        const ldouble t = (theta >= 0.0L ? 1.0L : -1.0L) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0L));
        const ldouble c = 1.0L / std::sqrt(t * t + 1.0L), sn = t * c;
        for (int k = 0; k < n; k++) {
          const ldouble akp = C[(size_t)k * n + p], akq = C[(size_t)k * n + q];
          C[(size_t)k * n + p] = c * akp - sn * akq;
          C[(size_t)k * n + q] = sn * akp + c * akq;
        }
        for (int k = 0; k < n; k++) {
          const ldouble apk = C[(size_t)p * n + k], aqk = C[(size_t)q * n + k];
          C[(size_t)p * n + k] = c * apk - sn * aqk;
          C[(size_t)q * n + k] = sn * apk + c * aqk;
        }
        C[(size_t)p * n + q] = C[(size_t)q * n + p] = 0.0L;
        for (int k = 0; k < n; k++) {
          const ldouble wkp = W[(size_t)k * n + p], wkq = W[(size_t)k * n + q];
          W[(size_t)k * n + p] = c * wkp - sn * wkq;
          W[(size_t)k * n + q] = sn * wkp + c * wkq;
        }
      }
  }
  return false;
}

/*!
 * Generalized symmetric eigenproblem K v = e N v, N positive definite, by
 * Cholesky reduction: N = L L^T, C = L^-1 K L^-T = W diag(e) W^T, V = L^-T W,
 * in extended precision: the formal energies and widths it gives are rounded
 * to double once, so the channel route sees them as accurately as the level
 * route sees the level parameters, however close two levels are.  K and N are
 * n x n row-major and are overwritten.  Returns false if N is not (safely)
 * positive definite.
 */
bool GeneralizedEigen(std::vector<ldouble> &K, std::vector<ldouble> &N, int n, std::vector<ldouble> &e,
                      std::vector<ldouble> &V, double &minPivot) {
  // Cholesky, lower triangle of N in place.  N = a^T a for the eigenvectors a
  // of Brune's transformation, so it is positive definite whenever the
  // parameters have a standard (real) equivalent; its smallest pivot measures
  // how far from that the parameters are, and how much the transformation
  // amplifies round-off.
  minPivot = 1.0e300;
  for (int j = 0; j < n; j++) {
    ldouble d = N[(size_t)j * n + j];
    for (int k = 0; k < j; k++) d -= N[(size_t)j * n + k] * N[(size_t)j * n + k];
    minPivot = std::min(minPivot, (double)d);
    if (!(d > kMinOverlapPivot) || !std::isfinite((double)d)) return false;
    d = std::sqrt(d);
    N[(size_t)j * n + j] = d;
    for (int i = j + 1; i < n; i++) {
      ldouble s = N[(size_t)i * n + j];
      for (int k = 0; k < j; k++) s -= N[(size_t)i * n + k] * N[(size_t)j * n + k];
      N[(size_t)i * n + j] = s / d;
    }
  }
  auto Lm = [&](int i, int j) { return N[(size_t)i * n + j]; };
  // X = L^-1 K (forward substitution, column by column), stored in K.
  for (int c = 0; c < n; c++)
    for (int i = 0; i < n; i++) {
      ldouble s = K[(size_t)i * n + c];
      for (int k = 0; k < i; k++) s -= Lm(i, k) * K[(size_t)k * n + c];
      K[(size_t)i * n + c] = s / Lm(i, i);
    }
  // C = X L^-T = (L^-1 X^T)^T: row r of C solves L y = (row r of X).
  std::vector<ldouble> C((size_t)n * n);
  for (int r = 0; r < n; r++)
    for (int i = 0; i < n; i++) {
      ldouble s = K[(size_t)r * n + i];
      for (int k = 0; k < i; k++) s -= Lm(i, k) * C[(size_t)r * n + k];
      C[(size_t)r * n + i] = s / Lm(i, i);
    }
  for (int r = 0; r < n; r++)
    for (int c = r + 1; c < n; c++) {
      const ldouble s = 0.5L * (C[(size_t)r * n + c] + C[(size_t)c * n + r]);
      C[(size_t)r * n + c] = C[(size_t)c * n + r] = s;
    }
  std::vector<ldouble> W;
  if (!Jacobi(C, n, e, W)) return false;
  // Ascending eigenvalues.
  std::vector<int> order(n);
  for (int k = 0; k < n; k++) order[k] = k;
  std::stable_sort(order.begin(), order.end(), [&e](int a, int b) { return e[a] < e[b]; });
  std::vector<ldouble> es(n);
  V.assign((size_t)n * n, 0.0L);
  for (int kk = 0; kk < n; kk++) {
    const int k = order[kk];
    es[kk] = e[k];
    // v = L^-T w (back substitution).
    for (int i = n - 1; i >= 0; i--) {
      ldouble s = W[(size_t)i * n + k];
      for (int j = i + 1; j < n; j++) s -= Lm(j, i) * V[(size_t)j * n + kk];
      V[(size_t)i * n + kk] = s / Lm(i, i);
    }
  }
  e.swap(es);
  return true;
}

}  // namespace

bool Update(JGroup *jg, const std::vector<int> &levels, const Config &configure, FormalLevels &f) {
  const bool park = (configure.paramMask & Config::USE_PARK_FORMALISM);
  const bool brune = !park && (configure.paramMask & Config::USE_BRUNE_FORMALISM);
  const bool rmc = (configure.paramMask & Config::USE_RMC_FORMALISM);
  const int mode = park ? 2 : brune ? 1 : 0;
  const int n = (int)levels.size();
  const int nch = jg->NumChannels();

  std::vector<char> particle(nch);
  for (int c = 0; c < nch; c++) particle[c] = (jg->GetChannel(c + 1)->GetRadType() == 'P');

  if (mode != 0) {
    // The Brune and Park representations cost an eigen-decomposition: rebuild
    // only when a parameter they are made of changed (ALevel::Stamp).
    std::vector<std::uint64_t> &key = f.scratch_;
    key.resize(4 + 2 * (size_t)n);
    key[0] = mode;
    key[1] = rmc;
    key[2] = n;
    key[3] = nch;
    for (int i = 0; i < n; i++) {
      ALevel *level = jg->GetLevel(levels[i]);
      key[4 + 2 * i] = reinterpret_cast<std::uintptr_t>(level);
      key[5 + 2 * i] = level->Stamp();
    }
    if (f.mode == mode && key == f.key) return false;
    f.key = key;
  }

  f.mode = mode;
  f.n = n;
  f.nch = nch;
  f.status = 0;
  f.lamShift.assign(nch, 0.0);
  f.gcut.assign((size_t)n * nch, 0.0);
  std::vector<double> E(n);
  for (int i = 0; i < n; i++) {
    ALevel *level = jg->GetLevel(levels[i]);
    E[i] = level->GetFitE();
    for (int c = 0; c < nch; c++) f.gcut[(size_t)i * nch + c] = Cut(level->GetFitGamma(c + 1));
  }
  f.act.clear();
  f.pos.assign(nch, -1);
  for (int c = 0; c < nch; c++) {
    bool any = false;
    for (int i = 0; i < n && !any; i++) any = (f.gcut[(size_t)i * nch + c] != 0.0);
    if (any) {
      f.pos[c] = (int)f.act.size();
      f.act.push_back(c);
    }
  }
  const int m = (int)f.act.size();
  f.usable = (n > 0 && m < n);
  if (!f.usable) f.status = 1;
  f.complexE = false;
  f.adjointOK = true;

  f.pencil = false;
  if (mode == 0) {
    f.identity = true;
    f.V.clear();
    f.e.resize(n);
    for (int i = 0; i < n; i++) {
      complex ei(E[i], 0.0);
      if (rmc)
        for (int c = 0; c < nch; c++) {
          const char t = jg->GetChannel(c + 1)->GetRadType();
          const double g = f.gcut[(size_t)i * nch + c];
          if ((t == 'M' || t == 'E') && g != 0.0) ei -= complex(0.0, g * g);
        }
      if (std::imag(ei) != 0.0) f.complexE = true;
      f.e[i] = ei;
    }
    f.emax = 0.0;
    for (int i = 0; i < n; i++) f.emax = std::max(f.emax, std::abs(f.e[i]));
    f.gamma.resize((size_t)n * m);
    for (int i = 0; i < n; i++)
      for (int a = 0; a < m; a++) f.gamma[(size_t)i * m + a] = f.gcut[(size_t)i * nch + f.act[a]];
    return true;
  }

  // Brune and Park: K and N as the level-matrix fill builds them.
  f.identity = false;
  if (rmc || !f.usable) {
    // RMC widths make K complex; the transformation is not needed when the
    // level route is taken anyway.
    if (rmc) f.status = 2;
    f.usable = false;
    f.adjointOK = false;
    return true;
  }
  std::vector<ldouble> K((size_t)n * n, 0.0L), N((size_t)n * n, 0.0L);
  std::vector<double> S((size_t)n * nch, 0.0);
  std::vector<double> parkNorm(n, 1.0);
  for (int i = 0; i < n; i++) {
    ALevel *level = jg->GetLevel(levels[i]);
    if (park) parkNorm[i] = level->GetParkNorm();
    for (int c = 0; c < nch; c++)
      if (particle[c]) S[(size_t)i * nch + c] = level->GetShiftFunction(c + 1);
  }
  // Reference shifts: Lambda_c = L_c + B_c - Sref_c, K -= Sref_c G_c G_c^T.
  std::vector<double> sref(nch, 0.0);
  for (int c = 0; c < nch; c++) {
    if (!particle[c]) continue;
    double sum = 0.0;
    int cnt = 0;
    for (int i = 0; i < n; i++)
      if (f.gcut[(size_t)i * nch + c] != 0.0) {
        sum += S[(size_t)i * nch + c];
        cnt++;
      }
    if (cnt > 0 && std::isfinite(sum)) sref[c] = sum / cnt;
    f.lamShift[c] = jg->GetChannel(c + 1)->GetBoundaryCondition() - sref[c];
  }
  bool finite = true;
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++) {
      ldouble kk = 0.0L, nn = 0.0L;
      if (park) {
        // Park (2021) Eqs. (9), (21)-(22), (28).
        ldouble overlap = 0.0L, shift = 0.0L;
        for (int c = 0; c < nch; c++) {
          if (!particle[c]) continue;
          const ldouble g = (ldouble)f.gcut[(size_t)i * nch + c] * f.gcut[(size_t)k * nch + c];
          if (g == 0.0) continue;
          shift += g * S[(size_t)i * nch + c];
          if (i != k) overlap -= g * ((ldouble)S[(size_t)i * nch + c] - S[(size_t)k * nch + c]) / ((ldouble)E[i] - E[k]);
        }
        if (i == k) overlap = parkNorm[i];
        nn = overlap;
        kk = E[i] * overlap + shift;
      } else {
        // Brune (2002): the shift-function term s0 + E s1.
        if (i == k) {
          kk = E[i];
          nn = 1.0L;
        }
        for (int c = 0; c < nch; c++) {
          if (!particle[c]) continue;
          const ldouble g = (ldouble)f.gcut[(size_t)i * nch + c] * f.gcut[(size_t)k * nch + c];
          if (g == 0.0) continue;
          const ldouble Si = S[(size_t)i * nch + c], Sk = S[(size_t)k * nch + c];
          if (i == k)
            kk += g * Si;
          else {
            const ldouble den = (ldouble)E[i] - E[k];
            kk += g * (Sk * E[i] - Si * E[k]) / den;
            nn -= g * (Si - Sk) / den;
          }
        }
      }
      for (int c = 0; c < nch; c++)
        if (sref[c] != 0.0) kk -= (ldouble)sref[c] * f.gcut[(size_t)i * nch + c] * f.gcut[(size_t)k * nch + c];
      if (!std::isfinite(kk) || !std::isfinite(nn)) finite = false;
      K[(size_t)i * n + k] = kk;
      N[(size_t)i * n + k] = nn;
    }
  // K of the Park level matrix is symmetric only up to round-off.
  for (int i = 0; i < n; i++)
    for (int k = i + 1; k < n; k++) {
      const ldouble s = 0.5L * (K[(size_t)i * n + k] + K[(size_t)k * n + i]);
      K[(size_t)i * n + k] = K[(size_t)k * n + i] = s;
      const ldouble t = 0.5L * (N[(size_t)i * n + k] + N[(size_t)k * n + i]);
      N[(size_t)i * n + k] = N[(size_t)k * n + i] = t;
    }
  std::vector<ldouble> e, V;
  double minPivot = 0.0;
  f.pencil = false;
  std::vector<ldouble> Kd, Nd;
  if (finite) {
    Kd = K;
    Nd = N;
  }
  if (!finite || !GeneralizedEigen(K, N, n, e, V, minPivot)) {
    if (Debug())
      std::cerr << "[channel-space] J-group " << jg->GetJ() << (jg->GetPi() < 0 ? "-" : "+") << " (" << n
                << " levels, " << m << " channels): " << (finite ? "overlap matrix N not positive definite, pivot "
                                                                 : "non-finite K or N (equal level energies?)")
                << (finite ? std::to_string(minPivot) : std::string())
                << (finite ? "; pencil route" : "; level-space route") << std::endl;
    f.adjointOK = false;
    if (!finite) {
      f.usable = false;
      f.status = 2;
      return true;
    }
    // The pencil K - E N, solved at each energy.
    f.pencil = true;
    f.K.swap(Kd);
    f.N.swap(Nd);
    f.V.clear();
    f.e.clear();
    f.emax = 0.0;
    f.gamma.resize((size_t)n * m);
    for (int i = 0; i < n; i++)
      for (int a = 0; a < m; a++) f.gamma[(size_t)i * m + a] = f.gcut[(size_t)i * nch + f.act[a]];
    return true;
  }
  f.K.clear();
  f.N.clear();
  f.e.resize(n);
  f.emax = 0.0;
  for (int i = 0; i < n; i++) {
    f.e[i] = (double)e[i];
    f.emax = std::max(f.emax, std::abs(f.e[i]));
  }
  f.V.resize((size_t)n * n);
  for (size_t k = 0; k < V.size(); k++) f.V[k] = (double)V[k];
  if (Debug()) {
    double vmax = 0.0;
    for (double v : f.V) vmax = std::max(vmax, std::fabs(v));
    std::cerr << "[channel-space] J " << jg->GetJ() << (jg->GetPi() < 0 ? "-" : "+") << " n " << n << ": N min pivot "
              << minPivot << ", e " << (double)e.front() << " .. " << (double)e.back() << ", max|V| " << vmax << std::endl;
  }
  // Gamma = V^T G over the active channels, in extended precision.
  f.gamma.assign((size_t)n * m, 0.0);
  {
    std::vector<ldouble> acc(m);
    for (int k = 0; k < n; k++) {
      std::fill(acc.begin(), acc.end(), 0.0L);
      for (int r = 0; r < n; r++) {
        const ldouble v = V[(size_t)r * n + k];
        if (v == 0.0L) continue;
        const double *g = &f.gcut[(size_t)r * nch];
        for (int a = 0; a < m; a++) acc[a] += v * g[f.act[a]];
      }
      for (int a = 0; a < m; a++) f.gamma[(size_t)k * m + a] = (double)acc[a];
    }
  }
  // The eigen-decomposition adjoint divides by e_l - e_k: formal levels that
  // both couple to a channel must not be (nearly) degenerate.
  const double emax = f.emax;
  std::vector<char> coupled(n, 0);
  for (int i = 0; i < n; i++)
    for (int a = 0; a < m; a++)
      if (f.gamma[(size_t)i * m + a] != 0.0) coupled[i] = 1;
  for (int i = 0; i < n && f.adjointOK; i++)
    for (int k = i + 1; k < n; k++)
      if (coupled[i] && coupled[k] && std::fabs((double)(e[i] - e[k])) < 1.0e-7 * (1.0 + emax)) {
        f.adjointOK = false;
        break;
      }
  return true;
}

}  // namespace ChannelSpace
