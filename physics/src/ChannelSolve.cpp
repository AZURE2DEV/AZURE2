#include "azure2/physics/ChannelSolve.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

// The arithmetic here was moved from the engine (AMatrixFunc, ChannelSpace) with its order of operations: every
// result is the same to the last digit.  A change of an expression changes the output files.

namespace azure2 {
namespace physics {
namespace channel_space {

namespace {
typedef std::complex<double> complex;
typedef std::vector<complex> vector_c;
using std::size_t;
}  // namespace

double lu_factor(complex *a, int *piv, int m) {
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

void lu_solve(const complex *lu, const int *piv, int m, complex *b) {
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

double hessenberg_lu(double *a, char *swapped, int n) {
  double pmin = 0.0, pmax = 0.0;
  for (int k = 0; k < n; k++) {
    double *rowK = &a[(size_t)k * n];
    swapped[k] = 0;
    if (k + 1 < n) {
      double *rowN = &a[(size_t)(k + 1) * n];
      if (std::fabs(rowN[k]) > std::fabs(rowK[k])) {
        for (int j = k; j < n; j++) std::swap(rowK[j], rowN[j]);
        swapped[k] = 1;
      }
    }
    const double mod = std::fabs(rowK[k]);
    if (!(mod > 0.0) || !std::isfinite(mod)) return 0.0;
    if (k == 0) pmin = pmax = mod;
    pmin = std::min(pmin, mod);
    pmax = std::max(pmax, mod);
    if (k + 1 < n) {
      double *rowN = &a[(size_t)(k + 1) * n];
      const double l = rowN[k] / rowK[k];
      rowN[k] = l;
      if (l != 0.0)
        for (int j = k + 1; j < n; j++) rowN[j] -= l * rowK[j];
    }
  }
  return n > 0 ? pmin / pmax : 1.0;
}

void hessenberg_solve(const double *lu, const char *swapped, int n, double *b, int nrhs) {
  for (int k = 0; k + 1 < n; k++) {
    double *bk = &b[(size_t)k * nrhs], *bn = &b[(size_t)(k + 1) * nrhs];
    if (swapped[k])
      for (int r = 0; r < nrhs; r++) std::swap(bk[r], bn[r]);
    const double l = lu[(size_t)(k + 1) * n + k];
    if (l != 0.0)
      for (int r = 0; r < nrhs; r++) bn[r] -= l * bk[r];
  }
  for (int i = n - 1; i >= 0; i--) {
    const double *row = &lu[(size_t)i * n];
    double *bi = &b[(size_t)i * nrhs];
    for (int j = i + 1; j < n; j++) {
      const double u = row[j];
      if (u == 0.0) continue;
      const double *bj = &b[(size_t)j * nrhs];
      for (int r = 0; r < nrhs; r++) bi[r] -= u * bj[r];
    }
    const double inv = 1.0 / row[i];
    for (int r = 0; r < nrhs; r++) bi[r] *= inv;
  }
}

double lu_factor_real(double *a, int *piv, int m) {
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

void lu_solve_real(const double *lu, const int *piv, int m, double *b, int nrhs) {
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

// Refused near a formal pole, where (K - E N) is nearly singular: the estimated relative error of R (round-off
// times the condition of the pencil) times the size of R Lambda must stay below kMaxPoleError.
bool solve_pencil(const FormalBlock &f, ChannelSolve &st, double inEnergy) {
  const int n = f.n, m = (int)f.act.size();
  // H - E T, formed in extended precision: H and E T are much larger than
  // their difference near a formal pole.  Upper Hessenberg.
  st.P.assign((size_t)n * n, 0.0);
  const long double e = inEnergy;
  double pmax = 0.0;
  for (int i = 0; i < n; i++)
    for (int k = std::max(0, i - 1); k < n; k++) {
      const double v = (double)(f.H[(size_t)i * n + k] - e * f.T[(size_t)i * n + k]);
      st.P[(size_t)i * n + k] = v;
      pmax = std::max(pmax, std::fabs(v));
    }
  st.hswap.resize(n);
  const double ratio = hessenberg_lu(st.P.data(), st.hswap.data(), n);
  if (!(ratio > 0.0)) {
    return false;
  }
  // X = (H - E T)^-1 Qt G;  R = (Z^T G)^T X  (= G^T (K - E N)^-1 G).
  st.X.assign(f.gq.begin(), f.gq.end());
  hessenberg_solve(st.P.data(), st.hswap.data(), n, st.X.data(), m);
  st.R.assign((size_t)m * m, complex(0.0, 0.0));
  double rl = 0.0, x2max = 0.0;
  for (int a = 0; a < m; a++) {
    double x2 = 0.0;
    for (int k = 0; k < n; k++) x2 += st.X[(size_t)k * m + a] * st.X[(size_t)k * m + a];
    x2max = std::max(x2max, x2 * std::abs(st.lamAct[a]));
    for (int b = a; b < m; b++) {
      double sab = 0.0, sba = 0.0;
      for (int k = 0; k < n; k++) {
        sab += f.gz[(size_t)k * m + a] * st.X[(size_t)k * m + b];
        sba += f.gz[(size_t)k * m + b] * st.X[(size_t)k * m + a];
      }
      const double r = 0.5 * (sab + sba);
      st.R[(size_t)a * m + b] = st.R[(size_t)b * m + a] = r;
      rl = std::max(rl, std::fabs(r) * std::max(std::abs(st.lamAct[a]), std::abs(st.lamAct[b])));
    }
  }
  // Round-off of R Lambda: the pencil is known to a few ulps of its entries,
  // dP, and dR = -G^T P^-1 dP P^-1 G, so the column of channel a of R Lambda
  // is uncertain by about ulps eps |P| |P^-1 G_a|^2 |Lambda_a| -- the analogue
  // of eta on the formal route (P^-1 G = Z X).
  const double err = kEnergyUlps * std::numeric_limits<double>::epsilon() * pmax * x2max;
  if (!(rl <= kMaxPoleRatio) || !(err <= kMaxPoleError)) {
    return false;
  }
  return true;
}

bool solve_formal(const FormalBlock &f, ChannelSolve &st, double inEnergy) {
  const int n = f.n, m = (int)f.act.size();
  st.dinv.resize(n);
  st.isDefl.assign(n, 0);
  // Levels near a pole are deflated: kept as unknowns of the bordered system
  // (SolveChannelSpace) rather than summed into R.  Level k contributes
  // w_k / (e_k - E) to R, w_k = Gamma_k Gamma_k^T, which the kernel turns into
  // terms of size s_k = |Gamma_k|^2_Lambda / |e_k - E|.  Near a formal pole
  // two things would go wrong that the level matrix does not suffer from: the
  // large terms cancel in (1 - Lambda R)^-1 (kappa_k = s_k), and the absolute
  // error of e_k - E (round-off of E and of the formal energy e_k, a few ulps
  // of the energies involved) becomes a relative error of s_k
  // (eta_k = s_k de / |e_k - E|).
  const double de = kEnergyUlps * std::numeric_limits<double>::epsilon() * (f.emax + std::fabs(inEnergy));
  for (int k = 0; k < n; k++) {
    const complex d = (f.e[k] - inEnergy) + f.elo[k];
    const double *g = &f.gamma[(size_t)k * m];
    double w = 0.0;
    for (int a = 0; a < m; a++) w += g[a] * g[a] * std::abs(st.lamAct[a]);
    st.dinv[k] = (d == complex(0.0, 0.0)) ? complex(0.0, 0.0) : 1.0 / d;
    const double inv = std::abs(st.dinv[k]);
    const double kappa = w * inv;
    if (w > 0.0 && (d == complex(0.0, 0.0) || !(kappa <= kDeflateRatio) ||
                    !(kappa * de * inv <= kMaxPoleError))) {
      // A deflated level is solved with the channels; the round-off of that
      // solve moves e_k - E by a few ulps of the (unit-scale) bordered system,
      // which matters within the level's width w_k: relative error
      // ~ eps w_k / (|e_k - E|^2 + w_k^2).  Narrower than that, the level
      // matrix is taken.
      const double ad = std::abs(d);
      if (!(kEnergyUlps * std::numeric_limits<double>::epsilon() * w / (ad * ad + w * w) <=
            kMaxPoleError)) {
        return false;
      }
      st.isDefl[k] = 1;
      st.defl.push_back(k);
    }
  }
  if ((int)st.defl.size() > kMaxDeflated) {
    return false;
  }
  // R = Gamma^T D^-1 Gamma, symmetric; real unless the formal energies are complex.
  st.R.assign((size_t)m * m, complex(0.0, 0.0));
  if (!f.complexE) {
    std::vector<double> &rr = st.rreal;
    rr.assign((size_t)m * m, 0.0);
    for (int k = 0; k < n; k++) {
      if (st.isDefl[k]) continue;
      const double w = std::real(st.dinv[k]);
      const double *g = &f.gamma[(size_t)k * m];
      for (int a = 0; a < m; a++) {
        const double ga = g[a] * w;
        if (ga == 0.0) continue;
        double *row = &rr[(size_t)a * m];
        for (int b = a; b < m; b++) row[b] += ga * g[b];
      }
    }
    for (int a = 0; a < m; a++)
      for (int b = a; b < m; b++) st.R[(size_t)a * m + b] = st.R[(size_t)b * m + a] = rr[(size_t)a * m + b];
  } else {
    for (int k = 0; k < n; k++) {
      if (st.isDefl[k]) continue;
      const complex w = st.dinv[k];
      const double *g = &f.gamma[(size_t)k * m];
      for (int a = 0; a < m; a++) {
        if (g[a] == 0.0) continue;
        const complex ga = g[a] * w;
        complex *row = &st.R[(size_t)a * m];
        for (int b = a; b < m; b++) row[b] += ga * g[b];
      }
    }
    for (int a = 0; a < m; a++)
      for (int b = a + 1; b < m; b++) st.R[(size_t)b * m + a] = st.R[(size_t)a * m + b];
  }
  return true;
}

Route solve_channel_space(const FormalBlock &f, const std::complex<double> *lo, double inEnergy, ChannelSolve &st) {
  const int nch = f.nch, m = (int)f.act.size();
  st.inverseValid = false;
  st.lam.resize(nch);
  for (int c = 0; c < nch; c++) {
    complex l = lo[c];
    l += f.lamShift[c];
    st.lam[c] = l;
  }
  st.lamAct.resize(m);
  for (int a = 0; a < m; a++) st.lamAct[a] = st.lam[f.act[a]];
  st.defl.clear();
  if (f.pencil) {
    if (!solve_pencil(f, st, inEnergy)) return Route::kPole;
  } else if (!solve_formal(f, st, inEnergy)) return Route::kPole;
  // The bordered system of the deflated levels p and the channels:
  //   S = [ D_p        -Gamma_p    ]
  //       [ -Lambda Gamma_p^T   1 - Lambda R ]
  // (just 1 - Lambda R without deflated levels).
  const int kd = (int)st.defl.size(), ns = kd + m;
  st.lu.assign((size_t)ns * ns, complex(0.0, 0.0));
  for (int p = 0; p < kd; p++) {
    const int lv = st.defl[p];
    st.lu[(size_t)p * ns + p] = (f.e[lv] - inEnergy) + f.elo[lv];
    for (int a = 0; a < m; a++) {
      const double g = f.gamma[(size_t)lv * m + a];
      st.lu[(size_t)p * ns + kd + a] = -g;
      st.lu[(size_t)(kd + a) * ns + p] = -st.lamAct[a] * g;
    }
  }
  for (int a = 0; a < m; a++)
    for (int b = 0; b < m; b++)
      st.lu[(size_t)(kd + a) * ns + kd + b] =
          (a == b ? complex(1.0, 0.0) : complex(0.0, 0.0)) - st.lamAct[a] * st.R[(size_t)a * m + b];
  st.piv.resize(ns);
  if (ns > 0 && !(lu_factor(st.lu.data(), st.piv.data(), ns) >= kMinPivotRatio))
    return Route::kPivot;
  st.F.resize((size_t)m * m);
  st.fcol.assign(m, 0);
  st.usolve.resize(nch);
  st.uvalid.assign(nch, 0);
  return Route::kTaken;
}

void solve_bordered(const FormalBlock &f, ChannelSolve &st, complex *wp, complex *wz) {
  const int kd = (int)st.defl.size(), m = (int)f.act.size(), ns = kd + m;
  if (ns == 0) return;
  vector_c &b = st.w3b;
  b.resize(ns);
  for (int p = 0; p < kd; p++) b[p] = wp[p];
  for (int a = 0; a < m; a++) b[kd + a] = wz[a];
  lu_solve(st.lu.data(), st.piv.data(), ns, b.data());
  for (int p = 0; p < kd; p++) wp[p] = b[p];
  for (int a = 0; a < m; a++) wz[a] = b[kd + a];
}

// With the deflated levels p, F_b = Gamma_p^T y_p + R_b + R z from the bordered system with right-hand side
// (Gamma_pb, Lambda R_b) (without deflated levels (1 - R Lambda)^-1 R_b).
void compute_f_column(const FormalBlock &f, ChannelSolve &st, int b) {
  const int m = (int)f.act.size(), kd = (int)st.defl.size();
  complex *col = &st.F[(size_t)b * m];
  vector_c &yp = st.wp, &z = st.wz;
  yp.resize(kd);
  z.resize(m);
  for (int p = 0; p < kd; p++) yp[p] = f.gamma[(size_t)st.defl[p] * m + b];
  for (int a = 0; a < m; a++) z[a] = st.lamAct[a] * st.R[(size_t)a * m + b];
  solve_bordered(f, st, yp.data(), z.data());
  for (int a = 0; a < m; a++) {
    complex s = st.R[(size_t)a * m + b];
    for (int c = 0; c < m; c++) s += st.R[(size_t)a * m + c] * z[c];
    for (int p = 0; p < kd; p++) s += f.gamma[(size_t)st.defl[p] * m + a] * yp[p];
    col[a] = s;
  }
  st.fcol[b] = 1;
}

// With t = Gamma_r^T D_r^-1 w_r the bordered system with right-hand side (w_p, Lambda t) gives y_p and z, and
// y_r = D_r^-1 (w_r + Gamma_r z).
void formal_apply(const FormalBlock &f, ChannelSolve &st, const complex *w, complex *y) {
  const int n = f.n, m = (int)f.act.size(), kd = (int)st.defl.size();
  vector_c &yp = st.wp, &z = st.wz;
  yp.resize(kd);
  z.assign(m, complex(0.0, 0.0));
  for (int p = 0; p < kd; p++) yp[p] = w[st.defl[p]];
  for (int k = 0; k < n; k++) {
    if (st.isDefl[k]) continue;
    const complex t = st.dinv[k] * w[k];
    if (t == complex(0.0, 0.0)) continue;
    const double *g = &f.gamma[(size_t)k * m];
    for (int a = 0; a < m; a++) z[a] += g[a] * t;
  }
  for (int a = 0; a < m; a++) z[a] *= st.lamAct[a];
  solve_bordered(f, st, yp.data(), z.data());
  for (int k = 0; k < n; k++) {
    if (st.isDefl[k]) continue;
    const double *g = &f.gamma[(size_t)k * m];
    complex s = w[k];
    for (int a = 0; a < m; a++) s += g[a] * z[a];
    y[k] = st.dinv[k] * s;
  }
  for (int p = 0; p < kd; p++) y[st.defl[p]] = yp[p];
}

// V A_f V^T u on the formal route (formal_apply); P^-1 u + P^-1 G z with (1 - Lambda R) z = Lambda G^T P^-1 u on
// the pencil route (P = K - E N).
void apply_a(const FormalBlock &f, ChannelSolve &st, const complex *u, complex *out) {
  const int n = f.n, m = (int)f.act.size();
  vector_c &w = st.w1, &q = st.w2;
  if (f.pencil) {
    // P = Qt^T (H - E T) Z^T:  P^-1 u = Z (H - E T)^-1 Qt u,  P^-1 G = Z X,
    // G^T P^-1 u = X^T Z^T u.
    std::vector<double> &re = st.rreal;
    re.assign((size_t)n * 4, 0.0);
    double *qu = re.data(), *zu = re.data() + 2 * (size_t)n;  // n x 2 each: Qt u and Z^T u
    for (int i = 0; i < n; i++) {
      const double *qt = &f.Qt[(size_t)i * n];
      double qr = 0.0, qi = 0.0, zr = 0.0, zi = 0.0;
      for (int k = 0; k < n; k++) {
        qr += qt[k] * std::real(u[k]);
        qi += qt[k] * std::imag(u[k]);
        const double zz = f.Z[(size_t)k * n + i];
        zr += zz * std::real(u[k]);
        zi += zz * std::imag(u[k]);
      }
      qu[2 * i] = qr;
      qu[2 * i + 1] = qi;
      zu[2 * i] = zr;
      zu[2 * i + 1] = zi;
    }
    hessenberg_solve(st.P.data(), st.hswap.data(), n, qu, 2);
    q.resize(m);
    for (int a = 0; a < m; a++) {
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += st.X[(size_t)k * m + a] * complex(zu[2 * k], zu[2 * k + 1]);
      q[a] = st.lamAct[a] * s;
    }
    solve_bordered(f, st, nullptr, q.data());
    vector_c &y = st.w1;
    y.resize(n);
    for (int k = 0; k < n; k++) {
      complex s(qu[2 * k], qu[2 * k + 1]);
      for (int a = 0; a < m; a++) s += st.X[(size_t)k * m + a] * q[a];
      y[k] = s;
    }
    for (int r = 0; r < n; r++) {
      const double *zr = &f.Z[(size_t)r * n];
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += zr[k] * y[k];
      out[r] = s;
    }
    return;
  }
  w.assign(n, complex(0.0, 0.0));
  if (f.identity)
    for (int k = 0; k < n; k++) w[k] = u[k];
  else
    for (int r = 0; r < n; r++) {
      const complex ur = u[r];
      if (ur == complex(0.0, 0.0)) continue;
      const double *vr = &f.V[(size_t)r * n];
      for (int k = 0; k < n; k++) w[k] += vr[k] * ur;
    }
  vector_c &y = st.w4;
  y.resize(n);
  formal_apply(f, st, w.data(), y.data());
  if (f.identity)
    for (int k = 0; k < n; k++) out[k] = y[k];
  else
    for (int r = 0; r < n; r++) {
      const double *vr = &f.V[(size_t)r * n];
      complex s(0.0, 0.0);
      for (int k = 0; k < n; k++) s += vr[k] * y[k];
      out[r] = s;
    }
}

complex channel_bilinear(const FormalBlock &f, ChannelSolve &st, const complex *u, const complex *v) {
  const int n = f.n;
  vector_c av(n);
  apply_a(f, st, v, av.data());
  complex s(0.0, 0.0);
  for (int k = 0; k < n; k++) s += u[k] * av[k];
  return s;
}

}  // namespace channel_space
}  // namespace physics
}  // namespace azure2
