#include "ThmDwVertex.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

/*
 * Numerics (docs/source/theory/thm_implementation.rst, "Distorted-wave
 * entrance vertex").
 *
 * With r = r_xA along z and u = r_sx in the xz plane, R_s = alpha r z^ + u and
 * R_a = r z^ + beta u (polar angles theta_s, theta_a, azimuth 0), the reduced
 * amplitude
 *
 *   h^l_{Ls La}(r) = sqrt(4pi/(2l+1)) 2pi Int u^2 du phi(u) Int dcos(theta_u)
 *                    f_Ls(R_s) f_La(R_a) K^l_{Ls La}(theta_s, theta_a),
 *   K = sum_M <Ls M La -M|l 0> Y_LsM(theta_s, 0) Y_La-M(theta_a, 0),
 *
 * f_L(R) = u_L(kR)/R, does not depend on the directions of k_aA and k_sF.  K is
 * the zero component of the bipolar harmonic {Y_Ls(R_s^) x Y_La(R_a^)}_l;
 * rotating R_s^ onto z leaves only M_s = 0, so
 *
 *   K = sqrt((2Ls+1)/(2l+1)) sum_nu (-1)^nu <Ls 0 La nu|l nu> Y_l nu(theta_s, 0)
 *       Y_La nu(theta_a - theta_s, 0),
 *
 * |nu| <= l: O(l) per (l, Ls, La) instead of O(min(Ls, La)).  The r derivative
 * at fixed u: dR_s/dr = alpha cos theta_s, dtheta_s/dr = -alpha sin theta_s/R_s,
 * dR_a/dr = cos theta_a, dtheta_a/dr = -sin theta_a/R_a.
 *
 * Quadrature: u from rmin to rmin + 26/kappa in 16-point Gauss-Legendre panels
 * of at most 4 radians of (k_sF + beta k_aA) u and 5 fm; cos theta_u on
 * 1.2 Lmax + 20 Gauss-Legendre nodes.  Partial waves L <= k R_max + 14 in each
 * channel (R_max = a + beta u_max, alpha a + u_max).  Only f_Ls depends on the
 * energy: the geometry, K and f_La are evaluated once per point for a batch of
 * energies.  Waves: ThmDistortion::Wave (Numerov from the origin, normalized
 * to F + T H^+, matched to COUL), or the Riccati-Bessel functions for a plane
 * wave, on a uniform grid read by six-point Lagrange interpolation (value and
 * derivative).
 */

namespace {

const double kGridStep = 0.02;   // MeV
const double kTail = 26.0;       // u_max = rmin + kTail/kappa
const int kGL = 16;              // nodes per u panel
const double kPanelPhase = 4.0;  // radians per u panel
const double kPanelMax = 5.0;    // fm
const int kLMargin = 14;
const int kLCap = 200;
const int kBatch = 8;            // energies per pass over the points
const int kSumBlocks = 64;        // blocks of the sum over the points (fixed order, any thread count)
const size_t kSumBytes = (size_t)32 << 20;  // at most this for their buffers

inline int Idx(int L, int M) { return L * (L + 1) / 2 + M; }

/*
 * Y_LM(theta, 0) for 0 <= M <= min(L, mMax), L <= lMax (Condon-Shortley phase)
 * and d/dtheta, by the standard stable recurrences; theta may be negative
 * (Y(-theta, 0) = Y(theta, pi)).  Stored at Idx(L, M); entries with M > mMax
 * are left as they are.
 */
void YTable(int lMax, int mMax, double theta, std::vector<double> &y, std::vector<double> &dy) {
  const int size = (lMax + 1) * (lMax + 2) / 2;
  if ((int)y.size() < size) y.assign(size, 0.0);
  if ((int)dy.size() < size) dy.assign(size, 0.0);
  const double x = std::cos(theta), s = std::sin(theta);
  const int top = std::min(mMax + 1, lMax);  // one order more for the derivative
  double pmm = 1.0 / std::sqrt(4.0 * M_PI);
  for (int M = 0; M <= top; M++) {
    if (M > 0) pmm = -std::sqrt((2.0 * M + 1.0) / (2.0 * M)) * s * pmm;
    y[Idx(M, M)] = pmm;
    if (M < lMax) {
      double p0 = pmm, p1 = std::sqrt(2.0 * M + 3.0) * x * pmm;
      y[Idx(M + 1, M)] = p1;
      for (int L = M + 2; L <= lMax; L++) {
        double a = std::sqrt((4.0 * L * L - 1.0) / ((double)L * L - (double)M * M));
        double b = std::sqrt(((L - 1.0) * (L - 1.0) - (double)M * M) / (4.0 * (L - 1.0) * (L - 1.0) - 1.0));
        double p2 = a * (x * p1 - b * p0);
        y[Idx(L, M)] = p2;
        p0 = p1;
        p1 = p2;
      }
    }
  }
  // dY_LM/dtheta = [sqrt((L-M)(L+M+1)) Y_L,M+1 - sqrt((L+M)(L-M+1)) Y_L,M-1]/2, Y_L,-1 = -Y_L,1.
  for (int L = 0; L <= lMax; L++)
    for (int M = 0; M <= std::min(L, mMax); M++) {
      double up = M + 1 <= L ? std::sqrt((double)(L - M) * (L + M + 1)) * y[Idx(L, M + 1)] : 0.0;
      double down = 0.0;
      if (M > 0)
        down = std::sqrt((double)(L + M) * (L - M + 1)) * y[Idx(L, M - 1)];
      else if (L > 0)
        down = -std::sqrt((double)L * (L + 1)) * y[Idx(L, 1)];
      dy[Idx(L, M)] = 0.5 * (up - down);
    }
}

// Six-point Lagrange weights (value and d/dx) at x on a uniform grid of step h
// with n nodes: the first node j (0 <= j <= n - 6).
struct Lagrange6 {
  int j = 0;
  double w[6], dw[6];
  void Set(double x, double h, int n) {
    double t = x / h;
    j = (int)std::floor(t) - 2;
    if (j < 0) j = 0;
    if (j > n - 6) j = n - 6;
    double s = t - j;
    static const double den[6] = {-120.0, 24.0, -12.0, 12.0, -24.0, 120.0};
    double d[6];
    for (int i = 0; i < 6; i++) d[i] = s - i;
    for (int i = 0; i < 6; i++) {
      double p = 1.0, dp = 0.0;
      for (int m = 0; m < 6; m++) {
        if (m == i) continue;
        dp = dp * d[m] + p;  // derivative of the running product
        p *= d[m];
      }
      w[i] = p / den[i];
      dw[i] = dp / den[i] / h;
    }
  }
  void Apply(const complex *u, complex &v, complex &dv) const {
    v = dv = complex(0.0, 0.0);
    for (int i = 0; i < 6; i++) {
      v += w[i] * u[j + i];
      dv += dw[i] * u[j + i];
    }
  }
};

// e^{i sigma_L}, L = 0..lMax (eta = 0: 1).
std::vector<complex> CoulombPhases(double eta, int lMax) {
  std::vector<complex> out(lMax + 1, complex(1.0, 0.0));
  if (eta == 0.0) return out;
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e(1.0, eta, &lnr, &arg);
  double sigma = arg.val;
  for (int l = 0; l <= lMax; l++) {
    if (l > 0) sigma += std::atan(eta / l);
    out[l] = std::polar(1.0, sigma);
  }
  return out;
}

// u_L(kR) for L = 0..lMax on j step, j < n: Riccati-Bessel for a plane wave,
// else ThmDistortion::Wave.  False (and the L) if a wave cannot be normalized.
bool WaveTable(const ThmDistortion &d, const ThmDistortion::Channel &c, int lMax, double step, int n,
               std::vector<std::vector<complex>> &u, int &bad) {
  u.assign(lMax + 1, std::vector<complex>());
  if (c.kind == ThmDistortion::Channel::PLANE) {
    for (int L = 0; L <= lMax; L++) u[L].assign(n, complex(0.0, 0.0));
    std::vector<double> jl(lMax + 1);
    for (int i = 1; i < n; i++) {
      double x = c.k * i * step;
      gsl_sf_bessel_jl_steed_array(lMax, x, jl.data());
      for (int L = 0; L <= lMax; L++) u[L][i] = x * jl[L];
    }
    return true;
  }
  for (int L = 0; L <= lMax; L++)
    if (!d.Wave(c, L, step, n, u[L])) {
      bad = L;
      return false;
    }
  return true;
}

std::string Number(double x) {
  std::ostringstream s;
  s.precision(6);
  s << x;
  return s.str();
}

}  // namespace

double ThmCG(int j1, int m1, int j2, int m2, int J, int M) {
  if (m1 + m2 != M || J < std::abs(j1 - j2) || J > j1 + j2 || std::abs(m1) > j1 || std::abs(m2) > j2 ||
      std::abs(M) > J)
    return 0.0;
  auto lf = [](int n) { return lgammal((long double)n + 1.0L); };
  long double pre = 0.5L * (logl(2.0L * J + 1.0L) + lf(J + j1 - j2) + lf(J - j1 + j2) + lf(j1 + j2 - J) -
                            lf(j1 + j2 + J + 1) + lf(J + M) + lf(J - M) + lf(j1 - m1) + lf(j1 + m1) + lf(j2 - m2) +
                            lf(j2 + m2));
  long double sum = 0.0L;
  int kMin = std::max(0, std::max(j2 - J - m1, j1 - J + m2));
  int kMax = std::min(j1 + j2 - J, std::min(j1 - m1, j2 + m2));
  for (int k = kMin; k <= kMax; k++) {
    long double t = expl(pre - lf(k) - lf(j1 + j2 - J - k) - lf(j1 - m1 - k) - lf(j2 + m2 - k) -
                         lf(J - j2 + m1 + k) - lf(J - j1 - m2 + k));
    sum += (k % 2) ? -t : t;
  }
  return (double)sum;
}

int ThmDwVertex::LIndex(int l) const {
  for (size_t i = 0; i < lvals.size(); i++)
    if (lvals[i] == l) return (int)i;
  return -1;
}

void ThmDwVertex::Factor(const double g[4], complex a[2], complex d[2]) {
  const double g11 = std::max(g[0], 0.0), g22 = std::max(g[1], 0.0);
  const complex g12(g[2], g[3]);
  a[0] = a[1] = d[0] = d[1] = complex(0.0, 0.0);
  if (g11 >= g22) {
    if (!(g11 > 0.0)) return;
    double r = std::sqrt(g11);
    a[0] = r;
    d[0] = g12 / r;
    d[1] = std::sqrt(std::max(0.0, g22 - std::norm(g12) / g11));
  } else {
    // Pivoted: G' = (g22, conj g12; g12, g11) for c' = (-1, B - 1).
    double r = std::sqrt(g22);
    d[0] = r;
    a[0] = std::conj(g12) / r;
    a[1] = std::sqrt(std::max(0.0, g11 - std::norm(g12) / g22));
  }
}

double ThmDwVertex::Vertex2(const double g[4], complex B) {
  const complex c1 = B - 1.0;
  const complex g12(g[2], g[3]);
  // c^+ G c with c = (c1, -1): |c1|^2 g11 + g22 - 2 Re(conj(c1) g12).
  return std::norm(c1) * g[0] + g[1] - 2.0 * std::real(std::conj(c1) * g12);
}

double ThmDwVertex::Vertex2Factored(const double g[4], complex B) {
  complex a[2], d[2];
  Factor(g, a, d);
  return std::norm(a[0] * (B - 1.0) - d[0]) + std::norm(a[1] * (B - 1.0) - d[1]);
}

std::string ThmDwVertex::Build(const ThmExperiment &x, const ThmDistortion::Kinematics &k, double a,
                               const std::vector<int> &ls, double eLo, double eHi,
                               const std::vector<double> &points) {
  const auto t0 = std::chrono::steady_clock::now();
  experiment = x.name;
  radius = a;
  lvals = ls;
  std::sort(lvals.begin(), lvals.end());
  lvals.erase(std::unique(lvals.begin(), lvals.end()), lvals.end());
  if (lvals.empty()) return "the entrance pair has no channel";
  if (!(a > 0.0)) return "the channel radius of the entrance pair is not positive";
  const int nl = (int)lvals.size();
  const int lTop = lvals.back();

  gridLo = eLo;
  gridStep = kGridStep;
  nE = std::max(4, (int)std::ceil((eHi - eLo) / kGridStep - 1.0e-9) + 1);
  if (!points.empty()) {
    dist.dataLo = *std::min_element(points.begin(), points.end());
    dist.dataHi = *std::max_element(points.begin(), points.end());
  }
  std::string why = dist.Setup(x, k, eLo);
  if (!why.empty()) return why;
  // Every data point must be reachable, as for R(E) (EData): a lab angle the
  // spectator cannot reach would otherwise be clamped (SpectatorCos) to a
  // different direction.  Grid energies beyond the data may be out of reach.
  for (double e : points) {
    std::string w = dist.CheckEnergy(e);
    if (!w.empty()) return w;
  }
  for (int e = 0; e < nE; e++) {
    std::string w = dist.CheckEnergy(eLo + e * gridStep);
    if (!w.empty() && dist.angleKind != ThmDistortion::LAB) return w;
  }
  alpha = k.mA / (k.mx + k.mA);
  beta = k.ms / k.ma;
  const double ka = dist.aa.k;
  const double kappa = dist.kappa;
  const double rmin = dist.rmin;
  uMax = rmin + kTail / kappa;

  // The acceptance: the accepted spectator directions (spectatorAngles=
  // and/or the |p_s| cut of a ps window; ThmDistortion::AngleNodes).
  angles = dist.angWindow;
  nNodes = 1;
  nAng = 0;
  if (angles) {
    // The model is linear in G (sum over the Cholesky components of
    // |X a_k + Y d_k|^2 = |X|^2 G11 + |Y|^2 G22 + 2 Re(X* Y G12)), so the
    // average over the directions is one vertex with the averaged G: one node.
    dist.SetAngleSlots(eHi);
    nAng = dist.angSlots * dist.angNodes;
  }

  // Partial waves and the u, cos theta_u quadrature.
  std::vector<double> ksE(nE), esfE(nE);
  double ksMax = 0.0;
  for (int e = 0; e < nE; e++) {
    esfE[e] = dist.EsF(eLo + e * gridStep);
    ksE[e] = std::sqrt(2.0 * dist.sf.mu * esfE[e]) / hbarc;
    ksMax = std::max(ksMax, ksE[e]);
  }
  const double rA = a + beta * uMax, rS = alpha * a + uMax;
  laMax = (int)std::ceil(ka * rA) + kLMargin;
  lsMax = (int)std::ceil(ksMax * rS) + kLMargin;
  if (laMax > kLCap || lsMax > kLCap)
    return "the partial-wave sums would need L up to " + Number(std::max(laMax, lsMax)) + " (kappa too small)";
  std::vector<double> uNode, uWt;
  {
    double width = std::min(kPanelMax, kPanelPhase / (ksMax + beta * ka));
    int panels = std::max(2, (int)std::ceil((uMax - rmin) / width));
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(kGL);
    for (int p = 0; p < panels; p++) {
      double lo = rmin + (uMax - rmin) * p / panels, hi = rmin + (uMax - rmin) * (p + 1) / panels;
      for (int i = 0; i < kGL; i++) {
        double xi, wi;
        gsl_integration_glfixed_point(lo, hi, i, &xi, &wi, t);
        uNode.push_back(xi);
        uWt.push_back(wi);
      }
    }
    gsl_integration_glfixed_table_free(t);
  }
  cNodes = (int)std::ceil(1.2 * std::max(laMax, lsMax)) + 20;
  std::vector<double> cNode(cNodes), cWt(cNodes);
  {
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(cNodes);
    for (int i = 0; i < cNodes; i++) gsl_integration_glfixed_point(-1.0, 1.0, i, &cNode[i], &cWt[i], t);
    gsl_integration_glfixed_table_free(t);
  }
  uNodes = (int)uNode.size();
  std::vector<double> phiU(uNodes);
  for (int i = 0; i < uNodes; i++) phiU[i] = dist.Phi(uNode[i]);

  // Triples (l, Ls, La): |Ls - La| <= l <= Ls + La, Ls + La + l even; the
  // coefficients sqrt((2Ls+1)/(2l+1)) (-1)^nu <Ls 0 La nu|l nu> (x2 for nu > 0).
  struct Triple {
    int li, l, ls, la;
    std::vector<double> c;   // nu = 0..min(l, La)
    std::vector<double> cm;  // <Ls m La 0|l m> sqrt((2La+1)/4pi), m = 0..min(l, Ls)
  };
  std::vector<Triple> triples;
  for (int li = 0; li < nl; li++) {
    int l = lvals[li];
    for (int Ls = 0; Ls <= lsMax; Ls++)
      for (int La = std::abs(Ls - l); La <= std::min(laMax, Ls + l); La++) {
        if ((Ls + La + l) % 2) continue;
        Triple t{li, l, Ls, La, {}, {}};
        double pre = std::sqrt((2.0 * Ls + 1.0) / (2.0 * l + 1.0));
        for (int nu = 0; nu <= std::min(l, La); nu++)
          t.c.push_back(pre * (nu % 2 ? -1.0 : 1.0) * (nu ? 2.0 : 1.0) * ThmCG(Ls, 0, La, nu, l, nu));
        for (int m = 0; m <= std::min(l, Ls); m++)
          t.cm.push_back(ThmCG(Ls, m, La, 0, l, m) * std::sqrt((2.0 * La + 1.0) / (4.0 * M_PI)));
        triples.push_back(t);
      }
  }
  const int nT = (int)triples.size();

  // The a + A waves, once.
  auto tableStep = [](const ThmDistortion::Channel &c, double e) {
    double kl = std::sqrt(2.0 * c.mu * (std::max(e, 0.0) + ThmDistortion::Depth(c))) / hbarc;
    return std::min(0.02, 0.1 / std::max(kl, 1.0e-6));
  };
  const double stepA = tableStep(dist.aa, dist.eAA);
  const int nA = (int)std::ceil(rA / stepA) + 8;
  std::vector<std::vector<complex>> uA;
  int bad = -1;
  if (!WaveTable(dist, dist.aa, laMax, stepA, nA, uA, bad))
    return "the a + A wave l = " + Number(bad) + " could not be normalized";
  std::vector<complex> sigA = CoulombPhases(dist.aa.eta, laMax);
  const double stepS = tableStep(dist.sf, *std::max_element(esfE.begin(), esfE.end()));
  const int nS = (int)std::ceil(rS / stepS) + 8;

  // h and dh/dr at r = a: hh[(t nE + e) 2 + {0,1}].
  std::vector<complex> hh((size_t)nT * nE * 2, complex(0.0, 0.0));
  const int nPts = uNodes * cNodes;
  for (int e0 = 0; e0 < nE; e0 += kBatch) {
    const int nb = std::min(kBatch, nE - e0);
    std::vector<std::vector<std::vector<complex>>> uS(nb);
    std::vector<std::string> errors(nb);
#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nb; b++) {
      ThmDistortion::Channel c = dist.SfAt(eLo + (e0 + b) * gridStep);
      int badL = -1;
      if (!WaveTable(dist, c, lsMax, stepS, nS, uS[b], badL))
        errors[b] = "the s + F wave l = " + Number(badL) + " at E = " + Number(eLo + (e0 + b) * gridStep) +
                    " MeV could not be normalized";
    }
    for (const std::string &err : errors)
      if (!err.empty()) return err;
    // The sum over the points in a fixed order, whatever the number of
    // threads and the schedule: up to kSumBlocks blocks of consecutive points
    // (fewer if their buffers would pass kSumBytes), each summed in point
    // order into its own buffer, and the buffers added in block order.
    const size_t blockSize = (size_t)nT * nb * 2;
    const int nBlocks = (int)std::max<size_t>(
        1, std::min<size_t>({(size_t)nPts, (size_t)kSumBlocks, kSumBytes / (blockSize * sizeof(complex))}));
    std::vector<complex> blockSums((size_t)nBlocks * blockSize, complex(0.0, 0.0));
#pragma omp parallel
    {
      std::vector<double> yl, dyl, ya, dya;
      std::vector<complex> fa(laMax + 1), dfa(laMax + 1), fs((size_t)nb * (lsMax + 1)), dfs((size_t)nb * (lsMax + 1));
      Lagrange6 lagA, lagS;
#pragma omp for schedule(dynamic, 1)
      for (int blk = 0; blk < nBlocks; blk++) {
        complex *acc = &blockSums[(size_t)blk * blockSize];
        const int ipEnd = (int)((long long)nPts * (blk + 1) / nBlocks);
        for (int ip = (int)((long long)nPts * blk / nBlocks); ip < ipEnd; ip++) {
          const int iu = ip / cNodes, ic = ip % cNodes;
          const double u = uNode[iu], ct = cNode[ic], st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
          const double wgt = 2.0 * M_PI * u * u * phiU[iu] * uWt[iu] * cWt[ic];
          if (wgt == 0.0) continue;
          const double xs = u * st, zs = alpha * a + u * ct, xa = beta * u * st, za = a + beta * u * ct;
          const double Rs = std::hypot(xs, zs), Ra = std::hypot(xa, za);
          if (!(Rs > 0.0) || !(Ra > 0.0)) continue;
          const double ths = std::atan2(xs, zs), tha = std::atan2(xa, za);
          const double cs = zs / Rs, ss = xs / Rs, ca = za / Ra, sa = xa / Ra;
          YTable(lTop, lTop, ths, yl, dyl);
          YTable(laMax, lTop, tha - ths, ya, dya);
          lagA.Set(Ra, stepA, nA);
          for (int L = 0; L <= laMax; L++) {
            complex v, dv;
            lagA.Apply(uA[L].data(), v, dv);
            fa[L] = v / Ra;
            dfa[L] = dv / Ra - v / (Ra * Ra);
          }
          lagS.Set(Rs, stepS, nS);
          for (int b = 0; b < nb; b++)
            for (int L = 0; L <= lsMax; L++) {
              complex v, dv;
              lagS.Apply(uS[b][L].data(), v, dv);
              fs[b * (lsMax + 1) + L] = v / Rs;
              dfs[b * (lsMax + 1) + L] = dv / Rs - v / (Rs * Rs);
            }
          const double dthsdr = -alpha * ss / Rs, dthadr = -sa / Ra;
          for (int t = 0; t < nT; t++) {
            const Triple &tr = triples[t];
            double K = 0.0, Ks = 0.0, Ka = 0.0;  // K and dK/dtheta_s, dK/dtheta_a
            for (int nu = 0; nu < (int)tr.c.size(); nu++) {
              const double yln = yl[Idx(tr.l, nu)], dyln = dyl[Idx(tr.l, nu)];
              const double yan = ya[Idx(tr.la, nu)], dyan = dya[Idx(tr.la, nu)];
              K += tr.c[nu] * yln * yan;
              Ka += tr.c[nu] * yln * dyan;
              Ks += tr.c[nu] * (dyln * yan - yln * dyan);
            }
            const complex fA = fa[tr.la], dfA = dfa[tr.la];
            const complex A1 = wgt * fA * K;
            const complex A2 = wgt * alpha * cs * fA * K;
            const complex A3 = wgt * (fA * (Ks * dthsdr + Ka * dthadr) + ca * dfA * K);
            complex *out = &acc[(size_t)t * nb * 2];
            for (int b = 0; b < nb; b++) {
              const complex f = fs[b * (lsMax + 1) + tr.ls], df = dfs[b * (lsMax + 1) + tr.ls];
              out[2 * b] += A1 * f;
              out[2 * b + 1] += A2 * df + A3 * f;
            }
          }
        }
      }
    }
    for (int blk = 0; blk < nBlocks; blk++) {
      const complex *acc = &blockSums[(size_t)blk * blockSize];
      for (int t = 0; t < nT; t++)
        for (int b = 0; b < nb; b++) {
          hh[((size_t)t * nE + e0 + b) * 2] += acc[((size_t)t * nb + b) * 2];
          hh[((size_t)t * nE + e0 + b) * 2 + 1] += acc[((size_t)t * nb + b) * 2 + 1];
        }
    }
  }

  // Per energy and node: s_m, d_m and the Gram matrix.
  w.assign((size_t)nE * nNodes, 0.0);
  q.assign((size_t)nE * nNodes, 0.0);
  G.assign((size_t)nE * nNodes * nl * 4, 0.0);
  Gd.assign((size_t)nE * nl * 4, 0.0);
  aw.assign((size_t)nE * nAng, 0.0);
  aq.assign((size_t)nE * nAng, 0.0);
  ath.assign((size_t)nE * nAng, 0.0);
  qd.assign(nE, 0.0);
  pd.assign(nE, 0.0);
  thd.assign(nE, 0.0);
  valid.assign(nE, 1);
  std::vector<double> ys, dys;
  for (int e = 0; e < nE; e++) {
    const double ks = ksE[e];
    const double etas = dist.sf.kind == ThmDistortion::Channel::PLANE
                            ? 0.0
                            : k.Zs * (k.Zx + k.ZA) * fstruc * dist.sf.mu / (hbarc * ks);
    std::vector<complex> sigS = CoulombPhases(etas, lsMax);
    // Reduced amplitudes (4pi)^2/(ks ka) (-i)^Ls i^La e^{i(sigma+sigma)} sqrt(4pi/(2l+1)) h.
    std::vector<complex> H(nT), dH(nT);
    for (int t = 0; t < nT; t++) {
      const Triple &tr = triples[t];
      complex ph = sigS[tr.ls] * sigA[tr.la];
      int n = ((tr.la - tr.ls) % 4 + 4) % 4;  // i^(La - Ls)
      static const complex ipow[4] = {complex(1, 0), complex(0, 1), complex(-1, 0), complex(0, -1)};
      ph *= ipow[n] * (16.0 * M_PI * M_PI / (ks * ka)) * std::sqrt(4.0 * M_PI / (2.0 * tr.l + 1.0));
      H[t] = ph * hh[((size_t)t * nE + e) * 2];
      dH[t] = ph * hh[((size_t)t * nE + e) * 2 + 1];
    }
    // The Gram entries at an angle theta_sa between k_sF and k_aA, for q (fm^-1).
    auto gram = [&](double xsa, double qq, double *g) {
      double phit = 0.0;
      for (int i = 0; i < uNodes; i++) {
        double z = qq * uNode[i];
        double j0 = z < 1.0e-4 ? 1.0 - z * z / 6.0 : std::sin(z) / z;
        phit += uWt[i] * uNode[i] * uNode[i] * phiU[i] * j0;
      }
      phit *= 4.0 * M_PI;
      const double th = std::acos(std::max(-1.0, std::min(1.0, xsa)));
      YTable(lsMax, lTop, th, ys, dys);
      for (int li = 0; li < nl; li++) {
        const int l = lvals[li];
        double g11 = 0.0, g22 = 0.0;
        complex g12(0.0, 0.0);
        for (int m = 0; m <= l; m++) {
          complex sm(0.0, 0.0), dm(0.0, 0.0);
          for (int t = 0; t < nT; t++) {
            const Triple &tr = triples[t];
            if (tr.li != li || m > tr.ls) continue;
            double c = tr.cm[m] * ys[Idx(tr.ls, m)];
            if (c == 0.0) continue;
            sm += c * H[t];
            dm += c * dH[t];
          }
          sm /= 4.0 * M_PI * phit;
          dm *= a / (4.0 * M_PI * phit);
          double mult = m ? 2.0 : 1.0;
          g11 += mult * std::norm(sm);
          g22 += mult * std::norm(dm);
          g12 += mult * std::conj(sm) * dm;
        }
        const double f = 4.0 * M_PI / (2.0 * l + 1.0);
        g[li * 4 + 0] = f * g11;
        g[li * 4 + 1] = f * g22;
        g[li * 4 + 2] = f * g12.real();
        g[li * 4 + 3] = f * g12.imag();
      }
      return phit;
    };
    // The spectatorAngle direction.
    {
      double th = 0.0;
      double xsa = dist.SpectatorCos(ks, &th);
      const double kb = beta * ka;
      qd[e] = std::sqrt(std::max(0.0, ks * ks + kb * kb - 2.0 * ks * kb * xsa));
      pd[e] = std::sqrt(std::max(0.0, ka * ka + alpha * alpha * ks * ks - 2.0 * ka * alpha * ks * xsa));
      thd[e] = std::acos(std::max(-1.0, std::min(1.0, xsa))) * 180.0 / M_PI;
      double phit = gram(xsa, qd[e], &Gd[(size_t)e * nl * 4]);
      if (!(std::fabs(phit) > 0.0) || !std::isfinite(phit))
        return "the plane-wave source phi~(q) vanishes at E = " + Number(eLo + e * gridStep) + " MeV";
    }
    if (angles) {
      // The accepted directions: event weight d cos(theta_cm) x acceptance x
      // |phi~(q)|^2 (the data are divided by |phi|^2 averaged over the same
      // events, and G is normalized by phi~ of its own node).
      std::vector<ThmDistortion::AngleNode> an;
      if (!dist.AngleNodes(ks, an)) {
        valid[e] = 0;
        continue;
      }
      double total = 0.0;
      std::vector<double> gk((size_t)nAng * nl * 4, 0.0);
      for (int i = 0; i < nAng; i++) {
        const ThmDistortion::AngleNode &n = an[i];
        double *ai = &aw[(size_t)e * nAng + i];
        aq[(size_t)e * nAng + i] = n.q * hbarc;
        ath[(size_t)e * nAng + i] = n.theta;
        if (!(n.w > 0.0)) continue;  // an empty branch
        double phit = gram(n.x, n.q, &gk[(size_t)i * nl * 4]);
        if (!(std::fabs(phit) > 0.0) || !std::isfinite(phit))
          return "the plane-wave source phi~(q) vanishes at E = " + Number(eLo + e * gridStep) + " MeV, q = " +
                 Number(n.q * hbarc) + " MeV/c";
        *ai = n.w * phit * phit;
        total += *ai;
      }
      if (!(total > 0.0) || !std::isfinite(total)) {
        valid[e] = 0;
        continue;
      }
      double qm = 0.0;
      for (int i = 0; i < nAng; i++) {
        const double wi = aw[(size_t)e * nAng + i] /= total;
        qm += wi * aq[(size_t)e * nAng + i];
        for (int c = 0; c < nl * 4; c++) G[(size_t)e * nl * 4 + c] += wi * gk[(size_t)i * nl * 4 + c];
      }
      w[e] = 1.0;
      q[e] = qm;
      continue;
    }
    w[e] = 1.0;
    q[e] = qd[e] * hbarc;
    std::copy(&Gd[(size_t)e * nl * 4], &Gd[(size_t)e * nl * 4] + nl * 4, &G[(size_t)e * nl * 4]);
  }

  // Every data point needs a reachable window; grid energies without one take
  // the nearest valid entry (they are beyond the data).
  for (double e : points) {
    double t = (e - gridLo) / gridStep;
    int i0 = std::max(0, std::min(nE - 1, (int)std::floor(t))), i1 = std::min(nE - 1, i0 + 1);
    if (!valid[i0] || !valid[i1]) {
      std::string why = dist.CheckWindow(e);
      return why.empty() ? "at E = " + Number(e) + " MeV the spectator-direction window is out of reach next to it"
                         : why;
    }
  }
  for (int e = 0; e < nE; e++) {
    if (valid[e]) continue;
    int best = -1;
    for (int f = 0; f < nE; f++)
      if (valid[f] && (best < 0 || std::abs(f - e) < std::abs(best - e))) best = f;
    if (best < 0) return "the spectator-direction window is out of reach at every energy";
    for (int i = 0; i < nAng; i++) {
      aw[(size_t)e * nAng + i] = aw[(size_t)best * nAng + i];
      aq[(size_t)e * nAng + i] = aq[(size_t)best * nAng + i];
      ath[(size_t)e * nAng + i] = ath[(size_t)best * nAng + i];
    }
    for (int i = 0; i < nNodes; i++) {
      w[(size_t)e * nNodes + i] = w[(size_t)best * nNodes + i];
      q[(size_t)e * nNodes + i] = q[(size_t)best * nNodes + i];
      for (int c = 0; c < nl * 4; c++)
        G[((size_t)e * nNodes + i) * nl * 4 + c] = G[((size_t)best * nNodes + i) * nl * 4 + c];
    }
  }

  std::ostringstream d;
  d.precision(6);
  d << "vertexModel=dw: surface term of the prior-form DWBA; a + A: " << dist.ChannelText(0)
    << "; s + F: " << dist.ChannelText(1) << "; bound state " << (dist.yukawa ? "yukawa" : "whittaker")
    << (rmin > 0.0 ? ", r >= " + Number(rmin) + " fm" : "") << "; "
    << (angles ? "spectator directions " + dist.AngleText() + " (" + Number(dist.angNodes) +
                     " nodes in cos theta_cm" + (dist.angSlots > 1 ? " per branch" : "") +
                     (dist.qCut && !dist.angAll ? ", |p_s| cut by the ps window" : "") +
                     "; weight d cos theta_cm x acceptance x |phi~(q)|^2)"
               : std::string("spectator angle ") +
                     (dist.angleKind == ThmDistortion::QF
                          ? "qf (k_sF along k_aA)"
                          : (dist.angleKind == ThmDistortion::LAB ? "lab " : "cm ") + Number(dist.angle) + " deg"));
  description = d.str();
  buildSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return "";
}

bool ThmDwVertex::Reached(double energy) const {
  if (nE == 0) return false;
  const double t = (energy - gridLo) / gridStep;
  if (t < 0.0 || t > nE - 1.0) return false;
  const int i0 = std::max(0, std::min(nE - 1, (int)std::floor(t))), i1 = std::min(nE - 1, i0 + 1);
  return valid.empty() || (valid[i0] && valid[i1]);
}

void ThmDwVertex::AngleNodesAt(double energy, std::vector<double> &weight, std::vector<double> &qk,
                               std::vector<double> &theta) const {
  weight.clear();
  qk.clear();
  theta.clear();
  if (!angles || nE == 0) return;
  const int e = std::max(0, std::min(nE - 1, (int)std::lround((energy - gridLo) / gridStep)));
  for (int i = 0; i < nAng; i++) {
    weight.push_back(aw[(size_t)e * nAng + i]);
    qk.push_back(aq[(size_t)e * nAng + i]);
    theta.push_back(ath[(size_t)e * nAng + i]);
  }
}

void ThmDwVertex::Interpolate(double energy, std::vector<double> &weight, std::vector<double> &qk,
                              std::vector<double> &g, std::vector<double> &gd, double &qDelta, double &pDelta,
                              bool *outside) const {
  const int nl = (int)lvals.size();
  double t = (energy - gridLo) / gridStep;
  bool out = t < 0.0 || t > nE - 1.0;
  if (outside) *outside = out;
  t = std::max(0.0, std::min(nE - 1.0, t));
  int j = std::min(std::max((int)std::floor(t), 1), nE - 3);
  double s = t - j;
  const double c[4] = {-s * (s - 1.0) * (s - 2.0) / 6.0, (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0,
                       -(s + 1.0) * s * (s - 2.0) / 2.0, (s + 1.0) * s * (s - 1.0) / 6.0};
  auto lag = [&](const std::vector<double> &v, size_t stride, size_t off) {
    double r = 0.0;
    for (int i = 0; i < 4; i++) r += c[i] * v[(size_t)(j - 1 + i) * stride + off];
    return r;
  };
  weight.assign(nNodes, 0.0);
  qk.assign(nNodes, 0.0);
  g.assign((size_t)nNodes * nl * 4, 0.0);
  gd.assign((size_t)nl * 4, 0.0);
  double total = 0.0;
  for (int k = 0; k < nNodes; k++) {
    weight[k] = std::max(0.0, lag(w, nNodes, k));
    total += weight[k];
    qk[k] = lag(q, nNodes, k);
    for (int m = 0; m < nl * 4; m++) g[(size_t)k * nl * 4 + m] = lag(G, (size_t)nNodes * nl * 4, (size_t)k * nl * 4 + m);
  }
  if (total > 0.0)
    for (double &v : weight) v /= total;
  for (int m = 0; m < nl * 4; m++) gd[m] = lag(Gd, (size_t)nl * 4, m);
  qDelta = lag(qd, 1, 0);
  pDelta = lag(pd, 1, 0);
}

void ThmDwVertex::Evaluate(double energy, At &out) const {
  std::vector<double> qk, g, gd;
  double qq, pp;
  Interpolate(energy, out.weight, qk, g, gd, qq, pp, &out.outside);
  const int nl = (int)lvals.size();
  out.nodes = nNodes;
  out.nl = nl;
  out.a.assign((size_t)nNodes * nl * 2, complex(0.0, 0.0));
  out.d.assign((size_t)nNodes * nl * 2, complex(0.0, 0.0));
  for (int k = 0; k < nNodes; k++)
    for (int li = 0; li < nl; li++) Factor(&g[((size_t)k * nl + li) * 4], &out.a[((size_t)k * nl + li) * 2], &out.d[((size_t)k * nl + li) * 2]);
}
