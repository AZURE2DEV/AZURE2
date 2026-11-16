// The linear profile of the objective's profiled norms and THM experiments (LinearProfile.h).  Moved unchanged from
// src/ThmExperiment.cpp (SolveThmProfile, ThmProfileDerivative, InvertSPD, ThmProfile's members) and from
// Segment::ProfileNormChiSquared (ProfiledNormChi2's arithmetic), with their order of operations.
#include "azure2/inference/LinearProfile.h"

#include <algorithm>
#include <cmath>

namespace azure2 {
namespace inference {

double ProfiledNormChi2(const std::vector<double> &m, const std::vector<double> &d, const std::vector<double> &e,
                        double &norm) {
  double Smm = 0.0, Smd = 0.0, Sdd = 0.0;
  for (size_t i = 0; i < m.size(); i++) {
    if (e[i] == 0.0) continue;
    double w = 1.0 / (e[i] * e[i]);
    Smm += m[i] * m[i] * w;
    Smd += m[i] * d[i] * w;
    Sdd += d[i] * d[i] * w;
  }
  if (Smd > 0.0 && Smm > 0.0) {
    norm = Smm / Smd;
    double chi = Sdd - Smd * Smd / Smm;
    return chi > 0.0 ? chi : 0.0;
  }
  // Degenerate: the norm at 1, each point's ((m - d 1)/(e 1))^2 as DataResidual(f, J, y, dy, 1.0).Chi2() forms it.
  norm = 1.0;
  double chi = 0.0;
  for (size_t i = 0; i < m.size(); i++) {
    const double difference = m[i] - d[i] * 1.0, error = e[i] * 1.0;
    chi += error != 0.0 ? (difference * difference) / (error * error) : 0.0;
  }
  return chi;
}

namespace {

// Inverse of a symmetric positive-definite q x q matrix (q <= 4), by Cholesky
// on the Jacobi-scaled matrix (unit diagonal), so that the columns' units
// (the model, 1, E, E^2) do not matter.  False if it is not positive definite
// to working precision.
bool InvertSPD(const std::vector<double> &Gin, int q, std::vector<double> &inv) {
  std::vector<double> dscale(q), G(q * q), L(q * q, 0.0);
  for (int i = 0; i < q; i++) {
    if (!(Gin[i * q + i] > 0.0) || !std::isfinite(Gin[i * q + i])) return false;
    dscale[i] = 1.0 / std::sqrt(Gin[i * q + i]);
  }
  for (int i = 0; i < q; i++)
    for (int j = 0; j < q; j++) G[i * q + j] = Gin[i * q + j] * dscale[i] * dscale[j];
  for (int j = 0; j < q; j++) {
    double sum = G[j * q + j];
    for (int k = 0; k < j; k++) sum -= L[j * q + k] * L[j * q + k];
    if (!(sum > 1.0e-12)) return false;  // singular (collinear columns) to working precision
    L[j * q + j] = std::sqrt(sum);
    for (int i = j + 1; i < q; i++) {
      double s = G[i * q + j];
      for (int k = 0; k < j; k++) s -= L[i * q + k] * L[j * q + k];
      L[i * q + j] = s / L[j * q + j];
    }
  }
  // inv(G) = L^-T L^-1, column by column.
  inv.assign(q * q, 0.0);
  for (int c = 0; c < q; c++) {
    std::vector<double> y(q, 0.0), x(q, 0.0);
    for (int i = 0; i < q; i++) {
      double s = (i == c) ? 1.0 : 0.0;
      for (int k = 0; k < i; k++) s -= L[i * q + k] * y[k];
      y[i] = s / L[i * q + i];
    }
    for (int i = q - 1; i >= 0; i--) {
      double s = y[i];
      for (int k = i + 1; k < q; k++) s -= L[k * q + i] * x[k];
      x[i] = s / L[i * q + i];
    }
    for (int i = 0; i < q; i++) inv[i * q + c] = x[i] * dscale[i] * dscale[c];
  }
  return true;
}

// Column j of the design matrix at a point: the model, then 1, E, E^2.
inline double Column(int j, bool withModel, double m, double energy) {
  if (withModel) {
    if (j == 0) return m;
    j--;
  }
  return j == 0 ? 1.0 : j == 1 ? energy : energy * energy;
}

}  // namespace

double ThmProfile::Background(double energy) const {
  if (terms == 0) return 0.0;
  double b = a[0] + (terms > 1 ? a[1] * energy : 0.0) + (terms > 2 ? a[2] * energy * energy : 0.0);
  return b / s;
}

double ThmProfile::Residual(double m, double d, double e, double energy) const {
  if (e == 0.0) return 0.0;
  return (s * m + s * Background(energy) - d) / e;
}

void ThmProfile::Reported(double value[4], double covariance[16]) const {
  value[0] = 1.0 / s;
  for (int k = 0; k < 3; k++) value[k + 1] = k < terms ? a[k] / s : 0.0;
  for (int i = 0; i < 16; i++) covariance[i] = 0.0;
  // Rows of d(n, b0, b1, b2)/d(c), c the profiled coefficients.
  const int q = (scaleProfiled ? 1 : 0) + (backgroundProfiled ? terms : 0);
  if (q == 0 || (int)cov.size() != q * q) return;
  std::vector<double> T(4 * q, 0.0);
  int col = 0;
  if (scaleProfiled) {
    T[0 * q + col] = -1.0 / (s * s);
    for (int k = 0; k < terms; k++) T[(k + 1) * q + col] = -a[k] / (s * s);
    col++;
  }
  if (backgroundProfiled)
    for (int k = 0; k < terms; k++) T[(k + 1) * q + col + k] = 1.0 / s;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      double sum = 0.0;
      for (int u = 0; u < q; u++)
        for (int v = 0; v < q; v++) sum += T[i * q + u] * cov[u * q + v] * T[j * q + v];
      covariance[i * 4 + j] = sum;
    }
}

ThmProfile SolveThmProfile(const std::vector<double> &m, const std::vector<double> &d,
                           const std::vector<double> &e, const std::vector<double> &energy, int terms) {
  ThmProfile p;
  p.terms = terms;
  const size_t n = m.size();
  for (size_t i = 0; i < n; i++)
    if (e[i] != 0.0) p.points++;

  // Normal equations over the columns in use; `withModel` includes the model.
  auto solve = [&](bool withModel, std::vector<double> &c, std::vector<double> &cov) -> bool {
    const int q = (withModel ? 1 : 0) + terms;
    if (q == 0) return false;
    std::vector<double> G(q * q, 0.0), h(q, 0.0), col(q);
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double w = 1.0 / (e[i] * e[i]);
      // Without the model column the target is d - m (s fixed at 1).
      const double y = withModel ? d[i] : d[i] - m[i];
      for (int j = 0; j < q; j++) col[j] = Column(j, withModel, m[i], energy[i]);
      for (int j = 0; j < q; j++) {
        h[j] += col[j] * y * w;
        for (int k = 0; k < q; k++) G[j * q + k] += col[j] * col[k] * w;
      }
    }
    if (!InvertSPD(G, q, cov)) return false;
    c.assign(q, 0.0);
    for (int j = 0; j < q; j++)
      for (int k = 0; k < q; k++) c[j] += cov[j * q + k] * h[k];
    return true;
  };

  std::vector<double> c, cov;
  if (terms == 0) {
    // The model column alone: s = S_md/S_mm, as Segment::ProfileNormChiSquared.
    double Smm = 0.0, Smd = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double w = 1.0 / (e[i] * e[i]);
      Smm += m[i] * m[i] * w;
      Smd += m[i] * d[i] * w;
    }
    if (Smd > 0.0 && Smm > 0.0) {
      p.s = Smd / Smm;
      p.cov.assign(1, 1.0 / Smm);
    } else {
      p.scaleProfiled = false;
      p.s = 1.0;
      p.status = "degenerate: no positive model-data overlap, norm left at 1";
    }
    p.backgroundProfiled = false;
  } else if (solve(true, c, cov) && c[0] > 0.0) {
    p.s = c[0];
    for (int k = 0; k < terms; k++) p.a[k] = c[k + 1];
    p.cov = cov;
  } else {
    p.scaleProfiled = false;
    p.s = 1.0;
    if (solve(false, c, cov)) {
      for (int k = 0; k < terms; k++) p.a[k] = c[k];
      p.cov = cov;
      p.status = "degenerate: model and background cannot both be profiled (singular, or n <= 0); "
                 "norm left at 1, background profiled";
    } else {
      p.backgroundProfiled = false;
      p.status = "degenerate: too few points or collinear background terms; norm left at 1, no background";
    }
  }
  p.chi2 = 0.0;
  for (size_t i = 0; i < n; i++) {
    double r = p.Residual(m[i], d[i], e[i], energy[i]);
    p.chi2 += r * r;
  }
  return p;
}

void ThmProfileDerivative(const ThmProfile &p, const std::vector<double> &m, const std::vector<double> &d,
                          const std::vector<double> &e, const std::vector<double> &energy,
                          const std::vector<double> &Jm, int nCols, std::vector<double> &J,
                          std::vector<double> &G) {
  const size_t n = m.size();
  J.assign(n * nCols, 0.0);
  G.assign(n * nCols, 0.0);
  const bool withModel = p.scaleProfiled;
  const int q = (withModel ? 1 : 0) + (p.backgroundProfiled ? p.terms : 0);
  // Coefficients in the order of the profiled columns, and the fixed scale.
  std::vector<double> c;
  if (withModel) c.push_back(p.s);
  if (p.backgroundProfiled)
    for (int k = 0; k < p.terms; k++) c.push_back(p.a[k]);
  const double s = p.s;
  std::vector<double> rho(n, 0.0);  // y~ - A~ c = -r
  for (size_t i = 0; i < n; i++) rho[i] = -p.Residual(m[i], d[i], e[i], energy[i]);

  std::vector<double> v(q), dc(q);
  for (int col = 0; col < nCols; col++) {
    // v = (dA~)^T rho - A~^T (dA~) c, dA~ = [Jm/e, 0, ...]; (dA~) c = s Jm/e.
    std::fill(v.begin(), v.end(), 0.0);
    for (size_t i = 0; i < n; i++) {
      if (e[i] == 0.0) continue;
      const double jt = Jm[i * nCols + col] / e[i];
      if (withModel) v[0] += jt * rho[i];
      for (int j = 0; j < q; j++) v[j] -= Column(j, withModel, m[i], energy[i]) / e[i] * s * jt;
    }
    for (int j = 0; j < q; j++) {
      dc[j] = 0.0;
      for (int k = 0; k < q; k++) dc[j] += p.cov[j * q + k] * v[k];
    }
    for (size_t i = 0; i < n; i++) {
      // d f_i/dp = s Jm_i + sum_j A_ij dc_j (unwhitened), f the curve on the data scale.
      double df = s * Jm[i * nCols + col];
      for (int j = 0; j < q; j++) df += Column(j, withModel, m[i], energy[i]) * dc[j];
      if (e[i] != 0.0) J[i * nCols + col] = df / e[i];
      G[i * nCols + col] = df / s;
    }
  }
}

}  // namespace inference
}  // namespace azure2
