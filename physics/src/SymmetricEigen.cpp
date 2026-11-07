#include "azure2/physics/SymmetricEigen.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace azure2 {
namespace physics {

namespace {

// Ascending order of the eigenvalues, eigenvectors permuted with them, and each vector's sign fixed (largest
// component positive).
void SortAndFix(EigenResult &r, std::size_t n) {
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return r.values[a] < r.values[b]; });
  std::vector<double> values(n), vectors(n * n);
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t src = order[k];
    values[k] = r.values[src];
    std::size_t big = 0;
    for (std::size_t i = 1; i < n; i++)
      if (std::fabs(r.vectors[i * n + src]) > std::fabs(r.vectors[big * n + src])) big = i;
    const double sign = (r.vectors[big * n + src] < 0.0) ? -1.0 : 1.0;
    for (std::size_t i = 0; i < n; i++) vectors[i * n + k] = sign * r.vectors[i * n + src];
  }
  r.values.swap(values);
  r.vectors.swap(vectors);
}

}  // namespace

EigenResult symmetric_eigen(const std::vector<double> &Ain, std::size_t n) {
  if (Ain.size() != n * n) throw std::invalid_argument("symmetric_eigen: matrix size does not match n");
  EigenResult r;
  std::vector<double> A = Ain;
  for (double x : A)
    if (!std::isfinite(x)) return r;
  // Symmetrize (the callers build symmetric matrices; this removes any rounding asymmetry).
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) A[i * n + j] = A[j * n + i] = 0.5 * (A[i * n + j] + A[j * n + i]);
  std::vector<double> V(n * n, 0.0);
  for (std::size_t i = 0; i < n; i++) V[i * n + i] = 1.0;
  double norm = 0.0;
  for (double x : A) norm += x * x;
  norm = std::sqrt(norm);
  const int maxSweeps = 100;
  bool converged = (n <= 1);
  for (int sweep = 0; sweep < maxSweeps && !converged; sweep++) {
    double off = 0.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
    off = std::sqrt(off);
    if (off <= 1e-300 || off <= 1e-17 * norm) {
      converged = true;
      break;
    }
    for (std::size_t p = 0; p < n; p++)
      for (std::size_t q = p + 1; q < n; q++) {
        const double apq = A[p * n + q];
        if (apq == 0.0) continue;
        // Golub and Van Loan Algorithm 8.5.1 (sym.schur2): c, s that zero A[p][q].
        const double tau = (A[q * n + q] - A[p * n + p]) / (2.0 * apq);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        for (std::size_t k = 0; k < n; k++) {  // A <- A J (columns p, q)
          const double akp = A[k * n + p], akq = A[k * n + q];
          A[k * n + p] = c * akp - s * akq;
          A[k * n + q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; k++) {  // A <- J^T A (rows p, q)
          const double apk = A[p * n + k], aqk = A[q * n + k];
          A[p * n + k] = c * apk - s * aqk;
          A[q * n + k] = s * apk + c * aqk;
        }
        A[p * n + q] = A[q * n + p] = 0.0;
        for (std::size_t k = 0; k < n; k++) {  // V <- V J
          const double vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  if (!converged) {
    double off = 0.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
    converged = std::sqrt(off) <= 1e-300 || std::sqrt(off) <= 1e-17 * norm;
  }
  if (!converged) return r;
  r.values.resize(n);
  for (std::size_t i = 0; i < n; i++) r.values[i] = A[i * n + i];
  r.vectors = V;
  SortAndFix(r, n);
  r.ok = true;
  return r;
}

EigenResult generalized_symmetric_eigen(const std::vector<double> &N, const std::vector<double> &M, std::size_t n) {
  if (N.size() != n * n || M.size() != n * n)
    throw std::invalid_argument("generalized_symmetric_eigen: matrix size does not match n");
  EigenResult r;
  // Cholesky M = L L^T (lower L); a non-positive pivot means M is not positive definite.
  std::vector<double> L(n * n, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    double d = M[j * n + j];
    for (std::size_t k = 0; k < j; k++) d -= L[j * n + k] * L[j * n + k];
    if (!(d > 0.0)) return r;
    L[j * n + j] = std::sqrt(d);
    for (std::size_t i = j + 1; i < n; i++) {
      double s = 0.5 * (M[i * n + j] + M[j * n + i]);
      for (std::size_t k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = s / L[j * n + j];
    }
  }
  // C = L^-1 N L^-T: first Y = L^-1 N (forward substitution per column), then C = L^-1 Y^T (Y^T = N L^-T).
  std::vector<double> Y(n * n);
  for (std::size_t col = 0; col < n; col++)
    for (std::size_t i = 0; i < n; i++) {
      double s = N[i * n + col];
      for (std::size_t k = 0; k < i; k++) s -= L[i * n + k] * Y[k * n + col];
      Y[i * n + col] = s / L[i * n + i];
    }
  std::vector<double> C(n * n);
  for (std::size_t col = 0; col < n; col++)
    for (std::size_t i = 0; i < n; i++) {
      double s = Y[col * n + i];  // (Y^T)[i][col]
      for (std::size_t k = 0; k < i; k++) s -= L[i * n + k] * C[k * n + col];
      C[i * n + col] = s / L[i * n + i];
    }
  EigenResult e = symmetric_eigen(C, n);
  if (!e.ok) return r;
  // x = L^-T y (back substitution), then x^T M x = 1.
  r.values = e.values;
  r.vectors.assign(n * n, 0.0);
  for (std::size_t k = 0; k < n; k++) {
    for (std::size_t ii = n; ii-- > 0;) {
      double s = e.vectors[ii * n + k];
      for (std::size_t j = ii + 1; j < n; j++) s -= L[j * n + ii] * r.vectors[j * n + k];
      r.vectors[ii * n + k] = s / L[ii * n + ii];
    }
    double norm = 0.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) norm += r.vectors[i * n + k] * M[i * n + j] * r.vectors[j * n + k];
    if (!(norm > 0.0)) return EigenResult{};
    const double scale = 1.0 / std::sqrt(norm);
    for (std::size_t i = 0; i < n; i++) r.vectors[i * n + k] *= scale;
  }
  // The values are already ascending; fix the signs on x (not y) so the rule refers to the returned vectors.
  SortAndFix(r, n);
  for (double v : r.values)
    if (!std::isfinite(v)) return EigenResult{};
  r.ok = true;
  return r;
}

}  // namespace physics
}  // namespace azure2
