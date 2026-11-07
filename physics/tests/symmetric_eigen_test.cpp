// The eigensolvers of physics/ (SymmetricEigen.h), checked by the defining identities, not by another solver:
//  1. a closed form: the 2 x 2 symmetric matrix [[a, b], [b, d]] has e = (a + d)/2 -+ sqrt(((a - d)/2)^2 + b^2);
//  2. A v = e v, V^T V = 1, ascending order, sign rule (largest component positive), for pseudo-random symmetric
//     matrices of n = 1 .. 20 with eigenvalues spread over 1e-3 .. 1e3;
//  3. the generalized problem: N v = e M v, v^T M v = 1, v^T N v = e, for random N and positive-definite M;
//  4. refusal: M with a negative eigenvalue gives ok = false; a NaN input gives ok = false.
// A wrong rotation sign, a missing back-substitution or an unnormalized vector fails 2-3 by orders of magnitude.
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "azure2/physics/SymmetricEigen.h"

using namespace azure2::physics;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

// A fixed linear congruential generator, so the matrices are the same on every platform.
struct Rng {
  unsigned long long s = 12345;
  double next() {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(s >> 11) / 9007199254740992.0 * 2.0 - 1.0;
  }
};

std::vector<double> RandomSymmetric(Rng &r, std::size_t n, double scale) {
  std::vector<double> A(n * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j <= i; j++) A[i * n + j] = A[j * n + i] = scale * r.next();
  return A;
}

double MaxAbs(const std::vector<double> &A) {
  double m = 0.0;
  for (double x : A) m = std::fmax(m, std::fabs(x));
  return m;
}

}  // namespace

int main() {
  double worstResidual = 0.0, worstOrtho = 0.0, worstGen = 0.0, worstClosed = 0.0;
  {
    const double a = 1.3, b = -0.7, d = -2.1;
    const EigenResult r = symmetric_eigen({a, b, b, d}, 2);
    const double h = std::sqrt(0.25 * (a - d) * (a - d) + b * b);
    check("2x2 converged", r.ok);
    worstClosed = std::fmax(std::fabs(r.values[0] - (0.5 * (a + d) - h)), std::fabs(r.values[1] - (0.5 * (a + d) + h)));
    check("2x2 closed form", worstClosed < 1e-15);
  }
  Rng rng;
  for (std::size_t n = 1; n <= 20; n++)
    for (int rep = 0; rep < 5; rep++) {
      // A = Q diag(e) Q^T with e spread over decades: build from a random symmetric matrix's eigenvectors would use
      // the code under test, so use a random matrix plus a large graded diagonal instead.
      std::vector<double> A = RandomSymmetric(rng, n, 1.0);
      for (std::size_t i = 0; i < n; i++) A[i * n + i] += std::pow(10.0, -3.0 + 6.0 * (double)i / (double)(n > 1 ? n - 1 : 1));
      const EigenResult r = symmetric_eigen(A, n);
      check("converged n=" + std::to_string(n), r.ok);
      if (!r.ok) continue;
      const double scale = MaxAbs(A);
      for (std::size_t k = 0; k < n; k++) {
        if (k > 0) check("ascending", r.values[k] >= r.values[k - 1]);
        std::size_t big = 0;
        for (std::size_t i = 0; i < n; i++) {
          double s = 0.0;
          for (std::size_t j = 0; j < n; j++) s += A[i * n + j] * r.vectors[j * n + k];
          worstResidual = std::fmax(worstResidual, std::fabs(s - r.values[k] * r.vectors[i * n + k]) / scale);
          if (std::fabs(r.vectors[i * n + k]) > std::fabs(r.vectors[big * n + k])) big = i;
        }
        check("sign rule", r.vectors[big * n + k] > 0.0);
        for (std::size_t k2 = 0; k2 < n; k2++) {
          double s = 0.0;
          for (std::size_t i = 0; i < n; i++) s += r.vectors[i * n + k] * r.vectors[i * n + k2];
          worstOrtho = std::fmax(worstOrtho, std::fabs(s - (k == k2 ? 1.0 : 0.0)));
        }
      }
      // generalized: M = I + 0.3 R^T R / n (positive definite), N random symmetric
      std::vector<double> R = RandomSymmetric(rng, n, 1.0), M(n * n, 0.0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) {
          double s = 0.0;
          for (std::size_t k = 0; k < n; k++) s += R[k * n + i] * R[k * n + j];
          M[i * n + j] = (i == j ? 1.0 : 0.0) + 0.3 * s / (double)n;
        }
      const std::vector<double> N = RandomSymmetric(rng, n, 5.0);
      const EigenResult g = generalized_symmetric_eigen(N, M, n);
      check("generalized converged n=" + std::to_string(n), g.ok);
      if (!g.ok) continue;
      const double ns = MaxAbs(N) + MaxAbs(M);
      for (std::size_t k = 0; k < n; k++) {
        double vMv = 0.0, vNv = 0.0;
        for (std::size_t i = 0; i < n; i++) {
          double nv = 0.0, mv = 0.0;
          for (std::size_t j = 0; j < n; j++) {
            nv += N[i * n + j] * g.vectors[j * n + k];
            mv += M[i * n + j] * g.vectors[j * n + k];
          }
          worstGen = std::fmax(worstGen, std::fabs(nv - g.values[k] * mv) / ns);
          vMv += g.vectors[i * n + k] * mv;
          vNv += g.vectors[i * n + k] * nv;
        }
        worstGen = std::fmax(worstGen, std::fabs(vMv - 1.0));
        worstGen = std::fmax(worstGen, std::fabs(vNv - g.values[k]) / ns);
      }
    }
  check("A v = e v to 1e-14 of max|A|", worstResidual < 1e-14);
  check("V^T V = 1 to 1e-14", worstOrtho < 1e-14);
  check("N v = e M v, v^T M v = 1 to 1e-13", worstGen < 1e-13);
  {
    const EigenResult r = generalized_symmetric_eigen({1.0, 0.0, 0.0, 2.0}, {1.0, 2.0, 2.0, 1.0}, 2);  // M eigenvalues -1, 3
    check("indefinite M refused", !r.ok);
    const EigenResult s = symmetric_eigen({1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 1.0}, 2);
    check("NaN refused", !s.ok);
  }
  std::printf("  2x2 closed form %.2g; residual |Av - ev|/max|A| %.2g; |V^T V - 1| %.2g; generalized %.2g\n", worstClosed,
              worstResidual, worstOrtho, worstGen);
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
