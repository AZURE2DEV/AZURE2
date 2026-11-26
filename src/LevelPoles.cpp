#include "LevelPoles.h"

#include <algorithm>
#include <cmath>
#include <exception>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace azure2 {

std::vector<NarrowPole> NarrowPolesOf(const std::vector<PoleLevel> &levels) {
  std::vector<NarrowPole> poles;
  const std::size_t n = levels.size();
  if (n < 2) return poles;
  // W = a a^T over the channels; only levels that share a channel and overlap can interfere.
  std::vector<double> W(n * n, 0.0);
  bool overlap = false;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k <= i; k++) {
      double w = 0.0;
      const std::size_t nc = std::min(levels[i].amplitude.size(), levels[k].amplitude.size());
      for (std::size_t c = 0; c < nc; c++) w += levels[i].amplitude[c] * levels[k].amplitude[c];
      W[i * n + k] = W[k * n + i] = w;
      if (k < i && w != 0.0 && std::abs(levels[i].energy - levels[k].energy) < levels[i].width + levels[k].width)
        overlap = true;
    }
  if (!overlap) return poles;

  // H = A + i B, A = diag(E - mean), B = -W/2, a complex symmetric matrix.  Its eigenvalues are those of the real
  // matrix [[A, -B], [B, A]] with eigenvectors (x, -i x); the others are their complex conjugates.  W is positive
  // semi-definite, so the eigenvalues of H are the ones with a negative imaginary part.
  double mean = 0.0;
  for (const PoleLevel &l : levels) mean += l.energy / n;
  const std::size_t m = 2 * n;
  gsl_matrix *M = gsl_matrix_calloc(m, m);
  gsl_vector_complex *eval = gsl_vector_complex_alloc(m);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(m, m);
  gsl_eigen_nonsymmv_workspace *work = gsl_eigen_nonsymmv_alloc(m);
  for (std::size_t i = 0; i < n; i++) {
    gsl_matrix_set(M, i, i, levels[i].energy - mean);
    gsl_matrix_set(M, n + i, n + i, levels[i].energy - mean);
    for (std::size_t k = 0; k < n; k++) {
      gsl_matrix_set(M, i, n + k, 0.5 * W[i * n + k]);       // -B
      gsl_matrix_set(M, n + i, k, -0.5 * W[i * n + k]);      // B
    }
  }
  // The engine's GSL error handler throws (GSLException); the handler itself is the process's and is left alone here,
  // since this runs inside parallel regions.
  bool solved = false;
  try {
    solved = gsl_eigen_nonsymmv(M, eval, evec, work) == GSL_SUCCESS;
  } catch (const std::exception &) {
    solved = false;
  }
  if (solved) {
    for (std::size_t e = 0; e < m; e++) {
      const gsl_complex value = gsl_vector_complex_get(eval, e);
      const double width = -2.0 * GSL_IMAG(value);
      if (!(width > 0.0) || !std::isfinite(width) || !std::isfinite(GSL_REAL(value))) continue;
      // the pole's weight on each level, from the upper half of its eigenvector
      double total = 0.0;
      std::vector<double> weight(n);
      for (std::size_t i = 0; i < n; i++) {
        const gsl_complex x = gsl_matrix_complex_get(evec, i, e);
        weight[i] = GSL_REAL(x) * GSL_REAL(x) + GSL_IMAG(x) * GSL_IMAG(x);
        total += weight[i];
      }
      if (!(total > 0.0)) continue;
      double narrowest = 0.0;
      for (std::size_t i = 0; i < n; i++)
        if (weight[i] >= 0.1 * total && levels[i].width > 0.0 && (narrowest == 0.0 || levels[i].width < narrowest))
          narrowest = levels[i].width;
      if (narrowest > 0.0 && width < 0.5 * narrowest) poles.push_back({GSL_REAL(value) + mean, width});
    }
  }
  gsl_eigen_nonsymmv_free(work);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(eval);
  gsl_matrix_free(M);
  std::sort(poles.begin(), poles.end(), [](const NarrowPole &a, const NarrowPole &b) { return a.energy < b.energy; });
  return poles;
}

}  // namespace azure2
