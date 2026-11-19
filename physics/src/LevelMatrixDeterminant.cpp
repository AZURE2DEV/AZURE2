// The determinant of the inverse level matrix (LevelMatrixDeterminant.h).
#include "azure2/physics/LevelMatrixDeterminant.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace azure2 {
namespace physics {

std::complex<double> level_matrix_determinant(const std::vector<std::vector<std::complex<double>>> &M) {
  using complex = std::complex<double>;
  const std::size_t n = M.size();
  const complex nan(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
  std::vector<std::vector<complex>> a = M;
  for (const std::vector<complex> &row : a) {
    if (row.size() != n) return nan;
    for (const complex &x : row)
      if (!std::isfinite(x.real()) || !std::isfinite(x.imag())) return nan;
  }
  complex det(1.0, 0.0);
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::abs(a[i][k]) > std::abs(a[p][k])) p = i;
    if (a[p][k] == complex(0.0, 0.0)) return complex(0.0, 0.0);
    if (p != k) {
      std::swap(a[p], a[k]);
      det = -det;
    }
    det *= a[k][k];
    for (std::size_t i = k + 1; i < n; i++) {
      const complex f = a[i][k] / a[k][k];
      for (std::size_t j = k; j < n; j++) a[i][j] -= f * a[k][j];
    }
  }
  return det;
}

}  // namespace physics
}  // namespace azure2
