// level_matrix_determinant (LevelMatrixDeterminant.h) against closed forms: the 2 x 2 rule ad - bc, the 3 x 3 cofactor
// expansion (a matrix whose first pivot is zero, so the row exchange and its sign are exercised), a singular matrix,
// a wrong shape; and the inverse level matrix of one level and one channel from fill_level_matrix, whose determinant
// is the Breit-Wigner denominator E_1 - E - gamma^2 (S - B + i P) written out by hand.  A sign error in the row
// exchange, a missing pivot or a transposed element fails these at the 1e-14 level.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "azure2/physics/CollisionMatrix.h"
#include "azure2/physics/LevelMatrixDeterminant.h"

using complex = std::complex<double>;
using matrix = std::vector<std::vector<complex>>;
namespace phys = azure2::physics;

namespace {

int failures = 0;

void Check(const char *what, complex got, complex want, double tol = 1e-13) {
  const double scale = std::max(1.0, std::abs(want));
  const bool ok = std::abs(got - want) <= tol * scale;
  std::printf("%s %s: got (%.15g, %.15g), want (%.15g, %.15g)\n", ok ? "ok" : "FAIL", what, got.real(), got.imag(),
              want.real(), want.imag());
  if (!ok) failures++;
}

complex Det3(const matrix &m) {  // cofactor expansion along the first row
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

}  // namespace

int main() {
  const complex i(0.0, 1.0);
  Check("empty", phys::level_matrix_determinant({}), 1.0);
  const matrix two = {{1.0 + 2.0 * i, -0.5}, {3.0, 0.25 - i}};
  Check("2x2", phys::level_matrix_determinant(two), two[0][0] * two[1][1] - two[0][1] * two[1][0]);
  const matrix three = {{0.0, 2.0 - i, 0.5}, {1.5 + i, -1.0, 2.0 * i}, {0.3, 4.0, -2.0 + 0.5 * i}};
  Check("3x3, zero first pivot", phys::level_matrix_determinant(three), Det3(three));
  const matrix singular = {{1.0, 2.0}, {2.0, 4.0}};
  Check("singular", phys::level_matrix_determinant(singular), 0.0);
  const complex bad = phys::level_matrix_determinant({{1.0, 2.0}, {3.0}});
  const bool nan = std::isnan(bad.real()) && std::isnan(bad.imag());
  std::printf("%s ragged rows give NaN\n", nan ? "ok" : "FAIL");
  if (!nan) failures++;

  // One level, one channel (the standard form): A^-1 = E_1 - E - gamma^2 (S - B + i P).
  const double E1 = 2.365, gamma = 0.71, S = -0.42, B = -0.3, P = 0.13, E = 2.3;
  const complex lo(S - B, P);
  const char rad = 'P';
  phys::LevelMatrixTables t;
  t.n = 1, t.nc = 1, t.E = &E1, t.gamma = &gamma, t.lo = &lo, t.rad_type = &rad;
  const double zero = 0.0;
  t.shift = &zero, t.shift_d1 = &zero, t.shift_d2 = &zero, t.boundary = &B;
  phys::LevelMatrixSolver::matrix_c M;
  phys::fill_level_matrix(t, E, M);
  Check("one level, one channel", phys::level_matrix_determinant(M), E1 - E - gamma * gamma * lo);

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
