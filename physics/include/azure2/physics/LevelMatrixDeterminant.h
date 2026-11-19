#ifndef AZURE2_PHYSICS_LEVEL_MATRIX_DETERMINANT_H
#define AZURE2_PHYSICS_LEVEL_MATRIX_DETERMINANT_H

// The determinant of the inverse level matrix of one J^pi block at one energy (docs/upgrade/physics_tab.md section 5):
//
//   D(E) = det (A^-1)(E),   (A^-1)_lm = (E_l - E) delta_lm - sum_c gamma_lc gamma_mc (S_c(E) - B_c + i P_c)
//
// (the standard form; fill_level_matrix of CollisionMatrix.h has the Brune and Park forms, whose matrix is passed
// here the same way).  The poles of the collision matrix are the zeros of D continued to complex energy; on the real
// axis |D(E)| has a minimum near each resonance, and for one level and one channel D = E_1 - E - gamma^2 (S - B + i P)
// is the denominator of the Breit-Wigner form (Lane and Thomas, Rev. Mod. Phys. 30 (1958) 257, section XII).  Unit:
// MeV^n for n levels.  Diagnostic only: the engine's calculations never form it (they factor A^-1 and solve,
// LevelMatrixSolver.h).

#include <complex>
#include <vector>

namespace azure2 {
namespace physics {

/// det M of a square complex matrix (rows of equal length n), by Gaussian elimination with partial pivoting on a copy.
/// 1 for n = 0; 0 when a pivot is exactly zero; NaN when an element is not finite or the rows are not n long.
std::complex<double> level_matrix_determinant(const std::vector<std::vector<std::complex<double>>> &M);

}  // namespace physics
}  // namespace azure2

#endif
