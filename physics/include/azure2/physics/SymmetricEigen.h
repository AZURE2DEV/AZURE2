#ifndef AZURE2_PHYSICS_SYMMETRICEIGEN_H
#define AZURE2_PHYSICS_SYMMETRICEIGEN_H

#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {

/*!
 * Real symmetric eigenproblems of the small sizes that occur in the level-parameter transformations (one J^pi block:
 * a few to a few tens of levels).  physics/ may not use GSL, so these are written out:
 *
 *  - symmetric_eigen: A x = e x, A symmetric (n x n), by the cyclic Jacobi method (G. H. Golub and C. F. Van Loan,
 *    Matrix Computations, 4th ed. (2013), Sec. 8.5): rotations until every off-diagonal element is below
 *    1e-300 or below 1e-17 of the Frobenius norm; at most 100 sweeps (n <= 50 converges in < 12).
 *  - generalized_symmetric_eigen: N x = e M x, N symmetric, M symmetric positive definite, by the Cholesky
 *    reduction M = L L^T, C = L^-1 N L^-T, C y = e y, x = L^-T y (Golub and Van Loan Sec. 8.7.2), the route GSL's
 *    gsl_eigen_gensymmv takes.  The eigenvectors are normalized to x^T M x = 1.
 *
 * Matrices are row-major std::vector<double> of size n * n.  Eigenvalues are returned in ascending order; the
 * eigenvectors are the columns of the returned matrix (vectors[i * n + k] is component i of eigenvector k).
 *
 * The sign of an eigenvector is not fixed by the eigenproblem.  Here it is fixed so that the component of largest
 * magnitude (the first of equal ones) is positive, which makes the result reproducible; GSL's choice (used by the
 * engine's EigenFunc) is a different one, so the engine and these functions can differ by the sign of a vector.
 */
struct EigenResult {
  bool ok = false;               ///< false: not converged, M not positive definite, or a non-finite input
  std::vector<double> values;    ///< ascending
  std::vector<double> vectors;   ///< row-major, columns are eigenvectors
};

EigenResult symmetric_eigen(const std::vector<double> &A, std::size_t n);

EigenResult generalized_symmetric_eigen(const std::vector<double> &N, const std::vector<double> &M, std::size_t n);

}  // namespace physics
}  // namespace azure2

#endif
