#ifndef LEVELMATRIXSOLVER_H
#define LEVELMATRIXSOLVER_H

#include "Constants.h"
#include <vector>

/// A reusable dense complex LU factorization of the level matrix.

/*!
 * The level matrix of one J-group is small (a handful of levels) but is rebuilt
 * and solved for every energy point, so allocation and format conversion cost
 * more than the arithmetic does.  LevelMatrixSolver keeps its scratch storage
 * alive between calls and exposes the operations the physics actually needs:
 *
 *   - Solve()/SolveInto() apply \f$ M^{-1} \f$ to one right-hand side, which is
 *     all that is required for the bilinear forms
 *     \f$ \gamma_c^T M^{-1} \gamma_{c'} \f$ that build the U- and T-matrices;
 *   - Inverse() materializes the full \f$ M^{-1} \f$, but only lazily and only
 *     for the callers that genuinely need every element (the adjoint gradient).
 *
 * The factorization is a right-looking LU with partial pivoting on contiguous
 * row-major storage.  For the level counts that occur in practice this beats a
 * blocked library kernel; above roughly a dozen levels a blocked LU would win
 * on the factorization itself, though the bilinear-form path avoids enough work
 * that the crossover is not reached in normal use.
 */

class LevelMatrixSolver {
 public:
  LevelMatrixSolver() :
    n_(0),
    singular_(false),
    inverse_valid_(false) {};

  /*!
   * Factors the matrix \f$ A \f$, which must be square.  Returns false if a
   * pivot is zero or non-finite, in which case the factorization is unusable
   * and singular() is set.
   */
  bool Decompose(const matrix_c &A);

  /*!
   * Returns the dimension of the factored matrix.
   */
  int size() const { return n_; };
  /*!
   * Returns true if Decompose() encountered a zero or non-finite pivot.
   */
  bool singular() const { return singular_; };

  /*!
   * Replaces the length-size() vector \f$ b \f$ with \f$ M^{-1} b \f$.
   */
  void Solve(complex *b) const;

  /*!
   * Writes \f$ M^{-1} b \f$ into \f$ x \f$, resizing it as needed.
   */
  void SolveInto(const vector_c &b, vector_c &x) const;

  /*!
   * Returns \f$ u^T M^{-1} v \f$ using a single triangular solve.
   */
  complex Bilinear(const vector_c &u, const vector_c &v) const;

  /*!
   * Returns the full inverse, computed on the first call following Decompose()
   * and cached thereafter.
   */
  const matrix_c &Inverse() const;

 private:
  int n_;
  bool singular_;
  /// LU factors of the level matrix, row-major, n_ by n_.
  std::vector<complex> lu_;
  /// Row interchanges chosen by partial pivoting.
  std::vector<int> piv_;
  mutable bool inverse_valid_;
  mutable matrix_c inverse_;
  mutable vector_c work_;
};

#endif
