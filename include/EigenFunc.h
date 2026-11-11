#ifndef EIGENFUNC_H
#define EIGENFUNC_H

#include "Constants.h"

/// A function class to solve a eigenvalue problems.

/*!
 * The EigenFunc class is a function class to solve eigenfunction problems.
 * The class is used by the parameter trasfomation subroutines in AZURE2.
 */

class EigenFunc {
 public:
  EigenFunc(const matrix_r &);
  EigenFunc(const matrix_r &, const std::vector<vector_r> &);
  /*!
   * Returns a vector with the calculated eigenvalues.
   */
  const vector_r &eigenvalues() const { return eigenvalues_; };
  /*!
   * Returns a matrix with the calculated eigenvectors.
   */
  const matrix_r &eigenvectors() const { return eigenvectors_; };
  /*!
   * False when GSL reported a failure or a non-finite eigenvalue: for the matrix pair, when the second matrix is not
   * positive definite (gsl_eigen_gensymmv's Cholesky step).  The results are then not to be used.  (AZURE2's GSL
   * error handler ignores these codes, so before 2026-11-01 unset values were read as eigenvalues.)
   */
  bool ok() const { return ok_; };

 private:
  vector_r eigenvalues_;
  matrix_r eigenvectors_;
  bool ok_ = false;
};

#endif
