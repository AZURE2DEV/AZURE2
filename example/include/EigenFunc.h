#ifndef EIGENFUNC_H
#define EIGENFUNC_H

#include "Constants.h"

///A function class to solve a eigenvalue problems.

/*!
 * The EigenFunc class is a function class to solve eigenfunction problems.  
 * The class is used by the parameter trasfomation subroutines in AZURE2.  
 */

class EigenFunc {
 public:
  EigenFunc(const matrix_r&);
  EigenFunc(const matrix_r&, const std::vector<vector_r >&);
  /*!
   * Returns a vector with the calculated eigenvalues.
   */
  const vector_r& eigenvalues() const {return eigenvalues_;};
  /*!
   * Returns a matrix with the calculated eigenvectors.  
   */
  const matrix_r& eigenvectors() const {return eigenvectors_;};
 private:
  vector_r eigenvalues_;
  matrix_r eigenvectors_;
};
  
#endif
