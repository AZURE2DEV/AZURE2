#ifndef MATRIXINV_H
#define MATRIXINV_H

#include "Constants.h"

///A Function class to perform matrix inversion.

/*!
 * The MatrixInv class is a function class for matrix inversion using the GSL functions.
 */

class MatrixInv {
 public:
  MatrixInv(const matrix_c&);
  /*!
   * The function returns the inverse as calculated by the constructor.
   */
  const matrix_c& inverse() const {return inverse_;};
 private:
 matrix_c inverse_;
};

#endif
