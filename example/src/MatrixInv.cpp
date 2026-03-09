#include "MatrixInv.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <iostream>

/*!
 *  The MatrixInv constructor takes a complex matrix as an argument and stores the
 *  inverse in a private member variable accessable by the MatrixInv::inverse() function.
 */

MatrixInv::MatrixInv(const matrix_c &A) {
  inverse_.clear();

  gsl_complex x;
  gsl_matrix_complex * m = gsl_matrix_complex_alloc (A.size(), A.size());
  for(int i=0;i<A.size();i++) {
    for(int ii=0;ii<A.size();ii++) {
      GSL_SET_COMPLEX(&x,real(A[i][ii]),imag(A[i][ii]));
      gsl_matrix_complex_set (m,i,ii,x);
    }
  }
  int psign;
  gsl_permutation *p = gsl_permutation_alloc(A.size());
  gsl_matrix_complex * mi = gsl_matrix_complex_alloc (A.size(), A.size());
  gsl_linalg_complex_LU_decomp(m,p,&psign);
  gsl_linalg_complex_LU_invert(m,p,mi);

  vector_c AI_row;
  for(int i=0;i<A.size();i++) {
    inverse_.push_back(AI_row);
    for(int ii=0;ii<A.size();ii++) {
      x=gsl_matrix_complex_get(mi,i,ii);
      complex inv(GSL_REAL(x),GSL_IMAG(x));
      inverse_[i].push_back(inv);
    }
  }

  gsl_matrix_complex_free(m);
  gsl_matrix_complex_free(mi);
  gsl_permutation_free(p);
}
