#include "MatrixInv.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <iostream>

/*!
 *  The MatrixInv constructor takes a complex matrix as an argument and stores the
 *  inverse in a private member variable accessable by the MatrixInv::inverse() function.
 */

MatrixInv::MatrixInv(const matrix_c &A) {
    int n = A.size();
    inverse_.assign(n, std::vector<std::complex<double>>(n));

    gsl_matrix_complex *m  = gsl_matrix_complex_alloc(n, n);
    gsl_matrix_complex *mi = gsl_matrix_complex_alloc(n, n);
    gsl_permutation *p     = gsl_permutation_alloc(n);

    // Fill matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            gsl_matrix_complex_set(m, i, j,
                gsl_complex_rect(std::real(A[i][j]), std::imag(A[i][j])));
        }
    }

    // LU decomposition + inversion
    int signum;
    gsl_linalg_complex_LU_decomp(m, p, &signum);
    gsl_linalg_complex_LU_invert(m, p, mi);

    // Copy back to inverse_
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            gsl_complex z = gsl_matrix_complex_get(mi, i, j);
            inverse_[i][j] = {GSL_REAL(z), GSL_IMAG(z)};
        }
    }

    gsl_matrix_complex_free(m);
    gsl_matrix_complex_free(mi);
    gsl_permutation_free(p);
}
