#include "MatrixInv.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <stdexcept>

/*!
 *  The MatrixInv constructor takes a complex matrix as an argument and stores the
 *  inverse in a private member variable accessable by the MatrixInv::inverse() function.
 */

MatrixInv::MatrixInv(const matrix_c &A) {
    int n = A.size();
    
    // Handle empty or invalid matrices gracefully
    if(n == 0) {
        inverse_.clear();
        return;
    }
    
    inverse_.assign(n, std::vector<std::complex<double>>(n));

    gsl_matrix_complex *m = nullptr;
    gsl_matrix_complex *mi = nullptr;
    gsl_permutation *p = nullptr;
    
    try {
        // Safe allocation with error checking
        m = gsl_matrix_complex_alloc(n, n);
        if(!m) {
            throw std::runtime_error("Failed to allocate GSL matrix m");
        }
        
        mi = gsl_matrix_complex_alloc(n, n);
        if(!mi) {
            throw std::runtime_error("Failed to allocate GSL matrix mi");
        }
        
        p = gsl_permutation_alloc(n);
        if(!p) {
            throw std::runtime_error("Failed to allocate GSL permutation");
        }

        // Fill matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                gsl_matrix_complex_set(m, i, j,
                    gsl_complex_rect(std::real(A[i][j]), std::imag(A[i][j])));
            }
        }

        // LU decomposition + inversion
        int signum;
        int lu_result = gsl_linalg_complex_LU_decomp(m, p, &signum);
        if(lu_result != GSL_SUCCESS) {
            throw std::runtime_error("GSL LU decomposition failed");
        }
        
        int inv_result = gsl_linalg_complex_LU_invert(m, p, mi);
        if(inv_result != GSL_SUCCESS) {
            throw std::runtime_error("GSL matrix inversion failed");
        }

        // Copy back to inverse_
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                gsl_complex z = gsl_matrix_complex_get(mi, i, j);
                inverse_[i][j] = {GSL_REAL(z), GSL_IMAG(z)};
            }
        }
        
        // Safe cleanup
        gsl_matrix_complex_free(m);
        gsl_matrix_complex_free(mi);
        gsl_permutation_free(p);
        
    } catch(...) {
        // Exception-safe cleanup - only free non-null pointers
        if(m) gsl_matrix_complex_free(m);
        if(mi) gsl_matrix_complex_free(mi);
        if(p) gsl_permutation_free(p);
        throw;
    }
}
