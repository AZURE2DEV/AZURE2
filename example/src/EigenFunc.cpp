#include "EigenFunc.h"
#include <gsl/gsl_eigen.h>
#include <math.h>

/*!
 *  This constructor calculates the eigenvectors and eigenvalues for a real symmetric matrix.
 */

EigenFunc::EigenFunc(const matrix_r &A) {
  eigenvalues_.clear();
  eigenvectors_.clear();

  gsl_matrix * m = gsl_matrix_alloc(A.size(),A.size());
  for(int i=0;i<A.size();i++) {
    for(int ii=0;ii<A.size();ii++) {
      gsl_matrix_set (m,i,ii,A[i][ii]);
    }
  }

  gsl_vector * eval = gsl_vector_alloc (A.size());
  gsl_matrix * evec = gsl_matrix_alloc (A.size(), A.size());
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (A.size());

  gsl_eigen_symmv (m,eval,evec,w);
  gsl_eigen_symmv_sort (eval,evec,GSL_EIGEN_SORT_VAL_ASC);

  for(int i=0;i<A.size();i++) {
    eigenvalues_.push_back(gsl_vector_get(eval,i));
    vector_r tempRow;
    eigenvectors_.push_back(tempRow);
    for(int ii=0;ii<A.size();ii++) {
      eigenvectors_[i].push_back(gsl_matrix_get(evec,i,ii));
    }
  }
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  gsl_matrix_free (m);
  gsl_eigen_symmv_free (w);

}

/*!
 *  This constructor calulates the eigenvectors and eigenvalues for a real symmetric matrix pair.
 */

EigenFunc::EigenFunc(const matrix_r &A, const std::vector<vector_r > &B) {
  eigenvalues_.clear();
  eigenvectors_.clear();

  gsl_matrix * n = gsl_matrix_alloc(A.size(),A.size());
  gsl_matrix * m = gsl_matrix_alloc(B.size(),B.size());
  for(int i=0;i<A.size();i++) {
    for(int ii=0;ii<A.size();ii++) {
      gsl_matrix_set (n,i,ii,A[i][ii]);
      gsl_matrix_set (m,i,ii,B[i][ii]);
    }
  }
  
  gsl_vector * eval = gsl_vector_alloc (A.size());
  gsl_matrix * evec = gsl_matrix_alloc (A.size(), A.size());

  gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc (A.size());

  gsl_eigen_gensymmv(n,m,eval,evec,w);
  gsl_eigen_gensymmv_sort (eval,evec,GSL_EIGEN_SORT_VAL_ASC);

  for(int i=0;i<A.size();i++) {
    eigenvalues_.push_back(gsl_vector_get(eval,i));
    vector_r tempRow;
    eigenvectors_.push_back(tempRow);
    for(int ii=0;ii<A.size();ii++) {
      eigenvectors_[i].push_back(gsl_matrix_get(evec,i,ii));
    }
  }

  vector_r tempNormVec;
  for(int i=0;i<A.size();i++) tempNormVec.push_back(0.0);

  for(int i=0;i<A.size();i++) {
    for(int ii=0;ii<A.size();ii++) {
      double sum=0.0;
      for(int iii=0;iii<A.size();iii++) {
	sum+=B[ii][iii]*eigenvectors_[iii][i];
      }
      tempNormVec[ii]=sum;
    }
    double sum=0.0;
    for(int ii=0;ii<A.size();ii++) {
      sum+=tempNormVec[ii]*eigenvectors_[ii][i];
    }
    for(int ii=0;ii<A.size();ii++) {
      eigenvectors_[ii][i]=eigenvectors_[ii][i]/sqrt(sum);
    }
  }
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  gsl_matrix_free (m);
  gsl_matrix_free (n);
  gsl_eigen_gensymmv_free (w);
}
