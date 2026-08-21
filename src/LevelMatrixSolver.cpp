#include "LevelMatrixSolver.h"
#include <cmath>

/*!
 * Right-looking LU with partial pivoting.  The pivot is chosen on \f$ |z|^2 \f$,
 * which orders candidates identically to \f$ |z| \f$ while avoiding a square
 * root per comparison.
 */

bool LevelMatrixSolver::Decompose(const matrix_c &A) {
  n_ = (int)A.size();
  singular_ = false;
  inverse_valid_ = false;
  if (n_ == 0) return true;

  lu_.resize((size_t)n_ * n_);
  piv_.resize(n_);
  for (int i = 0; i < n_; i++) {
    const vector_c &src = A[i];
    complex *dest = &lu_[(size_t)i * n_];
    for (int j = 0; j < n_; j++) dest[j] = src[j];
  }

  for (int k = 0; k < n_; k++) {
    int p = k;
    double best = std::norm(lu_[(size_t)k * n_ + k]);
    for (int i = k + 1; i < n_; i++) {
      double value = std::norm(lu_[(size_t)i * n_ + k]);
      if (value > best) {
        best = value;
        p = i;
      }
    }
    piv_[k] = p;
    if (p != k) {
      complex *rowK = &lu_[(size_t)k * n_];
      complex *rowP = &lu_[(size_t)p * n_];
      for (int j = 0; j < n_; j++) std::swap(rowK[j], rowP[j]);
    }
    complex pivot = lu_[(size_t)k * n_ + k];
    if ((std::real(pivot) == 0. && std::imag(pivot) == 0.) ||
        !std::isfinite(std::real(pivot)) || !std::isfinite(std::imag(pivot))) {
      singular_ = true;
      return false;
    }
    complex pivotInv = 1.0 / pivot;
    const complex *rowK = &lu_[(size_t)k * n_];
    for (int i = k + 1; i < n_; i++) {
      complex *rowI = &lu_[(size_t)i * n_];
      complex factor = (rowI[k] *= pivotInv);
      if (std::real(factor) == 0. && std::imag(factor) == 0.) continue;
      for (int j = k + 1; j < n_; j++) rowI[j] -= factor * rowK[j];
    }
  }
  return true;
}

/*!
 * Applies the stored row interchanges, then forward and back substitution.
 */

void LevelMatrixSolver::Solve(complex *b) const {
  for (int k = 0; k < n_; k++)
    if (piv_[k] != k) std::swap(b[k], b[piv_[k]]);
  for (int i = 1; i < n_; i++) {
    const complex *row = &lu_[(size_t)i * n_];
    complex sum = b[i];
    for (int j = 0; j < i; j++) sum -= row[j] * b[j];
    b[i] = sum;
  }
  for (int i = n_ - 1; i >= 0; i--) {
    const complex *row = &lu_[(size_t)i * n_];
    complex sum = b[i];
    for (int j = i + 1; j < n_; j++) sum -= row[j] * b[j];
    b[i] = sum / row[i];
  }
}

void LevelMatrixSolver::SolveInto(const vector_c &b, vector_c &x) const {
  x.assign(b.begin(), b.end());
  if (n_ > 0) this->Solve(&x[0]);
}

/*!
 * The bilinear form needs only \f$ M^{-1} v \f$, so one triangular solve
 * replaces the double loop over levels that an explicit inverse would require.
 */

complex LevelMatrixSolver::Bilinear(const vector_c &u, const vector_c &v) const {
  if (n_ == 0) return complex(0.0, 0.0);
  work_.assign(v.begin(), v.end());
  this->Solve(&work_[0]);
  complex sum(0.0, 0.0);
  for (int i = 0; i < n_; i++) sum += u[i] * work_[i];
  return sum;
}

const matrix_c &LevelMatrixSolver::Inverse() const {
  if (!inverse_valid_) {
    inverse_.resize(n_);
    for (int i = 0; i < n_; i++) inverse_[i].assign(n_, complex(0.0, 0.0));
    work_.resize(n_);
    for (int c = 0; c < n_; c++) {
      for (int i = 0; i < n_; i++) work_[i] = complex(0.0, 0.0);
      work_[c] = complex(1.0, 0.0);
      this->Solve(&work_[0]);
      for (int i = 0; i < n_; i++) inverse_[i][c] = work_[i];
    }
    inverse_valid_ = true;
  }
  return inverse_;
}
