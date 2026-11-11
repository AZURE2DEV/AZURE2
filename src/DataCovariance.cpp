#include "InputFiles.h"
#include "DataCovariance.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace DataCovariance {

std::shared_ptr<File> Read(const std::string &path, std::string &err, bool cached) {
  azure2::experiment::CovarianceFileAccess access;
  access.open = [cached](const std::string &file) { return InputFiles::Open(file, cached); };
  return azure2::experiment::ReadCovarianceFile(path, access, err);
}

int Cholesky(const std::vector<double> &A, int n, std::vector<double> &L) {
  L.assign((size_t)n * n, 0.0);
  for (int j = 0; j < n; j++) {
    const double *Lj = &L[(size_t)j * n];
    double d = A[(size_t)j * n + j];
    for (int k = 0; k < j; k++) d -= Lj[k] * Lj[k];
    if (!(d > 0.0) || !std::isfinite(d)) return j;
    const double ljj = std::sqrt(d);
    L[(size_t)j * n + j] = ljj;
#pragma omp parallel for schedule(static) if (n - j > 256)
    for (int i = j + 1; i < n; i++) {
      double *Li = &L[(size_t)i * n];
      double s = A[(size_t)i * n + j];
      for (int k = 0; k < j; k++) s -= Li[k] * Lj[k];
      Li[j] = s / ljj;
    }
  }
  return -1;
}

namespace {

// Lower-triangular solves with a row-major m x m factor.
void Forward(const std::vector<double> &L, int m, double *y) {
  for (int i = 0; i < m; i++) {
    const double *Li = &L[(size_t)i * m];
    double s = y[i];
    for (int k = 0; k < i; k++) s -= Li[k] * y[k];
    y[i] = s / Li[i];
  }
}

void Backward(const std::vector<double> &L, int m, double *y) {
  for (int i = m - 1; i >= 0; i--) {
    const double *Li = &L[(size_t)i * m];
    y[i] /= Li[i];
    const double yi = y[i];
    for (int k = 0; k < i; k++) y[k] -= Li[k] * yi;
  }
}

int Find(std::vector<int> &parent, int a) {
  while (parent[a] != a) a = parent[a] = parent[parent[a]];
  return a;
}

}  // namespace

int Group::SetUp(const File &f) {
  // Member position of every file row (-1: not fitted).
  std::vector<int> pos(f.N, -1);
  for (int a = 0; a < n; a++) pos[rows[a]] = a;
  // Connected components of B on the members: blocks link their fitted rows.
  std::vector<int> parent(n);
  for (int a = 0; a < n; a++) parent[a] = a;
  for (const Block &b : f.blocks) {
    int first = -1;
    for (int r : b.rows) {
      const int a = pos[r];
      if (a < 0) continue;
      if (first < 0) {
        first = a;
        continue;
      }
      const int x = Find(parent, first), y = Find(parent, a);
      if (x != y) parent[std::max(x, y)] = std::min(x, y);
    }
  }
  std::vector<int> compOf(n, -1), where(n, -1);
  comps_.clear();
  for (int a = 0; a < n; a++) {
    const int root = Find(parent, a);
    if (compOf[root] < 0) {
      compOf[root] = (int)comps_.size();
      comps_.emplace_back();
    }
    compOf[a] = compOf[root];
    where[a] = (int)comps_[compOf[a]].idx.size();
    comps_[compOf[a]].idx.push_back(a);
  }
  // B' = D^-1 B D^-1 per component, and its Cholesky factor.
  std::vector<std::vector<double>> A(comps_.size());
  for (size_t c = 0; c < comps_.size(); c++) {
    const int m = (int)comps_[c].idx.size();
    A[c].assign((size_t)m * m, 0.0);
    for (int q = 0; q < m; q++) {
      const int a = comps_[c].idx[q];
      A[c][(size_t)q * m + q] = f.diag[rows[a]] / (sigma[a] * sigma[a]);
    }
  }
  for (const Block &b : f.blocks) {
    const int bm = (int)b.rows.size();
    for (int i = 0; i < bm; i++) {
      const int a = pos[b.rows[i]];
      if (a < 0) continue;
      std::vector<double> &Ac = A[compOf[a]];
      const int m = (int)comps_[compOf[a]].idx.size();
      for (int j = 0; j < bm; j++) {
        const int e = pos[b.rows[j]];
        if (e < 0) continue;
        Ac[(size_t)where[a] * m + where[e]] += b.values[(size_t)i * bm + j] / (sigma[a] * sigma[e]);
      }
    }
  }
  int bad = -1;
#pragma omp parallel for schedule(dynamic, 1)
  for (int c = 0; c < (int)comps_.size(); c++) {
    const int fail = Cholesky(A[c], (int)comps_[c].idx.size(), comps_[c].L);
    if (fail >= 0) {
#pragma omp critical
      bad = comps_[c].idx[fail];
    }
    std::vector<double>().swap(A[c]);
  }
  if (bad >= 0) return bad;
  // The rank-1 terms on the members.
  absolute.clear();
  for (const std::vector<double> &v : f.vectors) {
    std::vector<double> sub(n);
    bool any = false;
    for (int a = 0; a < n; a++) any |= (sub[a] = v[rows[a]]) != 0.0;
    if (any) absolute.push_back(std::move(sub));
  }
  rel.clear();
  for (const std::vector<double> &fr : f.rel) {
    std::vector<double> sub(n);
    bool any = false;
    for (int a = 0; a < n; a++) any |= (sub[a] = fr[rows[a]]) != 0.0;
    if (any) rel.push_back(std::move(sub));
  }
  return -1;
}

int Group::LargestComponent() const {
  size_t m = 0;
  for (const Component &c : comps_) m = std::max(m, c.idx.size());
  return (int)m;
}

void Group::SolveB(double *y) const {
  std::vector<double> t;
  for (const Component &c : comps_) {
    const int m = (int)c.idx.size();
    if (m == 1) {
      y[c.idx[0]] /= c.L[0];
      continue;
    }
    t.resize(m);
    for (int q = 0; q < m; q++) t[q] = y[c.idx[q]];
    Forward(c.L, m, t.data());
    for (int q = 0; q < m; q++) y[c.idx[q]] = t[q];
  }
}

void Group::SolveBT(double *y) const {
  std::vector<double> t;
  for (const Component &c : comps_) {
    const int m = (int)c.idx.size();
    if (m == 1) {
      y[c.idx[0]] /= c.L[0];
      continue;
    }
    t.resize(m);
    for (int q = 0; q < m; q++) t[q] = y[c.idx[q]];
    Backward(c.L, m, t.data());
    for (int q = 0; q < m; q++) y[c.idx[q]] = t[q];
  }
}

void Group::LowRank(double *y) const {
  if (k_ == 0) return;
  std::vector<double> p(k_, 0.0), q(k_, 0.0);
  for (int a = 0; a < n; a++) {
    const double *Va = &V_[(size_t)a * k_];
    for (int j = 0; j < k_; j++) p[j] += Va[j] * y[a];
  }
  for (int i = 0; i < k_; i++)
    for (int j = 0; j < k_; j++) q[i] += T_[(size_t)i * k_ + j] * p[j];
  for (int a = 0; a < n; a++) {
    const double *Va = &V_[(size_t)a * k_];
    double s = 0.0;
    for (int j = 0; j < k_; j++) s += Va[j] * q[j];
    y[a] += s;
  }
}

void Group::Whiten(double *r) const {
  SolveB(r);
  LowRank(r);
}

void Group::WhitenTranspose(double *w) const {
  LowRank(w);
  SolveBT(w);
}

void Group::Factor(const std::vector<double> &m) {
  // Columns u'_j = D^-1 u_j: the absolute terms, then the relvector terms on the model.
  k_ = (int)(absolute.size() + rel.size());
  V_.assign((size_t)n * k_, 0.0);
  std::vector<double> col(n);
  for (int j = 0; j < k_; j++) {
    for (int a = 0; a < n; a++) {
      const double u = j < (int)absolute.size() ? absolute[j][a] : rel[j - absolute.size()][a] * m[a];
      col[a] = u / sigma[a];
    }
    SolveB(col.data());
    for (int a = 0; a < n; a++) V_[(size_t)a * k_ + j] = col[a];
  }
  T_.assign((size_t)k_ * k_, 0.0);
  if (k_ > 0) {
    // V^T V = P diag(s) P^T;  T = P diag((1/sqrt(1+s) - 1)/s) P^T, so that
    // (I + V T V^T)^2 = (I + V V^T)^-1 (the factor is finite as s -> 0: -1/2).
    gsl_matrix *G = gsl_matrix_calloc(k_, k_);
    for (int a = 0; a < n; a++) {
      const double *Va = &V_[(size_t)a * k_];
      for (int i = 0; i < k_; i++)
        for (int j = 0; j <= i; j++) *gsl_matrix_ptr(G, i, j) += Va[i] * Va[j];
    }
    for (int i = 0; i < k_; i++)
      for (int j = 0; j < i; j++) gsl_matrix_set(G, j, i, gsl_matrix_get(G, i, j));
    gsl_vector *s = gsl_vector_alloc(k_);
    gsl_matrix *P = gsl_matrix_alloc(k_, k_);
    gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(k_);
    gsl_eigen_symmv(G, s, P, ws);
    gsl_eigen_symmv_free(ws);
    for (int e = 0; e < k_; e++) {
      const double x = std::max(0.0, gsl_vector_get(s, e));
      const double c = x < 1.0e-8 ? -0.5 + 0.375 * x : (1.0 / std::sqrt(1.0 + x) - 1.0) / x;
      for (int i = 0; i < k_; i++)
        for (int j = 0; j < k_; j++)
          T_[(size_t)i * k_ + j] += gsl_matrix_get(P, i, e) * c * gsl_matrix_get(P, j, e);
    }
    gsl_matrix_free(P);
    gsl_vector_free(s);
    gsl_matrix_free(G);
  }
  m_ = m;
  factored_ = true;
  WriteModel();
}

void Group::EnsureFactored(const std::vector<double> &r) {
  if (factored_) return;
  std::lock_guard<std::mutex> lk(mutex_);
  if (factored_) return;
  std::vector<double> m(n);
  for (int a = 0; a < n; a++) m[a] = data[a] + sigma[a] * r[a];
  Factor(m);
}

void Group::Refreeze() {
  std::lock_guard<std::mutex> lk(mutex_);
  if (HasRelative()) factored_ = false;
}

void Group::WriteModel() const {
  if (!HasRelative() || modelFile.empty()) return;
  std::ofstream out(modelFile.c_str(), std::ios::app);
  if (!out) return;
  out.precision(12);
  out << "# covariance " << path << ": model m frozen for the relvector terms (lab units of the data)\n"
      << "# row  segmentKey  labEnergy  angle  data  m  sigma\n";
  for (int a = 0; a < n; a++)
    out << rows[a] + 1 << " " << label[a] << " " << data[a] << " " << m_[a] << " " << sigma[a] << "\n";
}

}  // namespace DataCovariance
