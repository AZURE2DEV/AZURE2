#ifndef DATACOVARIANCE_H
#define DATACOVARIANCE_H

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

/*!
 * Data covariance matrices (segment keyword "covariance <file>").
 *
 * A covariance file names the data files it covers and their number of data rows,
 * then gives the matrix as a sum of terms, in the units of the files' value column:
 *
 *     file data/a.dat 700        # rows 1..700 are the data rows of data/a.dat
 *     file data/b.dat 300        # rows 701..1000
 *     block n  i1 .. in  c11 c12 .. cnn   # n x n symmetric block on rows i1..in (1-based)
 *     diagonal v1 .. vN                     # adds diag(v)
 *     vector   v1 .. vN                     # adds v v^T
 *     relvector f1 .. fN                    # adds (f m)(f m)^T, m = the model (see below)
 *     matrix   c11 .. cNN                   # adds a full N x N symmetric matrix
 *
 * Numbers are free-form (any whitespace, line breaks anywhere); '#' starts a comment.
 * A data row is a line of a data file that is neither blank nor a comment, counted from 1.
 * AZURE2 uses the rows of the points actually fitted (the submatrix: points cut away by
 * the energy or angle limits, or in inactive segments, simply drop out).
 *
 * The group's chi2 is r^T M^-1 r over AZURE2's usual normalized residuals
 * r_i = (fit_i - n d_i) / (n sigma_i), with M = D^-1 C D^-1 and D = diag(sigma_i) the
 * points' (lab) uncertainties, so normalizations and the lab-to-c.m. conversion act on
 * it exactly as on the error column.  The value does not depend on D (any D > 0 gives
 * v^T C^-1 v, v = model - data in the file's units); AZURE2 sets sigma_i = sqrt(C_ii)
 * when it reads the data, with the data value standing in for m in the relvector terms.
 *
 * The matrix is never formed.  C = B + U U^T: B = the diagonal, block and matrix terms
 * (block diagonal once the blocks are merged into connected components, each factored on
 * its own, B = L_B L_B^T), U = the k vector and relvector terms.  With V = L_B^-1 U and the
 * k x k matrix T = P diag((1/sqrt(1+s) - 1)/s) P^T from V^T V = P diag(s) P^T, the
 * whitening W = (I + V T V^T) L_B^-1 has W^T W = C^-1, so w = W v gives chi2 = w.w and
 * C^-1 v = W^T w, at O(n k) per evaluation.  B must be positive definite by itself: every
 * fitted row needs a diagonal or block variance of its own.
 *
 * relvector terms (multiplicative uncertainties, free of Peelle's pertinent puzzle) use
 * m = the model in the data's units, m_i = fit_i / (n conv_i) = d_i + sigma_i r_i, FROZEN
 * at the first evaluation of the run (the starting parameters) so that C stays constant
 * through a minimization; Refreeze() re-freezes it at the next evaluation (iterative
 * refits).  The m used is written to covarianceModel.out.
 */
namespace DataCovariance {

/// Canonical form of a path for matching data files (weakly_canonical, as a string).
std::string CanonicalPath(const std::string &path);

/// A symmetric block of the matrix on a set of (0-based) rows, m x m row-major.
struct Block {
  std::vector<int> rows;
  std::vector<double> values;
};

/// A parsed covariance file: its row layout and the terms.
struct File {
  std::string path;                                 ///< as written on the segment line
  std::vector<std::string> files;                   ///< data files, as written in the header
  std::vector<std::string> canonical;               ///< the same, CanonicalPath
  std::vector<int> counts;                          ///< data rows of each file
  std::vector<int> offsets;                         ///< 0-based row of each file's first record
  int N = 0;                                        ///< total rows
  std::vector<double> diag;                         ///< the diagonal terms (N)
  std::vector<Block> blocks;                        ///< the block and matrix terms
  std::vector<std::vector<double>> vectors;         ///< absolute rank-1 terms, each of size N
  std::vector<std::vector<double>> rel;             ///< relvector fractions, each of size N
  std::vector<double> diagB;                        ///< diagonal of B = diag + blocks (N)
  /// C_ii with \p model standing in for the model in the relvector terms.
  double Variance(int row, double model) const;
  /// 0-based row of record \p record (1-based) of data file \p canonicalFile; -1 if none.
  int RowOf(const std::string &canonicalFile, int record) const;
  /// Index of \p canonicalFile in the header, or -1.
  int FileIndex(const std::string &canonicalFile) const;
  // Bookkeeping for the setup report.
  int numSigmaMismatch = 0;                         ///< points whose error column differs > 1 % from sqrt(C_ii)
  double worstSigmaMismatch = 0.0;
  std::string firstSigmaMismatch;
};

/// Read a covariance file.  On failure returns null and writes the reason to \p err.
std::shared_ptr<File> Read(const std::string &path, std::string &err);

/// One covariance group: the fitted points of the segments naming one file, in row order,
/// and the whitening W of M = D^-1 C D^-1 (W^T W = M^-1).  Shared (one object) by every
/// copy of the EData, so a freeze of the model is seen by all of them.
struct Group {
  std::string path;                                 ///< the covariance file
  std::vector<std::pair<int, int>> members;         ///< (1-based segment, 0-based point), row order
  std::vector<int> rows;                            ///< the members' 0-based rows in the file
  int n = 0;
  int N = 0;                                        ///< rows in the file
  int numSegments = 0;
  std::vector<double> sigma;                        ///< D: the points' lab uncertainties
  std::vector<double> data;                         ///< the points' lab data values
  std::vector<std::string> label;                   ///< per member: "segKey E angle", for covarianceModel.out
  std::string modelFile;                            ///< where the frozen model is written ("" = nowhere)
  std::vector<std::vector<double>> absolute;        ///< absolute rank-1 terms on the members (lab units)
  std::vector<std::vector<double>> rel;             ///< relvector fractions on the members
  /// Build B' = D^-1 B D^-1 on the members from \p f and factor its components.  Returns -1
  /// on success, else the member (0-based) where B' is not positive definite.
  int SetUp(const File &f);
  int NumComponents() const { return (int)comps_.size(); }
  int LargestComponent() const;
  /// Build the low-rank part for the model m (lab units; ignored without relvector terms).
  void Factor(const std::vector<double> &m);
  /// Factor once, with m = data + sigma r (r = the members' normalized residuals now), if not
  /// yet done; thread-safe.
  void EnsureFactored(const std::vector<double> &r);
  /// Forget the frozen model: the next EnsureFactored freezes it again.
  void Refreeze();
  bool Factored() const { return factored_; }
  bool HasRelative() const { return !rel.empty(); }
  const std::vector<double> &Model() const { return m_; }
  /// w = W r in place (n values).
  void Whiten(double *r) const;
  /// z = W^T w in place, so that W^T W r = M^-1 r.
  void WhitenTranspose(double *w) const;

 private:
  struct Component {
    std::vector<int> idx;                           ///< members, ascending
    std::vector<double> L;                          ///< Cholesky factor of B' on them
  };
  std::vector<Component> comps_;
  int k_ = 0;                                       ///< rank-1 terms
  std::vector<double> V_;                           ///< n x k row-major, column j = L_B^-1 D^-1 u_j
  std::vector<double> T_;                           ///< k x k
  std::vector<double> m_;                           ///< the frozen model (lab units)
  std::atomic<bool> factored_{false};
  std::mutex mutex_;
  void SolveB(double *y) const;                     ///< y = L_B^-1 y
  void SolveBT(double *y) const;                    ///< y = L_B^-T y
  void LowRank(double *y) const;                    ///< y = (I + V T V^T) y
  void WriteModel() const;
};

/// Cholesky factor of the n x n symmetric positive-definite \p A (row-major) into \p L.
/// Returns -1 on success, else the (0-based) pivot that was not positive.
int Cholesky(const std::vector<double> &A, int n, std::vector<double> &L);

}  // namespace DataCovariance

#endif
