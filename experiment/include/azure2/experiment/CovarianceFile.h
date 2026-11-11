#ifndef AZURE2_EXPERIMENT_COVARIANCEFILE_H
#define AZURE2_EXPERIMENT_COVARIANCEFILE_H

#include <functional>
#include <istream>
#include <memory>
#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * A data covariance file (a data set's `covariance <file>`; the grammar and the meaning of the matrix are documented
 * with the engine's use of it, src/DataCovariance.h): the data files it covers with their numbers of data rows, and
 * the terms of the matrix in the units of the files' value column.  The one reader of the format: the engine's set-up
 * (DataCovariance::Read) and the build checks of the validation (Schema::Build) call it.
 */

/// A symmetric block of the matrix on a set of (0-based) rows, m x m row-major.
struct CovarianceBlock {
  std::vector<int> rows;
  std::vector<double> values;
};

/// A parsed covariance file: its row layout and the terms.
struct CovarianceFile {
  std::string path;                                 ///< as written on the data set
  std::vector<std::string> files;                   ///< data files, as written in the header
  std::vector<std::string> canonical;               ///< the same as opened, CanonicalPath
  std::vector<int> counts;                          ///< data rows of each file
  std::vector<int> offsets;                         ///< 0-based row of each file's first record
  int N = 0;                                        ///< total rows
  std::vector<double> diag;                         ///< the diagonal terms (N)
  std::vector<CovarianceBlock> blocks;              ///< the block and matrix terms
  std::vector<std::vector<double>> vectors;         ///< absolute rank-1 terms, each of size N
  std::vector<std::vector<double>> rel;             ///< relvector fractions, each of size N
  std::vector<double> diagB;                        ///< diagonal of B = diag + blocks (N)
  /// C_ii with \p model standing in for the model in the relvector terms (the units of the value column).
  double Variance(int row, double model) const;
  /// 0-based row of record \p record (1-based) of data file \p canonicalFile; -1 if none.
  int RowOf(const std::string &canonicalFile, int record) const;
  /// Index of \p canonicalFile in the header, or -1.
  int FileIndex(const std::string &canonicalFile) const;
  // Bookkeeping for the engine's set-up report.
  int numSigmaMismatch = 0;                         ///< points whose error column differs > 1 % from sqrt(C_ii)
  double worstSigmaMismatch = 0.0;
  std::string firstSigmaMismatch;
};

/// Canonical form of a path for matching data files (weakly_canonical, as a string).
std::string CanonicalPath(const std::string &path);

/// How the reader reaches the files: the covariance file and the data files its header names.
struct CovarianceFileAccess {
  /// The path a name is opened at (empty: as written, relative to the working directory, as the engine opens it).
  std::function<std::string(const std::string &)> resolve;
  /// Opens a resolved path (empty: std::ifstream); a null or failed stream is a file that cannot be read.
  std::function<std::unique_ptr<std::istream>(const std::string &)> open;
};

/// Read the covariance file \p path.  On failure returns null and writes the reason to \p err.
std::shared_ptr<CovarianceFile> ReadCovarianceFile(const std::string &path, const CovarianceFileAccess &access,
                                                   std::string &err);

}  // namespace experiment
}  // namespace azure2

#endif
