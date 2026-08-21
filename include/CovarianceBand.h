#ifndef COVARIANCEBAND_H
#define COVARIANCEBAND_H

#include "AZUREGrad.h"
#include <string>
#include <vector>
#include <map>

class CNuc;
class EData;
class EPoint;
class Config;

/// Fitted-parameter covariance, with each column tagged by its parameter
/// identity (cols[a]) so a covariance saved from a fit can be matched onto a
/// later run whose parameters are ordered or fixed differently.  M is dense
/// (nPacked x nPacked) in ParamIndexMap packed order.
struct BandCovariance {
  std::vector<ParamDesc> cols;  ///< empty when loaded from a plain covariance.dat
  std::vector<std::vector<double>> M;
  bool empty() const { return M.empty(); }
  int size() const { return (int)M.size(); }
};

/// Packed-column indices of a pmap that correspond to R-matrix parameters (level
/// energies and reduced widths) -- the only parameters a cross-section band is
/// sensitive to (dT/dn = dT/dshift = 0).  Normalizations and energy shifts are
/// excluded, so the band covariance and covariance.dat span just these columns.
std::vector<int> RMatrixPackedColumns(const ParamIndexMap &pmap);

/// Densify the R-matrix sub-block of Minuit's packed lower-triangular covariance
/// (Data()) into a BandCovariance, using pmap for the column identities.  Norm
/// and energy-shift columns are dropped.  Empty on size mismatch.
BandCovariance BuildBandCovarianceFromMinuit(const std::vector<double> &covData,
                                             const ParamIndexMap &pmap);

/// Write/read a BandCovariance to/from a plain text file: the bare N x N matrix,
/// one row per line, no header or column metadata.  False on I/O error.  A
/// loaded covariance therefore has no column identities (cov.cols is empty), so
/// its columns are taken to be in the current packed parameter order.
bool SaveBandCovariance(const std::string &path, const BandCovariance &cov);
/// Read a covariance matrix from file; false if it cannot be read.
bool LoadBandCovariance(const std::string &path, BandCovariance &cov);

/// Project a saved covariance onto pmap's R-matrix columns (RMatrixPackedColumns
/// order).  If the saved covariance carries column identities (same run) they are
/// matched by identity, leaving unmatched entries zero.  If it was loaded from
/// covariance.dat (no identities) its columns are assumed already in R-matrix
/// order and are copied through when the dimension matches the R-matrix count.
std::vector<std::vector<double>> RemapCovarianceToParamMap(const BandCovariance &saved,
                                                           const ParamIndexMap &pmap);

/// Inputs for the cross-section band: covariance M and per-point sensitivities
/// grad[point], both in R-matrix column order (RMatrixPackedColumns; norms and
/// energy shifts excluded).  dXS(g) = sqrt(g^T M g) (SAMMY Eq. IV E4.2); sum
/// gradient rows first to combine correlated points.
struct BandData {
  CNuc *compound;
  std::vector<std::vector<double>> M;
  std::map<EPoint *, std::vector<double>> grad;
  /// One-sigma cross-section uncertainty from a gradient vector, \f$\sqrt{g^T M g}\f$.
  double dXS(const std::vector<double> &g) const;
};

/// Build a BandData from savedCov for the current (best-fit-filled) compound/data.
/// False if the analytic sensitivities are unavailable for this model.
bool BuildBandData(CNuc *compound, EData *data, const Config &config,
                   const BandCovariance &savedCov, BandData &out);

#endif
