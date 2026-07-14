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
  std::vector<ParamDesc> cols;
  std::vector<std::vector<double> > M;
  bool empty() const { return cols.empty() || M.empty(); }
  int size() const { return (int)cols.size(); }
};

/// Densify Minuit's packed lower-triangular covariance (Data()) into a
/// BandCovariance, using pmap for the column identities.  Empty on size mismatch.
BandCovariance BuildBandCovarianceFromMinuit(const std::vector<double>& covData,
                                             const ParamIndexMap& pmap);

/// Write/read a BandCovariance to/from a small text file.  False on I/O error.
bool SaveBandCovariance(const std::string& path, const BandCovariance& cov);
bool LoadBandCovariance(const std::string& path, BandCovariance& cov);

/// Project a saved covariance onto pmap's packed columns (matched by identity).
/// Unmatched entries stay zero, which is correct for bands since dT/dn = 0.
std::vector<std::vector<double> > RemapCovarianceToParamMap(const BandCovariance& saved,
                                                            const ParamIndexMap& pmap);

/// Inputs for the cross-section band: covariance M and per-point sensitivities
/// grad[point], both in the current packed order.  dXS(g) = sqrt(g^T M g)
/// (SAMMY Eq. IV E4.2); sum gradient rows first to combine correlated points.
struct BandData {
  CNuc* compound;
  std::vector<std::vector<double> > M;
  std::map<EPoint*, std::vector<double> > grad;
  double dXS(const std::vector<double>& g) const;
};

/// Build a BandData from savedCov for the current (best-fit-filled) compound/data.
/// False if the analytic sensitivities are unavailable for this model.
bool BuildBandData(CNuc* compound, EData* data, const Config& config,
                   const BandCovariance& savedCov, BandData& out);

#endif
