#include "CovarianceBand.h"
#include "CNuc.h"
#include "EData.h"
#include "AZUREParams.h"
#include "Config.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <tuple>

namespace {

// Lower-triangular packed index used by ROOT::Minuit2::MnUserCovariance.
inline int TriIndex(int i, int j) {
  return (i <= j) ? (j * (j + 1)) / 2 + i : (i * (i + 1)) / 2 + j;
}

// A parameter's identity, used to match columns across runs.  R-matrix columns
// are keyed on (kind, jGroup, level, channel); norm/shift columns on
// (kind, segment) with the unused slots set to -1.
std::tuple<int,int,int,int> IdentityKey(const ParamDesc& d) {
  int kind = (int)d.kind;
  switch(d.kind) {
    case ParamKind::LevelEnergy: return std::make_tuple(kind, d.jGroup, d.level, -1);
    case ParamKind::Gamma:       return std::make_tuple(kind, d.jGroup, d.level, d.channel);
    case ParamKind::Norm:        return std::make_tuple(kind, d.segment, -1, -1);
    case ParamKind::EnergyShift: return std::make_tuple(kind, d.segment, -1, -1);
  }
  return std::make_tuple(kind, -1, -1, -1);
}

}  // namespace

BandCovariance BuildBandCovarianceFromMinuit(const std::vector<double>& covData,
                                             const ParamIndexMap& pmap) {
  BandCovariance cov;
  const int n = pmap.NumPacked();
  if(n <= 0) return cov;
  // Minuit stores the lower triangle: n(n+1)/2 entries for n variable params.
  if((int)covData.size() != n * (n + 1) / 2) return cov;

  cov.cols.resize(n);
  for(int a = 0; a < n; a++) cov.cols[a] = pmap.Desc(pmap.PackedToFull(a));

  cov.M.assign(n, std::vector<double>(n, 0.0));
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      cov.M[a][b] = covData[TriIndex(a, b)];
  return cov;
}

bool SaveBandCovariance(const std::string& path, const BandCovariance& cov) {
  std::ofstream out(path.c_str());
  if(!out) return false;
  const int n = cov.size();
  out << "# AZURE2 parameter covariance for cross-section bands\n";
  out << "# columns: kind jGroup level channel segment\n";
  out << n << "\n";
  for(int a = 0; a < n; a++) {
    const ParamDesc& d = cov.cols[a];
    out << (int)d.kind << " " << d.jGroup << " " << d.level << " "
        << d.channel << " " << d.segment << "\n";
  }
  out.precision(17);
  out << std::scientific;
  for(int a = 0; a < n; a++) {
    for(int b = 0; b < n; b++) out << cov.M[a][b] << (b + 1 < n ? " " : "");
    out << "\n";
  }
  return (bool)out;
}

bool LoadBandCovariance(const std::string& path, BandCovariance& cov) {
  std::ifstream in(path.c_str());
  if(!in) return false;
  cov.cols.clear();
  cov.M.clear();

  // Skip comment lines starting with '#'.
  std::string line;
  int n = -1;
  while(std::getline(in, line)) {
    if(line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    if(iss >> n) break;
  }
  if(n <= 0) return false;

  cov.cols.resize(n);
  for(int a = 0; a < n; a++) {
    int kind, j, level, channel, segment;
    if(!(in >> kind >> j >> level >> channel >> segment)) return false;
    ParamDesc d;
    d.kind = (ParamKind)kind;
    d.jGroup = j; d.level = level; d.channel = channel; d.segment = segment;
    d.fixed = false;
    cov.cols[a] = d;
  }

  cov.M.assign(n, std::vector<double>(n, 0.0));
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      if(!(in >> cov.M[a][b])) return false;
  return true;
}

std::vector<std::vector<double> > RemapCovarianceToParamMap(const BandCovariance& saved,
                                                            const ParamIndexMap& pmap) {
  const int n = pmap.NumPacked();
  std::vector<std::vector<double> > M(n, std::vector<double>(n, 0.0));
  if(saved.empty() || n <= 0) return M;

  // Map each stored column's identity to its stored index.
  std::map<std::tuple<int,int,int,int>, int> savedIndex;
  for(int a = 0; a < saved.size(); a++) savedIndex[IdentityKey(saved.cols[a])] = a;

  // For each current packed column, find the matching stored column (if any).
  std::vector<int> curToSaved(n, -1);
  for(int a = 0; a < n; a++) {
    auto it = savedIndex.find(IdentityKey(pmap.Desc(pmap.PackedToFull(a))));
    if(it != savedIndex.end()) curToSaved[a] = it->second;
  }

  for(int a = 0; a < n; a++) {
    if(curToSaved[a] < 0) continue;
    for(int b = 0; b < n; b++) {
      if(curToSaved[b] < 0) continue;
      M[a][b] = saved.M[curToSaved[a]][curToSaved[b]];
    }
  }
  return M;
}

double BandData::dXS(const std::vector<double>& g) const {
  const int n = (int)M.size();
  if((int)g.size() != n) return 0.0;
  // v = g^T M g = sum_a g_a * (sum_b M_ab g_b).  Skip zero rows (norm/shift
  // columns and any parameter the point is insensitive to) for speed.
  double v = 0.0;
  for(int a = 0; a < n; a++) {
    if(g[a] == 0.0) continue;
    double Mg = 0.0;
    const std::vector<double>& Ma = M[a];
    for(int b = 0; b < n; b++) Mg += Ma[b] * g[b];
    v += g[a] * Mg;
  }
  return (v > 0.0) ? std::sqrt(v) : 0.0;
}

bool BuildBandData(CNuc* compound, EData* data, const Config& config,
                   const BandCovariance& savedCov, BandData& out) {
  out.compound = compound;
  out.M.clear();
  out.grad.clear();
  if(savedCov.empty()) return false;

  // Parameter-index map for the current run, using the same fixed mask (in
  // Minuit order) as the fitter -- mirrors AZURECalc::ResidualJacobian.
  AZUREParams tp;
  compound->FillMnParams(tp.GetMinuitParams(), &config);
  data->FillMnParams(tp.GetMinuitParams());
  const int nMn = tp.GetMinuitParams().Params().size();
  std::vector<bool> fixed(nMn);
  for(int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();

  ParamIndexMap pmap = BuildParamIndexMap(compound, data, fixed);

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(config.paramMask & Config::USE_BRUNE_FORMALISM) {
    compound->CalcShiftFunctions(config);
    shiftDeriv = BuildShiftDerivTable(compound, config);
    sdp = &shiftDeriv;
  }

  out.M = RemapCovarianceToParamMap(savedCov, pmap);
  return ComputeModelGradients(compound, data, config, pmap, sdp, out.grad);
}
