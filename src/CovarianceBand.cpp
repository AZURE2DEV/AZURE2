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
// (kind, segment), THM coherent-background ones on (kind, index), with the
// unused slots set to -1.
std::tuple<int, int, int, int> IdentityKey(const ParamDesc &d) {
  int kind = (int)d.kind;
  switch (d.kind) {
    case ParamKind::LevelEnergy: return std::make_tuple(kind, d.jGroup, d.level, -1);
    case ParamKind::Gamma: return std::make_tuple(kind, d.jGroup, d.level, d.channel);
    case ParamKind::Norm: return std::make_tuple(kind, d.segment, -1, -1);
    case ParamKind::EnergyShift: return std::make_tuple(kind, d.segment, -1, -1);
    case ParamKind::ThmCoherent: return std::make_tuple(kind, d.channel, -1, -1);  // channel = its index
  }
  return std::make_tuple(kind, -1, -1, -1);
}

}  // namespace

std::vector<int> BandPackedColumns(const ParamIndexMap &pmap) {
  std::vector<int> cols;
  const int n = pmap.NumPacked();
  for (int a = 0; a < n; a++) {
    ParamKind k = pmap.Desc(pmap.PackedToFull(a)).kind;
    if (k == ParamKind::LevelEnergy || k == ParamKind::Gamma || k == ParamKind::ThmCoherent) cols.push_back(a);
  }
  return cols;
}

BandCovariance BuildBandCovarianceFromMinuit(const std::vector<double> &covData,
                                             const ParamIndexMap &pmap) {
  BandCovariance cov;
  const int n = pmap.NumPacked();
  if (n <= 0) return cov;
  // Minuit stores the lower triangle: n(n+1)/2 entries for n variable params.
  if ((int)covData.size() != n * (n + 1) / 2) return cov;

  // Keep only the band columns (R-matrix, and THM coherent background if any):
  // the band is insensitive to norms and
  // energy shifts, so their rows/columns never contribute to dXS.  Dropping them
  // makes the saved matrix and the free-parameter count reflect the R-matrix
  // parameters alone (stable across fit and extrapolation runs).
  const std::vector<int> rc = BandPackedColumns(pmap);
  const int m = (int)rc.size();
  cov.cols.resize(m);
  for (int a = 0; a < m; a++) cov.cols[a] = pmap.Desc(pmap.PackedToFull(rc[a]));

  cov.M.assign(m, std::vector<double>(m, 0.0));
  for (int a = 0; a < m; a++)
    for (int b = 0; b < m; b++)
      cov.M[a][b] = covData[TriIndex(rc[a], rc[b])];
  return cov;
}

bool SaveBandCovariance(const std::string &path, const BandCovariance &cov) {
  std::ofstream out(path.c_str());
  if (!out) return false;
  // Bare N x N matrix, one row per line -- no header, comments, or column
  // metadata, so the file is trivially readable by external tools.
  const int n = cov.size();
  out.precision(17);
  out << std::scientific;
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) out << cov.M[a][b] << (b + 1 < n ? " " : "");
    out << "\n";
  }
  return (bool)out;
}

bool LoadBandCovariance(const std::string &path, BandCovariance &cov) {
  std::ifstream in(path.c_str());
  if (!in) return false;
  cov.cols.clear();  // plain matrix carries no column identities
  cov.M.clear();

  // Read every non-empty line as a row of doubles.  The file must be square.
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::vector<double> row;
    double v;
    while (iss >> v) row.push_back(v);
    if (!row.empty()) cov.M.push_back(row);
  }
  const int n = (int)cov.M.size();
  if (n <= 0) return false;
  for (int a = 0; a < n; a++)
    if ((int)cov.M[a].size() != n) {
      cov.M.clear();
      return false;
    }
  return true;
}

std::vector<std::vector<double>> RemapCovarianceToParamMap(const BandCovariance &saved,
                                                           const ParamIndexMap &pmap) {
  // Target columns are the current run's free band parameters (BandPackedColumns).
  const std::vector<int> rc = BandPackedColumns(pmap);
  const int m = (int)rc.size();
  std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0));
  if (saved.empty() || m <= 0) return M;

  // A covariance loaded from covariance.dat has no column identities, so its
  // columns are taken to be in BandPackedColumns order.  Copy it through when the
  // dimension matches the R-matrix count; a mismatch is reported and rejected by
  // BuildBandData, so here it simply leaves M zero.
  if (saved.cols.empty()) {
    if (saved.size() == m) M = saved.M;
    return M;
  }

  // Map each stored column's identity to its stored index.
  std::map<std::tuple<int, int, int, int>, int> savedIndex;
  for (int a = 0; a < saved.size(); a++) savedIndex[IdentityKey(saved.cols[a])] = a;

  // For each current R-matrix column, find the matching stored column (if any).
  std::vector<int> curToSaved(m, -1);
  for (int a = 0; a < m; a++) {
    auto it = savedIndex.find(IdentityKey(pmap.Desc(pmap.PackedToFull(rc[a]))));
    if (it != savedIndex.end()) curToSaved[a] = it->second;
  }

  for (int a = 0; a < m; a++) {
    if (curToSaved[a] < 0) continue;
    for (int b = 0; b < m; b++) {
      if (curToSaved[b] < 0) continue;
      M[a][b] = saved.M[curToSaved[a]][curToSaved[b]];
    }
  }
  return M;
}

double BandData::dXS(const std::vector<double> &g) const {
  const int n = (int)M.size();
  if ((int)g.size() != n) return 0.0;
  // v = g^T M g = sum_a g_a * (sum_b M_ab g_b).  Skip zero rows (norm/shift
  // columns and any parameter the point is insensitive to) for speed.
  double v = 0.0;
  for (int a = 0; a < n; a++) {
    if (g[a] == 0.0) continue;
    double Mg = 0.0;
    const std::vector<double> &Ma = M[a];
    for (int b = 0; b < n; b++) Mg += Ma[b] * g[b];
    v += g[a] * Mg;
  }
  return (v > 0.0) ? std::sqrt(v) : 0.0;
}

bool BuildBandData(CNuc *compound, EData *data, const Config &config,
                   const BandCovariance &savedCov, BandData &out,
                   const std::vector<bool> *fixedMask) {
  out.compound = compound;
  out.M.clear();
  out.grad.clear();
  if (savedCov.empty()) return false;

  // Parameter-index map for the current run, using the same fixed mask (in
  // Minuit order) as the fitter: the caller's (param.par included), else the
  // .azr's -- mirrors AZURECalc::ResidualJacobian.
  AZUREParams tp;
  compound->FillMnParams(tp.GetMinuitParams(), &config);
  data->FillMnParams(tp.GetMinuitParams());
  const int nMn = tp.GetMinuitParams().Params().size();
  std::vector<bool> fixed(nMn);
  for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();
  if (fixedMask && (int)fixedMask->size() == nMn) fixed = *fixedMask;

  ParamIndexMap pmap = BuildParamIndexMap(compound, data, fixed);

  // The band spans only the free R-matrix parameters (level energies and reduced
  // widths) and THM coherent-background ones; norms and energy shifts have zero
  // sensitivity and are excluded.
  const std::vector<int> rc = BandPackedColumns(pmap);
  const int m = (int)rc.size();

  // A covariance loaded from covariance.dat has no column identities; its
  // dimension must then equal the number of free band parameters.
  if (savedCov.cols.empty() && savedCov.size() != m) {
    const bool coherent = data->NumThmCoherentParams() > 0;
    config.outStream << "Data covariance size mismatch: covariance.dat has "
                     << savedCov.size() << " rows/columns but the model has " << m
                     << (coherent ? " free R-matrix and THM coherent-background parameters (level energies, "
                                    "reduced widths and cbkg values; "
                                  : " free R-matrix parameters (level energies and reduced widths; ")
                     << "normalizations and energy shifts are excluded). "
                        "Skipping the uncertainty band."
                     << std::endl;
    return false;
  }

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (config.paramMask & Config::USE_BRUNE_FORMALISM) {
    compound->CalcShiftFunctions(config);
    shiftDeriv = BuildShiftDerivTable(compound, config);
    sdp = &shiftDeriv;
  }

  out.M = RemapCovarianceToParamMap(savedCov, pmap);  // m x m, BandPackedColumns order

  // ComputeModelGradients yields full packed rows; reduce each to the band
  // columns so the sensitivities line up with the covariance.
  std::map<EPoint *, vector_r> fullGrad;
  if (!ComputeModelGradients(compound, data, config, pmap, sdp, fullGrad, /*skipTHM=*/true))
    return false;

  // THM (HOES) points: the adjoint differentiates the T-matrix observable, not
  // the HOES one, so their rows come from ComputeTHMRows (central differences of
  // the HOES model).  What the output file shows for a THM point is the model m
  // next to the data scaled by the profiled norm, d n*; the parameters move
  // both.  The band is therefore that of the model as it lies against the data,
  // q = s(p) m(p) / s(p0), s = 1/n*, whose gradient at the best fit is
  //   dq/dp = J_m + m (ds/dp)/s
  // (ComputeTHMRows' J_m and ds/dp; a fixed THM norm has ds/dp = 0 and the band
  // is that of m).  The overall scale of the THM model, arbitrary in THM, thus
  // gets no band: a change of the parameters that only rescales m leaves q put.
  {
    const vector_r full = tp.GetMinuitParams().Params();
    std::vector<THMRows> thm = ComputeTHMRows(compound, data, config, full, pmap);
    const int nCols = pmap.NumPacked();
    for (const THMRows &tr : thm) {
      ESegment *seg = data->GetSegment(tr.segment);
      size_t i = 0;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        EPoint *pt = seg->GetPoint(pid + 1);
        if (!pt) continue;
        vector_r &row = fullGrad[pt];
        row.assign(nCols, 0.0);
        for (int c = 0; c < nCols; c++) {
          double g = tr.Jm[i * nCols + c];
          if (tr.profiled && tr.s != 0.0) g += tr.m[i] * tr.ds[c] / tr.s;
          // A THM experiment: through its shared norm and background.
          if (!tr.G.empty()) g = tr.G[i * nCols + c];
          row[c] = g;
        }
        i++;
      }
    }
  }
  for (std::map<EPoint *, vector_r>::const_iterator it = fullGrad.begin();
       it != fullGrad.end(); ++it) {
    const vector_r &f = it->second;
    vector_r g(m, 0.0);
    for (int k = 0; k < m; k++)
      if (rc[k] < (int)f.size()) g[k] = f[rc[k]];
    out.grad[it->first] = g;
  }
  return true;
}
