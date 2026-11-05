#include "PriorCovariance.h"
#include "Config.h"
#include "InputFiles.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>

namespace {

// strtod over the whole token ("nan" / "inf" are refused: a prior has finite numbers).
bool Number(const std::string &token, double &v) {
  if (token.empty()) return false;
  char *end = nullptr;
  v = std::strtod(token.c_str(), &end);
  return end && *end == '\0' && std::isfinite(v);
}

// In-place Cholesky factor of the symmetric n x n row-major A: A = L L^T, L in the lower triangle (upper zeroed).
// False when a pivot is not positive (or below 1e-13 of its diagonal entry: numerically singular); bad = that row.
bool Cholesky(std::vector<double> &A, int n, int &bad) {
  for (int j = 0; j < n; j++) {
    double d = A[(size_t)j * n + j];
    const double diag = d;
    for (int k = 0; k < j; k++) d -= A[(size_t)j * n + k] * A[(size_t)j * n + k];
    if (!(d > 0.0) || d <= 1.0e-13 * std::fabs(diag)) {
      bad = j;
      return false;
    }
    const double ljj = std::sqrt(d);
    A[(size_t)j * n + j] = ljj;
    for (int i = j + 1; i < n; i++) {
      double s = A[(size_t)i * n + j];
      for (int k = 0; k < j; k++) s -= A[(size_t)i * n + k] * A[(size_t)j * n + k];
      A[(size_t)i * n + j] = s / ljj;
    }
    for (int i = 0; i < j; i++) A[(size_t)i * n + j] = 0.0;
  }
  return true;
}

std::string Trim(const std::string &s) {
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

}  // namespace

std::string PriorCovariance::ResolvePath(const Config &configure) {
  const std::string &p = configure.priorCovariance;
  if (p.empty() || p[0] == '/') return p;
  const size_t slash = configure.configfile.find_last_of('/');
  if (slash == std::string::npos) return p;
  return configure.configfile.substr(0, slash + 1) + p;
}

bool PriorCovariance::Read(const std::string &path, File &out, std::string &err, bool cached) {
  out = File();
  std::unique_ptr<std::istream> inStream = InputFiles::Open(path, cached);
  std::istream &in = *inStream;
  if (!in) {
    err = "cannot open the prior covariance file " + path;
    return false;
  }
  std::vector<std::vector<std::string>> rows;  // non-comment lines, as tokens
  std::vector<int> lineOf;
  std::string raw;
  int lineNumber = 0;
  while (std::getline(in, raw)) {
    lineNumber++;
    const std::string line = Trim(raw);
    if (line.empty()) continue;
    if (line[0] == '#') {
      // "#parametrization <0|1|2>" (param.sav's first line, "#parametrization 1.000000e+00 0.000000e+00" there):
      // only before the names line.
      std::istringstream f(line.substr(1));
      std::string key, value;
      f >> key >> value;
      if (key == "parametrization" && rows.empty()) {
        double v;
        if (!Number(value, v) || (v != 0.0 && v != 1.0 && v != 2.0)) {
          err = path + " line " + std::to_string(lineNumber) + ": the #parametrization tag must be 0, 1 or 2";
          return false;
        }
        out.basis = (int)v;
      }
      continue;
    }
    std::istringstream f(line.substr(0, line.find('#')));
    std::vector<std::string> tokens;
    std::string t;
    while (f >> t) tokens.push_back(t);
    if (tokens.empty()) continue;
    rows.push_back(tokens);
    lineOf.push_back(lineNumber);
  }
  if (rows.size() < 2) {
    err = path + ": expected a line of parameter names, a line of means and the covariance rows";
    return false;
  }
  out.names = rows[0];
  const int n = (int)out.names.size();
  std::map<std::string, int> seen;
  for (const std::string &name : out.names)
    if (seen[name]++) {
      err = path + ": parameter " + name + " is listed twice";
      return false;
    }
  if ((int)rows[1].size() != n) {
    err = path + ": size mismatch, " + std::to_string(n) + " names but " + std::to_string(rows[1].size()) +
          " means (line " + std::to_string(lineOf[1]) + ")";
    return false;
  }
  out.mean.resize(n);
  for (int i = 0; i < n; i++)
    if (!Number(rows[1][i], out.mean[i])) {
      err = path + " line " + std::to_string(lineOf[1]) + ": '" + rows[1][i] + "' is not a number";
      return false;
    }
  if ((int)rows.size() - 2 != n) {
    err = path + ": size mismatch, " + std::to_string(n) + " parameters but " + std::to_string(rows.size() - 2) +
          " covariance rows";
    return false;
  }
  out.cov.assign((size_t)n * n, 0.0);
  for (int i = 0; i < n; i++) {
    const std::vector<std::string> &row = rows[i + 2];
    if ((int)row.size() != n) {
      err = path + " line " + std::to_string(lineOf[i + 2]) + ": size mismatch, covariance row " +
            std::to_string(i + 1) + " has " + std::to_string(row.size()) + " entries, expected " + std::to_string(n);
      return false;
    }
    for (int j = 0; j < n; j++)
      if (!Number(row[j], out.cov[(size_t)i * n + j])) {
        err = path + " line " + std::to_string(lineOf[i + 2]) + ": '" + row[j] + "' is not a number";
        return false;
      }
  }
  // Symmetric to 1e-10 relative to the scale of the entry, sqrt(C_ii C_jj); then exactly symmetric.
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      const double a = out.cov[(size_t)i * n + j], b = out.cov[(size_t)j * n + i];
      const double scale = std::max(std::sqrt(std::fabs(out.cov[(size_t)i * n + i] * out.cov[(size_t)j * n + j])),
                                    std::max(std::fabs(a), std::fabs(b)));
      if (std::fabs(a - b) > 1.0e-10 * scale) {
        std::ostringstream m;
        m.precision(10);
        m << path << ": the covariance is not symmetric: C(" << out.names[i] << ", " << out.names[j] << ") = " << a
          << " but C(" << out.names[j] << ", " << out.names[i] << ") = " << b;
        err = m.str();
        return false;
      }
      out.cov[(size_t)i * n + j] = out.cov[(size_t)j * n + i] = 0.5 * (a + b);
    }
  std::vector<double> L = out.cov;
  int bad = -1;
  if (!Cholesky(L, n, bad)) {
    err = path + ": the covariance is not positive definite (Cholesky factorization fails at " + out.names[bad] + ")";
    return false;
  }
  return true;
}

bool PriorCovariance::Setup(const File &file, const std::vector<std::string> &paramNames,
                            const std::vector<bool> &fixed, int runBasis, std::string &err,
                            std::vector<std::string> &notes) {
  static const char *basisName[] = {"0 (standard)", "1 (Brune)", "2 (Park)"};
  if (file.basis >= 0 && file.basis != runBasis) {
    const std::string runName = (runBasis >= 0 && runBasis <= 2) ? basisName[runBasis] : "?";
    err = "the prior covariance file " + path_ + " holds parametrization " + basisName[file.basis] +
          " amplitudes and this run uses parametrization " + runName +
          "; it is not converted.  Rerun the fit the prior came from in this run's parametrization, or run this one "
          "in the file's";
    return false;
  }
  std::map<std::string, int> index;
  for (int i = 0; i < (int)paramNames.size(); i++) index.emplace(paramNames[i], i);
  const int m = (int)file.names.size();
  std::vector<int> keep;
  fullIndex_.clear();
  names_.clear();
  marginalized_.clear();
  mu_.clear();
  for (int i = 0; i < m; i++) {
    auto it = index.find(file.names[i]);
    if (it == index.end()) {
      err = "parameter " + file.names[i] + " of the prior covariance file " + path_ +
            " is not a parameter of this model";
      return false;
    }
    if (it->second < (int)fixed.size() && fixed[it->second]) {
      marginalized_.push_back(file.names[i]);
      continue;
    }
    keep.push_back(i);
    fullIndex_.push_back(it->second);
    names_.push_back(file.names[i]);
    mu_.push_back(file.mean[i]);
  }
  n_ = (int)keep.size();
  if (!marginalized_.empty()) {
    std::string list;
    for (const std::string &s : marginalized_) list += (list.empty() ? "" : ", ") + s;
    notes.push_back("Note: prior covariance " + path_ + ": " + list + (marginalized_.size() == 1 ? " is" : " are") +
                    " fixed in this run and marginalized (row and column dropped: the prior on the other " +
                    std::to_string(n_) + " parameter" + (n_ == 1 ? "" : "s") + " is their marginal).");
  }
  // C restricted to the kept parameters (the marginal), C = L L^T, L^-1 by forward substitution, C^-1 = L^-T L^-1.
  std::vector<double> L((size_t)n_ * n_);
  for (int a = 0; a < n_; a++)
    for (int b = 0; b < n_; b++) L[(size_t)a * n_ + b] = file.cov[(size_t)keep[a] * m + keep[b]];
  int bad = -1;
  if (!Cholesky(L, n_, bad)) {  // cannot happen for a sub-block of a positive definite matrix
    err = "the marginal prior covariance is not positive definite (at " + names_[bad] + ")";
    return false;
  }
  linv_.assign((size_t)n_ * n_, 0.0);
  double logDet = 0.0;
  for (int i = 0; i < n_; i++) logDet += 2.0 * std::log(L[(size_t)i * n_ + i]);
  for (int c = 0; c < n_; c++) {  // column c of L^-1: solve L x = e_c
    for (int i = c; i < n_; i++) {
      double s = (i == c) ? 1.0 : 0.0;
      for (int k = c; k < i; k++) s -= L[(size_t)i * n_ + k] * linv_[(size_t)k * n_ + c];
      linv_[(size_t)i * n_ + c] = s / L[(size_t)i * n_ + i];
    }
  }
  cinv_.assign((size_t)n_ * n_, 0.0);
  for (int a = 0; a < n_; a++)
    for (int b = 0; b < n_; b++) {
      double s = 0.0;
      for (int k = std::max(a, b); k < n_; k++) s += linv_[(size_t)k * n_ + a] * linv_[(size_t)k * n_ + b];
      cinv_[(size_t)a * n_ + b] = s;
    }
  logNorm_ = -0.5 * (n_ * std::log(2.0 * M_PI) + logDet);
  return true;
}

int PriorCovariance::Load(const Config &configure, const std::vector<std::string> &paramNames,
                          const std::vector<bool> &fixed, int runBasis, std::shared_ptr<const PriorCovariance> &out) {
  out.reset();
  if (configure.priorCovariance.empty()) return 0;
  auto prior = std::make_shared<PriorCovariance>();
  prior->path_ = ResolvePath(configure);
  File file;
  std::string err;
  std::vector<std::string> notes;
  if (!Read(prior->path_, file, err, configure.cacheInputFiles) ||
      !prior->Setup(file, paramNames, fixed, runBasis, err, notes)) {
    configure.outStream << "**ERROR: priorCovariance: " << err << "." << std::endl;
    return -1;
  }
  for (const std::string &note : notes) configure.Note(note);
  configure.Note("Note: prior covariance " + prior->path_ + ": a multivariate Gaussian prior on " +
                 std::to_string(prior->Size()) + " parameter" + (prior->Size() == 1 ? "" : "s") +
                 " (Total-Prior-Covariance-Chi-Squared in chiSquared.out).");
  out = prior;
  return 1;
}

void PriorCovariance::Residuals(const std::vector<double> &full, std::vector<double> &r) const {
  r.assign(n_, 0.0);
  std::vector<double> d(n_);
  for (int a = 0; a < n_; a++) d[a] = (fullIndex_[a] < (int)full.size() ? full[fullIndex_[a]] : mu_[a]) - mu_[a];
  for (int i = 0; i < n_; i++) {
    double s = 0.0;
    for (int k = 0; k <= i; k++) s += linv_[(size_t)i * n_ + k] * d[k];
    r[i] = s;
  }
}

double PriorCovariance::Chi2(const std::vector<double> &full) const {
  std::vector<double> r;
  Residuals(full, r);
  double c = 0.0;
  for (double v : r) c += v * v;
  return c;
}

void PriorCovariance::AddGradient(const std::vector<double> &full, std::vector<double> &grad) const {
  std::vector<double> d(n_);
  for (int a = 0; a < n_; a++) d[a] = (fullIndex_[a] < (int)full.size() ? full[fullIndex_[a]] : mu_[a]) - mu_[a];
  for (int a = 0; a < n_; a++) {
    if (fullIndex_[a] >= (int)grad.size()) continue;
    double s = 0.0;
    for (int b = 0; b < n_; b++) s += cinv_[(size_t)a * n_ + b] * d[b];
    grad[fullIndex_[a]] += 2.0 * s;
  }
}

void PriorCovariance::AddCurvature(std::vector<double> &g2) const {
  for (int a = 0; a < n_; a++)
    if (fullIndex_[a] < (int)g2.size()) g2[fullIndex_[a]] += 2.0 * cinv_[(size_t)a * n_ + a];
}
