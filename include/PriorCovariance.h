#ifndef PRIORCOVARIANCE_H
#define PRIORCOVARIANCE_H

#include <memory>
#include <string>
#include <vector>

class Config;

/*!
 * A multivariate Gaussian prior on a set of fit parameters (SAMMY's "PUP", a Bayesian sequential update): a fit
 * starts from an earlier evaluation's parameters mu and their covariance C, and
 *
 *     chi2 += (p - mu)^T C^-1 (p - mu)
 *
 * over the listed parameters, on their values as the fit's parameter vector holds them (the param.sav space:
 * level energies in MeV, reduced-width amplitudes of the run's parametrization, norms, shifts, ...).
 *
 * The file (the <config> key priorCovariance=<path>, relative to the .azr's directory) is plain text with '#'
 * comments: an optional "#parametrization <0|1|2>" line (param.sav's tag), then the parameter NAMES on one line, their
 * MEANS on the next, and the n rows of the symmetric n x n COVARIANCE in the same order.
 *
 * Setup() matches the names onto the run's parameter vector.  A listed parameter that is fixed in the run is
 * marginalized (its row and column dropped: the marginal of a Gaussian is the Gaussian of the remaining block), the
 * rest is factored once, C = L L^T, and every evaluation uses the whitened residuals r = L^-1 (p - mu), r.r being the
 * chi2 term, 2 C^-1 (p - mu) = 2 L^-T r its gradient and L^-1 the Jacobian of the n least-squares rows.
 */
class PriorCovariance {
 public:
  /// The file as read: names, means, the dense n x n covariance (row-major) and the basis tag (-1 = none).
  struct File {
    std::vector<std::string> names;
    std::vector<double> mean;
    std::vector<double> cov;
    int basis = -1;
  };

  /// Parse the file.  False with err set (bad numbers, a row of the wrong length, too few or too many rows, a size
  /// mismatch between names and means, a duplicate name, a non-symmetric matrix at 1e-10 relative to
  /// sqrt(C_ii C_jj), or a matrix that is not positive definite).
  /// cached: read through InputFiles' cache (Config::cacheInputFiles, a pyazr session's rebuilds).
  static bool Read(const std::string &path, File &out, std::string &err, bool cached = false);
  /// The key's path resolved against the .azr's directory (an absolute path is kept); empty without the key.
  static std::string ResolvePath(const Config &configure);
  /*!
   * Match a read file onto a run.  paramNames / fixed: the run's whole parameter vector (Minuit order) and its fixed
   * flags; runBasis: AZUREParams::BasisForConfig.  Refused (false, err set) when the file's basis tag differs from
   * the run's or a name is not a parameter of the model.  Fixed listed parameters are marginalized; notes receives
   * the Note lines to print.
   */
  bool Setup(const File &file, const std::vector<std::string> &paramNames, const std::vector<bool> &fixed,
             int runBasis, std::string &err, std::vector<std::string> &notes);
  /*!
   * Read and set up from the run's <config> key, printing the Notes and, on a refusal, the ERROR line to
   * configure.outStream.  Returns 0 without the key (out reset), 1 with a prior set up (out), -1 on a refusal.
   */
  static int Load(const Config &configure, const std::vector<std::string> &paramNames,
                  const std::vector<bool> &fixed, int runBasis, std::shared_ptr<const PriorCovariance> &out);

  /// The resolved path of the file.
  const std::string &Path() const { return path_; }
  /// Number of parameters the prior acts on (after marginalization; may be 0).
  int Size() const { return n_; }
  /// Full-vector (Minuit) index of each prior parameter, in the file's order.
  const std::vector<int> &FullIndex() const { return fullIndex_; }
  const std::vector<std::string> &Names() const { return names_; }
  const std::vector<double> &Mean() const { return mu_; }
  /// The listed parameters that are fixed in this run (marginalized out).
  const std::vector<std::string> &Marginalized() const { return marginalized_; }

  /// Whitened residuals r = L^-1 (p - mu), p the full parameter vector.
  void Residuals(const std::vector<double> &full, std::vector<double> &r) const;
  /// The chi2 term (p - mu)^T C^-1 (p - mu) = r.r.
  double Chi2(const std::vector<double> &full) const;
  /// grad[full index] += 2 C^-1 (p - mu) (indices beyond grad's size are skipped).
  void AddGradient(const std::vector<double> &full, std::vector<double> &grad) const;
  /// g2[full index] += 2 (C^-1)_ii, the diagonal curvature (Minuit's G2).
  void AddCurvature(std::vector<double> &g2) const;
  /// (L^-1)_ij, the Jacobian of residual row i with respect to prior parameter j (zero for j > i).
  double Linv(int i, int j) const { return linv_[(size_t)i * n_ + j]; }
  /// (C^-1)_ij.
  double Cinv(int i, int j) const { return cinv_[(size_t)i * n_ + j]; }
  /*!
   * Append the n whitened rows r = L^-1 (p - mu) to a residual vector and, when jac is given, their Jacobian rows to
   * a row-major Jacobian with nCols columns: column colOfFull(full index) of row i holds (L^-1)_ij (the rows are
   * linear in p).  colOfFull returns -1 for a parameter without a column.
   */
  template <class ColOf>
  void AppendRows(const std::vector<double> &full, ColOf colOfFull, int nCols, std::vector<double> &res,
                  std::vector<double> *jac) const {
    std::vector<double> r;
    Residuals(full, r);
    for (int i = 0; i < n_; i++) {
      res.push_back(r[i]);
      if (!jac) continue;
      const size_t base = jac->size();
      jac->resize(base + (size_t)nCols, 0.0);
      for (int j = 0; j <= i; j++) {
        const int a = colOfFull(fullIndex_[j]);
        if (a >= 0 && a < nCols) (*jac)[base + (size_t)a] = Linv(i, j);
      }
    }
  }
  /// log of the normalization of the Gaussian, -0.5 log det(2 pi C) (MCMC log prior).
  double LogNormalization() const { return logNorm_; }

 private:
  std::string path_;
  int n_ = 0;
  std::vector<int> fullIndex_;
  std::vector<std::string> names_, marginalized_;
  std::vector<double> mu_, linv_, cinv_;
  double logNorm_ = 0.0;
};

#endif
