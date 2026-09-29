#ifndef THM_EXPERIMENT_H
#define THM_EXPERIMENT_H

#include <string>
#include <vector>

/*!
 * THM experiments: the `experiment[<name>] key=value ...` lines of the <thm>
 * block (docs/source/theory/thm_implementation.rst, "THM experiments").
 *
 * An experiment groups THM data segments measured together.  All of its
 * segments share ONE profiled normalization, and optionally a smooth
 * background b(E) = b0 + b1 E + b2 E^2 (model units, E the c.m. energy of the
 * THM entrance pair) added to the folded HOES model.  The norm and the
 * background are linear parameters, eliminated in closed form by weighted
 * linear least squares (SolveThmProfile) instead of being fitted.
 */

/// A nucleus of the three-body reaction: the built-in table (AME2020 nuclear
/// masses) or the explicit form "Z,A,mass" (mass in u, nuclear).
struct ThmNuclide {
  std::string name;
  int Z = 0;
  int A = 0;
  double mass = 0.0;  ///< nuclear mass, u
  /// Fills `out` from a table name or "Z,A,mass"; "" on success, else the reason.
  static std::string Parse(const std::string &text, ThmNuclide &out);
  /// The table entry with this Z and A, or null.
  static const ThmNuclide *Find(int Z, int A);
};

/// One experiment[<name>] record, merged over the lines that name it.
struct ThmExperiment {
  std::string name;
  /// <segmentsData> line numbers (segment keys, counting inactive lines), ascending.
  std::vector<int> segments;
  /// Number of background terms: 0 none, 1 const, 2 linear, 3 quadratic.
  int backgroundTerms = 0;
  /// Kinematics of the three-body reaction (optional; all four or none).
  bool hasKinematics = false;
  ThmNuclide beam, target, spectator;
  double beamEnergy = 0.0;  ///< lab MeV
  /// Coulomb line-shape factor N_C of the spectator (`lineshape=on`,
  /// ThmLineshape.h); needs the kinematics.  Off by default.
  bool lineshape = false;
  /// Keys given so far (a key may not be repeated).
  std::vector<std::string> keys;
  static const char *BackgroundName(int terms);
};

/*!
 * Parses one `experiment[<name>] key=value ...` line (comment already
 * stripped, trimmed) into `experiments`, merging with an earlier line of the
 * same name.  Returns "" or what is wrong.  Keys: segments (required),
 * background, beam, target, spectator, Ebeam, lineshape (on|off); ps, theta
 * and distortion are reserved and refused ("not implemented yet").
 */
std::string ParseThmExperimentLine(const std::string &line, std::vector<ThmExperiment> &experiments);

/// Checks the complete set once the block is read: segments given, kinematics
/// all-or-none, lineshape=on only with them, no segment in two experiments.  "" or what is wrong.
std::string CheckThmExperiments(const std::vector<ThmExperiment> &experiments);

/*!
 * The closed-form profile of one experiment's linear parameters.
 *
 * Points i (every point of every segment of the experiment) carry the model
 * m_i (folded HOES), the data d_i, the error e_i and the c.m. energy E_i.  The
 * fitted curve on the data scale is
 *   f_i = s m_i + a_0 + a_1 E_i + a_2 E_i^2,
 * s = 1/n the scale on the model (n the AZURE2 norm, which multiplies the data)
 * and a_k = s b_k.  c = (s, a_0, ...) minimizes chi^2 = sum ((f_i - d_i)/e_i)^2:
 *   c = G^-1 h,  G = A~^T A~,  h = A~^T y~,  A~ = A/e,  y~ = d/e,
 * and cov(c) = G^-1.  Points with e_i = 0 carry no weight (as everywhere).
 *
 * Degenerate cases: when G is singular, or s <= 0 (no positive overlap of
 * model and data: the n = 1/s it implies is meaningless), the scale is not
 * profiled -- s = 1 (n = 1, as ESegment::ProfileNormChiSquared does) -- and
 * only the background is, if its own normal matrix is regular; otherwise
 * nothing is (s = 1, b = 0).  `status` says which.
 */
struct ThmProfile {
  int terms = 0;               ///< background terms (0..3)
  bool scaleProfiled = true;   ///< s is a profiled parameter (else fixed at 1)
  bool backgroundProfiled = true;  ///< the a_k are profiled (false only if singular)
  double s = 1.0;
  double a[3] = {0.0, 0.0, 0.0};
  /// Covariance of the profiled coefficients, q x q row-major, q = number of
  /// profiled columns, in the order (s if profiled, a_0, ...).
  std::vector<double> cov;
  double chi2 = 0.0;
  int points = 0;              ///< points with e_i != 0
  std::string status = "profiled";
  double Norm() const { return 1.0 / s; }
  /// b(E) in model units: (a_0 + a_1 E + a_2 E^2)/s.
  double Background(double energy) const;
  /// Residual (f - d)/e (0 if e = 0).
  double Residual(double m, double d, double e, double energy) const;
  /// (norm, b0, b1, b2) and their 4 x 4 covariance (row-major), from cov by
  /// n = 1/s, b_k = a_k/s; zero rows for what is not profiled.
  void Reported(double value[4], double covariance[16]) const;
};

/// The profile of m, d, e, E (all of equal length) with `terms` background terms.
ThmProfile SolveThmProfile(const std::vector<double> &m, const std::vector<double> &d,
                           const std::vector<double> &e, const std::vector<double> &energy, int terms);

/*!
 * Total derivatives through the profile.  `Jm` is d m_i/d p_c (n x nCols,
 * row-major).  Returns in `J` the derivative of the standardized residual
 * r_i = (f_i - d_i)/e_i with the linear parameters at their optimum c*(p),
 * and in `G` the derivative of the curve as the output files show it next to
 * data scaled by n*(p0): q_i = f_i(p) / s(p0), i.e. G = (s Jm + A dc/dp)/s.
 *
 * With A~(p) = A(p)/e, rho = y~ - A~c (minus the residual) and
 * dA~/dp = [Jm/e, 0, ...] (only the model column depends on p), the
 * derivative of the least-squares solution is exact (Golub & Pereyra 1973):
 *   dc/dp = G^-1 ( (dA~/dp)^T rho - A~^T (dA~/dp) c ),
 *   dr/dp = (dA~/dp) c + A~ dc/dp.
 * The second term of dc/dp alone gives dr/dp = P_perp (dA~/dp) c, the
 * projection onto the orthogonal complement of the columns of A~ (Kaufman's
 * approximation, exact at zero residual); the first carries the curvature
 * that the residual adds.  With the model column alone this is the
 * ds/dp = (sum d J/e^2 - 2 s sum m J/e^2)/S_mm of ComputeTHMRows.  A column
 * that is not profiled has no row in dc/dp.
 */
void ThmProfileDerivative(const ThmProfile &profile, const std::vector<double> &m,
                          const std::vector<double> &d, const std::vector<double> &e,
                          const std::vector<double> &energy, const std::vector<double> &Jm, int nCols,
                          std::vector<double> &J, std::vector<double> &G);

/// What the output file and pyazr report for one experiment.
struct ThmExperimentReport {
  std::string name;
  std::vector<int> segments;   ///< segment keys in use
  std::string background;      ///< none | const | linear | quadratic
  int points = 0;
  double chi2 = 0.0;
  double value[4] = {1.0, 0.0, 0.0, 0.0};  ///< norm, b0, b1, b2
  double covariance[16] = {0.0};
  std::string status;
};

#endif
