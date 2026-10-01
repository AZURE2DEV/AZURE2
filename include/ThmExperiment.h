#ifndef THM_EXPERIMENT_H
#define THM_EXPERIMENT_H

#include <memory>
#include <string>
#include <vector>

struct ThmWeightTable;

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
  /*!
   * Spectator-momentum window (`ps=`, ThmLineshape.h ThmSpectatorWindow):
   * the HOES model at E is the average of the cross section over the
   * spectator momentum p_s in [psMin, psMax] (MeV/c) with the weight of the
   * accepted events, |phi(p_s)|^2 p_s^2 dp_s.  DELTA (default) is the
   * quasi-free point p_s = 0 (or the scalar spectatorEnergy).  Needs the
   * kinematics.
   */
  enum PsKind { PS_DELTA, PS_HULTHEN, PS_GAUSS, PS_TABLE };
  PsKind psKind = PS_DELTA;
  double psA = 0.2317, psB = 1.202;  ///< Hulthen a, b (fm^-1; Tribble 2014 eq. 4.4)
  double psFwhm = 0.0;               ///< gauss: FWHM of |phi(p)|^2 (MeV/c)
  double psMin = 0.0, psMax = 0.0;   ///< window (MeV/c); a table: its first and last p
  std::string psTable;               ///< table: the file as written
  std::vector<double> psTableP, psTableW;  ///< table rows (Config::ReadThmBlock loads them)
  int psNodes = 16;                  ///< Gauss-Legendre nodes in p_s (`psNodes=`)
  /*!
   * Distortion factor R(E) multiplying the model of every segment
   * (`distortion=`, ThmDistortion.h): none (default), coulomb (point-Coulomb
   * waves in a + A and s + F), optical (per channel `opticalAA=`,
   * `opticalSF=`: plane, coulomb or a Woods-Saxon potential) or a table
   * w(E) (`table:<file>`, the weight[k] format).  coulomb and optical need the
   * kinematics.
   */
  enum DistortionKind { DIST_NONE, DIST_COULOMB, DIST_OPTICAL, DIST_TABLE };
  DistortionKind distortion = DIST_NONE;
  std::string distortionTable;  ///< table: the file as written
  std::shared_ptr<const ThmWeightTable> distortionWeights;  ///< table rows (Config::ReadThmBlock)
  /// opticalAA= / opticalSF=: 0 plane, 1 coulomb, 2 Woods-Saxon with p[10] =
  /// V,R,a,W,RW,aW,WD,RD,aD,RC (MeV, fm), 3 a global optical potential
  /// (ThmOptical.h: `<name>` or `<name>:extrapolate`, index `global`).
  struct Optical {
    int kind = 1;
    double p[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int global = -1;           ///< kind 3: index in ThmGlobalOpticals()
    bool extrapolate = false;  ///< kind 3: allowed outside its validity range (warned)
  };
  Optical opticalAA, opticalSF;
  /// spectatorAngle=: 0 qf (k_sF along k_aA, default), 1 lab degrees, 2 cm:degrees.
  int angleKind = 0;
  double angle = 0.0;
  bool hasDistortionRef = false;
  double distortionRef = 0.0;  ///< MeV (default: the middle of the data range)
  bool distortionRatioPW = true;  ///< distortionRatio=dwpw (default) or dw
  bool boundYukawa = false;       ///< boundState=yukawa (default whittaker)
  double boundRmin = 0.0;         ///< boundState=...:rmin (fm)
  /*!
   * Entrance vertex (`vertexModel=`, ThmDwVertex.h): pw (default) the
   * plane-wave M_l of the first section, dw the surface term of the prior-form
   * DWBA with the experiment's distorted waves (distortion=coulomb|optical),
   * which then replaces R(E) (not applied).
   */
  bool vertexDW = false;
  /*!
   * Angular window of the exit pair (`theta=`, ThmAngular.h): the model of
   * every segment is the HOES dsigma/dOmega averaged over theta_cm in
   * [thetaMin, thetaMax] (degrees, c.m. of the exit pair relative to p_xA)
   * instead of the angle-integrated cross section.  `theta=all` (default)
   * keeps the latter.
   */
  bool hasTheta = false;
  double thetaMin = 0.0, thetaMax = 180.0;
  /// Keys given so far (a key may not be repeated).
  std::vector<std::string> keys;
  static const char *BackgroundName(int terms);
};

/*!
 * Parses one `experiment[<name>] key=value ...` line (comment already
 * stripped, trimmed) into `experiments`, merging with an earlier line of the
 * same name.  Returns "" or what is wrong.  Keys: segments (required),
 * background, beam, target, spectator, Ebeam, lineshape (on|off), ps
 * (delta | hulthen:pmin-pmax | hulthen:a,b:pmin-pmax | gauss:FWHM:pmin-pmax |
 * table:file), psNodes, distortion (none | coulomb | optical | table:file),
 * opticalAA, opticalSF, spectatorAngle, distortionRef, distortionRatio,
 * boundState, theta (all | thmin-thmax, degrees, 0 <= thmin <= thmax <= 180),
 * vertexModel (pw | dw).
 */
std::string ParseThmExperimentLine(const std::string &line, std::vector<ThmExperiment> &experiments);

/// Checks the complete set once the block is read: segments given, kinematics
/// all-or-none, lineshape=on and a ps window only with them, psNodes only with a
/// window, no segment in two experiments, distortion keys consistent (and
/// distortionRatio=dw not with a ps window, which would weight the model by the
/// momentum distribution twice).  "" or what is wrong.
std::string CheckThmExperiments(const std::vector<ThmExperiment> &experiments);

/*!
 * The rules that keep the Coulomb treatment of a THM experiment consistent
 * with the global options of the <thm> block
 * (docs/source/theory/thm_implementation.rst, "Coulomb effects: what each
 * option contains").  Refused ("" or the reason, prefixed by the experiment):
 *  - coulombIntegral=1 with vertexModel=dw: the DW vertex is the surface term;
 *    its external prior term, whose plane-wave limit is C_l, is not computed.
 *  - coulombIntegral=1 with the computed distortion factor R(E) and a
 *    distorted a + A wave (distortion=coulomb, or optical with opticalAA not
 *    plane): the a + A Coulomb wave of R already contains the x-A Coulomb
 *    interaction that C_l adds (Z_a = Z_x + Z_s; the prior operator outside
 *    the radius, V_xA + V_sA - U_aA, vanishes in the zero range of R), so C_l
 *    would count it twice.
 * Warned (appended to `warnings` if given, without the "WARNING: " prefix):
 *  - coulombIntegral=1 with distortion=table: if the table is a DWBA with a
 *    distorted a + A wave, the same double counting.
 */
std::string CheckThmCoulombConsistency(const std::vector<ThmExperiment> &experiments, bool coulombIntegral,
                                       std::vector<std::string> *warnings = nullptr);

/// Reads a `ps=table:` file: two columns, p_s (MeV/c, >= 0, strictly
/// increasing) and the event weight w >= 0 per unit p_s; '#' starts a comment;
/// at least two rows and some positive weight.  "" or what is wrong.
std::string ReadThmPsTable(const std::string &path, std::vector<double> &p, std::vector<double> &w);

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
