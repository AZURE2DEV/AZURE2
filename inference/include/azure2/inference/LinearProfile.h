#ifndef AZURE2_INFERENCE_LINEARPROFILE_H
#define AZURE2_INFERENCE_LINEARPROFILE_H

#include <string>
#include <vector>

namespace azure2 {
namespace inference {

/*!
 * The parameters of the objective that enter the model linearly and are profiled in closed form rather than fitted
 * (T1a and T1b of the objective; docs/upgrade/objective_arrays.md section 2): a profiled norm of one segment, and
 * the norm and polynomial background of a THM experiment.  Arithmetic on arrays: which points enter (the segments of
 * an experiment, the points with a model record) and the model in the frame of the data are the engine's gathers
 * (Segment::ProfileNormChiSquared, WorkingCopy::ProfileThmGroup, src/ThmNorm.cpp).
 *
 * Every array holds one entry per point, all of equal length: the model m_i and the datum d_i with its error e_i in
 * the frame and unit of the data (b, b/sr, ...), and, for the background, the c.m. energy E_i (MeV).  Points with
 * e_i = 0 carry no weight (as everywhere in the objective).
 */

/*!
 * The profiled norm of one segment (T1a).  With the norm n multiplying the data (residual m - n d, error n e), the
 * chi^2 = sum ((m - n d)/(n e))^2 is minimized at n* = S_mm / S_md, giving chi^2_min = S_dd - S_md^2 / S_mm, where
 * S_ab = sum a_i b_i / e_i^2; a negative chi^2_min (rounding) is returned as 0.  Degenerate (S_md <= 0 or S_mm <= 0:
 * no positive model-data overlap): n = 1 and the unscaled chi^2 = sum ((m - d)/e)^2, formed point by point as the
 * objective's residual with the norm 1 forms it (ModelPass.h DataResidual(f, J, y, dy, 1).Chi2()), so the fit still
 * sees a finite, honest value.  Returns chi^2; `norm` receives n* (dimensionless; 1 if degenerate).
 */
double ProfiledNormChi2(const std::vector<double> &m, const std::vector<double> &d, const std::vector<double> &e,
                        double &norm);

/*!
 * The closed-form profile of one experiment's linear parameters (T1b).
 *
 * Points i (every point of every segment of the experiment) carry the model
 * m_i (folded HOES), the data d_i, the error e_i and the c.m. energy E_i.  The
 * fitted curve on the data scale is
 *   f_i = s m_i + a_0 + a_1 E_i + a_2 E_i^2,
 * s = 1/n the scale on the model (n the AZURE2 norm, which multiplies the data)
 * and a_k = s b_k.  c = (s, a_0, ...) minimizes chi^2 = sum ((f_i - d_i)/e_i)^2:
 *   c = G^-1 h,  G = A~^T A~,  h = A~^T y~,  A~ = A/e,  y~ = d/e,
 * and cov(c) = G^-1.  Points with e_i = 0 carry no weight (as everywhere).
 * Units: s dimensionless, a_k in the unit of the data times MeV^-k.
 *
 * Degenerate cases: when G is singular, or s <= 0 (no positive overlap of
 * model and data: the n = 1/s it implies is meaningless), the scale is not
 * profiled -- s = 1 (n = 1, as ProfiledNormChi2 does) -- and
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

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_LINEARPROFILE_H
