#ifndef AZURE2_THMPROFILE_H
#define AZURE2_THMPROFILE_H

#include <string>
#include <vector>

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

#endif
