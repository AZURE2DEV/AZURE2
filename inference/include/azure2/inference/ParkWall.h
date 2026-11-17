#ifndef AZURE2_INFERENCE_PARKWALL_H
#define AZURE2_INFERENCE_PARKWALL_H

#include <cmath>
#include <vector>

namespace azure2 {
namespace inference {

/*!
 * Park's wall J > 0 (T6 of the objective; docs/upgrade/objective_arrays.md section 2) as arithmetic on arrays.
 *
 * In the Park formalism each R-matrix level carries the overlap J = 1 - sum_c gamma_c^2 S'_c (dimensionless; gamma_c
 * the reduced-width amplitude of particle channel c in MeV^1/2, S'_c = dS_c/dE at the level's energy in MeV^-1).  A
 * level with J <= 0 has no standard R-matrix counterpart, and the objective holds the fit off that region with the
 * wall sum over levels with J < 0 of (J/scale)^2.  A level is on the wall when J < 0: a J >= 0 and a NaN J are not
 * (OnParkWall, the one rule of the value, the gradient and the rows).
 *
 * The derivatives, with S''_c = d^2S_c/dE^2 (MeV^-2) and the level energy E (MeV):
 *   dJ/dgamma_c = -2 gamma_c S'_c,    dJ/dE = -sum_c gamma_c^2 S''_c.
 * The gradient and the least-squares rows are written in two forms, each with its own rounding, kept as the paths
 * had them before the arrays: the gradient as factor * dJ, factor = 2 J/scale^2 (ParkWallGradient); a row as
 * r = J/scale with dJ/scale (ParkWallRow).  Their Gauss-Newton product 2 r dr agrees with the gradient to rounding,
 * not bit for bit.
 *
 * Which level is in the R matrix and which channel is a particle channel is the engine's gather (ParkLevelsOf,
 * include/Objective.h): this header holds no rule of eligibility.
 */
struct ParkLevel {
  int group = 0;   ///< the engine's J group (1-based), for the caller's map of the derivatives
  int level = 0;   ///< the level in its J group (1-based)
  double J = 0.0;  ///< 1 - sum_c gamma_c^2 S'_c
  std::vector<int> channels;  ///< the particle channels (1-based, the J group's numbering); empty in a value gather
  std::vector<double> gamma;  ///< gamma_c, MeV^1/2, per entry of channels
  std::vector<double> dS;     ///< S'_c, MeV^-1
  std::vector<double> d2S;    ///< S''_c, MeV^-2
};

/// The level is on the wall: J < 0 (false for J >= 0 and for a NaN J).
inline bool OnParkWall(double J) { return J < 0.0; }

/// The wall from 0, in list order: the sum over the levels on the wall of (J/scale)^2.
inline double ParkWallChi2(const std::vector<ParkLevel> &levels, double scale) {
  double chi2 = 0.0;
  for (const ParkLevel &l : levels)
    if (OnParkWall(l.J)) chi2 += std::pow(l.J / scale, 2.0);
  return chi2;
}

/// The gradient's form: d(J/scale)^2/dgamma_c = factor (-2 gamma_c S'_c) into dGamma (resized to the channels) and
/// d/dE = factor (-sum_c gamma_c^2 S''_c) into dE, factor = 2 J/scale^2.  False, with zeros, off the wall.
inline bool ParkWallGradient(const ParkLevel &l, double scale, std::vector<double> &dGamma, double &dE) {
  dGamma.assign(l.channels.size(), 0.0);
  dE = 0.0;
  if (!OnParkWall(l.J)) return false;
  const double s2 = scale * scale;
  const double factor = 2.0 * l.J / s2;
  double sum = 0.0;
  for (size_t c = 0; c < l.channels.size(); c++) {
    const double g = l.gamma[c];
    dGamma[c] = factor * (-2.0 * g * l.dS[c]);
    sum -= g * g * l.d2S[c];
  }
  dE = factor * sum;
  return true;
}

/// The row's form: r = J/scale, dr/dgamma_c = -2 gamma_c S'_c/scale into dGamma (resized), dr/dE =
/// (-sum_c gamma_c^2 S''_c)/scale into dE.  Off the wall r = 0, zeros, and false (the row is kept, with no Jacobian).
inline bool ParkWallRow(const ParkLevel &l, double scale, double &r, std::vector<double> &dGamma, double &dE) {
  dGamma.assign(l.channels.size(), 0.0);
  dE = 0.0;
  r = 0.0;
  if (!OnParkWall(l.J)) return false;
  double sum = 0.0;
  for (size_t c = 0; c < l.channels.size(); c++) {
    const double g = l.gamma[c];
    dGamma[c] = -2.0 * g * l.dS[c] / scale;
    sum -= g * g * l.d2S[c];
  }
  dE = sum / scale;
  r = l.J / scale;
  return true;
}

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_PARKWALL_H
