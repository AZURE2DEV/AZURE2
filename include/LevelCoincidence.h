#ifndef LEVEL_COINCIDENCE_H
#define LEVEL_COINCIDENCE_H

#include <cmath>

/*!
 * Two levels of a J-group closer than this (MeV) are coincident: the divided
 * differences of their shift functions in the Brune and Park level matrices,
 * (S_a - S_b)/(E_a - E_b), are replaced by their analytic limits in S'(E) and
 * S''(E).  At this gap a linear expansion of S is exact to about 1e-9 of S';
 * the difference quotient is still free of round-off noise well above it.
 */
constexpr double kLevelCoincidenceGap = 1.0e-9;

inline bool LevelsCoincide(double energyA, double energyB) {
  return std::fabs(energyA - energyB) < kLevelCoincidenceGap;
}

#endif
