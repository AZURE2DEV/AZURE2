#ifndef AZURE2_PHYSICS_LEVEL_PAIR_SHIFT_H
#define AZURE2_PHYSICS_LEVEL_PAIR_SHIFT_H

// The shift function of one channel between two levels of one J^pi, as the Brune (2002) and Park (2021) level
// matrices and the Brune -> formal transformation hold it (docs/source/developer/level_matrix.rst).  Header only:
// the level-matrix fills call it per level pair and channel in the evaluation loop, in double and (the transformed
// route of ChannelSpace) in long double.

#include <cmath>

namespace azure2 {
namespace physics {

/*!
 * Two levels of one J^pi closer than this (MeV) are coincident: LevelPairShift replaces the divided differences by
 * their Taylor forms about the pair.  Above it the divided difference loses about 2 eps |S| / gap to round-off
 * (4e-10 |S| here); below it the Taylor forms are off by gap^2 S'''/12 in the slope and gap S'''/6 in its
 * derivatives, both below the 1e-8 to which S' itself is known (Richardson differences, ChannelFunctions.h).
 */
constexpr double kCoincidentLevelGap_MeV = 1e-6;

/// Are two level energies (MeV) closer than kCoincidentLevelGap_MeV?
template <class Real>
bool levels_coincide(Real Ea_MeV, Real Eb_MeV) {
  return std::fabs(Ea_MeV - Eb_MeV) < Real(kCoincidentLevelGap_MeV);
}

/// The shift function S_c and its first and second energy derivatives at a level energy E (MeV; S' in MeV^-1, S'' in
/// MeV^-2).  S' enters coincident levels (and Park's matrix), S'' only the analytic gradients of coincident levels.
template <class Real>
struct ShiftAtLevel {
  Real E = 0, S = 0, dS = 0, d2S = 0;
};

/*!
 * The chord of S_c between levels a and b: D(E) = s0 + s1 E through (E_a, S_a) and (E_b, S_b),
 *   s1   = (S_a - S_b) / (E_a - E_b)                          the slope (Park's overlap, Brune's N; MeV^-1),
 *   s0   = (S_b E_a - S_a E_b) / (E_a - E_b)                  D at E = 0 (Brune's K, the transformation's N),
 *   D(E) = (S_a (E - E_b) - S_b (E - E_a)) / (E_a - E_b)      D at the channel energy E (Brune's level matrix),
 * and their derivatives by the level energies, which the analytic gradients need:
 *   ds1/dE_a = (S'_a - s1) / (E_a - E_b),       ds1/dE_b = (s1 - S'_b) / (E_a - E_b),
 *   ds0/dE_a = (S_b - S'_a E_b - s0) / (E_a - E_b),
 *   dD/dE_a  = (S'_a (E - E_b) + S_b - D) / (E_a - E_b),   dD/dE_b = (D - S_a - S'_b (E - E_a)) / (E_a - E_b).
 *
 * Every quantity is a fraction over the one denominator `gap`: the members give the numerators, and the caller
 * divides by gap (or multiplies by 1/gap, inverse_level_gap) where its sum needs it.  A derivative's numerator takes
 * the quotient it contains (s1, s0 or D) as the caller computed it.
 *
 * Coincident levels (levels_coincide): gap is 1 and the numerators are the Taylor forms about the pair, with the
 * means S'_m = (S'_a + S'_b) / 2 and S''_m = (S''_a + S''_b) / 2 (the mean of the slopes at the two ends is the
 * slope of the chord to second order in the gap):
 *   s1 = S'_m,   s0 = S_a - S'_m E_a,   D(E) = S_a + S'_m (E - E_a),
 *   ds1/dE_a = ds1/dE_b = S''_m / 2,   ds0/dE_a = (S'_a - S'_b) / 2 - S''_m E_a / 2,
 *   dD/dE_a = (S'_a - S'_b) / 2 + S''_m (E - E_a) / 2,   dD/dE_b = S''_m (E - E_a) / 2.
 * At E_a = E_b these are the limits of the divided differences.
 */
template <class Real>
struct LevelPairShift {
  ShiftAtLevel<Real> a, b;
  Real gap;         ///< E_a - E_b (MeV); 1 for coincident levels
  bool coincident;  ///< levels_coincide(E_a, E_b)

  Real slope() const { return coincident ? MeanSlope() : a.S - b.S; }
  Real intercept() const { return coincident ? a.S - MeanSlope() * a.E : b.S * a.E - a.S * b.E; }
  Real chord(Real E) const { return coincident ? a.S + MeanSlope() * (E - a.E) : a.S * (E - b.E) - b.S * (E - a.E); }
  Real slope_by_Ea(Real s1) const { return coincident ? HalfMeanCurvature() : a.dS - s1; }
  Real slope_by_Eb(Real s1) const { return coincident ? HalfMeanCurvature() : -b.dS + s1; }
  Real intercept_by_Ea(Real s0) const {
    return coincident ? HalfSlopeDifference() - HalfMeanCurvature() * a.E : b.S - a.dS * b.E - s0;
  }
  Real chord_by_Ea(Real E, Real D) const {
    return coincident ? HalfSlopeDifference() + HalfMeanCurvature() * (E - a.E) : a.dS * (E - b.E) + b.S - D;
  }
  Real chord_by_Eb(Real E, Real D) const {
    return coincident ? HalfMeanCurvature() * (E - a.E) : -a.S - b.dS * (E - a.E) + D;
  }

 private:
  Real MeanSlope() const { return Real(0.5) * (a.dS + b.dS); }
  Real HalfMeanCurvature() const { return Real(0.25) * (a.d2S + b.d2S); }
  Real HalfSlopeDifference() const { return Real(0.5) * (a.dS - b.dS); }
};

/// The denominator of the fractions of LevelPairShift: E_a - E_b, or 1 for coincident levels.
template <class Real>
Real level_gap(Real Ea_MeV, Real Eb_MeV) {
  return levels_coincide(Ea_MeV, Eb_MeV) ? Real(1) : Ea_MeV - Eb_MeV;
}

/// 1 / level_gap, for the callers that multiply.
inline double inverse_level_gap(double Ea_MeV, double Eb_MeV) { return 1.0 / level_gap(Ea_MeV, Eb_MeV); }

/// The chord of S_c between levels a and b (a != b; LevelPairShift).
template <class Real>
LevelPairShift<Real> level_pair_shift(const ShiftAtLevel<Real> &a, const ShiftAtLevel<Real> &b) {
  return LevelPairShift<Real>{a, b, level_gap(a.E, b.E), levels_coincide(a.E, b.E)};
}

}  // namespace physics
}  // namespace azure2

#endif
