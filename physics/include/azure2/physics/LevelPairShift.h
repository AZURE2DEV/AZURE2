#ifndef AZURE2_PHYSICS_LEVEL_PAIR_SHIFT_H
#define AZURE2_PHYSICS_LEVEL_PAIR_SHIFT_H

// The shift function of one channel between two levels of one J^pi, as the Brune (2002) and Park (2021) level
// matrices and the Brune -> formal transformation hold it (docs/source/developer/level_matrix.rst).  Header only:
// the level-matrix fills call it per level pair and channel in the evaluation loop, in double and (the transformed
// route of ChannelSpace) in long double.

namespace azure2 {
namespace physics {

/// The shift function S_c, its first and second energy derivatives at a level energy E (MeV; S' in MeV^-1, S'' in
/// MeV^-2).  S'' enters only the analytic gradients.
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
 */
template <class Real>
struct LevelPairShift {
  ShiftAtLevel<Real> a, b;
  Real gap;  ///< E_a - E_b (MeV)

  Real slope() const { return a.S - b.S; }
  Real intercept() const { return b.S * a.E - a.S * b.E; }
  Real chord(Real E) const { return a.S * (E - b.E) - b.S * (E - a.E); }
  Real slope_by_Ea(Real s1) const { return a.dS - s1; }
  Real slope_by_Eb(Real s1) const { return -b.dS + s1; }
  Real intercept_by_Ea(Real s0) const { return b.S - a.dS * b.E - s0; }
  Real chord_by_Ea(Real E, Real D) const { return a.dS * (E - b.E) + b.S - D; }
  Real chord_by_Eb(Real E, Real D) const { return -a.S - b.dS * (E - a.E) + D; }
};

/// The denominator of the fractions of LevelPairShift: E_a - E_b.
template <class Real>
Real level_gap(Real Ea_MeV, Real Eb_MeV) {
  return Ea_MeV - Eb_MeV;
}

/// 1 / level_gap, for the callers that multiply.
inline double inverse_level_gap(double Ea_MeV, double Eb_MeV) { return 1.0 / level_gap(Ea_MeV, Eb_MeV); }

/// The chord of S_c between levels a and b (a != b; LevelPairShift).
template <class Real>
LevelPairShift<Real> level_pair_shift(const ShiftAtLevel<Real> &a, const ShiftAtLevel<Real> &b) {
  return LevelPairShift<Real>{a, b, level_gap(a.E, b.E)};
}

}  // namespace physics
}  // namespace azure2

#endif
