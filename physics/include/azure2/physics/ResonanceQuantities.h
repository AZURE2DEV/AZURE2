#ifndef AZURE2_PHYSICS_RESONANCEQUANTITIES_H
#define AZURE2_PHYSICS_RESONANCEQUANTITIES_H

#include <vector>

namespace azure2 {
namespace physics {

/*!
 * The quantities of an isolated resonance that a model datum compares with a literature value
 * (docs/upgrade/data_kinds.md section 4.7): the resonance strength, the branching ratio and the total width, as
 * functions of the observed partial widths of one level at its observed energy E_r, with their partial derivatives.
 * Widths in any one energy unit (the results carry it); spins as 2J (integers).
 *
 * The widths are the observed ones, Gamma_c = 2 P_c(E_r) gamma~_c^2 / (1 + sum_c' gamma~_c'^2 dS_c'/dE) in Brune's
 * parameters (level_parameters.md section 2); the caller sums them over the channels of a pair (all l and s of that
 * pair) and over all channels for Gamma.
 */

/// omega = (2J + 1)(1 + delta_aA) / ((2 j_a + 1)(2 j_A + 1)) (data_kinds.md eq. 8): two_J = 2J of the level, two_ja
/// and two_jA = 2j of the two entrance particles; identical: the entrance particles are identical (delta_aA = 1).
double statistical_factor(int two_J, int two_ja, int two_jA, bool identical);

/// A value with its partial derivatives with respect to the arguments of the function that made it.
struct ResonanceStrength {
  double value = 0.0;    ///< omega Gamma_a Gamma_b / Gamma
  double d_entrance = 0.0;  ///< d/dGamma_a = omega Gamma_b / Gamma
  double d_exit = 0.0;      ///< d/dGamma_b = omega Gamma_a / Gamma
  double d_total = 0.0;     ///< d/dGamma  = -omega Gamma_a Gamma_b / Gamma^2
};

/*!
 * omega gamma = omega Gamma_a Gamma_b / Gamma (data_kinds.md eq. 8), with Gamma_a and Gamma_b the observed widths to
 * the entrance and exit pairs and Gamma the observed total width, treated as three independent arguments: a caller
 * whose Gamma is the sum of the partial widths adds d_total to d_entrance and d_exit (the chain rule).  Gamma = 0
 * gives 0 and derivatives 0 (a level without width has no strength).
 */
ResonanceStrength resonance_strength(double omega, double Gamma_a, double Gamma_b, double Gamma);

struct BranchingRatio {
  double value = 0.0;    ///< Gamma_b / Gamma
  double d_exit = 0.0;   ///< d/dGamma_b = 1 / Gamma
  double d_total = 0.0;  ///< d/dGamma  = -Gamma_b / Gamma^2
};

/// Gamma_b / Gamma, with Gamma_b and Gamma independent (as resonance_strength).  Gamma = 0 gives 0 and derivatives 0.
BranchingRatio branching_ratio(double Gamma_b, double Gamma);

/// Gamma = sum_c |Gamma_c| of observed partial widths, signed as <levels> gives them (the sign is that of the
/// amplitude); d Gamma / d Gamma_c = sign(Gamma_c) (0 for Gamma_c = 0), into *d_widths when it is not null.
double total_width(const std::vector<double> &widths, std::vector<double> *d_widths = nullptr);

}  // namespace physics
}  // namespace azure2

#endif
