#ifndef AZURE2_EXPERIMENT_LEGENDREEXPANSION_H
#define AZURE2_EXPERIMENT_LEGENDREEXPANSION_H

#include <string>

namespace azure2 {
namespace experiment {

/*!
 * The normalisations of the c.m. Legendre coefficients of an angular distribution (quantity legendre, field
 * 'expansion'; docs/upgrade/data_kinds.md section 4.4, EXFOR's LEG, LEG/RS, LEG/RSL, LEG/2L2, LEG/4PI).  With
 *
 *   dsigma/dOmega_cm(theta) = sum_L A_L P_L(cos theta)   (A_L in b/sr, A_0 = sigma / 4 pi),
 *
 * the coefficient of order L in each expansion is
 *
 *   a_over_a0  A_L / A_0                    (1; EXFOR LEG/RS)
 *   absolute   A_L                          (b/sr; EXFOR LEG)
 *   rsl        A_L / ((2L + 1) A_0)         (1; EXFOR LEG/RSL)
 *   2l2        2 A_L / (2L + 1)             (b/sr; EXFOR LEG/2L2: dsigma/dOmega = sum_L (2L+1)/2 a_L P_L)
 *   4pi        4 pi A_L                     (b; EXFOR LEG/4PI: dsigma/dOmega = (1/4 pi) sum_L a_L P_L)
 *
 * The engine computes A_L / A_0 or A_L (the expansion's base, Normalised); the rest is the factor of
 * LegendreExpansionFactor, which depends on L only.
 */
enum class LegendreExpansion { AOverA0, Absolute, Rsl, TwoL2, FourPi };

/// The expansion of a word of the field 'expansion' (a_over_a0, absolute, rsl, 2l2, 4pi); AOverA0 for any other word.
LegendreExpansion LegendreExpansionOf(const std::string &word);
/// The word of the field 'expansion' for an expansion.
const char *LegendreExpansionName(LegendreExpansion e);
/// True when the expansion divides by A_0 (a_over_a0, rsl): the engine's base is A_L / A_0, else A_L.
bool LegendreExpansionNormalised(LegendreExpansion e);
/// The coefficient of order L in words, e.g. "A_L / ((2L+1) A_0)", A_L of dsigma/dOmega_cm = sum_L A_L P_L.
const char *LegendreExpansionFormula(LegendreExpansion e);
/// The factor on the base of order L (>= 0): 1 (a_over_a0, absolute), 1/(2L+1) (rsl), 2/(2L+1) (2l2), 4 pi (4pi).
double LegendreExpansionFactor(LegendreExpansion e, int order);

}  // namespace experiment
}  // namespace azure2

#endif
