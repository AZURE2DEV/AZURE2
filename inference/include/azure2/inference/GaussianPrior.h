#ifndef AZURE2_INFERENCE_GAUSSIANPRIOR_H
#define AZURE2_INFERENCE_GAUSSIANPRIOR_H

namespace azure2 {
namespace inference {

/*!
 * One Gaussian prior on one parameter x with centre x0 and width sigma > 0, all three in the unit of the parameter's
 * kind (ParamKindTable(): MeV for a level energy or an energy shift, eV for U_e, atoms/cm^2 for an areal density, ...):
 *
 *     chi2 = ((x - x0) / sigma)^2,    d chi2/dx = 2 (x - x0) / sigma^2,    d^2 chi2/dx^2 = 2 / sigma^2.
 *
 * The one statement of the formula for the objective's diagonal priors (normalizations, energy shifts, screening,
 * density, Doppler temperature, thickness, apparatus parameters, abundances and their sum, the nuisance priors),
 * their gradients, Minuit's G2 and the chiSquared.out writer.  Which parameters carry such a prior, and whether one is replaced by the prior covariance,
 * is the engine's gather (include/DiagonalPriors.h; the list: DiagonalPriors.h).  sigma = 0 is "no prior": the caller does not call.
 *
 * The value is the square of the pull r = (x - x0)/sigma, formed as r * r (one division, one product).
 */
inline double GaussianPriorChi2(double x, double centre, double sigma) {
  const double r = (x - centre) / sigma;
  return r * r;
}

/// d chi2/dx of GaussianPriorChi2.
inline double GaussianPriorGradient(double x, double centre, double sigma) { return 2.0 * (x - centre) / (sigma * sigma); }

/// d^2 chi2/dx^2 of GaussianPriorChi2 (it does not depend on x).
inline double GaussianPriorCurvature(double sigma) { return 2.0 / (sigma * sigma); }

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_GAUSSIANPRIOR_H
