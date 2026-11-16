#ifndef AZURE2_INFERENCE_POSTERIORSUPPORT_H
#define AZURE2_INFERENCE_POSTERIORSUPPORT_H

#include <vector>

namespace azure2 {
namespace inference {

/*!
 * The support of the posterior in the parameters (docs/upgrade/objective.md section 6, "The support"): the fit vectors
 * whose free components lie inside the limits the fit gives them, of every parameter kind.  The one rule of the
 * ensemble sampler (both sampling bases: it tests the fit vector it evaluates the objective at) and of any other
 * log-probability of the objective.
 *
 * Plain data: a list of limits by full index, each in the unit of its parameter's kind (ParamKindTable()).  The
 * engine makes it from the fit's parameters (PosteriorSupportOfFit, PosteriorSupportOfParameters: the engine's
 * PosteriorSupport.h).
 */
struct PosteriorSupport {
  /// The limits of one free parameter: its index in the full vector, lower and upper (+-infinity where it has none).
  struct Limit {
    int index;
    double lower, upper;
  };
  std::vector<Limit> limits;

  /// Whether the full fit vector `fitVector` lies inside every limit (a vector too short for a limit is outside, and
  /// so is a component that is not a number).
  bool Contains(const std::vector<double> &fitVector) const;
};

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_POSTERIORSUPPORT_H
