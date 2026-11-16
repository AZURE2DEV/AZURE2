// The support of the posterior: the fit's limits of every parameter kind (azure2/inference/PosteriorSupport.h).
#include "azure2/inference/PosteriorSupport.h"

namespace azure2 {
namespace inference {

bool PosteriorSupport::Contains(const std::vector<double> &fitVector) const {
  for (const Limit &l : limits)
    if (l.index >= (int)fitVector.size() || !(fitVector[l.index] >= l.lower && fitVector[l.index] <= l.upper)) return false;
  return true;
}

}  // namespace inference
}  // namespace azure2
