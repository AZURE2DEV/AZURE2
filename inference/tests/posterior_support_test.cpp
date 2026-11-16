// inference_posterior_support: the rule of the posterior's support (azure2/inference/PosteriorSupport.h): a fit
// vector is inside when every limited component lies in its closed interval.  A limit taken as open, a missing
// component taken as inside or a NaN taken as inside fails here.
#include "azure2/inference/PosteriorSupport.h"

#include <limits>
#include <vector>

#include "check.h"

using azure2::inference::PosteriorSupport;

int main() {
  const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
  PosteriorSupport support;
  check("no limit: every vector is inside", support.Contains({}) && support.Contains({1.0, -5.0}));
  // Component 0 in [0, +inf) (a normalization), component 2 in [2.36, 2.38] MeV (a level energy); 1 unlimited.
  support.limits = {{0, 0.0, inf}, {2, 2.36, 2.38}};
  check("inside", support.Contains({1.0, -7.0, 2.37}));
  check("on a limit is inside", support.Contains({0.0, 0.0, 2.36}) && support.Contains({0.0, 0.0, 2.38}));
  check("below a lower limit", !support.Contains({-1e-12, 0.0, 2.37}));
  check("above an upper limit", !support.Contains({1.0, 0.0, 2.3800001}));
  check("an unlimited component is free", support.Contains({1.0, 1e300, 2.37}) && support.Contains({1.0, -inf, 2.37}));
  check("a vector too short for a limit is outside", !support.Contains({1.0, 0.0}) && !support.Contains({}));
  check("a component that is not a number is outside", !support.Contains({nan, 0.0, 2.37}) && !support.Contains({1.0, 0.0, nan}));
  check("an infinite component against an infinite limit", support.Contains({inf, 0.0, 2.37}));
  return result();
}
