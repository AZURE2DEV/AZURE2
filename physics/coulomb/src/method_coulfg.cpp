// Method::coulfg: A. R. Barnett's COULFG (coulfg/, C translation from LLNL), which returns F and G only.  The
// derivatives follow from F, G at l + 1 by the recurrence (Abramowitz and Stegun 14.2.1)
//   u'_l = [(l + 1)/rho + eta/(l + 1)] u_l - sqrt(1 + eta^2/(l + 1)^2) u_{l+1},   u = F or G,
// so COULFG is asked for l and l + 1.
#include "coulfg2.h"

#include <cmath>
#include <vector>

#include "methods.h"

namespace azure2 {
namespace physics {
namespace coulomb {
namespace detail {

Status coulfg(int l, double eta, double rho, CoulombWaves &out) {
  // COULFG indexes its arrays by l itself, from l to l + 1.
  std::vector<double> fc(l + 2, 0.0), gc(l + 2, 0.0);
  int m1 = 0;
  if (Coulomb_FG(rho, eta, (double)l, (double)(l + 1), fc.data(), gc.data(), &m1) != 0) return Status::failed;
  const double lp1 = l + 1.0;
  const double s = lp1 / rho + eta / lp1;
  const double r = std::hypot(1.0, eta / lp1);
  out.F = fc[l];
  out.G = gc[l];
  out.Fp = s * fc[l] - r * fc[l + 1];
  out.Gp = s * gc[l] - r * gc[l + 1];
  return Status::ok;
}

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
