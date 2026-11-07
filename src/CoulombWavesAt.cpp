// The engine's one call of the Coulomb functions (physics/coulomb, docs/upgrade/coulomb.md).  In a file of its own
// so that a user of CoulombTable alone (azure2-coulomb-table, its tests) does not link CoulFunc.
#include <sstream>
#include <stdexcept>

#include "CoulFunc.h"
#include "azure2/physics/coulomb/Coulomb.h"

CoulWaves CoulombWavesAt(int l, double eta, double rho) {
  const azure2::physics::coulomb::CoulombWaves w = azure2::physics::coulomb::coulomb(l, eta, rho);
  if (w.status != azure2::physics::coulomb::Status::ok) {
    std::ostringstream msg;
    msg << "Coulomb functions failed at l = " << l << ", eta = " << eta << ", rho = " << rho << " (physics/coulomb, "
        << azure2::physics::coulomb::method_name(azure2::physics::coulomb::Method::coulsm) << ")";
    throw std::runtime_error(msg.str());
  }
  return CoulWaves{w.F, w.Fp, w.G, w.Gp};
}

