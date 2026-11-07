// One adaptor per algorithm, all with one signature (docs/upgrade/coulomb.md section 1.3).  Each fills F, Fp, G, Gp
// of out for one l and returns ok or failed; the arguments are already checked, and sigma and the NaN rule are set
// by Coulomb.cpp.
#ifndef AZURE2_PHYSICS_COULOMB_METHODS_H
#define AZURE2_PHYSICS_COULOMB_METHODS_H

#include "azure2/physics/coulomb/Coulomb.h"

namespace azure2 {
namespace physics {
namespace coulomb {
namespace detail {

Status coulsm(int l, double eta, double rho, CoulombWaves &out);
Status coul(int l, double eta, double rho, CoulombWaves &out);
Status coulfg(int l, double eta, double rho, CoulombWaves &out);
Status coulcc(int l, double eta, double rho, CoulombWaves &out);

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2

#endif
