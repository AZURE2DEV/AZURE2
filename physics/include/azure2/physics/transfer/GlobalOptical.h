#ifndef AZURE2_PHYSICS_TRANSFER_GLOBALOPTICAL_H
#define AZURE2_PHYSICS_TRANSFER_GLOBALOPTICAL_H

#include <string>
#include <vector>

#include "azure2/physics/transfer/OpticalPotential.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The global optical sets: one table, one row per published parametrisation (docs/upgrade/transfer.md, section "The
 * engine's THM observable on this library").  A row gives, for a light projectile (Z_p, A_p) on a target (Z_t, A_t)
 * at the projectile's LAB energy E (MeV), the OpticalPotential of Eq. (2.1) with the TargetMass radius convention
 * (R = r0 A_t^{1/3}, the convention of every set below) and the uniform-sphere Coulomb radius parameter of the set.
 * The surface term of every set is the derivative form normalised as -4 a_D W_D df/dr, which is Eq. (2.1)'s, so
 * W_D maps one to one.  A negative imaginary depth of a formula (Liang's volume term below 23 MeV) is set to 0, and a
 * term of depth 0 carries no geometry.
 *
 * Central terms only: the spin-orbit terms of the published sets are not in the rows (the THM weight is spinless;
 * a row gains its spin-orbit term with its primary source or an independent code's table and a test).
 *
 * A new set is one row of GlobalOpticalSets() (src/transfer/GlobalOptical.cpp) and one test in
 * physics/tests/transfer_global_optical_test.cpp (the recipe "Add a global optical set" of CLAUDE.md).
 */

/// A light ion a set describes, by charge and mass number.
struct LightIon {
  int z = 0;
  int a = 0;
};

/// One row of the table.
struct GlobalOpticalSet {
  const char *name;                ///< the name users write (opticalAA=<name>)
  const char *projectiles;         ///< the light ions in words ("d", "n p", ...)
  std::vector<LightIon> describes; ///< the same as (Z, A)
  const char *reference;           ///< the publication
  double a_min, a_max;             ///< target mass numbers the set was fitted to
  double e_min_MeV, e_max_MeV;     ///< projectile lab energies the set was fitted to (MeV)
  /// The potential for (z_p, a_p) on (z_t, a_t) at the projectile lab energy (MeV), evaluated wherever asked (the
  /// validity range is the caller's business).
  OpticalPotential (*at)(int z_p, int a_p, int z_t, int a_t, double e_lab_MeV);
};

/// The table, in the order the documentation lists the sets (constant).
const std::vector<GlobalOpticalSet> &GlobalOpticalSets();
/// The index of a set by name, or -1.
int GlobalOpticalIndex(const std::string &name);
/// The names, separated by ", ".
std::string GlobalOpticalNames();
/// Whether the set is a potential for the projectile (z_p, a_p).
bool Describes(const GlobalOpticalSet &set, int z_p, int a_p);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
