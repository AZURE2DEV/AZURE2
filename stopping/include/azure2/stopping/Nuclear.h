#ifndef AZURE2_STOPPING_NUCLEAR_H
#define AZURE2_STOPPING_NUCLEAR_H

#include <array>

#include "azure2/stopping/Material.h"
#include "azure2/stopping/Tables.h"

namespace azure2 {
namespace stopping {

/*!
 * SRIM-2013's nuclear stopping (docs/upgrade/stopping.md section 5), ported from LunaSRIM 1.0.0
 * (lunasrim/stopping.py reduced_energy, nuclear_stopping): the ZBL reduced energy, and for every pair (Z1, Z2) the
 * fit of SNUC03.dat (the pair-specific coefficients SRIM's tables use), with ZBL's high-energy limit above reduced
 * energy 30:
 *
 *   eps = E 32.53 M2 / (Z1 Z2 (M1 + M2) (Z1^0.23 + Z2^0.23))                          (E in keV, lab; masses in u)
 *   s_n = 0.5 ln(1 + c1 eps) / (eps + c2 eps^c3 + c4 eps^0.5)   for eps <= 30;   ln(eps) / (2 eps) above
 *   S_n = s_n 8.462 Z1 Z2 M1 / ((M1 + M2) (Z1^0.23 + Z2^0.23))                      eV/(1e15 atoms/cm^2)
 *
 * The fit and the limit do not meet at eps = 30 (S_n jumps there by up to 18 %, for H on H): an integral over the
 * energy needs a break point there (StoppingBreakpoints in Stopping.h).  M2 is the caller's; Stopping() passes SRIM's
 * isotope-averaged atomic weight of the element (SCOEF03.DAT column 4).
 */
constexpr double kNuclearFitMaxReducedEnergy = 30.0;  ///< the SNUC03.dat fit up to here, ln(eps)/(2 eps) above
/// Provisional, not a quoted number: SRIM's published accuracy is for its total stopping (Electronic.h); no source
/// for the nuclear part alone was found (stopping.md section 7).
constexpr double kNuclearRelUncertainty = 0.10;

/// The ZBL reduced energy at lab energy \p energy_keV (lunasrim reduced_energy).
double ZblReducedEnergy(int z1, double m1_u, int z2, double m2_u, double energy_keV);
/// The reduced nuclear stopping s_n(eps) with SNUC03.dat coefficients \p c.
double SrimReducedNuclearStopping(const std::array<double, 4> &c, double reduced_energy);
/// S_n of \p projectile on element \p z2 of mass \p m2_u at lab energy \p energy_MeV, eV/(1e15 atoms/cm^2)
/// (lunasrim nuclear_stopping).
double NuclearStopping(const SrimData &data, const Projectile &projectile, int z2, double m2_u, double energy_MeV);
/// The same with SrimTables().
double NuclearStopping(const Projectile &projectile, int z2, double m2_u, double energy_MeV);
/// The lab energy (MeV) at which the reduced energy is 30, where S_n changes form.
double NuclearFormChange_MeV(const Projectile &projectile, int z2, double m2_u);

}  // namespace stopping
}  // namespace azure2

#endif
