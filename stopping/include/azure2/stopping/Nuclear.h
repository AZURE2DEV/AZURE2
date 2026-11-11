#ifndef AZURE2_STOPPING_NUCLEAR_H
#define AZURE2_STOPPING_NUCLEAR_H

#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

/*!
 * Nuclear stopping with the ZBL universal interatomic potential: J. F. Ziegler, J. P. Biersack and U. Littmark,
 * "The Stopping and Range of Ions in Solids" (Pergamon, 1985), the universal screening length
 * a_U = 0.8854 a0 / (Z1^0.23 + Z2^0.23) and the analytic fit of the reduced nuclear stopping
 * (docs/upgrade/stopping.md section 5):
 *
 *   epsilon = 32.53 M2 E / (Z1 Z2 (M1 + M2) (Z1^0.23 + Z2^0.23))          (E in keV, lab; masses in u)
 *   s_n(epsilon) = ln(1 + 1.1383 epsilon) / (2 (epsilon + 0.01321 epsilon^0.21226 + 0.19593 epsilon^0.5)),
 *                  epsilon <= 30;  ln(epsilon) / (2 epsilon) above
 *   S_n(E) = 8.462 Z1 Z2 M1 s_n(epsilon) / ((M1 + M2) (Z1^0.23 + Z2^0.23))   eV/(1e15 atoms/cm^2)
 */
constexpr double kNuclearRelUncertainty = 0.10;  ///< assigned (stopping.md section 7)

/// The reduced energy epsilon of the ZBL formula at lab energy \p energy_keV.
double ZblReducedEnergy(int z1, double m1_u, int z2, double m2_u, double energy_keV);
/// The reduced nuclear stopping s_n(epsilon) of the ZBL fit.
double ZblReducedStopping(double reduced_energy);
/// S_n of \p projectile on element \p z2 (mass \p m2_u) at lab energy \p energy_MeV, eV/(1e15 atoms/cm^2).
double NuclearStopping(const Projectile &projectile, int z2, double m2_u, double energy_MeV);

}  // namespace stopping
}  // namespace azure2

#endif
