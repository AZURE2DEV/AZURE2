#ifndef AZURE2_STOPPING_STOPPING_H
#define AZURE2_STOPPING_STOPPING_H

#include <string>
#include <vector>

#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

/// The citation every use of these stopping powers owes (stopping/data/srim/README.md).
constexpr const char *kSrimCitation =
    "J. F. Ziegler, M. D. Ziegler and J. P. Biersack, SRIM - The stopping and range of ions in matter (2010), "
    "Nucl. Instrum. Methods Phys. Res. B 268 (2010) 1818, doi:10.1016/j.nimb.2010.02.091 (www.srim.org)";
/// What the code is: SRIM-2013's stopping as re-implemented by LunaSRIM, ported to C++.
constexpr const char *kSrimImplementation =
    "SRIM-2013 stopping powers, ported to C++ from LunaSRIM 1.0.0 (J. Skowronski, github.com/skowrons94/LunaSRIM, MIT)";

/// The stopping cross section of a material, per atom of the material (eV/(1e15 atoms/cm^2)).
struct StoppingValue {
  double electronic = 0.0;
  double nuclear = 0.0;
  double total = 0.0;            ///< electronic + nuclear (a user override is all in total and electronic)
  double rel_uncertainty = 0.0;  ///< relative 1-sigma of total (stopping.md section 7)
  bool from_user = false;        ///< a UserStopping of the material gave it
};

/*!
 * epsilon(E) of \p material for \p projectile at lab energy \p energy_MeV (> 0).  A UserStopping of the material for
 * this projectile replaces the evaluation.  Otherwise SRIM's Bragg sum (stopping.md section 6):
 *   S_e = sum_i f_i S_e,i x FadeOut(PhaseFactor_i x C),   S_n = sum_i f_i S_n,i(M_i),
 * f_i the atom fractions, S_e,i ElectronicStopping in the material's phase, C the material's compound correction,
 * M_i SRIM's atomic weight of element i.  The compound correction enters once, as in SRIM's stopping tables, TRIM
 * and lunasrim.straggling.energy_loss (lunasrim.stopping.compound_stopping, which follows the list mode of SRIM's SR
 * module, applies it a second time below 10 MeV; stopping.md section 6).  The uncertainty: the electronic parts are
 * fully correlated (one evaluation), the nuclear part is added in quadrature.  Throws std::invalid_argument for a
 * material without constituents or with a compound correction <= 0.
 */
StoppingValue Stopping(const Material &material, const Projectile &projectile, double energy_MeV);

/// The energies (MeV, ascending) where epsilon(E) of the material changes form: the branch boundaries of the
/// electronic stopping (2 keV/u, the end of the light-ion fit and 1.2 times it, 30 MeV/u for heavy ions; continuous,
/// with a kink) and the reduced energy 30 of each element's nuclear stopping (a jump).  The integrals of
/// EnergyLoss.h split there.  Empty for a material with user stopping for this projectile.
std::vector<double> StoppingBreakpoints(const Material &material, const Projectile &projectile);

/// The UserStopping of \p material for \p projectile, or nullptr if there is none.
const UserStopping *FindUserStopping(const Material &material, const Projectile &projectile);

/// A text that says where epsilon of this material and projectile comes from: SRIM-2013 with its citation, LunaSRIM
/// as the implementation, the data source, the phase and compound correction used, per element the SCOEF03.DAT row,
/// mean excitation energy, mass, phase factor and the uncertainties; or the user source.
std::string StoppingProvenance(const Material &material, const Projectile &projectile);

}  // namespace stopping
}  // namespace azure2

#endif
