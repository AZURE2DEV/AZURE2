#ifndef AZURE2_STOPPING_STOPPING_H
#define AZURE2_STOPPING_STOPPING_H

#include <string>

#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

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
 * this projectile replaces the evaluation.  Otherwise Bragg additivity over the constituents,
 * epsilon = sum_i f_i (S_e,i + S_n,i), times the compound correction factor (1 in this branch).  The uncertainty of
 * a compound treats the constituents' uncertainties as fully correlated (their sum, weighted), and adds the nuclear
 * part's in quadrature.
 */
StoppingValue Stopping(const Material &material, const Projectile &projectile, double energy_MeV);

/*!
 * The core-and-bond correction of compound stopping (Ziegler and Manoyan, Nucl. Instrum. Methods B 35 (1988) 215),
 * the factor epsilon_compound / epsilon_Bragg.  Not implemented in this branch: it returns 1 for every material, and
 * the place is here so that the next branch changes this one function (docs/upgrade/stopping.md section 6).
 */
double CompoundCorrectionFactor(const Material &material, const Projectile &projectile, double energy_MeV);

/// The UserStopping of \p material for \p projectile, or nullptr if there is none.
const UserStopping *FindUserStopping(const Material &material, const Projectile &projectile);

/// A text that says where epsilon of this material and projectile comes from: the user source, or per element the
/// table rows (with their origin line), I, the Bethe corrections and the nuclear stopping, and the uncertainties.
std::string StoppingProvenance(const Material &material, const Projectile &projectile);

}  // namespace stopping
}  // namespace azure2

#endif
