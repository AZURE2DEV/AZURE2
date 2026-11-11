#ifndef AZURE2_STOPPING_STRAGGLING_H
#define AZURE2_STOPPING_STRAGGLING_H

#include <functional>
#include <string>
#include <vector>

#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

/*!
 * Energy-loss straggling (docs/upgrade/stopping.md section 10).  Variances in MeV^2, areal densities in atoms/cm^2,
 * energies in MeV (lab).  The specific variance dOmega^2/dN is per atom of the material (MeV^2 cm^2).
 *
 * Two models, as lunasrim/straggling.py has them:
 *   kTrim (the default): TRIM's straggling, W = w^2/12 per unit path with w = min(k E, w_Bohr), where
 *     w_Bohr = 177 Z1 (Z2 N)^1/2 eV/A^1/2 x ((1/beta^2 - 1/2) / (1/beta^2 - 1)) (Bohr with TRIM's relativistic factor;
 *     177^2/12 is TRIM's rounding of 4 pi e^4), k^2 = 3.2 x 12 N (Z1 + Z2)^(8/3) / (1000 M1) with E in eV (the
 *     low-energy limit; it acts only below a few hundred eV per u), N in atoms/A^3, and Z2 TRIM's mean target charge,
 *     an integer (TrimMeanTargetCharge).  Per atom: dOmega^2/dN = W / N, which does not depend on N.
 *   kBohr: N. Bohr, Phil. Mag. 30 (1915) 581: dOmega^2/dN = 4 pi z1^2 e^4 Z2 with e^2 = alpha hbar c and the exact
 *     atom-fraction-weighted Z2.
 *
 * Propagation (B. Tschalaer, Nucl. Instrum. Methods 61 (1968) 141; 64 (1968) 237), the solution of
 * dOmega^2/dx = W(E) - 2 S'(E) Omega^2 that lunasrim.straggling.energy_loss integrates step by step:
 *   sigma^2(x) = S(E_out)^2 [ sigma_0^2 / S(E_0)^2 + int_{E_out}^{E_0} (dOmega^2/dx)(E) / S(E)^3 dE ],
 * sigma_0 the beam's energy spread at the entrance.  With S constant it is sigma_0^2 + Omega^2 (the beam spread added
 * in quadrature); for a thin layer it is sigma_0^2 + dOmega^2/dx N.  The mean energy uses the total stopping; the
 * spread only the electronic straggling (nuclear collisions give a low-energy tail that a mean and a width do not
 * describe).
 */
enum class StragglingModel { kTrim, kBohr };

/// "trim", "bohr".
std::string StragglingModelName(StragglingModel model);
/// The inverse; throws std::invalid_argument for any other word.
StragglingModel StragglingModelFromName(const std::string &name);

/// dOmega_B^2/dN (MeV^2 cm^2): 4 pi z1^2 e^4 z2.
double BohrSpecificVariance(int z1, double z2);
/// The atom-fraction-weighted Z2 of a material (Bragg additivity of the Bohr variance).
double MeanTargetCharge(const Material &material);
/// TRIM's mean target charge: Z = round(f_i Z_i + Z) over the constituents in order, rounding half to even as
/// TRIM.exe (Visual Basic's CInt) and lunasrim do.
int TrimMeanTargetCharge(const Material &material);
/// TRIM's dOmega^2/dN (MeV^2 cm^2) of \p projectile at \p energy_MeV in \p material (see above).
double TrimSpecificVariance(const Material &material, const Projectile &projectile, double energy_MeV);

/// dOmega^2/dN (MeV^2 cm^2) of \p material for \p projectile at \p energy_MeV under \p model.
double SpecificStragglingVariance(StragglingModel model, const Material &material, const Projectile &projectile,
                                  double energy_MeV);

/// The widths after a layer (MeV, 1 sigma).
struct StragglingWidths {
  double energy_out_MeV = 0.0;
  double straggling = 0.0;  ///< from the layer alone, propagated (sigma_0 = 0)
  double beam = 0.0;        ///< the entrance spread sigma_0 carried through the layer: sigma_0 S(E_out)/S(E_0)
  double total = 0.0;       ///< sqrt(straggling^2 + beam^2)
};
/// Tschalaer's propagation for any epsilon(E) (eV/(1e15 atoms/cm^2)) and dOmega^2/dN(E) (MeV^2 cm^2); the integrals
/// split at \p breaks_MeV (EnergyLoss.h).
StragglingWidths PropagateStraggling(const StoppingFunction &epsilon, const std::function<double(double)> &specific_variance,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV,
                                     const std::vector<double> &breaks_MeV = {});
/// The same on a material (total stopping, its break points) with a straggling model.
StragglingWidths PropagateStraggling(StragglingModel model, const Material &material, const Projectile &projectile,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV);

}  // namespace stopping
}  // namespace azure2

#endif
