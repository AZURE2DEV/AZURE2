#ifndef AZURE2_STOPPING_ELECTRONIC_H
#define AZURE2_STOPPING_ELECTRONIC_H

#include <array>

#include "azure2/stopping/Material.h"
#include "azure2/stopping/Tables.h"

namespace azure2 {
namespace stopping {

/*!
 * SRIM-2013's electronic stopping (docs/upgrade/stopping.md section 4), ported from LunaSRIM 1.0.0
 * (lunasrim/stopping.py; the function each one comes from is named at it).  Inside, SRIM's variable is the energy per
 * mass unit e = E [keV] / M1 [u] in keV/u, M1 the caller's projectile mass; stopping cross sections are in
 * eV/(1e15 atoms/cm^2) per atom of the target element.  The public function ElectronicStopping takes the library's
 * units (MeV, lab).  Arithmetic is std::pow, std::log, std::exp: LunaSRIM's vbmath (Visual Basic's x87 arithmetic) is
 * not reproduced, so the agreement with LunaSRIM is to rounding (stopping.md section 9).
 *
 *   H (Z1 = 1):  SRIM's 10-parameter fit from 2 keV/u to 1000 keV/u, linear in e to Bethe at 1.2 x 1000 keV/u, Bethe
 *                above; below 2 keV/u proportional to the velocity.
 *   He, Li:      the same with their own coefficients, fit to 2000 keV/u, Bethe above 2400 keV/u.
 *   Z1 >= 4:     helium at the same e times the ion polynomial (of Z1) and the target polynomial (of Z2) in ln e, and
 *                Z1^2 / 4, for 2 .. 30000 keV/u; velocity-proportional below; above 30 MeV/u scaled with helium's
 *                Bethe value and SRIM's high-energy effective charge.
 *   Phase:       a gas target uses the gas mean excitation energy and no second shell term in Bethe, and the
 *                low-energy phase factor 0.92 (or 1/0.92) applies when the phase differs from the element's natural
 *                one (not for noble gases), faded out with velocity together with the compound correction.
 */
constexpr double kSrimVelocityLimit_keV_u = 2.0;           ///< velocity-proportional below (all ions)
constexpr double kProtonFitMax_keV_u = 1000.0;             ///< end of the proton fit
constexpr double kHeliumLithiumFitMax_keV_u = 2000.0;      ///< end of the helium and lithium fits
constexpr double kBetheJoinFactor = 1.2;                   ///< light ions: Bethe alone above 1.2 x the end of the fit
constexpr double kHeavyIonScalingMax_keV_u = 30000.0;      ///< heavy ions: helium scaling up to here
constexpr double kSrimMaxEnergy_keV_u = 1.0e7;             ///< refused above (SRIM's stopping is 0 there)
constexpr double kSrimPhaseFactor = 0.92;                  ///< lunasrim/stopping.py PHASE_FACTOR

/// SRIM's 10-parameter fit at \p e_keV_u (lunasrim _fit).
double SrimFit(const std::array<double, 10> &c, double e_keV_u);
/// Bethe-Bloch with SRIM's shell, Barkas and Bloch terms for ion charge \p z1 in \p target (lunasrim bethe).
double SrimBethe(const SrimElement &target, int z1, double e_keV_u, Phase phase);
/// H, He, Li: the fit \p c up to \p fit_max_keV_u, linear join to Bethe, velocity-proportional below 2 keV/u
/// (lunasrim _light_ion).
double SrimLightIon(const SrimElement &target, const std::array<double, 10> &c, int z1, double e_keV_u,
                    double fit_max_keV_u, Phase phase);
/// The high-energy effective-charge factor of a heavy ion with parameters \p hi (lunasrim _hi_charge).
double SrimHighEnergyCharge(const std::array<double, 5> &hi, double e_keV_u);
/// Ion \p z1 in element \p z2 without the phase factor (lunasrim electronic_stopping_element).
double SrimElementStopping(const SrimData &data, int z1, int z2, double e_keV_u, Phase phase);
/// The low-energy phase factor: 1 for noble gases and for an element in its natural phase, 0.92 for a natural gas
/// as a solid, 1/0.92 for a natural solid as a gas (lunasrim phase_factor).
double SrimPhaseFactor(int z2, Phase phase);
/// \p factor applied fully at low velocity and faded out at high velocity (around 1.2 MeV/u for H, 400 keV/u for
/// heavier ions): 1 + (factor - 1) / (1 + exp(1.48 (sqrt(e/25) - c))), c = 7 for H, 4 otherwise
/// (lunasrim fade_correction).
double SrimFadeOut(double factor, double e_keV_u, int z1);

/// An electronic stopping value and its relative uncertainty.
struct StoppingWithUncertainty {
  double value = 0.0;            ///< eV/(1e15 atoms/cm^2)
  double rel_uncertainty = 0.0;
};

/*!
 * The electronic stopping of \p projectile in element \p z2 at lab energy \p energy_MeV in a target of phase \p phase
 * whose compound correction is \p compound_correction (1 for Bragg's rule): SrimElementStopping times
 * SrimFadeOut(SrimPhaseFactor x compound_correction) (lunasrim electronic_stopping).  The uncertainty is SRIM's
 * published accuracy for the projectile (ElectronicRelUncertainty).  Throws std::domain_error for E <= 0 or above
 * kSrimMaxEnergy_keV_u, std::out_of_range for Z1 or Z2 outside 1..92, std::invalid_argument for a compound
 * correction <= 0 (LunaSRIM reads 0 as 1, the convention of SRIM's SR.IN; here it must be said).
 */
StoppingWithUncertainty ElectronicStopping(const SrimData &data, const Projectile &projectile, int z2, double energy_MeV,
                                           Phase phase, double compound_correction);
/// The same with SrimTables().
StoppingWithUncertainty ElectronicStopping(const Projectile &projectile, int z2, double energy_MeV, Phase phase,
                                           double compound_correction);

/*!
 * SRIM's published accuracy of its electronic stopping: the average deviation of the experimental stopping powers
 * from SRIM-2010's (the stopping of SRIM-2013), in the table "Accuracy of SRIM Stopping Power Calculations" of
 * www.srim.org/SRIM/SRIM2011.htm (read 2026-11-11): H ions 4.0 % (about 8300 data points), He ions 3.9 % (6500),
 * Li ions 4.8 % (1400), Be-U ions 5.8 % (9000).  Used as the relative 1-sigma uncertainty.
 */
double ElectronicRelUncertainty(int z1);

}  // namespace stopping
}  // namespace azure2

#endif
