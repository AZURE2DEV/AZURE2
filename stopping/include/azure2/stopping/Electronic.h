#ifndef AZURE2_STOPPING_ELECTRONIC_H
#define AZURE2_STOPPING_ELECTRONIC_H

#include "azure2/stopping/Material.h"
#include "azure2/stopping/Tables.h"

namespace azure2 {
namespace stopping {

/*!
 * Electronic stopping of an ion in one element (docs/upgrade/stopping.md section 4).  All values in
 * eV/(1e15 atoms/cm^2).
 *
 * Protons and helium-4: the low-energy fit of the coefficient table (Andersen-Ziegler 1977 / Ziegler 1977 forms as
 * adopted in ICRU Report 49, 1993) up to T1 = fit_max of the row (1 MeV/u), Bethe's formula with shell, Barkas and
 * Bloch corrections from T2 = kBetheFrom_MeV_u (3 MeV/u), and linear interpolation in log S - log T between them.
 * d and t use the proton row at the same velocity, 3He the helium row at the same velocity.  Heavier ions:
 * the proton value at the same velocity times the squared effective charge (Pierce and Blann 1968), with a larger
 * stated uncertainty.
 */
constexpr double kBetheFrom_MeV_u = 3.0;          ///< Bethe alone from this energy per nucleon (MeV/u)
constexpr double kBetheRelUncertainty = 0.02;     ///< assigned to the Bethe region (stopping.md section 7)
constexpr double kHeavyIonRelUncertainty = 0.20;  ///< effective-charge scaling, v >= v0 Z1^(2/3)
constexpr double kSlowHeavyIonRelUncertainty = 0.50;  ///< effective-charge scaling below that velocity

/// The proton fit (form AZ_p) at proton energy \p t_keV: velocity-proportional below 10 keV.
double ProtonFitStopping(const CoefficientRow &row, double t_keV);
/// The helium fit (form Z_He) at helium-4 energy \p t_MeV: velocity-proportional below 1 keV.
double HeliumFitStopping(const CoefficientRow &row, double t_MeV);

/// The terms of the stopping number L = L0 - C/Z2 + z1 L1 + z1^2 L2 (stopping.md section 4.2).
struct BetheTerms {
  double l0 = 0.0;      ///< ln(2 m c^2 beta^2 gamma^2 / I) - beta^2
  double shell = 0.0;   ///< C/Z2, Barkas and Berger (1964), held at eta = 0.13 below
  double barkas = 0.0;  ///< L1, Ashley, Ritchie and Brandt (1972)
  double bloch = 0.0;   ///< L2, Bloch (1933)
};
/// The terms at \p beta2 = (v/c)^2 for an ion of charge \p z1 on element \p z2 with mean excitation energy I (eV).
BetheTerms BetheStoppingTerms(int z1, int z2, double mean_excitation_eV, double beta2);
/// Bethe's stopping cross section K z1^2 Z2 / beta^2 x L, K = 4 pi r_e^2 m_e c^2 (eV/(1e15 atoms/cm^2)).
double BetheStopping(int z1, int z2, double mean_excitation_eV, double beta2);
/// Bloch's correction L2 = -y^2 sum_n 1/(n (n^2 + y^2)), y = z1 alpha / beta.
double BlochCorrection(double y);
/// Barkas and Berger's shell correction C (not divided by Z2) at eta = beta gamma >= 0.13, I in eV.
double ShellCorrectionBarkasBerger(double mean_excitation_eV, double eta);
/// Ashley-Ritchie-Brandt Barkas term L1 = F(b / x^1/2) / (Z2^1/2 x^3/2), x = beta^2 / (alpha^2 Z2).
double BarkasCorrection(int z2, double beta2);

/// The equilibrium charge fraction q/z1 = 1 - exp(-0.95 v_r), v_r = v / (v0 z1^(2/3)) (Pierce and Blann 1968).
double EffectiveChargeFraction(int z1, double beta2);

/// An electronic stopping value and its relative uncertainty.
struct StoppingWithUncertainty {
  double value = 0.0;            ///< eV/(1e15 atoms/cm^2)
  double rel_uncertainty = 0.0;
};
/// The electronic stopping of \p projectile in element \p z2 at lab energy \p energy_MeV.
StoppingWithUncertainty ElectronicStopping(const Projectile &projectile, int z2, double energy_MeV);

}  // namespace stopping
}  // namespace azure2

#endif
