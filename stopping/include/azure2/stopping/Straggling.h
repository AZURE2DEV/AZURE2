#ifndef AZURE2_STOPPING_STRAGGLING_H
#define AZURE2_STOPPING_STRAGGLING_H

#include <array>
#include <string>
#include <vector>

#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

/*!
 * Energy-loss straggling (docs/upgrade/stopping.md section 10).  Variances in MeV^2, areal densities in atoms/cm^2,
 * energies in MeV (lab).
 *
 * Bohr (N. Bohr, Phil. Mag. 30 (1915) 581): dOmega_B^2/dN = 4 pi z1^2 e^4 Z2 per atom, e^2 = alpha hbar c; for a
 * material the atom-fraction-weighted Z2.
 *
 * Chu (W. K. Chu, Phys. Rev. A 13 (1976) 2057): Omega^2 = H(E/M1, Z2) Omega_B^2, H from Chu's Hartree-Fock-Slater
 * calculation.  H is data, read from a ChuTable (E/M1 in keV/u x Z2, linear interpolation in E/M1); no table is
 * shipped in this branch (stopping.md section 10.2).
 *
 * Charge-state term (Q. Yang, D. J. O'Connor, Z. Wang, Nucl. Instrum. Methods B 61 (1991) 149):
 *   Omega^2 / Omega_B^2 = gamma^2 H + (Delta Omega / Omega_B)^2,
 *   hydrogen ions (gamma = 1):  (Delta)^2 = B1 G / ((E - B2)^2 + G^2),  G = B3 (1 - exp(-B4 E)),  E in MeV/u;
 *   helium and heavier:  (Delta)^2 = Z1^(4/3) / Z2^(1/3) x C1 G / ((e - C2)^2 + G^2),  G = C3 (1 - exp(-C4 e)),
 *                        e = E / (Z1^(3/2) Z2^(1/2)), E in MeV/u.
 * The constants B1..B4, C1..C4 are data (YangCoefficients); none are shipped in this branch.
 *
 * Propagation (B. Tschalaer, Nucl. Instrum. Methods 61 (1968) 141; 64 (1968) 237): with S(E) the stopping and x the
 * areal density travelled,
 *   sigma^2(x) = S(E(x))^2 [ sigma_0^2 / S(E_0)^2 + int_0^x (dOmega^2/dx') / S(E(x'))^2 dx' ]
 *             = S(E_out)^2 [ sigma_0^2 / S(E_0)^2 + int_{E_out}^{E_0} (dOmega^2/dx)(E) / S(E)^3 dE ],
 * where sigma_0 is the beam's energy spread at the entrance.  With S constant it is sigma_0^2 + Omega^2 (the beam
 * spread added in quadrature); for a thin layer it is sigma_0^2 + dOmega^2/dx N.
 */

/// dOmega_B^2/dN (MeV^2 cm^2): 4 pi z1^2 e^4 z2.
double BohrSpecificVariance(int z1, double z2);
/// The atom-fraction-weighted Z2 of a material (Bragg additivity of the Bohr variance).
double MeanTargetCharge(const Material &material);

/// Chu's factor H(E/M1, Z2) as a table: for each Z2 a row of H on the common grid energy_keV_u (ascending).
struct ChuTable {
  std::vector<double> energy_keV_u;
  std::vector<int> z2;
  std::vector<std::vector<double>> h;  ///< h[i][k]: Z2 = z2[i], E/M1 = energy_keV_u[k]
  std::string source;
};
/// H at \p energy_keV_u for target Z2: linear in E/M1, 1 above the grid (Bohr's limit), the first value below it.
/// Throws std::out_of_range if the table has no row for z2.
double ChuFactor(const ChuTable &table, int z2, double energy_keV_u);

/// The constants of Yang, O'Connor and Wang (1991): b for hydrogen ions, c for helium and heavier.
struct YangCoefficients {
  std::array<double, 4> b{};
  std::array<double, 4> c{};
  std::string source;
};
/// (Delta Omega / Omega_B)^2 of the charge-state term for projectile charge \p z1 on target \p z2 at \p energy_MeV_u.
double YangChargeStateTerm(const YangCoefficients &yang, int z1, int z2, double energy_MeV_u);

/// Which straggling to use, with its data (the caller owns the tables; null pointers where not needed).
enum class StragglingKind { kBohr, kChu, kChuYang };
struct StragglingModel {
  StragglingKind kind = StragglingKind::kBohr;
  const ChuTable *chu = nullptr;          ///< needed by kChu and kChuYang
  const YangCoefficients *yang = nullptr; ///< needed by kChuYang
};

/// dOmega^2/dN (MeV^2 cm^2) of \p material for \p projectile at \p energy_MeV under \p model (constituents added
/// with their atom fractions).  Throws std::invalid_argument if the model lacks its table.
double SpecificStragglingVariance(const StragglingModel &model, const Material &material, const Projectile &projectile,
                                  double energy_MeV);

/// The widths after a layer (MeV, 1 sigma).
struct StragglingWidths {
  double energy_out_MeV = 0.0;
  double straggling = 0.0;  ///< from the layer alone, propagated (sigma_0 = 0)
  double beam = 0.0;        ///< the entrance spread sigma_0 carried through the layer: sigma_0 S(E_out)/S(E_0)
  double total = 0.0;       ///< sqrt(straggling^2 + beam^2)
};
/// Tschalaer's propagation for any epsilon(E) (eV/(1e15 atoms/cm^2)) and dOmega^2/dN(E) (MeV^2 cm^2).
StragglingWidths PropagateStraggling(const StoppingFunction &epsilon, const std::function<double(double)> &specific_variance,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV);
/// The same on a material with a straggling model.
StragglingWidths PropagateStraggling(const StragglingModel &model, const Material &material, const Projectile &projectile,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV);

}  // namespace stopping
}  // namespace azure2

#endif
