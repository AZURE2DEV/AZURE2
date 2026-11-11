#ifndef AZURE2_STOPPING_MATERIAL_H
#define AZURE2_STOPPING_MATERIAL_H

#include <functional>
#include <string>
#include <vector>

#include "azure2/stopping/Tables.h"

namespace azure2 {
namespace stopping {

/*!
 * Units of the whole module, stated once (docs/upgrade/stopping.md section 2):
 *   energies              MeV, kinetic, laboratory frame, of the projectile;
 *   stopping cross sections  eV/(1e15 atoms/cm^2), per atom of the material (all atoms counted);
 *   areal densities       atoms/cm^2 (all atoms of the material);
 *   straggling variances  MeV^2.
 * An energy loss is Delta E [MeV] = epsilon [eV/(1e15 atoms/cm^2)] x N [atoms/cm^2] x 1e-21.
 */
constexpr double kStoppingToMeVcm2 = 1.0e-21;  ///< eV/(1e15 atoms/cm^2) -> MeV cm^2

/// A projectile: charge number and atomic mass (u).
struct Projectile {
  int z = 1;
  double mass_u = 1.00782503;
  std::string label = "p";
};

/// The light projectiles, atomic masses of AME2020 (u).
Projectile Proton();
Projectile Deuteron();
Projectile Triton();
Projectile Helium3();
Projectile Alpha();
/// Any other ion: charge number \p z >= 1, atomic mass \p mass_u (u).
Projectile Ion(int z, double mass_u, const std::string &label);

/// beta^2 = (v/c)^2 of a projectile of kinetic energy \p energy_MeV (relativistic).
double Beta2(const Projectile &projectile, double energy_MeV);

/// One element of a material and its share of the atoms.
struct Constituent {
  int z = 0;
  double atom_fraction = 0.0;  ///< normalised so that the fractions of a material sum to 1
};

/*!
 * A measured or user-given stopping cross section that replaces the evaluation for one projectile in one material
 * (stopping.md section 6): a table (energies ascending, MeV; values in eV/(1e15 atoms/cm^2), interpolated linearly
 * in log E - log epsilon, refused outside its range) or a formula epsilon(E_MeV).  Exactly one of the two is given.
 */
struct UserStopping {
  int projectile_z = 1;
  double projectile_mass_u = 1.00782503;    ///< matched to the projectile within 0.01 u
  std::vector<double> energy_MeV;           ///< table: lab energies (MeV), ascending
  std::vector<double> stopping;             ///< table: eV/(1e15 atoms/cm^2), per atom of the material
  std::function<double(double)> formula;    ///< formula: epsilon(E [MeV]) in eV/(1e15 atoms/cm^2)
  double rel_uncertainty = 0.0;             ///< stated relative 1-sigma uncertainty
  std::string provenance;                   ///< where the numbers come from (shown by StoppingProvenance)
};

/*!
 * A target material: elements with atomic fractions, phase, optional mass density (g/cm^3; 0 = not given, needed
 * only to convert to or from lengths), and user stopping that overrides the evaluation for given projectiles.
 */
struct Material {
  std::string name;
  std::vector<Constituent> constituents;
  Phase phase = Phase::kSolid;
  double mass_density_g_cm3 = 0.0;
  std::vector<UserStopping> user_stopping;
};

/// Parse a chemical formula ("Ta2O5", "TiN", "H2O", "C", "Si3N4", "CH4"; element symbols as physics/'s
/// ElementSymbol, counts optional and possibly decimal: "Ti0.5Al0.5N") into normalised atomic fractions.  Repeated
/// elements are summed.  Throws std::invalid_argument for an unknown symbol, an element above Z = 92 or a bad count.
std::vector<Constituent> ParseFormula(const std::string &formula);
/// A material from a formula; the name is the formula.
Material MaterialFromFormula(const std::string &formula, Phase phase, double mass_density_g_cm3);
/// A material from explicit fractions (normalised here; each z in 1..92, fractions > 0).
Material MaterialFromFractions(const std::string &name, const std::vector<Constituent> &constituents, Phase phase,
                               double mass_density_g_cm3);
/// The atom-fraction-weighted standard atomic weight (u): mass per atom of the material.
double MeanAtomicWeight(const Material &material);
/// The atom fraction of element \p z in the material (0 if absent): divide a per-atom stopping cross section by it
/// to get one per atom of that element.
double AtomFraction(const Material &material, int z);

/*!
 * Read a two-column table into a UserStopping.  The first non-empty line that is not a comment names the units:
 *   "E[MeV] eps[eV/1e15atoms/cm2]"  (energy unit MeV or keV; stopping unit eV/1e15atoms/cm2, MeV*cm2/g or keV/um)
 * then one "energy value" pair per line; '#' starts a comment.  MeV*cm2/g is converted with the material's mean
 * atomic weight; keV/um needs the material's mass density.  Throws std::invalid_argument for an unknown unit, a
 * missing density, fewer than two rows, energies not ascending or a value <= 0.
 */
UserStopping ReadStoppingTable(const std::string &text, const Material &material, const Projectile &projectile,
                               double rel_uncertainty, const std::string &provenance);

}  // namespace stopping
}  // namespace azure2

#endif
