#ifndef AZURE2_STOPPING_TABLES_H
#define AZURE2_STOPPING_TABLES_H

#include <array>
#include <string>
#include <vector>

namespace azure2 {
namespace stopping {

/*!
 * The data tables of the stopping module (docs/upgrade/stopping.md sections 3 and 4).  Both are text files in
 * stopping/data/, compiled into the library at build time; every row carries its source.
 *
 *   elements.tsv                 Z = 1..92: standard atomic weight (u), mean excitation energy I (eV), phase.
 *   electronic_coefficients.tsv  per element and projectile (p, He): the low-energy fit coefficients A1..A5, the
 *                                energy per nucleon up to which they are used alone, the assigned relative
 *                                uncertainty, the source and a note.
 *
 * ElementTable() and CoefficientTable() parse the compiled text once and return a const reference (the constant
 * tables of docs/upgrade/ARCHITECTURE.md section 9).  The parsers are public so that a test or a user can read a
 * table of the same format from any text.
 */

/// The state of a material, or of an element at room temperature and pressure.
enum class Phase { kSolid, kLiquid, kGas };

/// "solid", "liquid", "gas".
std::string PhaseName(Phase phase);
/// The inverse of PhaseName; throws std::invalid_argument for any other word.
Phase PhaseFromName(const std::string &name);

/// One row of elements.tsv.
struct ElementRecord {
  int z = 0;
  std::string symbol;
  double atomic_weight_u = 0.0;       ///< standard atomic weight A_r (u)
  double atomic_weight_unc_u = 0.0;   ///< its stated interval (u)
  std::string atomic_weight_source;
  double mean_excitation_eV = 0.0;    ///< I (eV), ICRU 37 / ICRU 49
  Phase phase = Phase::kSolid;        ///< state at room temperature and pressure
};

/// The two fit forms of the low-energy electronic stopping (stopping.md section 4.1).
enum class FitForm {
  kAndersenZieglerProton,  ///< "AZ_p": protons, T in keV
  kZieglerHelium           ///< "Z_He": helium-4, T in MeV
};

/// One row of electronic_coefficients.tsv.
struct CoefficientRow {
  int z = 0;                       ///< target element
  int projectile_z = 0;            ///< 1 (p) or 2 (He)
  FitForm form = FitForm::kAndersenZieglerProton;
  std::array<double, 5> a{};       ///< A1..A5 as published
  double fit_max_MeV_u = 0.0;      ///< energy per nucleon up to which the row is used alone (MeV/u)
  double rel_uncertainty = 0.0;    ///< assigned relative 1-sigma uncertainty of the row
  std::string source;
  std::string note;                ///< "-" when none
  std::string origin;              ///< "<file>:<line>" the row was read from
};

/// Parse elements.tsv text (comment lines start with '#', one header line, tab-separated columns).  Throws
/// std::invalid_argument naming \p origin and the line for a malformed row.
std::vector<ElementRecord> ParseElementTable(const std::string &text, const std::string &origin);
/// Parse electronic_coefficients.tsv text; throws std::invalid_argument as above.
std::vector<CoefficientRow> ParseCoefficientTable(const std::string &text, const std::string &origin);

/// The compiled elements.tsv, index z - 1.
const std::vector<ElementRecord> &ElementTable();
/// The compiled electronic_coefficients.tsv.
const std::vector<CoefficientRow> &CoefficientTable();

/// The record of element z; throws std::out_of_range outside 1..92.
const ElementRecord &Element(int z);
/// The coefficient row of target z for projectile charge 1 or 2; throws std::out_of_range if there is none.
const CoefficientRow &Coefficients(int z, int projectile_z);

}  // namespace stopping
}  // namespace azure2

#endif
