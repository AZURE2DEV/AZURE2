#ifndef AZURE2_STOPPING_TABLES_H
#define AZURE2_STOPPING_TABLES_H

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace azure2 {
namespace stopping {

/*!
 * SRIM-2013's data files (docs/upgrade/stopping.md section 3; stopping/data/srim/README.md): SCOEF03.DAT, one row
 * per element Z = 1..92 with the electronic stopping coefficients and the element's masses, densities and mean
 * excitation energies, and SNUC03.dat, the coefficients of the nuclear stopping fit for every pair (Z1, Z2).  They
 * are J. F. Ziegler's (www.srim.org), copied unchanged from LunaSRIM 1.0.0 (github.com/skowrons94/LunaSRIM), whose
 * reader (lunasrim/data_files.py) documents the column layout this parser follows.
 *
 * SrimTables() is the constant table of docs/upgrade/ARCHITECTURE.md section 9, parsed once: from the directory named
 * by the environment variable AZURE2_SRIM_DATA_DIR if it is set (in the process's environment record,
 * azure2/support/Environment.h), otherwise from the copy compiled into the library (when the build found
 * stopping/data/srim/; otherwise it refuses and says how to get the files).  ReadSrimData() reads a directory the
 * caller names into a table the caller owns; SrimTablesIn() is the table of a directory the caller names (the engine
 * names the run's, Config::environment.srimDataDirectory).
 */

/// The state of a material: SRIM distinguishes solid and gas targets (mean excitation energy, shell correction,
/// low-energy phase factor).  A liquid is a solid for SRIM.
enum class Phase { kSolid, kGas };

/// "solid", "gas".
std::string PhaseName(Phase phase);
/// The inverse of PhaseName; throws std::invalid_argument for any other word.
Phase PhaseFromName(const std::string &name);

/// One row of SCOEF03.DAT (columns numbered as lunasrim/data_files.py numbers them, 1-based after the two names).
struct SrimElement {
  int z = 0;                              ///< column 1
  std::string symbol;                     ///< the first quoted name, without padding
  std::string name;                       ///< the second quoted name, without padding
  double most_abundant_mass_u = 0.0;      ///< column 3: mass of the most abundant isotope (u)
  double atomic_weight_u = 0.0;           ///< column 4: isotope-averaged atomic weight (u)
  double solid_density_g_cm3 = 0.0;       ///< column 5
  double gas_density_g_cm3 = 0.0;         ///< column 9
  std::array<double, 10> proton{};        ///< columns 11-20: proton stopping fit
  std::array<double, 10> helium{};        ///< columns 21-30: helium stopping fit
  std::array<double, 10> lithium{};       ///< columns 31-40: lithium stopping fit
  std::array<double, 13> ion_poly{};      ///< columns 41-53: heavy-ion scaling polynomial (this element as the ion)
  std::array<double, 13> target_poly{};   ///< columns 54-66: heavy-ion scaling polynomial (this element as target)
  double mean_excitation_solid_eV = 0.0;  ///< column 67
  std::array<double, 2> mean_excitation_gas_eV{};  ///< columns 68-69: gas, for Z1 = 1 and for Z1 > 1
  std::array<double, 11> shell_c{};       ///< columns 70-80: shell correction polynomial in ln(E [MeV/u])
  std::array<double, 8> shell_d{};        ///< columns 81-88: second shell polynomial (solids only)
  std::array<double, 5> high_energy_charge{};  ///< columns 89-93: high-energy effective charge (this element as ion)
  std::string origin;                     ///< "<file>:<line>" the row was read from
};

/// SRIM's data: the 92 element rows and the nuclear stopping coefficients of every pair.
struct SrimData {
  std::vector<SrimElement> elements;            ///< index z - 1
  std::vector<std::array<double, 4>> nuclear;   ///< index (z1 - 1) * 92 + (z2 - 1), symmetric in z1, z2
  std::string source;                           ///< where the files were read from
};

constexpr int kSrimMaxZ = 92;

/// Parse the text of SCOEF03.DAT and SNUC03.dat.  Throws std::invalid_argument naming \p source and the line for a
/// malformed row, a missing element or a missing pair.
SrimData ParseSrimData(const std::string &scoef_text, const std::string &snuc_text, const std::string &source);
/// Read SCOEF03.DAT and SNUC03.dat from \p directory (file names matched without regard to case, as LunaSRIM does).
/// Throws std::runtime_error if a file is missing.
SrimData ReadSrimData(const std::string &directory);

/// The environment variable that points SrimTables() to a directory with the two files.
constexpr const char *kSrimDataDirVariable = "AZURE2_SRIM_DATA_DIR";
/// The data every function without an explicit SrimData uses (see above).  Throws std::runtime_error if neither the
/// variable nor a compiled-in copy gives the files.
const SrimData &SrimTables();
/// The tables of \p directory, read now and owned by the caller, or the compiled-in copy (the process's, parsed once)
/// when \p directory is empty.  Throws as ReadSrimData, or as SrimTables() when the copy is missing.
std::shared_ptr<const SrimData> SrimTablesIn(const std::string &directory);

/// The row of element \p z of \p data; throws std::out_of_range outside 1..92.
const SrimElement &ElementData(const SrimData &data, int z);
/// The same in SrimTables().
const SrimElement &ElementData(int z);
/// The four SNUC03.dat coefficients of the pair (z1, z2); throws std::out_of_range outside 1..92.
const std::array<double, 4> &NuclearCoefficients(const SrimData &data, int z1, int z2);

}  // namespace stopping
}  // namespace azure2

#endif
