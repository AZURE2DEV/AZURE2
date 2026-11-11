// SRIM-2013's data files (docs/upgrade/stopping.md section 3): parsed completely, the columns where LunaSRIM's reader
// puts them, read from a directory exactly as from the compiled copy, and the parser refuses malformed text.  Run a
// second time (ctest stopping_srim_directory) with AZURE2_SRIM_DATA_DIR set: then SrimTables() must be the directory's.
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "azure2/physics/Elements.h"
#include "azure2/stopping/Tables.h"
#include "check.h"

using namespace azure2::stopping;

namespace {
template <std::size_t N>
bool Same(const std::array<double, N> &a, const std::array<double, N> &b) {
  for (std::size_t k = 0; k < N; ++k) {
    if (a[k] != b[k]) return false;
  }
  return true;
}
}  // namespace

int main(int argc, char **argv) {
  const std::string dir = argc > 2 ? argv[2] : "../data/srim";
  const char *env = std::getenv(kSrimDataDirVariable);
  const SrimData &d = SrimTables();
  if (env != nullptr && env[0] != '\0') {
    check::True(std::string("with ") + kSrimDataDirVariable + " set, SrimTables() reads that directory", d.source == env);
  } else {
    check::True("without the variable, SrimTables() is the compiled copy", d.source.find("compiled in") != std::string::npos);
  }
  check::True("92 elements, 92 x 92 nuclear pairs", d.elements.size() == 92 && d.nuclear.size() == 92 * 92);
  bool symbols = true, positive = true, symmetric = true;
  for (int z = 1; z <= 92; ++z) {
    const SrimElement &e = ElementData(d, z);
    symbols = symbols && e.z == z && e.symbol == azure2::physics::ElementSymbol(z);
    positive = positive && e.atomic_weight_u > 0 && e.mean_excitation_solid_eV > 0 && e.most_abundant_mass_u > 0;
    for (int z2 = 1; z2 <= 92; ++z2) symmetric = symmetric && NuclearCoefficients(d, z, z2) == NuclearCoefficients(d, z2, z);
  }
  check::True("SCOEF03.DAT: Z = 1..92 in order, symbols are physics/'s ElementSymbol(Z)", symbols);
  check::True("SCOEF03.DAT: masses and I > 0", positive);
  check::True("SNUC03.dat: coefficients symmetric in (Z1, Z2)", symmetric);
  // Columns, values read off SCOEF03.DAT's first line (H) and SNUC03.dat's first line (1 1) by eye.
  const SrimElement &h = ElementData(d, 1);
  check::Rel("H: column 3, most abundant isotope 1.0078 u", h.most_abundant_mass_u, 1.0078, 0.0);
  check::Rel("H: column 4, atomic weight 1.0079 u", h.atomic_weight_u, 1.0079, 0.0);
  check::Rel("H: column 11, first proton coefficient", h.proton[0], 0.0137722064, 0.0);
  check::Rel("H: column 18, eighth proton coefficient", h.proton[7], 0.251198273, 0.0);
  check::Rel("H: column 21, first helium coefficient", h.helium[0], 1.62392840E-2, 0.0);
  check::Rel("H: column 67, I solid 19.4 eV", h.mean_excitation_solid_eV, 19.4, 0.0);
  check::Rel("H: column 69, I gas (Z1 > 1) 19.7 eV", h.mean_excitation_gas_eV[1], 19.7, 0.0);
  check::Rel("H: column 70, first shell_c coefficient", h.shell_c[0], 2.70659733E-2, 0.0);
  check::Rel("H: column 88, last shell_d coefficient", h.shell_d[7], 1.11332707E-9, 0.0);
  check::Rel("SNUC03.dat (1, 1): first coefficient", NuclearCoefficients(d, 1, 1)[0], 3.31122449, 0.0);
  check::Rel("SNUC03.dat (92, 92): last coefficient", NuclearCoefficients(d, 92, 92)[3], 0.243177191, 0.0);
  check::True("row origins name the file and line", h.origin.find("SCOEF03.DAT:1") != std::string::npos);

  // The directory reader gives the same numbers as the compiled copy (or as the variable's directory).
  const SrimData r = ReadSrimData(dir);
  bool same = r.nuclear == d.nuclear;
  for (int z = 1; z <= 92; ++z) {
    const SrimElement &a = r.elements[z - 1], &b = d.elements[z - 1];
    same = same && a.symbol == b.symbol && a.name == b.name && a.atomic_weight_u == b.atomic_weight_u &&
           a.most_abundant_mass_u == b.most_abundant_mass_u && Same(a.proton, b.proton) && Same(a.helium, b.helium) &&
           Same(a.lithium, b.lithium) && Same(a.ion_poly, b.ion_poly) && Same(a.target_poly, b.target_poly) &&
           a.mean_excitation_solid_eV == b.mean_excitation_solid_eV &&
           Same(a.mean_excitation_gas_eV, b.mean_excitation_gas_eV) && Same(a.shell_c, b.shell_c) &&
           Same(a.shell_d, b.shell_d) && Same(a.high_energy_charge, b.high_energy_charge);
  }
  check::True("ReadSrimData(" + dir + ") equals SrimTables() in every number", same);
  check::Throws<std::runtime_error>("a directory without the files is refused", [] { ReadSrimData("/nonexistent-srim-dir"); });

  check::Throws<std::out_of_range>("ElementData(0) refused", [] { ElementData(0); });
  check::Throws<std::out_of_range>("NuclearCoefficients(93, 1) refused", [&] { NuclearCoefficients(d, 93, 1); });
  const std::string h_line = "\" H\" \"Hydrogen\" 1 1 1.0078 1.0079" + std::string(" 0.5") + "\n";
  check::Throws<std::invalid_argument>("a short SCOEF03.DAT row is refused", [&] { ParseSrimData(h_line, "", "t"); });
  check::Throws<std::invalid_argument>("a SNUC03.dat row with five fields is refused", [] {
    ParseSrimData("", " 1  1 3.3 0.005 0.45\n", "t");
  });
  check::Throws<std::invalid_argument>("an incomplete SCOEF03.DAT is refused (no row for Z = 1)", [] { ParseSrimData("", "", "t"); });
  check::Throws<std::invalid_argument>("phase 'liquid' is refused (SRIM knows solid and gas)", [] { PhaseFromName("liquid"); });
  return check::Done();
}
