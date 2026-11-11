#include "azure2/stopping/Tables.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "TableText.h"

namespace azure2 {
namespace stopping {

namespace {

[[noreturn]] void Fail(const std::string &where, const std::string &what) {
  throw std::invalid_argument(where + ": " + what);
}

double Number(const std::string &text, const std::string &where) {
  char *end = nullptr;
  const double v = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0') Fail(where, "'" + text + "' is not a number");
  return v;
}

std::string Trim(const std::string &s) {
  const auto a = s.find_first_not_of(' ');
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(' ') - a + 1);
}

/// A line of SCOEF03.DAT: two quoted names, then numbers separated by blanks.
std::vector<std::string> ScoefTokens(const std::string &line, const std::string &where) {
  std::vector<std::string> tokens;
  std::size_t i = 0;
  while (i < line.size()) {
    if (std::isspace(static_cast<unsigned char>(line[i]))) {
      ++i;
    } else if (line[i] == '"') {
      const auto close = line.find('"', i + 1);
      if (close == std::string::npos) Fail(where, "unterminated quoted name");
      tokens.push_back(line.substr(i + 1, close - i - 1));
      i = close + 1;
    } else {
      const auto start = i;
      while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
      tokens.push_back(line.substr(start, i - start));
    }
  }
  return tokens;
}

template <std::size_t N>
void Columns(std::array<double, N> &out, const std::vector<double> &col, int first) {
  for (std::size_t k = 0; k < N; ++k) out[k] = col[static_cast<std::size_t>(first) + k];
}

/// The text of a file, carriage returns kept (the parsers drop them).
std::string FileText(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("stopping: cannot read " + path.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

/// \p name in \p directory, or a file whose name differs only in case.
std::filesystem::path FindFile(const std::string &directory, const std::string &name) {
  const std::filesystem::path exact = std::filesystem::path(directory) / name;
  if (std::filesystem::exists(exact)) return exact;
  std::string lower = name;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
  std::error_code error;
  for (const auto &entry : std::filesystem::directory_iterator(directory, error)) {
    std::string candidate = entry.path().filename().string();
    std::transform(candidate.begin(), candidate.end(), candidate.begin(), [](unsigned char c) { return std::tolower(c); });
    if (candidate == lower) return entry.path();
  }
  throw std::runtime_error("stopping: no " + name + " in '" + directory + "' (SRIM-2013's data file; see " +
                           "stopping/data/srim/README.md)");
}

}  // namespace

std::string PhaseName(Phase phase) { return phase == Phase::kGas ? "gas" : "solid"; }

Phase PhaseFromName(const std::string &name) {
  if (name == "solid") return Phase::kSolid;
  if (name == "gas") return Phase::kGas;
  throw std::invalid_argument("phase '" + name + "': expected solid or gas");
}

SrimData ParseSrimData(const std::string &scoef_text, const std::string &snuc_text, const std::string &source) {
  SrimData data;
  data.source = source;
  data.elements.resize(kSrimMaxZ);
  std::vector<bool> seen(kSrimMaxZ, false);
  {
    std::istringstream in(scoef_text);
    std::string line;
    int number = 0;
    while (std::getline(in, line)) {
      ++number;
      if (!line.empty() && line.back() == '\r') line.pop_back();
      const std::string where = source + "/SCOEF03.DAT:" + std::to_string(number);
      const auto tokens = ScoefTokens(line, where);
      if (tokens.empty()) continue;
      if (tokens.size() != 95) Fail(where, std::to_string(tokens.size()) + " fields, 95 expected (2 names, 93 numbers)");
      std::vector<double> col(94);  // col[1..93], as lunasrim/data_files.py numbers them
      for (int k = 1; k <= 93; ++k) col[k] = Number(tokens[k + 1], where);
      const int z = static_cast<int>(col[1]);
      if (z != col[1] || z < 1 || z > kSrimMaxZ) Fail(where, "Z outside 1..92");
      if (seen[z - 1]) Fail(where, "element Z = " + std::to_string(z) + " given twice");
      seen[z - 1] = true;
      SrimElement &e = data.elements[z - 1];
      e.z = z;
      e.symbol = Trim(tokens[0]);
      e.name = Trim(tokens[1]);
      e.most_abundant_mass_u = col[3];
      e.atomic_weight_u = col[4];
      e.solid_density_g_cm3 = col[5];
      e.gas_density_g_cm3 = col[9];
      Columns(e.proton, col, 11);
      Columns(e.helium, col, 21);
      Columns(e.lithium, col, 31);
      Columns(e.ion_poly, col, 41);
      Columns(e.target_poly, col, 54);
      e.mean_excitation_solid_eV = col[67];
      Columns(e.mean_excitation_gas_eV, col, 68);
      Columns(e.shell_c, col, 70);
      Columns(e.shell_d, col, 81);
      Columns(e.high_energy_charge, col, 89);
      e.origin = where;
      if (!(e.atomic_weight_u > 0.0) || !(e.mean_excitation_solid_eV > 0.0)) Fail(where, "mass and I must be > 0");
    }
  }
  for (int z = 1; z <= kSrimMaxZ; ++z) {
    if (!seen[z - 1]) Fail(source + "/SCOEF03.DAT", "no row for Z = " + std::to_string(z));
  }
  data.nuclear.assign(kSrimMaxZ * kSrimMaxZ, std::array<double, 4>{});
  std::vector<bool> pair(kSrimMaxZ * kSrimMaxZ, false);
  {
    std::istringstream in(snuc_text);
    std::string line;
    int number = 0;
    while (std::getline(in, line)) {
      ++number;
      std::istringstream fields(line);
      std::vector<std::string> f;
      std::string token;
      while (fields >> token) f.push_back(token);
      if (f.empty()) continue;
      const std::string where = source + "/SNUC03.dat:" + std::to_string(number);
      if (f.size() != 6) Fail(where, std::to_string(f.size()) + " fields, 6 expected (Z1 Z2 and 4 coefficients)");
      const double z1d = Number(f[0], where), z2d = Number(f[1], where);
      const int z1 = static_cast<int>(z1d), z2 = static_cast<int>(z2d);
      if (z1 != z1d || z2 != z2d || z1 < 1 || z2 < 1 || z1 > kSrimMaxZ || z2 > kSrimMaxZ) Fail(where, "Z outside 1..92");
      std::array<double, 4> c{};
      for (int k = 0; k < 4; ++k) c[k] = Number(f[2 + k], where);
      for (const auto &[a, b] : {std::pair<int, int>{z1, z2}, std::pair<int, int>{z2, z1}}) {
        data.nuclear[(a - 1) * kSrimMaxZ + (b - 1)] = c;
        pair[(a - 1) * kSrimMaxZ + (b - 1)] = true;
      }
    }
  }
  for (int k = 0; k < kSrimMaxZ * kSrimMaxZ; ++k) {
    if (!pair[k]) {
      Fail(source + "/SNUC03.dat", "no row for the pair (" + std::to_string(k / kSrimMaxZ + 1) + ", " +
                                       std::to_string(k % kSrimMaxZ + 1) + ")");
    }
  }
  return data;
}

SrimData ReadSrimData(const std::string &directory) {
  return ParseSrimData(FileText(FindFile(directory, "SCOEF03.DAT")), FileText(FindFile(directory, "SNUC03.dat")),
                       directory);
}

const SrimData &SrimTables() {
  static const SrimData data = [] {
    const char *directory = std::getenv(kSrimDataDirVariable);
    if (directory != nullptr && directory[0] != '\0') return ReadSrimData(directory);
    const std::string scoef = EmbeddedScoefText();
    if (scoef.empty()) {
      throw std::runtime_error(std::string("stopping: SRIM-2013's data files are not compiled into this build; set ") +
                               kSrimDataDirVariable + " to a directory holding SCOEF03.DAT and SNUC03.dat (e.g. " +
                               "the data folder of LunaSRIM; see stopping/data/srim/README.md)");
    }
    return ParseSrimData(scoef, EmbeddedSnucText(), "stopping/data/srim (compiled in)");
  }();
  return data;
}

const SrimElement &ElementData(const SrimData &data, int z) {
  if (z < 1 || z > static_cast<int>(data.elements.size())) {
    throw std::out_of_range("stopping: element Z = " + std::to_string(z) + " outside 1..92");
  }
  return data.elements[z - 1];
}

const SrimElement &ElementData(int z) { return ElementData(SrimTables(), z); }

const std::array<double, 4> &NuclearCoefficients(const SrimData &data, int z1, int z2) {
  if (z1 < 1 || z2 < 1 || z1 > kSrimMaxZ || z2 > kSrimMaxZ) {
    throw std::out_of_range("stopping: nuclear stopping pair (" + std::to_string(z1) + ", " + std::to_string(z2) +
                            ") outside 1..92");
  }
  return data.nuclear[(z1 - 1) * kSrimMaxZ + (z2 - 1)];
}

}  // namespace stopping
}  // namespace azure2
