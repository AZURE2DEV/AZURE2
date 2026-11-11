#include "azure2/stopping/Material.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <stdexcept>

#include "azure2/physics/Constants.h"
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace stopping {

namespace {
constexpr double kAvogadro = 6.02214076e23;  ///< 1/mol, exact in the 2019 SI
}

Projectile Proton() { return {1, 1.00782503223, "p"}; }
Projectile Deuteron() { return {1, 2.01410177812, "d"}; }
Projectile Triton() { return {1, 3.01604928132, "t"}; }
Projectile Helium3() { return {2, 3.01602932197, "3He"}; }
Projectile Alpha() { return {2, 4.00260325413, "a"}; }

Projectile Ion(int z, double mass_u, const std::string &label) {
  if (z < 1 || !(mass_u > 0.0)) throw std::invalid_argument("stopping: ion '" + label + "' needs z >= 1 and a mass > 0");
  return {z, mass_u, label};
}

double Beta2(const Projectile &projectile, double energy_MeV) {
  const double gamma = 1.0 + energy_MeV / (projectile.mass_u * physics::kAtomicMassUnit_MeV);
  return 1.0 - 1.0 / (gamma * gamma);
}

std::vector<Constituent> ParseFormula(const std::string &formula) {
  std::map<int, double> counts;
  std::vector<int> order;
  std::size_t i = 0;
  const auto fail = [&](const std::string &why) {
    throw std::invalid_argument("stopping: formula '" + formula + "': " + why);
  };
  while (i < formula.size()) {
    if (!std::isupper(static_cast<unsigned char>(formula[i]))) fail("expected an element symbol at position " + std::to_string(i));
    std::string symbol(1, formula[i++]);
    while (i < formula.size() && std::islower(static_cast<unsigned char>(formula[i]))) symbol += formula[i++];
    const int z = physics::ElementZ(symbol);
    if (z < 1) fail("unknown element '" + symbol + "'");
    if (z > 92) fail("element '" + symbol + "' above Z = 92 has no stopping data");
    std::size_t start = i;
    while (i < formula.size() && (std::isdigit(static_cast<unsigned char>(formula[i])) || formula[i] == '.')) ++i;
    double count = 1.0;
    if (i > start) {
      const std::string text = formula.substr(start, i - start);
      char *end = nullptr;
      count = std::strtod(text.c_str(), &end);
      if (*end != '\0' || !(count > 0.0)) fail("bad count '" + text + "'");
    }
    if (!counts.count(z)) order.push_back(z);
    counts[z] += count;
  }
  if (counts.empty()) fail("empty");
  double total = 0.0;
  for (const auto &kv : counts) total += kv.second;
  std::vector<Constituent> out;
  for (int z : order) out.push_back({z, counts[z] / total});
  return out;
}

Material MaterialFromFractions(const std::string &name, const std::vector<Constituent> &constituents, Phase phase,
                               double mass_density_g_cm3) {
  if (constituents.empty()) throw std::invalid_argument("stopping: material '" + name + "' has no constituents");
  double total = 0.0;
  for (const auto &c : constituents) {
    if (c.z < 1 || c.z > 92) throw std::invalid_argument("stopping: material '" + name + "': Z = " + std::to_string(c.z) + " outside 1..92");
    if (!(c.atom_fraction > 0.0)) throw std::invalid_argument("stopping: material '" + name + "': fractions must be > 0");
    total += c.atom_fraction;
  }
  if (mass_density_g_cm3 < 0.0) throw std::invalid_argument("stopping: material '" + name + "': negative density");
  Material m;
  m.name = name;
  m.phase = phase;
  m.mass_density_g_cm3 = mass_density_g_cm3;
  for (const auto &c : constituents) m.constituents.push_back({c.z, c.atom_fraction / total});
  return m;
}

Material MaterialFromFormula(const std::string &formula, Phase phase, double mass_density_g_cm3) {
  return MaterialFromFractions(formula, ParseFormula(formula), phase, mass_density_g_cm3);
}

double MeanAtomicWeight(const Material &material) {
  double a = 0.0;
  for (const auto &c : material.constituents) a += c.atom_fraction * Element(c.z).atomic_weight_u;
  return a;
}

double AtomFraction(const Material &material, int z) {
  double f = 0.0;
  for (const auto &c : material.constituents) {
    if (c.z == z) f += c.atom_fraction;
  }
  return f;
}

UserStopping ReadStoppingTable(const std::string &text, const Material &material, const Projectile &projectile,
                               double rel_uncertainty, const std::string &provenance) {
  const auto fail = [&](const std::string &why) {
    throw std::invalid_argument("stopping table (" + provenance + "): " + why);
  };
  std::istringstream in(text);
  std::string line;
  double energy_factor = 0.0;  // to MeV
  double value_factor = 0.0;   // to eV/(1e15 atoms/cm^2)
  UserStopping user;
  user.projectile_z = projectile.z;
  user.projectile_mass_u = projectile.mass_u;
  user.rel_uncertainty = rel_uncertainty;
  user.provenance = provenance;
  while (std::getline(in, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream fields(line);
    std::string a, b, extra;
    if (!(fields >> a)) continue;
    if (!(fields >> b) || (fields >> extra)) fail("expected two columns in '" + line + "'");
    if (energy_factor == 0.0) {
      if (a == "E[MeV]") energy_factor = 1.0;
      else if (a == "E[keV]") energy_factor = 1.0e-3;
      else fail("first column unit '" + a + "': expected E[MeV] or E[keV]");
      if (b == "eps[eV/1e15atoms/cm2]") {
        value_factor = 1.0;
      } else if (b == "S[MeV*cm2/g]") {
        // MeV cm^2/g x (g per atom) = MeV cm^2 per atom; / 1e-21 to eV/(1e15 atoms/cm^2)
        value_factor = MeanAtomicWeight(material) / kAvogadro / kStoppingToMeVcm2;
      } else if (b == "S[keV/um]") {
        if (!(material.mass_density_g_cm3 > 0.0)) fail("keV/um needs the mass density of material '" + material.name + "'");
        const double atoms_per_cm3 = material.mass_density_g_cm3 * kAvogadro / MeanAtomicWeight(material);
        // keV/um = 1e-3 MeV / 1e-4 cm = 10 MeV/cm; / (atoms/cm^3) = MeV cm^2 per atom
        value_factor = 10.0 / atoms_per_cm3 / kStoppingToMeVcm2;
      } else {
        fail("second column unit '" + b + "': expected eps[eV/1e15atoms/cm2], S[MeV*cm2/g] or S[keV/um]");
      }
      continue;
    }
    char *end = nullptr;
    const double e = std::strtod(a.c_str(), &end);
    if (*end != '\0') fail("'" + a + "' is not a number");
    const double v = std::strtod(b.c_str(), &end);
    if (*end != '\0') fail("'" + b + "' is not a number");
    if (!(e > 0.0) || !(v > 0.0)) fail("energies and values must be > 0");
    if (!user.energy_MeV.empty() && !(e * energy_factor > user.energy_MeV.back())) fail("energies must ascend");
    user.energy_MeV.push_back(e * energy_factor);
    user.stopping.push_back(v * value_factor);
  }
  if (energy_factor == 0.0) fail("no unit line");
  if (user.energy_MeV.size() < 2) fail("fewer than two rows");
  return user;
}

}  // namespace stopping
}  // namespace azure2
