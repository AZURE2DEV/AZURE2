// The unit table of format 3 (Units.h; docs/upgrade/project_format3.md section 7.2).
#include "azure2/experiment/Units.h"

#include <cctype>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace experiment {

const std::vector<Dimension> &Dimensions() {
  // {dimension, its base unit (the engine's, factor 1), what it measures}
  static const std::vector<Dimension> table = {
      {"energy", "MeV", "energies, widths, U_e, energy scales"},
      {"sqrt_energy", "MeV^1/2", "the coefficient b of an energy shift b sqrt(E), and of straggling c sqrt(dE)"},
      {"angle", "deg", "detection angles"},
      {"cross_section", "b", "angle-integrated cross sections"},
      {"diff_cross_section", "b/sr", "differential cross sections"},
      {"sfactor", "MeV*b", "astrophysical S factors"},
      {"dimensionless", "1", "ratios, normalizations, abundances, efficiencies"},
      {"areal_density", "atoms/cm2", "areal densities of a target's active nucleus"},
      {"stopping", "MeV*cm2/atom", "stopping cross sections per active atom"},
      {"temperature", "K", "target temperatures (free-gas Doppler broadening)"},
      {"length", "fm", "scattering lengths, flight paths"},
      {"time", "ns", "time-of-flight resolution"},
      {"anc", "fm^-1/2", "asymptotic normalization coefficients"},
      {"yield", "1/particle", "reactions per incident beam particle"},
      {"yield_per_sr", "1/particle/sr", "reactions per incident particle per steradian"},
      {"per_charge", "1/C", "reactions per unit of beam charge (needs the charge state q)"},
      {"counts_per_charge", "counts/C", "detected counts per unit of beam charge (needs q and the efficiency)"},
      {"arbitrary", "arbitrary", "arbitrary units: the normalization is free (no prior) or profiled"},
  };
  return table;
}

const std::vector<Unit> &Units() {
  // {token, dimension, factor to the dimension's base unit (Dimensions() above), aliases if any}
  static const std::vector<Unit> table = {
      {"MeV", "energy", 1.0},
      {"keV", "energy", 1e-3},
      {"eV", "energy", 1e-6},
      {"meV", "energy", 1e-9},
      {"ueV", "energy", 1e-12, {"\xC2\xB5" "eV"}},
      {"MeV^1/2", "sqrt_energy", 1.0},
      {"keV^1/2", "sqrt_energy", std::sqrt(1e-3)},
      {"deg", "angle", 1.0},
      {"rad", "angle", 180.0 / 3.14159265358979323846},
      {"b", "cross_section", 1.0},
      {"mb", "cross_section", 1e-3},
      {"ub", "cross_section", 1e-6, {"\xC2\xB5" "b"}},
      {"nb", "cross_section", 1e-9},
      {"b/sr", "diff_cross_section", 1.0},
      {"mb/sr", "diff_cross_section", 1e-3},
      {"ub/sr", "diff_cross_section", 1e-6, {"\xC2\xB5" "b/sr"}},
      {"nb/sr", "diff_cross_section", 1e-9},
      {"MeV*b", "sfactor", 1.0},
      {"keV*b", "sfactor", 1e-3},
      {"eV*b", "sfactor", 1e-6},
      {"1", "dimensionless", 1.0},
      {"%", "dimensionless", 1e-2},
      {"atoms/cm2", "areal_density", 1.0},
      {"atoms/b", "areal_density", 1e24},
      {"MeV*cm2/atom", "stopping", 1.0},
      {"keV*cm2/atom", "stopping", 1e-3},
      {"eV*cm2/atom", "stopping", 1e-6},
      {"K", "temperature", 1.0},
      {"fm", "length", 1.0},
      {"m", "length", 1e15},
      {"ns", "time", 1.0},
      {"ps", "time", 1e-3},
      {"us", "time", 1e3, {"\xC2\xB5" "s"}},
      {"fm^-1/2", "anc", 1.0},
      {"1/particle", "yield", 1.0},
      {"1/particle/sr", "yield_per_sr", 1.0},
      {"1/uC", "per_charge", 1e6, {"1/\xC2\xB5" "C"}},
      {"1/mC", "per_charge", 1e3},
      {"1/C", "per_charge", 1.0},
      {"counts/uC", "counts_per_charge", 1e6, {"counts/\xC2\xB5" "C"}},
      {"counts/mC", "counts_per_charge", 1e3},
      {"counts/C", "counts_per_charge", 1.0},
      {"arbitrary", "arbitrary", 1.0},
  };
  return table;
}

const Unit *FindUnit(const std::string &token) {
  for (const Unit &u : Units()) {
    if (token == u.token) return &u;
    for (const std::string &a : u.aliases)
      if (token == a) return &u;
  }
  return nullptr;
}

const Unit *FindUnitIgnoringCase(const std::string &token) {
  auto lower = [](std::string s) {
    for (char &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
  };
  const std::string t = lower(token);
  const Unit *found = nullptr;
  for (const Unit &u : Units()) {
    if (lower(u.token) != t) continue;
    if (found) return nullptr;  // meV and MeV: no fix when the case decides
    found = &u;
  }
  return found;
}

std::vector<std::string> UnitsOfDimension(const std::string &dimension) {
  std::vector<std::string> out;
  for (const Unit &u : Units())
    if (dimension == u.dimension) out.push_back(u.token);
  return out;
}

const Dimension *FindDimension(const std::string &name) {
  for (const Dimension &d : Dimensions())
    if (name == d.name) return &d;
  return nullptr;
}

double Rescale(double value, double factor) {
  if (factor == 1.0) return value;
  if (factor >= 1.0) {
    const double n = std::round(factor);
    return std::fabs(factor - n) <= 1e-12 * n ? value * n : value * factor;
  }
  const double inverse = 1.0 / factor, n = std::round(inverse);
  return std::fabs(inverse - n) <= 1e-12 * n ? value / n : value * factor;
}

double FactorBetween(const Unit &from, const Unit &to) { return from.factor / to.factor; }

double ToBase(double value, const Unit *unit) { return unit ? Rescale(value, unit->factor) : value; }

double ToUnit(double value, const Unit *from, const std::string &to) {
  const Unit *target = FindUnit(to);
  if (!target || (from && std::string(from->dimension) != target->dimension))
    throw std::invalid_argument("experiment::ToUnit: " + std::string(from ? from->token : "(none)") + " -> " + to);
  return from ? Rescale(value, FactorBetween(*from, *target)) : value;
}

double ParticlesPerUnitCharge(const Unit &unit, int chargeState) {
  return 1.0 / (unit.factor * chargeState * physics::kElementaryCharge_C);
}

}  // namespace experiment
}  // namespace azure2
