#include "azure2/stopping/Stopping.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Nuclear.h"

namespace azure2 {
namespace stopping {

namespace {

/// A user table at energy_MeV: linear in log E - log epsilon; refused outside the table.
double UserTableValue(const UserStopping &user, double energy_MeV) {
  const auto &e = user.energy_MeV;
  if (energy_MeV < e.front() || energy_MeV > e.back()) {
    std::ostringstream o;
    o << "stopping: " << energy_MeV << " MeV is outside the user table (" << e.front() << " - " << e.back()
      << " MeV) of " << user.provenance << "; the table is not extrapolated";
    throw std::out_of_range(o.str());
  }
  const auto hi = std::upper_bound(e.begin(), e.end(), energy_MeV);
  const std::size_t k = std::min<std::size_t>(static_cast<std::size_t>(hi - e.begin()), e.size() - 1);
  const double x = std::log(energy_MeV / e[k - 1]) / std::log(e[k] / e[k - 1]);
  return std::exp(std::log(user.stopping[k - 1]) + x * std::log(user.stopping[k] / user.stopping[k - 1]));
}

}  // namespace

const UserStopping *FindUserStopping(const Material &material, const Projectile &projectile) {
  for (const auto &u : material.user_stopping) {
    if (u.projectile_z == projectile.z && std::fabs(u.projectile_mass_u - projectile.mass_u) < 0.01) return &u;
  }
  return nullptr;
}

StoppingValue Stopping(const Material &material, const Projectile &projectile, double energy_MeV) {
  if (!(energy_MeV > 0.0)) throw std::domain_error("stopping: energy must be > 0");
  StoppingValue v;
  if (const UserStopping *user = FindUserStopping(material, projectile)) {
    const bool table = !user->energy_MeV.empty();
    if (table == static_cast<bool>(user->formula)) {
      throw std::invalid_argument("stopping: user stopping '" + user->provenance + "' must give a table or a formula, not both or neither");
    }
    v.total = table ? UserTableValue(*user, energy_MeV) : user->formula(energy_MeV);
    if (!(v.total > 0.0)) throw std::domain_error("stopping: user stopping '" + user->provenance + "' is not positive");
    v.electronic = v.total;
    v.rel_uncertainty = user->rel_uncertainty;
    v.from_user = true;
    return v;
  }
  if (material.constituents.empty()) throw std::invalid_argument("stopping: material '" + material.name + "' has no constituents");
  if (!(material.compound_correction > 0.0)) {
    throw std::invalid_argument("stopping: material '" + material.name + "': the compound correction must be > 0");
  }
  const SrimData &data = SrimTables();
  double electronic_unc = 0.0;
  for (const auto &c : material.constituents) {
    const StoppingWithUncertainty e =
        ElectronicStopping(data, projectile, c.z, energy_MeV, material.phase, material.compound_correction);
    const double n = NuclearStopping(data, projectile, c.z, ElementData(data, c.z).atomic_weight_u, energy_MeV);
    v.electronic += c.atom_fraction * e.value;
    v.nuclear += c.atom_fraction * n;
    electronic_unc += c.atom_fraction * e.value * e.rel_uncertainty;  // fully correlated
  }
  v.total = v.electronic + v.nuclear;
  v.rel_uncertainty = std::hypot(electronic_unc, kNuclearRelUncertainty * v.nuclear) / v.total;
  return v;
}

std::vector<double> StoppingBreakpoints(const Material &material, const Projectile &projectile) {
  std::vector<double> breaks;
  if (FindUserStopping(material, projectile) != nullptr) return breaks;
  const double per_u = projectile.mass_u * 1.0e-3;  // keV/u -> MeV
  breaks.push_back(kSrimVelocityLimit_keV_u * per_u);
  if (projectile.z <= 3) {
    const double fit_max = projectile.z == 1 ? kProtonFitMax_keV_u : kHeliumLithiumFitMax_keV_u;
    breaks.push_back(fit_max * per_u);
    breaks.push_back(fit_max * kBetheJoinFactor * per_u);
  } else {
    // helium's own branches at the same e, inside the heavy-ion scaling, and the end of the scaling
    breaks.push_back(kHeliumLithiumFitMax_keV_u * per_u);
    breaks.push_back(kHeliumLithiumFitMax_keV_u * kBetheJoinFactor * per_u);
    breaks.push_back(kHeavyIonScalingMax_keV_u * per_u);
  }
  for (const auto &c : material.constituents) {
    breaks.push_back(NuclearFormChange_MeV(projectile, c.z, ElementData(c.z).atomic_weight_u));
  }
  std::sort(breaks.begin(), breaks.end());
  breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
  return breaks;
}

std::string StoppingProvenance(const Material &material, const Projectile &projectile) {
  std::ostringstream o;
  o << "stopping of " << projectile.label << " (Z1 = " << projectile.z << ", " << projectile.mass_u << " u) in "
    << material.name << " (phase " << PhaseName(material.phase) << ", compound correction "
    << material.compound_correction << "), eV/(1e15 atoms/cm^2) per atom of the material\n";
  if (const UserStopping *user = FindUserStopping(material, projectile)) {
    o << "  user " << (user->energy_MeV.empty() ? "formula" : "table") << ": " << user->provenance
      << ", relative uncertainty " << user->rel_uncertainty << "\n";
    return o.str();
  }
  const SrimData &data = SrimTables();
  o << "  " << kSrimImplementation << "\n  cite SRIM: " << kSrimCitation << "\n  data: " << data.source
    << " (SCOEF03.DAT, SNUC03.dat of SRIM-2013)\n";
  o << "  Bragg's rule; ";
  if (material.compound_correction == 1.0) {
    o << "no compound correction (C = 1)\n";
  } else {
    o << "compound correction C = " << material.compound_correction
      << " on the electronic stopping, faded out with velocity as SRIM does (applied once)\n";
  }
  o << "  electronic: " << (projectile.z == 1 ? "SRIM's H fit to 1000 keV/u, Bethe above 1200 keV/u"
                            : projectile.z <= 3 ? "SRIM's He/Li fit to 2000 keV/u, Bethe above 2400 keV/u"
                                                : "SRIM's He stopping scaled to the ion (2 keV/u - 30 MeV/u), "
                                                  "high-energy effective charge above")
    << ", velocity-proportional below 2 keV/u; uncertainty " << ElectronicRelUncertainty(projectile.z)
    << " (SRIM's published accuracy for this ion, www.srim.org/SRIM/SRIM2011.htm)\n";
  for (const auto &c : material.constituents) {
    const SrimElement &el = ElementData(data, c.z);
    o << "  " << el.symbol << " (atom fraction " << c.atom_fraction << "): row " << el.origin << ", I = "
      << (material.phase == Phase::kGas ? el.mean_excitation_gas_eV[projectile.z == 1 ? 0 : 1] : el.mean_excitation_solid_eV)
      << " eV (" << PhaseName(material.phase) << "), M = " << el.atomic_weight_u << " u, phase factor "
      << SrimPhaseFactor(c.z, material.phase) << (IsNaturalGas(c.z) ? " (natural gas)" : " (natural solid)")
      << "; nuclear: SNUC03.dat pair (" << projectile.z << ", " << c.z << "), reduced energy 30 at "
      << NuclearFormChange_MeV(projectile, c.z, el.atomic_weight_u) << " MeV\n";
  }
  o << "  nuclear uncertainty " << kNuclearRelUncertainty << " (provisional, stopping.md section 7)\n";
  return o.str();
}

}  // namespace stopping
}  // namespace azure2
