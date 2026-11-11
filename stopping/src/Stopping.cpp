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

double CompoundCorrectionFactor(const Material &, const Projectile &, double) { return 1.0; }

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
  double electronic_unc = 0.0;
  for (const auto &c : material.constituents) {
    const StoppingWithUncertainty e = ElectronicStopping(projectile, c.z, energy_MeV);
    const double n = NuclearStopping(projectile, c.z, Element(c.z).atomic_weight_u, energy_MeV);
    v.electronic += c.atom_fraction * e.value;
    v.nuclear += c.atom_fraction * n;
    electronic_unc += c.atom_fraction * e.value * e.rel_uncertainty;  // fully correlated
  }
  const double factor = CompoundCorrectionFactor(material, projectile, energy_MeV);
  v.electronic *= factor;
  v.nuclear *= factor;
  electronic_unc *= factor;
  v.total = v.electronic + v.nuclear;
  v.rel_uncertainty = std::hypot(electronic_unc, kNuclearRelUncertainty * v.nuclear) / v.total;
  return v;
}

std::string StoppingProvenance(const Material &material, const Projectile &projectile) {
  std::ostringstream o;
  o << "stopping of " << projectile.label << " (Z1 = " << projectile.z << ", " << projectile.mass_u << " u) in "
    << material.name << " (" << PhaseName(material.phase) << "), eV/(1e15 atoms/cm^2) per atom of the material\n";
  if (const UserStopping *user = FindUserStopping(material, projectile)) {
    o << "  user " << (user->energy_MeV.empty() ? "formula" : "table") << ": " << user->provenance
      << ", relative uncertainty " << user->rel_uncertainty << "\n";
    return o.str();
  }
  o << "  Bragg additivity; compound correction: none (factor 1, docs/upgrade/stopping.md section 6)\n";
  for (const auto &c : material.constituents) {
    const ElementRecord &el = Element(c.z);
    o << "  " << el.symbol << " (atom fraction " << c.atom_fraction << "): I = " << el.mean_excitation_eV
      << " eV (ICRU 37/49), A_r = " << el.atomic_weight_u << " u (" << el.atomic_weight_source << ")\n";
    const int table_z = projectile.z <= 2 ? projectile.z : 1;
    const CoefficientRow &row = Coefficients(c.z, table_z);
    o << "    electronic: " << row.source << ", row " << row.origin << " (used alone to " << row.fit_max_MeV_u
      << " MeV/u, uncertainty " << row.rel_uncertainty << ")";
    if (row.note != "-") o << "; note: " << row.note;
    o << "\n    above " << kBetheFrom_MeV_u << " MeV/u: Bethe with Barkas-Berger shell, Ashley-Ritchie-Brandt Barkas and "
      << "Bloch corrections (uncertainty " << kBetheRelUncertainty << "), log-log interpolation in between\n";
    if (projectile.z > 2) {
      o << "    heavy ion: proton value at the same velocity x (Pierce-Blann effective charge)^2, uncertainty "
        << kHeavyIonRelUncertainty << " (" << kSlowHeavyIonRelUncertainty << " below v = v0 Z1^(2/3))\n";
    }
    if (el.phase != material.phase) {
      o << "    phase: the data are for the element as " << PhaseName(el.phase) << ", the material is "
        << PhaseName(material.phase) << " (no phase correction in this branch)\n";
    }
  }
  o << "  nuclear: ZBL universal (Ziegler, Biersack, Littmark 1985), uncertainty " << kNuclearRelUncertainty << "\n";
  return o.str();
}

}  // namespace stopping
}  // namespace azure2
