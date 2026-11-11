#include "azure2/stopping/Nuclear.h"

#include <cmath>
#include <stdexcept>

namespace azure2 {
namespace stopping {

namespace {
double Screening(int z1, int z2) { return std::pow(static_cast<double>(z1), 0.23) + std::pow(static_cast<double>(z2), 0.23); }
}  // namespace

double ZblReducedEnergy(int z1, double m1_u, int z2, double m2_u, double energy_keV) {
  return energy_keV * 32.53 * m2_u / (static_cast<double>(z1) * z2 * (m1_u + m2_u) * Screening(z1, z2));
}

double SrimReducedNuclearStopping(const std::array<double, 4> &c, double eps) {
  if (!(eps > 0.0)) throw std::domain_error("stopping: reduced energy must be > 0");
  if (eps > kNuclearFitMaxReducedEnergy) return std::log(eps) / (2.0 * eps);
  // log1p where LunaSRIM writes Log(1 + c1 eps): the same to 1e-16 / (c1 eps) relative (below 1e-12 for eps >= 1e-4),
  // and free of the rounding noise of 1 + x at eps < 1e-12, which the range integral down to E = 0 cannot converge on
  return 0.5 * std::log1p(c[0] * eps) / (eps + c[1] * std::pow(eps, c[2]) + c[3] * std::pow(eps, 0.5));
}

double NuclearStopping(const SrimData &data, const Projectile &projectile, int z2, double m2_u, double energy_MeV) {
  if (!(energy_MeV > 0.0)) throw std::domain_error("stopping: energy must be > 0");
  const int z1 = projectile.z;
  const double m1 = projectile.mass_u;
  const double eps = ZblReducedEnergy(z1, m1, z2, m2_u, energy_MeV * 1.0e3);
  const double sn = SrimReducedNuclearStopping(NuclearCoefficients(data, z1, z2), eps);
  return sn * 8.462 * z1 * z2 * m1 / ((m1 + m2_u) * Screening(z1, z2));
}

double NuclearStopping(const Projectile &projectile, int z2, double m2_u, double energy_MeV) {
  return NuclearStopping(SrimTables(), projectile, z2, m2_u, energy_MeV);
}

double NuclearFormChange_MeV(const Projectile &projectile, int z2, double m2_u) {
  return kNuclearFitMaxReducedEnergy / ZblReducedEnergy(projectile.z, projectile.mass_u, z2, m2_u, 1.0) * 1.0e-3;
}

}  // namespace stopping
}  // namespace azure2
