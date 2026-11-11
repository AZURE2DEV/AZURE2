#include "azure2/stopping/Nuclear.h"

#include <cmath>
#include <stdexcept>

namespace azure2 {
namespace stopping {

double ZblReducedEnergy(int z1, double m1_u, int z2, double m2_u, double energy_keV) {
  const double screening = std::pow(static_cast<double>(z1), 0.23) + std::pow(static_cast<double>(z2), 0.23);
  return 32.53 * m2_u * energy_keV / (z1 * z2 * (m1_u + m2_u) * screening);
}

double ZblReducedStopping(double reduced_energy) {
  const double e = reduced_energy;
  if (!(e > 0.0)) throw std::domain_error("stopping: ZBL reduced energy must be > 0");
  if (e > 30.0) return std::log(e) / (2.0 * e);
  return std::log1p(1.1383 * e) / (2.0 * (e + 0.01321 * std::pow(e, 0.21226) + 0.19593 * std::sqrt(e)));
}

double NuclearStopping(const Projectile &projectile, int z2, double m2_u, double energy_MeV) {
  if (!(energy_MeV > 0.0)) throw std::domain_error("stopping: energy must be > 0");
  const int z1 = projectile.z;
  const double m1 = projectile.mass_u;
  const double screening = std::pow(static_cast<double>(z1), 0.23) + std::pow(static_cast<double>(z2), 0.23);
  const double eps = ZblReducedEnergy(z1, m1, z2, m2_u, energy_MeV * 1.0e3);
  return 8.462 * z1 * z2 * m1 * ZblReducedStopping(eps) / ((m1 + m2_u) * screening);
}

}  // namespace stopping
}  // namespace azure2
