// The astrophysical S factor (AstrophysicalFactor.h).
#include "azure2/physics/AstrophysicalFactor.h"

#include <cmath>
#include <limits>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace physics {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

double sfactor_conversion(const ChannelKinematics &kin, double E_MeV) {
  if (!(E_MeV > 0.0) || kin.z1 * kin.z2 <= 0) return std::numeric_limits<double>::quiet_NaN();
  ChannelKinematics nonRelativistic = kin;
  nonRelativistic.relativistic = false;
  return E_MeV * std::exp(2.0 * kPi * eta(nonRelativistic, E_MeV));
}

double modified_sfactor_conversion(const ChannelKinematics &kin, double E_MeV, double g_perMeV) {
  const double s = sfactor_conversion(kin, E_MeV);
  return g_perMeV == 0.0 ? s : s * std::exp(g_perMeV * E_MeV);
}

}  // namespace physics
}  // namespace azure2
