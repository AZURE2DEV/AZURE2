#include "azure2/stopping/Straggling.h"

#include <cmath>
#include <stdexcept>

#include "azure2/physics/Constants.h"
#include "azure2/stopping/Stopping.h"

namespace azure2 {
namespace stopping {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kEv2A2ToMeV2cm2 = 1.0e-28;  ///< eV^2 A^2 -> MeV^2 cm^2
}  // namespace

std::string StragglingModelName(StragglingModel model) { return model == StragglingModel::kBohr ? "bohr" : "trim"; }

StragglingModel StragglingModelFromName(const std::string &name) {
  if (name == "trim") return StragglingModel::kTrim;
  if (name == "bohr") return StragglingModel::kBohr;
  throw std::invalid_argument("straggling model '" + name + "': expected trim or bohr");
}

double BohrSpecificVariance(int z1, double z2) {
  const double e2_MeV_cm = physics::kFineStructure * physics::kHbarC_MeV_fm * 1.0e-13;
  return 4.0 * kPi * z1 * z1 * z2 * e2_MeV_cm * e2_MeV_cm;
}

double MeanTargetCharge(const Material &material) {
  double z = 0.0;
  for (const auto &c : material.constituents) z += c.atom_fraction * c.z;
  return z;
}

int TrimMeanTargetCharge(const Material &material) {
  // lunasrim/straggling.py energy_loss: zz = int(round(st / tot * z + zz)); std::nearbyint rounds half to even in the
  // default rounding mode, as Python's round does
  double z = 0.0;
  for (const auto &c : material.constituents) z = std::nearbyint(c.atom_fraction * c.z + z);
  return static_cast<int>(z);
}

double TrimSpecificVariance(const Material &material, const Projectile &projectile, double energy_MeV) {
  // lunasrim/straggling.py w_trim with the atom density N = 1 atom/A^3 (W / N does not depend on N)
  const int zz = TrimMeanTargetCharge(material);
  const int z1 = projectile.z;
  const double bohr = 177.0 * z1 * std::sqrt(static_cast<double>(zz));
  const double k_low = std::sqrt(3.2 * 12.0 / (projectile.mass_u * 1000) * std::pow(static_cast<double>(zz + z1), 8.0 / 3.0));
  const double e = energy_MeV * 1.0e6;  // eV
  double w = k_low * e;
  if (w > bohr) {
    const double q = 1 - 1 / std::pow(1 + e / (projectile.mass_u * 1000 * 931494), 2);
    w = bohr * ((1 / q - 0.5) / (1 / q - 1));
  }
  return w * w / 12.0 * kEv2A2ToMeV2cm2;
}

double SpecificStragglingVariance(StragglingModel model, const Material &material, const Projectile &projectile,
                                  double energy_MeV) {
  if (model == StragglingModel::kBohr) return BohrSpecificVariance(projectile.z, MeanTargetCharge(material));
  return TrimSpecificVariance(material, projectile, energy_MeV);
}

StragglingWidths PropagateStraggling(const StoppingFunction &epsilon, const std::function<double(double)> &specific_variance,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV,
                                     const std::vector<double> &breaks_MeV) {
  if (!(beam_sigma_MeV >= 0.0)) throw std::invalid_argument("stopping: beam spread must be >= 0");
  StragglingWidths w;
  w.energy_out_MeV = EnergyAfter(epsilon, energy_in_MeV, areal_density, breaks_MeV);
  // S in MeV cm^2 per atom
  const auto s = [&epsilon](double e) { return epsilon(e) * kStoppingToMeVcm2; };
  const double s_in = s(energy_in_MeV);
  const double s_out = s(w.energy_out_MeV);
  const double integral = IntegrateInSqrtE(
      [&](double e) {
        const double se = s(e);
        return specific_variance(e) / (se * se * se);
      },
      w.energy_out_MeV, energy_in_MeV, kQuadratureRelTol, breaks_MeV);
  w.straggling = s_out * std::sqrt(integral);
  w.beam = beam_sigma_MeV * s_out / s_in;
  w.total = std::hypot(w.straggling, w.beam);
  return w;
}

StragglingWidths PropagateStraggling(StragglingModel model, const Material &material, const Projectile &projectile,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV) {
  return PropagateStraggling(
      MaterialStopping(material, projectile),
      [&](double e) { return SpecificStragglingVariance(model, material, projectile, e); }, energy_in_MeV, areal_density,
      beam_sigma_MeV, StoppingBreakpoints(material, projectile));
}

}  // namespace stopping
}  // namespace azure2
