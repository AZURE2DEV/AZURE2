#include "azure2/stopping/Straggling.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/Constants.h"
#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Stopping.h"

namespace azure2 {
namespace stopping {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// The squared effective charge fraction gamma^2 Yang's formula multiplies H with: 1 for hydrogen ions; for helium
/// the Ziegler-Biersack-Littmark (1985) helium effective charge; for heavier ions the Pierce-Blann fraction squared.
double YangGamma2(int z1, int z2, double energy_MeV, double mass_u) {
  if (z1 == 1) return 1.0;
  if (z1 == 2) {
    static const double c[6] = {0.2865, 0.1266, -0.001429, 0.02402, -0.01135, 0.001475};
    const double e = std::max(0.0, std::log(energy_MeV / mass_u * 1.0e3));  // ln(E [keV/u])
    double x = c[0];
    double power = 1.0;
    for (int i = 1; i < 6; ++i) {
      power *= e;
      x += c[i] * power;
    }
    const double w0 = 7.6 - e;
    const double w = 1.0 + (0.007 + 0.00005 * z2) * std::exp(-w0 * w0);
    return (1.0 - std::exp(-x)) * w * w;
  }
  const double beta2 = 1.0 - 1.0 / std::pow(1.0 + energy_MeV / (mass_u * physics::kAtomicMassUnit_MeV), 2);
  const double q = EffectiveChargeFraction(z1, beta2);
  return q * q;
}

}  // namespace

double BohrSpecificVariance(int z1, double z2) {
  const double e2_MeV_cm = physics::kFineStructure * physics::kHbarC_MeV_fm * 1.0e-13;
  return 4.0 * kPi * z1 * z1 * z2 * e2_MeV_cm * e2_MeV_cm;
}

double MeanTargetCharge(const Material &material) {
  double z = 0.0;
  for (const auto &c : material.constituents) z += c.atom_fraction * c.z;
  return z;
}

double ChuFactor(const ChuTable &table, int z2, double energy_keV_u) {
  const auto it = std::find(table.z2.begin(), table.z2.end(), z2);
  if (it == table.z2.end()) throw std::out_of_range("stopping: Chu table (" + table.source + ") has no Z2 = " + std::to_string(z2));
  const auto &h = table.h[static_cast<std::size_t>(it - table.z2.begin())];
  const auto &e = table.energy_keV_u;
  if (e.empty() || h.size() != e.size()) throw std::invalid_argument("stopping: Chu table (" + table.source + ") is malformed");
  if (energy_keV_u >= e.back()) return 1.0;
  if (energy_keV_u <= e.front()) return h.front();
  const std::size_t k = static_cast<std::size_t>(std::upper_bound(e.begin(), e.end(), energy_keV_u) - e.begin());
  return h[k - 1] + (h[k] - h[k - 1]) * (energy_keV_u - e[k - 1]) / (e[k] - e[k - 1]);
}

double YangChargeStateTerm(const YangCoefficients &yang, int z1, int z2, double energy_MeV_u) {
  if (z1 == 1) {
    const auto &b = yang.b;
    const double g = b[2] * (1.0 - std::exp(-b[3] * energy_MeV_u));
    return b[0] * g / ((energy_MeV_u - b[1]) * (energy_MeV_u - b[1]) + g * g);
  }
  const auto &c = yang.c;
  const double e = energy_MeV_u / (std::pow(static_cast<double>(z1), 1.5) * std::sqrt(static_cast<double>(z2)));
  const double g = c[2] * (1.0 - std::exp(-c[3] * e));
  return std::pow(static_cast<double>(z1), 4.0 / 3.0) / std::cbrt(static_cast<double>(z2)) * c[0] * g /
         ((e - c[1]) * (e - c[1]) + g * g);
}

double SpecificStragglingVariance(const StragglingModel &model, const Material &material, const Projectile &projectile,
                                  double energy_MeV) {
  if (model.kind != StragglingKind::kBohr && model.chu == nullptr) {
    throw std::invalid_argument("stopping: Chu straggling needs a Chu table (none is shipped; docs/upgrade/stopping.md section 10)");
  }
  if (model.kind == StragglingKind::kChuYang && model.yang == nullptr) {
    throw std::invalid_argument("stopping: Yang straggling needs its coefficients (none are shipped; docs/upgrade/stopping.md section 10)");
  }
  double variance = 0.0;
  const double e_keV_u = energy_MeV / projectile.mass_u * 1.0e3;
  for (const auto &c : material.constituents) {
    double ratio = 1.0;
    if (model.kind != StragglingKind::kBohr) {
      ratio = ChuFactor(*model.chu, c.z, e_keV_u);
      if (model.kind == StragglingKind::kChuYang) {
        ratio = YangGamma2(projectile.z, c.z, energy_MeV, projectile.mass_u) * ratio +
                YangChargeStateTerm(*model.yang, projectile.z, c.z, e_keV_u * 1.0e-3);
      }
    }
    variance += c.atom_fraction * ratio * BohrSpecificVariance(projectile.z, c.z);
  }
  return variance;
}

StragglingWidths PropagateStraggling(const StoppingFunction &epsilon, const std::function<double(double)> &specific_variance,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV) {
  if (!(beam_sigma_MeV >= 0.0)) throw std::invalid_argument("stopping: beam spread must be >= 0");
  StragglingWidths w;
  w.energy_out_MeV = EnergyAfter(epsilon, energy_in_MeV, areal_density);
  // S in MeV cm^2 per atom
  const auto s = [&epsilon](double e) { return epsilon(e) * kStoppingToMeVcm2; };
  const double s_in = s(energy_in_MeV);
  const double s_out = s(w.energy_out_MeV);
  const double integral = IntegrateInSqrtE(
      [&](double e) {
        const double se = s(e);
        return specific_variance(e) / (se * se * se);
      },
      w.energy_out_MeV, energy_in_MeV, kQuadratureRelTol);
  w.straggling = s_out * std::sqrt(integral);
  w.beam = beam_sigma_MeV * s_out / s_in;
  w.total = std::hypot(w.straggling, w.beam);
  return w;
}

StragglingWidths PropagateStraggling(const StragglingModel &model, const Material &material, const Projectile &projectile,
                                     double energy_in_MeV, double areal_density, double beam_sigma_MeV) {
  return PropagateStraggling(
      MaterialStopping(material, projectile),
      [&](double e) { return SpecificStragglingVariance(model, material, projectile, e); }, energy_in_MeV, areal_density,
      beam_sigma_MeV);
}

}  // namespace stopping
}  // namespace azure2
