#include "azure2/stopping/Electronic.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace azure2 {
namespace stopping {

namespace {

/// sum_k c_k x^k in ascending order, each power by std::pow (lunasrim _poly; the same order of the sum).
template <std::size_t N>
double Poly(const std::array<double, N> &c, double x) {
  double s = 0.0;
  for (std::size_t k = 0; k < N; ++k) s += c[k] * std::pow(x, static_cast<double>(k));
  return s;
}

double Log10(double x) { return std::log(x) / std::log(10.0); }

}  // namespace

double SrimFit(const std::array<double, 10> &c, double e) {
  // lunasrim _fit: S_low = c1 e^c2 + c3 e^c4 + c9 e^0.5, S_high = c5 / e^c6 (ln(c7/e + c8 e) + c10 / e^0.2),
  // S = S_low S_high / (S_low + S_high)
  const double slow = c[0] * std::pow(e, c[1]) + c[2] * std::pow(e, c[3]) + c[8] * std::pow(e, 0.5);
  const double shigh = (c[4] / std::pow(e, c[5])) * (std::log(c[6] / e + c[7] * e) + c[9] / std::pow(e, 0.2));
  return slow * shigh / (slow + shigh);
}

double SrimBethe(const SrimElement &target, int z1, double e, Phase phase) {
  // lunasrim bethe: e in keV/u; 931494 keV/u, 1022000 eV = 2 m_e c^2, 1961.16 = 1 / (4 pi r_e^2 m_e c^2) in SRIM's
  // units, 137.036 = 1/alpha, as SRIM writes them.
  const bool gas = phase == Phase::kGas;
  const double z2 = target.z;
  const double beta2 = 1.0 - 1.0 / std::pow(1.0 + e / 931494.0, 2);
  const double lx = std::log(e / 1000.0);
  const double shell_d = gas ? 0.0 : Poly(target.shell_d, lx);
  const double shell_c = Poly(target.shell_c, lx);
  double mean_i = 0.0;
  double q = 0.0;
  if (gas) {
    mean_i = z1 == 1 ? target.mean_excitation_gas_eV[0] : target.mean_excitation_gas_eV[1];
    q = 1.35 / std::pow(e, 0.4);
  } else {
    mean_i = target.mean_excitation_solid_eV;
    q = 1.5 / std::pow(e, 0.4) + 45000.0 / (z2 * std::pow(e, 1.6));
  }
  const double barkas = (z1 * (0.001 * e) * q) / (0.001 * e + q);
  const double y = std::pow(z1 / 137.036, 2) / beta2;
  const double bloch = -y * (1.20206 - y * (1.042 + -0.8549 * y + 0.343 * y * y));
  const double l = (std::log(1022000.0 * beta2 / (1.0 - beta2)) - beta2 - shell_d - shell_c - std::log(mean_i) +
                    barkas + bloch);
  return (static_cast<double>(z1) * z1 * z2) / (beta2 * 1961.16) * l;
}

double SrimLightIon(const SrimElement &target, const std::array<double, 10> &c, int z1, double e, double fit_max,
                    Phase phase) {
  // lunasrim _light_ion
  const double emin = kSrimVelocityLimit_keV_u;
  const double ehigh = fit_max * kBetheJoinFactor;
  if (e > ehigh) return SrimBethe(target, z1, e, phase);
  const double s = SrimFit(c, std::min(std::max(e, emin), fit_max));
  if (e > fit_max) {
    const double sb = SrimBethe(target, z1, ehigh, phase);
    return s + (sb - s) * (e - fit_max) / (ehigh - fit_max);
  }
  if (e < emin) return s * std::pow(e / emin, 0.5);
  return s;
}

double SrimHighEnergyCharge(const std::array<double, 5> &hi, double e) {
  // lunasrim _hi_charge
  const double d = Log10(e) - hi[1];
  const double w = d >= 0 ? hi[3] : hi[2];
  const double x = std::exp(-d / w) + 1.0;
  return hi[0] + (1.0 - hi[0]) / x + hi[4] * std::exp(-std::fabs(1.3 + hi[1] - Log10(e)));
}

double SrimElementStopping(const SrimData &data, int z1, int z2, double e, Phase phase) {
  // lunasrim electronic_stopping_element
  const SrimElement &t = ElementData(data, z2);
  if (z1 == 1) return SrimLightIon(t, t.proton, 1, e, kProtonFitMax_keV_u, phase);
  if (z1 == 2) return SrimLightIon(t, t.helium, 2, e, kHeliumLithiumFitMax_keV_u, phase);
  if (z1 == 3) return SrimLightIon(t, t.lithium, 3, e, kHeliumLithiumFitMax_keV_u, phase);
  // heavy ions: helium stopping at the same e, scaled
  const double emin = kSrimVelocityLimit_keV_u;
  const double emax = kHeavyIonScalingMax_keV_u;
  const double ec = std::min(std::max(e, emin), emax);
  const double lne = std::log(ec);
  const SrimElement &ion = ElementData(data, z1);
  const double s = SrimLightIon(t, t.helium, 2, ec, kHeliumLithiumFitMax_keV_u, phase) * Poly(ion.ion_poly, lne) * z1 *
                   z1 / 4.0 * Poly(t.target_poly, lne);
  if (e < emin) return s * std::pow(e / emin, 0.5);
  if (e > emax) {
    const double r_hi = SrimHighEnergyCharge(ion.high_energy_charge, e) * SrimBethe(t, 2, e, phase);
    const double r_lo = SrimHighEnergyCharge(ion.high_energy_charge, emax) * SrimBethe(t, 2, emax, phase);
    return s * r_hi / r_lo;
  }
  return s;
}

double SrimPhaseFactor(int z2, Phase phase) {
  // lunasrim phase_factor
  for (int noble : {2, 10, 18, 36, 54, 86}) {
    if (z2 == noble) return 1.0;
  }
  const bool natural_gas = IsNaturalGas(z2);
  if (natural_gas == (phase == Phase::kGas)) return 1.0;
  return natural_gas ? kSrimPhaseFactor : 1.0 / kSrimPhaseFactor;
}

double SrimFadeOut(double factor, double e, int z1) {
  // lunasrim fade_correction (above x = 700 the fade is complete; std::exp would overflow near 709.78)
  const double c = z1 == 1 ? 7.0 : 4.0;
  const double x = 1.48 * (std::sqrt(e / 25.0) - c);
  if (x > 700) return 1.0;
  return 1.0 + (factor - 1.0) / (1.0 + std::exp(x));
}

double ElectronicRelUncertainty(int z1) {
  if (z1 == 1) return 0.040;
  if (z1 == 2) return 0.039;
  if (z1 == 3) return 0.048;
  return 0.058;
}

StoppingWithUncertainty ElectronicStopping(const SrimData &data, const Projectile &projectile, int z2, double energy_MeV,
                                           Phase phase, double compound_correction) {
  // lunasrim electronic_stopping
  const int z1 = projectile.z;
  if (z1 < 1 || z1 > kSrimMaxZ) throw std::out_of_range("stopping: projectile Z1 = " + std::to_string(z1) + " outside 1..92");
  const double e = energy_MeV * 1.0e3 / projectile.mass_u;
  if (!(energy_MeV > 0.0) || e > kSrimMaxEnergy_keV_u) {
    std::ostringstream o;
    o << "stopping: energy " << energy_MeV << " MeV of " << projectile.label << " outside SRIM's range (0, "
      << kSrimMaxEnergy_keV_u << " keV/u]";
    throw std::domain_error(o.str());
  }
  const double s = SrimElementStopping(data, z1, z2, e, phase);
  if (!(compound_correction > 0.0)) throw std::invalid_argument("stopping: the compound correction must be > 0");
  return {s * SrimFadeOut(SrimPhaseFactor(z2, phase) * compound_correction, e, z1), ElectronicRelUncertainty(z1)};
}

StoppingWithUncertainty ElectronicStopping(const Projectile &projectile, int z2, double energy_MeV, Phase phase,
                                           double compound_correction) {
  return ElectronicStopping(SrimTables(), projectile, z2, energy_MeV, phase, compound_correction);
}

}  // namespace stopping
}  // namespace azure2
