#include "azure2/stopping/Electronic.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace stopping {

namespace {

using physics::kAtomicMassUnit_MeV;
using physics::kElectronMass_MeV;
using physics::kFineStructure;
using physics::kHbarC_MeV_fm;

constexpr double kPi = 3.14159265358979323846;
constexpr double kProtonNuclearMass_u = 1.007276466621;  ///< the variable of the AZ_p form is proton keV / this
constexpr double kAlphaNuclearMass_u = 4.001506179127;   ///< the variable of the Z_He form is He-4 MeV at this mass

/// K = 4 pi r_e^2 m_e c^2 in eV/(1e15 atoms/cm^2): r_e = alpha hbar c / m_e c^2 (fm -> cm 1e-13), m_e c^2 in eV.
double BetheConstant() {
  const double re_cm = kFineStructure * kHbarC_MeV_fm / kElectronMass_MeV * 1.0e-13;
  return 4.0 * kPi * re_cm * re_cm * kElectronMass_MeV * 1.0e6 * 1.0e15;
}

/// beta^2 at kinetic energy per nucleon mass unit t_u (MeV per u): gamma = 1 + t_u / (u c^2).
double Beta2PerU(double t_MeV_u) {
  const double gamma = 1.0 + t_MeV_u / kAtomicMassUnit_MeV;
  return 1.0 - 1.0 / (gamma * gamma);
}

/// F(W) of Ashley, Ritchie and Brandt (Phys. Rev. B 5 (1972) 2393) for the Barkas term, the 47-point table of
/// Geant4's G4EmCorrections (Geant4 Software License 1.0); linear interpolation, F(0.02) below the table, F ~ 1/W
/// above it.
double BarkasF(double w) {
  static const std::array<std::array<double, 2>, 47> kTable = {{
      {0.02, 21.5}, {0.03, 20.0}, {0.04, 18.0}, {0.05, 15.6}, {0.06, 15.0}, {0.07, 14.0}, {0.08, 13.5}, {0.09, 13.0},
      {0.1, 12.2},  {0.2, 9.25},  {0.3, 7.0},   {0.4, 6.0},   {0.5, 4.5},   {0.6, 3.5},   {0.7, 3.0},   {0.8, 2.5},
      {0.9, 2.0},   {1.0, 1.7},   {1.2, 1.2},   {1.3, 1.0},   {1.4, 0.86},  {1.5, 0.7},   {1.6, 0.61},  {1.7, 0.52},
      {1.8, 0.5},   {1.9, 0.43},  {2.0, 0.42},  {2.1, 0.3},   {2.4, 0.2},   {3.0, 0.13},  {3.08, 0.1},  {3.1, 0.09},
      {3.3, 0.08},  {3.5, 0.07},  {3.8, 0.06},  {4.0, 0.051}, {4.1, 0.04},  {4.8, 0.03},  {5.0, 0.024}, {5.1, 0.02},
      {6.0, 0.013}, {6.5, 0.01},  {7.0, 0.009}, {7.1, 0.008}, {8.0, 0.006}, {9.0, 0.0032}, {10.0, 0.0025}}};
  if (w <= kTable.front()[0]) return kTable.front()[1];
  if (w >= kTable.back()[0]) return kTable.back()[1] * kTable.back()[0] / w;
  std::size_t k = 1;
  while (kTable[k][0] < w) ++k;
  const auto &lo = kTable[k - 1];
  const auto &hi = kTable[k];
  return lo[1] + (hi[1] - lo[1]) * (w - lo[0]) / (hi[0] - lo[0]);
}

/// The screening parameter b of the Barkas term per target (Geant4 G4EmCorrections, after ICRU Report 49).
double BarkasB(int z2) {
  if (z2 == 1) return 1.8;
  if (z2 == 2) return 0.6;
  if (z2 <= 10) return 1.8;
  if (z2 <= 17) return 1.4;
  if (z2 == 18) return 1.8;
  if (z2 <= 25) return 1.4;
  if (z2 <= 50) return 1.35;
  return 1.3;
}

/// The low-energy fit of the row at energy per nucleon mass unit t_u (MeV/u), at the same velocity.
double FitAt(const CoefficientRow &row, double t_MeV_u) {
  if (row.form == FitForm::kAndersenZieglerProton) return ProtonFitStopping(row, t_MeV_u * 1.0e3);
  return HeliumFitStopping(row, t_MeV_u * kAlphaNuclearMass_u);
}

/// Electronic stopping of hydrogen (z1 = 1) or helium (z1 = 2) ions at t_u (MeV/u, the same-velocity energy).
StoppingWithUncertainty LightIonStopping(int z1, int z2, double t_MeV_u) {
  const CoefficientRow &row = Coefficients(z2, z1);
  const double i_eV = Element(z2).mean_excitation_eV;
  const double t1 = row.fit_max_MeV_u;
  const double t2 = kBetheFrom_MeV_u;
  if (t_MeV_u <= t1) return {FitAt(row, t_MeV_u), row.rel_uncertainty};
  if (t_MeV_u >= t2) return {BetheStopping(z1, z2, i_eV, Beta2PerU(t_MeV_u)), kBetheRelUncertainty};
  // between the fit and Bethe: linear in log S - log T
  const double s1 = std::log(FitAt(row, t1));
  const double s2 = std::log(BetheStopping(z1, z2, i_eV, Beta2PerU(t2)));
  const double x = std::log(t_MeV_u / t1) / std::log(t2 / t1);
  return {std::exp(s1 + (s2 - s1) * x), std::max(row.rel_uncertainty, kBetheRelUncertainty)};
}

}  // namespace

double ProtonFitStopping(const CoefficientRow &row, double t_keV) {
  if (!(t_keV > 0.0)) throw std::domain_error("stopping: proton fit at T <= 0");
  double scale = 1.0;
  double t = t_keV;
  if (t < 10.0) {
    scale = std::sqrt(t / 10.0);
    t = 10.0;
  }
  const double s_low = row.a[1] * std::pow(t, 0.45);
  const double s_high = row.a[2] / t * std::log(1.0 + row.a[3] / t + row.a[4] * t);
  return scale * s_low * s_high / (s_low + s_high);
}

double HeliumFitStopping(const CoefficientRow &row, double t_MeV) {
  if (!(t_MeV > 0.0)) throw std::domain_error("stopping: helium fit at T <= 0");
  double scale = 1.0;
  double t = t_MeV;
  if (t < 1.0e-3) {
    scale = std::sqrt(t / 1.0e-3);
    t = 1.0e-3;
  }
  const double s_low = row.a[0] * std::pow(1000.0 * t, row.a[1]);
  const double s_high = row.a[2] / t * std::log(1.0 + row.a[3] / t + row.a[4] * t);
  return scale * s_low * s_high / (s_low + s_high);
}

double BlochCorrection(double y) {
  // -y^2 sum_{n>=1} 1/(n (n^2 + y^2)); the tail beyond N is y^2-independent to O(1/N^4): 1/(2 N^2)
  constexpr int kTerms = 2000;
  double sum = 0.0;
  for (int n = kTerms; n >= 1; --n) sum += 1.0 / (n * (static_cast<double>(n) * n + y * y));
  sum += 1.0 / (2.0 * kTerms * static_cast<double>(kTerms));
  return -y * y * sum;
}

double ShellCorrectionBarkasBerger(double mean_excitation_eV, double eta) {
  const double e2 = 1.0 / (eta * eta);
  const double e4 = e2 * e2;
  const double e6 = e4 * e2;
  const double i = mean_excitation_eV;
  return (0.422377 * e2 + 0.0304043 * e4 - 0.00038106 * e6) * 1.0e-6 * i * i +
         (3.858019 * e2 - 0.1667989 * e4 + 0.00157955 * e6) * 1.0e-9 * i * i * i;
}

double BarkasCorrection(int z2, double beta2) {
  const double x = beta2 / (kFineStructure * kFineStructure * z2);
  const double w = BarkasB(z2) / std::sqrt(x);
  return BarkasF(w) / (std::sqrt(static_cast<double>(z2)) * x * std::sqrt(x));
}

BetheTerms BetheStoppingTerms(int z1, int z2, double mean_excitation_eV, double beta2) {
  BetheTerms t;
  const double gamma2 = 1.0 / (1.0 - beta2);
  const double eta = std::sqrt(beta2 * gamma2);
  t.l0 = std::log(2.0 * kElectronMass_MeV * 1.0e6 * beta2 * gamma2 / mean_excitation_eV) - beta2;
  t.shell = ShellCorrectionBarkasBerger(mean_excitation_eV, std::max(eta, 0.13)) / z2;
  t.barkas = BarkasCorrection(z2, beta2);
  t.bloch = BlochCorrection(z1 * kFineStructure / std::sqrt(beta2));
  return t;
}

double BetheStopping(int z1, int z2, double mean_excitation_eV, double beta2) {
  const BetheTerms t = BetheStoppingTerms(z1, z2, mean_excitation_eV, beta2);
  const double l = t.l0 - t.shell + z1 * t.barkas + z1 * z1 * t.bloch;
  return BetheConstant() * z1 * z1 * z2 / beta2 * l;
}

double EffectiveChargeFraction(int z1, double beta2) {
  const double v_over_v0 = std::sqrt(beta2) / kFineStructure;
  const double vr = v_over_v0 / std::pow(static_cast<double>(z1), 2.0 / 3.0);
  return 1.0 - std::exp(-0.95 * vr);
}

StoppingWithUncertainty ElectronicStopping(const Projectile &projectile, int z2, double energy_MeV) {
  if (!(energy_MeV > 0.0)) throw std::domain_error("stopping: energy must be > 0");
  // the same-velocity energy per mass unit: (gamma - 1) u c^2 = E / M (M in u)
  const double t_MeV_u = energy_MeV / projectile.mass_u;
  if (projectile.z <= 2) return LightIonStopping(projectile.z, z2, t_MeV_u);
  // heavier ions: the proton value at the same velocity times (q z1)^2
  const double beta2 = Beta2(projectile, energy_MeV);
  const double q = EffectiveChargeFraction(projectile.z, beta2) * projectile.z;
  const StoppingWithUncertainty p = LightIonStopping(1, z2, t_MeV_u);
  const double vr = std::sqrt(beta2) / kFineStructure / std::pow(static_cast<double>(projectile.z), 2.0 / 3.0);
  return {q * q * p.value, vr >= 1.0 ? kHeavyIonRelUncertainty : kSlowHeavyIonRelUncertainty};
}

}  // namespace stopping
}  // namespace azure2
