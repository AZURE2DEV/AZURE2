// Width scales of a channel: the Wigner limit, the single-particle width, the photon penetrability, the Weisskopf
// estimate of a gamma width (Widths.h).
#include "azure2/physics/Widths.h"

#include <cmath>
#include <limits>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace physics {

double wigner_limit(double mu_u, double a_fm) {
  if (!(mu_u > 0.0) || !(a_fm > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double mu_MeV = mu_u * kAtomicMassUnit_MeV;  // MeV/c^2
  return 1.5 * std::pow(kHbarC_MeV_fm, 2.0) / (mu_MeV * std::pow(a_fm, 2.0));
}

double single_particle_width(const ChannelKinematics &kin, int l, double a_fm, double E_MeV, coulomb::Method method) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (!(E_MeV > 0.0) || !(a_fm > 0.0) || !(kin.mu_u > 0.0)) return nan;
  const double r = rho(kin, a_fm, E_MeV);
  const coulomb::CoulombWaves w = coulomb::coulomb(l, eta(kin, E_MeV), r, method);
  if (w.status != coulomb::Status::ok) return nan;
  return 2.0 * penetrability(w, r) * wigner_limit(kin.mu_u, a_fm);
}

double photon_penetrability(double E_gamma_MeV, int L) {
  return std::pow(std::fabs(E_gamma_MeV) / kHbarC_MeV_fm, 2.0 * L + 1);
}

double weisskopf_width_eV(bool magnetic, int L, double E_gamma_MeV, double mass_number) {
  // eV MeV^-(2L+1), L = 1..4: the two-digit textbook coefficients for r0 = 1.2 fm (Widths.h).
  constexpr double kElectric[4] = {6.8e-2, 4.9e-8, 2.3e-14, 6.8e-21};
  constexpr double kMagnetic[4] = {2.1e-2, 1.5e-8, 6.8e-15, 2.0e-21};
  if (L < 1 || L > 4 || !(E_gamma_MeV > 0.0) || !(mass_number > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double coefficient = magnetic ? kMagnetic[L - 1] : kElectric[L - 1];
  const int thirds = magnetic ? 2 * L - 2 : 2 * L;  // the power of A, in thirds: R^(2L) or R^(2L-2)
  // A whole power is written as one (A^2, not A^(6/3)), as the table has it.
  const double a_power = thirds == 0 ? 1.0 : std::pow(mass_number, thirds % 3 == 0 ? thirds / 3 : thirds / 3.0);
  return coefficient * a_power * std::pow(E_gamma_MeV, 2 * L + 1);
}

}  // namespace physics
}  // namespace azure2
