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
  if (L < 1 || L > 4 || !(E_gamma_MeV > 0.0) || !(mass_number > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  // The closed form (Widths.h; Bohr and Mottelson vol. 1, eqs. 3C-38, 3C-39), R = r0 A^(1/3).
  constexpr double kR0_fm = 1.2;
  double doubleFactorial = 1.0;  // (2L+1)!!
  for (int k = 2 * L + 1; k > 1; k -= 2) doubleFactorial *= k;
  const double R = kR0_fm * std::cbrt(mass_number);
  const double common = 2.0 * (L + 1) / (L * doubleFactorial * doubleFactorial) * std::pow(3.0 / (L + 3), 2) *
                        kFineStructure * kHbarC_MeV_fm * std::pow(E_gamma_MeV / kHbarC_MeV_fm, 2 * L + 1);  // MeV fm^-2L
  const double width_MeV =
      magnetic ? common * 40.0 * std::pow(kHbarC_MeV_fm / (2.0 * kProtonMass_MeV), 2) * std::pow(R, 2 * L - 2)
               : common * std::pow(R, 2 * L);
  return width_MeV * 1.0e6;
}

}  // namespace physics
}  // namespace azure2
