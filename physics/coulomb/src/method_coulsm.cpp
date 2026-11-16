// Method::coulsm: cwfcomp with the SAMMY-side patch (coulsm/), AZURE2's external/coul from 2025-09-05 until it was
// replaced by this directory.  The call is the one CoulFunc::operator() made, so the numbers are the same bits, apart
// from F and F' at eta = 0 since 2026-11-16.
// At eta = 0 the routine climbs in l from sin and cos for F and G alike; upward is F's unstable direction for rho < l,
// so F and F' of l >= 1 are taken from riccati_bessel instead (docs/upgrade/coulomb.md section 10).
// The third-party header first: it defines the M_* constants that <cmath> would otherwise define before it.
#include "cwfcomp_accurate.H"

#include <complex>
#include <vector>

#include "methods.h"

namespace azure2 {
namespace physics {
namespace coulomb {
namespace detail {

Status coulsm(int l, double eta, double rho, CoulombWaves &out) {
  ::coulsm::Coulomb_wave_functions_accurate cwf(true, std::complex<double>((double)l, 0.0), std::complex<double>(eta, 0.0));
  std::complex<double> F, dF, G, dG;
  cwf.compute(std::complex<double>(rho, 0.0), F, dF, G, dG);
  out.F = F.real();
  out.Fp = dF.real();
  out.G = G.real();
  out.Gp = dG.real();
  // eta = 0: the routine's G, G' (stable upward recurrence) and its l = 0 (sin, cos) are kept; F_l = rho j_l(rho)
  // and F_l' = F_{l-1} - (l/rho) F_l replace its F, F', which lose a factor l(2l+1)/rho^2 per step in l
  // (1.6e-10 at l = 1, 2.9e13 at l = 4 for rho = 7.5e-4).
  if (eta == 0.0 && l > 0) {
    std::vector<double> u(l + 1);
    riccati_bessel(l, rho, u.data());
    out.F = u[l];
    out.Fp = u[l - 1] - l / rho * u[l];
  }
  return Status::ok;
}

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
