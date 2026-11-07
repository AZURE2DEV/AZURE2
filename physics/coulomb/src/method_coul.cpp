// Method::coul: N. Michel's cwfcomp (coul/), as in AZURE2 before 2025-09.  One Coulomb_wave_functions object per
// call; the object owns all its working state.
// The third-party header first: it defines the M_* constants that <cmath> would otherwise define before it.
#include "cwfcomp.H"

#include <complex>

#include "methods.h"

namespace azure2 {
namespace physics {
namespace coulomb {
namespace detail {

Status coul(int l, double eta, double rho, CoulombWaves &out) {
  ::coul::Coulomb_wave_functions cwf(true, std::complex<double>((double)l, 0.0), std::complex<double>(eta, 0.0));
  const std::complex<double> z(rho, 0.0);
  std::complex<double> F, dF, G, dG;
  cwf.F_dF(z, F, dF);
  cwf.G_dG(z, G, dG);
  out.F = F.real();
  out.Fp = dF.real();
  out.G = G.real();
  out.Gp = dG.real();
  return Status::ok;
}

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
