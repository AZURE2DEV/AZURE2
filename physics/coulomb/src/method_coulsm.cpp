// Method::coulsm: cwfcomp with the SAMMY-side patch (coulsm/), AZURE2's external/coul from 2025-09-05 until it was
// replaced by this directory.  The call is the one CoulFunc::operator() made, so the numbers are the same bits.
// The third-party header first: it defines the M_* constants that <cmath> would otherwise define before it.
#include "cwfcomp_accurate.H"

#include <complex>

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
  return Status::ok;
}

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
