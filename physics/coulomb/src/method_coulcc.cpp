// Method::coulcc: COULCC (coulcc/, the C++17 translation that reproduces coulcc.f bit for bit), called with real rho,
// eta and lambda = l, F and G with their derivatives, no phase (Coulomb.cpp sets sigma for every method alike).
#include "coulcc.h"

#include <complex>

#include "methods.h"

namespace azure2 {
namespace physics {
namespace coulomb {
namespace detail {

Status coulcc(int l, double eta, double rho, CoulombWaves &out) {
  ::coulcc::Options opt;
  opt.kind = ::coulcc::Functions::CoulombNoPhase;
  const ::coulcc::Result r = ::coulcc::coulcc(std::complex<double>(rho, 0.0), std::complex<double>(eta, 0.0),
                                               std::complex<double>((double)l, 0.0), 1, opt);
  if (!r.ok() || r.F.empty() || r.G.empty() || r.Fp.empty() || r.Gp.empty()) return Status::failed;
  out.F = r.F[0].real();
  out.Fp = r.Fp[0].real();
  out.G = r.G[0].real();
  out.Gp = r.Gp[0].real();
  return Status::ok;
}

}  // namespace detail
}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
