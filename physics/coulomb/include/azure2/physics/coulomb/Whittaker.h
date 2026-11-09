#ifndef AZURE2_PHYSICS_COULOMB_WHITTAKER_H
#define AZURE2_PHYSICS_COULOMB_WHITTAKER_H

#include "azure2/physics/coulomb/Coulomb.h"

namespace azure2 {
namespace physics {
namespace coulomb {

/*!
 * The decaying solution of the bound-channel radial equation u'' - [1 + 2 eta_b/rho + l(l+1)/rho^2] u = 0 for
 * rho = kappa r > 0, eta_b = Z1 Z2 e^2 mu / (hbar^2 kappa) >= 0 and integer l >= 0 (docs/upgrade/coulomb.md
 * section 1.4): the Whittaker function
 *
 *     W_l(rho) = W_{-eta_b, l+1/2}(2 rho) = e^{-rho} (2 rho)^{l+1} U(l + 1 + eta_b, 2l + 2, 2 rho),
 *
 * normalized so that W_l(rho) ~ e^{-rho} (2 rho)^{-eta_b} for large rho (unit ANC).  Primes are d/drho.
 *
 * U is computed from its integral representation (DLMF 13.4.4) by a double-exponential quadrature centred on the
 * maximum of the integrand and summed in logarithms; no GSL.  Against mpmath (40 digits) the relative error of W is
 * below 1e-12 on l = 0-10, eta_b = 0-1000, rho = 1e-4-200 (physics_whittaker test).  Every function is pure: the
 * same arguments give the same bits in any thread and in any order of calls.
 *
 * Status: invalid_argument for l < 0 or l > 10000, eta_b < 0 or not finite, rho <= 0 or not finite; out_of_range
 * (whittaker only) when W or W' is outside the range of a double -- whittaker_scaled then still has the value.  On
 * anything but ok every number of the result is NaN.
 */
struct WhittakerWave {
  double W, Wp;  ///< W_l(rho), dW_l/drho
  Status status;
};

/// W = mantissa * 10^log10_scale, mantissa in [1, 10) (W > 0 for every valid argument); rho W'/W.
struct WhittakerScaled {
  double mantissa;
  long log10_scale;
  double rho_dlog;  ///< rho W'_l/W_l, the bound-channel shift function at rho = kappa a
  Status status;
};

WhittakerWave whittaker(int l, double eta_b, double rho) noexcept;
WhittakerScaled whittaker_scaled(int l, double eta_b, double rho) noexcept;

}  // namespace coulomb
}  // namespace physics
}  // namespace azure2

#endif
