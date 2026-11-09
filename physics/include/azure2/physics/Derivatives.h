#ifndef AZURE2_PHYSICS_DERIVATIVES_H
#define AZURE2_PHYSICS_DERIVATIVES_H

#include <functional>

namespace azure2 {
namespace physics {

/*!
 * Numerical derivatives of the channel functions (docs/upgrade/coulomb.md section 1.4, "Bound channel functions in
 * physics/").
 */

/*!
 * df/dx at x with step h by the algorithm of GSL's gsl_deriv_central (GSL 2.x, deriv/deriv.c, B. Gough; GPL-3, as
 * AZURE2): the 5-point central rule on x +- h, x +- h/2, and, when the rounding error estimate is below the truncation
 * error estimate, once more at h (round / (2 trunc))^(1/3), kept if its error is smaller and it agrees.  For the same f
 * it gives the same bits as GSL (tests/upgrade/derivatives_gsl_test.cpp).  abserr, if not null, gets GSL's estimate.
 */
double central_derivative(const std::function<double(double)> &f, double x, double h, double *abserr = nullptr);

/// Which side of x a one-sided stencil extends to.
enum class Side { Below, Above };

/*!
 * d^2f/dx^2 at x: the Richardson-extrapolated central second difference (4 D(h) - D(2h)) / 3 with step h, or, when
 * x is within 4h of the edge of f's domain (one_sided), the stencil (2 f0 - 5 f1 + 4 f2 - f3) / g^2 on x, x -+ g,
 * x -+ 2g, x -+ 3g toward `side`, as the engine's ShftFunc (Below) and CoulFunc (Above) have it (h = 1e-2, g = 1e-3).
 */
double second_derivative(const std::function<double(double)> &f, double x, double h, double g, bool one_sided,
                         Side side);

}  // namespace physics
}  // namespace azure2

#endif
