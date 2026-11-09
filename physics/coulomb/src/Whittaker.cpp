// The Whittaker function W_{-eta_b, l+1/2}(2 rho) without GSL (docs/upgrade/coulomb.md section 1.4).
#include "azure2/physics/coulomb/Whittaker.h"

#include <cmath>
#include <limits>

namespace azure2 {
namespace physics {
namespace coulomb {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kPi = 3.14159265358979323846;
const double kLn10 = 2.30258509299404568402;

/// ln Gamma(x) for x >= 1: x shifted to >= 15 by the recurrence, then Stirling's series to x^-11 (truncation below
/// 1e-17 relative there).  Written here because std::lgamma writes the global signgam (not thread-safe in glibc).
double LogGamma(double x) {
  double shift = 0.0, product = 1.0;
  while (x < 15.0) {
    product *= x;
    x += 1.0;
    if (product > 1e280) {  // keep the product finite (never reached for x >= 1, kept for safety)
      shift += std::log(product);
      product = 1.0;
    }
  }
  shift += std::log(product);
  const double r = 1.0 / x, r2 = r * r;
  // Bernoulli terms B_2k / (2k (2k-1) x^(2k-1)), k = 1..6
  const double series = r * (1.0 / 12.0 + r2 * (-1.0 / 360.0 + r2 * (1.0 / 1260.0 + r2 * (-1.0 / 1680.0 + r2 * (1.0 / 1188.0 + r2 * (-691.0 / 360360.0))))));
  return (x - 0.5) * std::log(x) - x + 0.91893853320467274178 + series - shift;
}

/*!
 * The integrals I_0 = int_0^inf f(t) dt and I_1 = int_0^inf t f(t) dt, f(t) = e^{-z t} t^{a-1} (1+t)^{b-a-1},
 * a >= 1, z > 0, so that U(a, b, z) = I_0 / Gamma(a) and U(a+1, b+1, z) = I_1 / Gamma(a+1) (DLMF 13.4.4): the second
 * integrand is the first times t, so one set of nodes serves both.  f has one maximum t0 > 0 (or decreases from
 * t = 0 when a = 1 and b <= 2 + z); with t = t0 exp(u), u = (pi/2) sinh s, both tails decay double-exponentially in s
 * and the trapezoidal rule converges geometrically in 1/h.  Terms are kept relative to the term at s = 0, so nothing
 * overflows (the tails are cut at e^-40 of it, in the logarithm without the factor cosh s).  The step is halved from 1/2, reusing the previous nodes; the error of the rule roughly squares at each
 * halving, so two levels that agree to 1e-9 leave the finer one good to rounding (physics_whittaker measures it).
 * Returns ln I_0 and I_1 / I_0.
 */
struct Integrals {
  double logI0, ratio;
};

Integrals LogIntegrals(double a, double b, double z) {
  const double c = b - a - 1.0;
  // g'(t) = 0 for g = -z t + (a-1) ln t + c ln(1+t):  z t^2 + (z - (a-1) - c) t - (a-1) = 0.
  const double qa = z, qb = z - (a - 1.0) - c, qc = -(a - 1.0);
  double t0 = 0.0;
  const double disc = qb * qb - 4.0 * qa * qc;
  if (disc >= 0.0) {
    // stable roots: q = -(qb + sign(qb) sqrt(disc)) / 2, roots q/qa and qc/q
    const double q = -0.5 * (qb + std::copysign(std::sqrt(disc), qb));
    const double r1 = q / qa, r2 = (q != 0.0) ? qc / q : 0.0;
    t0 = std::fmax(r1, r2);
  }
  if (!(t0 > 0.0) || !std::isfinite(t0)) t0 = 1.0 / z;  // decreasing from t = 0: centre on its decay length
  const double logT0 = std::log(t0);

  // At node s: t, cosh s, and ln of f(t) t (dt/ds = t (pi/2) cosh s; ln t = ln t0 + u; the factor (pi/2) cosh s is
  // applied outside the logarithm).  false past the range of a double.
  auto logTerm = [&](double s, double &t, double &ch, double &v) {
    const double e = std::exp(s), sh = 0.5 * (e - 1.0 / e);
    ch = 0.5 * (e + 1.0 / e);
    const double u = 0.5 * kPi * sh;
    if (!(u < 700.0 && u > -700.0)) return false;
    t = t0 * std::exp(u);
    if (!(t > 0.0) || !std::isfinite(t)) return false;
    v = -z * t + a * (logT0 + u) + c * std::log1p(t);
    return true;
  };
  double t = 0.0, ch = 1.0, ref = 0.0;
  if (!logTerm(0.0, t, ch, ref)) return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  // Sum over the nodes k h, both sides of 0 (k = 0, 1, 2, ... or the odd k = 1, 3, ... after a halving), while terms
  // exceed e^-40 of the reference.
  auto sumNodes = [&](double h, bool oddOnly, double &s0, double &s1) {
    s0 = s1 = 0.0;
    for (int side = 0; side < 2; side++) {
      const double sign = side == 0 ? 1.0 : -1.0;
      const long first = oddOnly ? 1 : (side == 0 ? 0 : 1), step = oddOnly ? 2 : 1;
      for (long k = first; k < 100000; k += step) {
        const double s = sign * k * h;
        double tk, chk, v;
        if (!logTerm(s, tk, chk, v)) break;
        const double w = std::exp(v - ref) * chk;
        s0 += w;
        s1 += w * tk;
        // past the bulk the tails fall monotonically: the first term below e^-40 ends that side
        if (v - ref < -40.0 && std::fabs(s) > 1.0) break;
      }
    }
  };
  double h = 0.5, sum0, sum1;
  sumNodes(h, false, sum0, sum1);
  double est0 = sum0 * h, est1 = sum1 * h;
  for (int level = 0; level < 10; level++) {
    h *= 0.5;
    double add0, add1;
    sumNodes(h, true, add0, add1);
    sum0 += add0;
    sum1 += add1;
    const double next0 = sum0 * h, next1 = sum1 * h;
    const bool converged = std::fabs(next0 - est0) <= 1e-9 * next0 && std::fabs(next1 - est1) <= 1e-9 * next1;
    est0 = next0;
    est1 = next1;
    if (converged) break;
  }
  return {ref + std::log(0.5 * kPi * est0), est1 / est0};
}

bool Valid(int l, double eta_b, double rho) {
  return l >= 0 && l <= 10000 && std::isfinite(eta_b) && eta_b >= 0.0 && std::isfinite(rho) && rho > 0.0;
}

/// ln W and rho W'/W (primes d/drho).
void LogWhittaker(int l, double eta_b, double rho, double &logW, double &rhoDlog) {
  const double a = l + 1.0 + eta_b, b = 2.0 * l + 2.0, z = 2.0 * rho;
  const Integrals in = LogIntegrals(a, b, z);
  logW = -rho + (l + 1.0) * std::log(z) - LogGamma(a) + in.logI0;
  // W = e^{-rho} (2rho)^{l+1} U(a,b,2rho), dU/dz = -a U(a+1,b+1,z), a U(a+1,b+1,z)/U(a,b,z) = I_1/I_0:
  // rho W'/W = -rho + (l+1) - 2 rho a U(a+1,b+1,2rho)/U(a,b,2rho) = -rho + (l+1) - 2 rho I_1/I_0
  rhoDlog = -rho + (l + 1.0) - 2.0 * rho * in.ratio;
}

}  // namespace

WhittakerScaled whittaker_scaled(int l, double eta_b, double rho) noexcept {
  if (!Valid(l, eta_b, rho)) return {kNaN, 0, kNaN, Status::invalid_argument};
  double logW, rhoDlog;
  LogWhittaker(l, eta_b, rho, logW, rhoDlog);
  if (!std::isfinite(logW) || !std::isfinite(rhoDlog)) return {kNaN, 0, kNaN, Status::failed};
  const double log10W = logW / kLn10;
  double intPart = std::floor(log10W);
  double mantissa = std::pow(10.0, log10W - intPart);
  if (mantissa >= 10.0) {  // rounding at the boundary
    mantissa /= 10.0;
    intPart += 1.0;
  }
  return {mantissa, (long)intPart, rhoDlog, Status::ok};
}

WhittakerWave whittaker(int l, double eta_b, double rho) noexcept {
  if (!Valid(l, eta_b, rho)) return {kNaN, kNaN, Status::invalid_argument};
  double logW, rhoDlog;
  LogWhittaker(l, eta_b, rho, logW, rhoDlog);
  if (!std::isfinite(logW) || !std::isfinite(rhoDlog)) return {kNaN, kNaN, Status::failed};
  const double W = std::exp(logW);
  const double Wp = W * rhoDlog / rho;
  if (!(W > 0.0) || !std::isfinite(W) || !std::isfinite(Wp) || (Wp == 0.0 && rhoDlog != 0.0) ||
      W < std::numeric_limits<double>::min())
    return {kNaN, kNaN, Status::out_of_range};
  return {W, Wp, Status::ok};
}

}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
