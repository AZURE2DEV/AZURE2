#include "ThmFunc.h"
#include "Constants.h"
#include "PPair.h"
#include "WhitFunc.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

//! forward finite-difference step for dj_l/drho (as in mrmpy's model.py)
static const double kThmFDStep = 1.0e-6;

double ThmSphericalBessel(int l, double x) {
  return gsl_sf_bessel_jl(l, x);
}

double ThmRho(double mu, double E, double B, double radius) {
  // p = sqrt(2 mu (E + B)) in MeV/c; rho = p r / (hbar c), dimensionless.
  return std::sqrt(2.0 * mu * (E + B)) * radius / hbarc;
}

void ThmBesselParts(int l, double mu, double E, double B, double radius,
                    double &jl, double &rhoDjl) {
  double rho = ThmRho(mu, E, B, radius);
  jl = gsl_sf_bessel_jl(l, rho);

  /* dj_l/drho by forward finite difference (matches mrmpy), then rho*dj_l/drho.

  Maybe would be better to switch to the analytic j_l'(x) = (l/x) j_l(x) - j_{l+1}(x)
  (first-order FD truncation ~1e-6 is the larger error for rho);
  on pause while we need to match mrmpy's behavior */

  rhoDjl = rho * (gsl_sf_bessel_jl(l, rho + kThmFDStep) - jl) / kThmFDStep;
}

double ThmFormFactor(int l, double b, double mu, double E, double B,
                     double radius) {
  double jl, rhoDjl;
  ThmBesselParts(l, mu, E, B, radius, jl, rhoDjl);
  return (b - 1.0) * jl - rhoDjl;
}

/* ------------------------------------------------------------------------
 * External Coulomb term of the THM entrance vertex
 *
 *   C_l(E) = 2 eta k Int_a^inf dr O_l(k r)/O_l(k a) j_l(p r)
 *
 * See include/ThmFunc.h for the physics and references.  Numerics, by regime
 * (q^2(r) = k^2 - 2 eta k/r - l(l+1)/r^2 is the local squared wave number,
 * S = Int_a^r_tp |q| dr the WKB barrier integral from a to the outer turning
 * point r_tp):
 *
 *  - Real axis, all regimes: composite 20-point Gauss-Legendre panels of width
 *    10/(p + max(|q|, k)) (at most a third of r), i.e. <= ~10 radians of
 *    phase or e-folds of decay per panel.  The integrand is analytic there.
 *    A panel needs O_l at 20 nodes; the Taylor segments supply them from one
 *    set of coefficients each, and O_l varies on the scale 1/k (not 1/p), so a
 *    segment serves several panels.
 *
 *  - E > 0, S < kBarrierSwitch (above or near the barrier): the real-axis
 *    part runs from a to R = rho_m / k, where the asymptotic series of H+
 *    (DLMF 33.11.1) converges to double precision.  O_l on [a, R] is obtained
 *    by integrating the Coulomb equation inwards from that series value at R
 *    in Taylor segments (RadialTaylor); inwards all solutions grow like G_l
 *    through the barrier, so this is stable, and no library Coulomb function is
 *    needed (useGSL has no effect).  The real-axis quadrature may stop early
 *    when |O(r)/O(a)| has fallen so far that the rest is negligible (restBound).
 *    The conditionally convergent tail R -> inf is done exactly:
 *    j_l = (h1_l + h2_l)/2 and each half is integrated along a ray into the
 *    complex r plane on which it decays, H+ h1 ~ exp(i(k+p)r) along
 *    r = R + i t and H+ h2 ~ exp(i(k-p)r) along r = R - i t (p > k).  By
 *    Cauchy and Jordan's lemma this is the value of the real integral (and its
 *    Abel limit).  H+ at complex r comes from the asymptotic series,
 *    theta = rho - eta ln 2rho - l pi/2 + sigma_l.
 *
 *  - E > 0, S >= kBarrierSwitch (deep sub-barrier, e.g. 12C+12C below
 *    ~0.3 MeV): G_l(ka) grows like e^S and overflows a double for S > ~700
 *    (the built-in cwfcomp routine then returns garbage).  But F_l/G_l ~ e^{-2S}
 *    and everything beyond the turning point is suppressed by e^{-S}, so
 *    O(r)/O(a) = G(r)/G(a) to double precision.  G is the solution that
 *    decays outwards, which is the dominant one when integrating INWARDS: it is
 *    integrated from r_c (where the barrier integral from a reaches
 *    kBarrierSwitch) down to a with the Taylor series, starting from the WKB
 *    log-derivative; the error of that starting value dies off like
 *    e^{-2 S(a, r_c)}.  Works for any eta.
 *
 *  - E < 0: O_l -> W_{-eta, l+1/2}(2 kappa r) from WhitFunc::Scaled (the
 *    exponent carried separately), which decays monotonically; integrated
 *    until the rest is negligible.
 *
 *  - E = 0 exactly: the E -> 0 limit from either side, the zero-energy decaying
 *    solution O(r)/O(a) -> sqrt(r/a) K_{2l+1}(sqrt(8 c r)) / K_{2l+1}(sqrt(8 c a)),
 *    c = eta k.  C_l(0) is real and C_l(E) is continuous through threshold.
 * ------------------------------------------------------------------------ */

namespace {

const int kGLOrder = 20;              // nodes per panel
const double kPanelPhase = 10.0;      // max radians / e-folds per panel
const double kRelTol = 1.0e-13;       // early-stop tolerance (relative)
const double kBarrierSwitch = 32.0;   // WKB barrier integral for the inward mode
const double kTaylorPhase = 8.0;      // max radians / e-folds of O_l per Taylor seed
const int kMaxPanels = 200000;        // runaway guard

struct GLRule {
  double x[kGLOrder];
  double w[kGLOrder];
  int last;  // index of the node nearest +1
};

const GLRule &glRule() {
  static const GLRule rule = [] {
    GLRule r;
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(kGLOrder);
    for (int i = 0; i < kGLOrder; i++) gsl_integration_glfixed_point(-1.0, 1.0, i, &r.x[i], &r.w[i], t);
    gsl_integration_glfixed_table_free(t);
    r.last = 0;
    for (int i = 1; i < kGLOrder; i++)
      if (r.x[i] > r.x[r.last]) r.last = i;
    return r;
  }();
  return rule;
}

double sphJl(int l, double x) {
  gsl_sf_result r;
  if (gsl_sf_bessel_jl_e(l, x, &r) != GSL_SUCCESS) return 0.0;
  return r.val;
}

/// Spherical Hankel h^(1)_l (sign = +1) or h^(2)_l (sign = -1) at complex x.
complex sphHankel(int l, complex x, int sign) {
  const complex is(0.0, (double)sign);
  complex sum(0.0, 0.0), term(1.0, 0.0);
  const complex inv2x = 1.0 / (2.0 * x);
  for (int m = 0; m <= l; m++) {
    sum += term;
    term *= is * (double)((l + m + 1) * (l - m)) / (double)(m + 1) * inv2x;
  }
  complex pre(1.0, 0.0);
  for (int m = 0; m <= l; m++) pre *= -is;  // (-+i)^(l+1)
  return pre * std::exp(is * x) / x * sum;
}

/// Asymptotic series of H+_l(eta, rho) = G + i F (DLMF 33.11.1) at complex rho.
/// Returns false if it does not reach double precision before diverging.
bool hPlusAsymptotic(int l, double eta, double sigma, complex rho, complex &h, complex *dh = nullptr) {
  const complex a(1.0 + l, eta), b(-(double)l, eta);
  const complex inv2irho = 1.0 / (complex(0.0, 2.0) * rho);
  complex term(1.0, 0.0), sum(1.0, 0.0), dsum(0.0, 0.0);  // dsum = rho dS/drho
  double lastRatio = 1.0e300;
  bool ok = false;
  for (int n = 0; n < 400; n++) {
    complex f = (a + (double)n) * (b + (double)n) / (double)(n + 1) * inv2irho;
    double ratio = std::abs(f);
    if (ratio > 1.0 && ratio > lastRatio) break;  // divergent tail reached
    lastRatio = ratio;
    term *= f;
    sum += term;
    dsum -= (double)(n + 1) * term;
    if (std::abs(term) < 1.0e-17 * std::abs(sum)) {
      ok = true;
      break;
    }
  }
  const complex theta = rho - eta * std::log(2.0 * rho) - 0.5 * M_PI * l + sigma;
  const complex e = std::exp(complex(0.0, 1.0) * theta);
  h = e * sum;
  // d/drho [e^{i theta} S] = e^{i theta} [i (1 - eta/rho) S + S'].
  if (dh) *dh = e * (complex(0.0, 1.0) * (1.0 - eta / rho) * sum + dsum / rho);
  return ok;
}

/// Local squared wave number (fm^-2); k2 = +k^2 (E > 0), -kappa^2 (E < 0) or 0.
double localQ2(double k2, double c, int l, double r) {
  return k2 - 2.0 * c / r - l * (l + 1.0) / (r * r);
}

/*!
 * A solution of the radial Coulomb equation r^2 u'' = (l(l+1) + 2 c r - k2 r^2) u
 * (c = eta k, k2 = k^2) near r0 as a Taylor series in x = r - r0:
 *
 *   r0^2 (n+2)(n+1) u_{n+2} = [l(l+1) + 2 c r0 - k2 r0^2 - n(n-1)] u_n
 *                            + (2 c - 2 k2 r0) u_{n-1} - k2 u_{n-2}
 *                            - 2 r0 (n+1) n u_{n+1}.
 *
 * Used only for |x| <= r0/3 (the singularity is at r = 0), where it converges
 * geometrically.  The coefficients are real, so real and imaginary parts (G
 * and F) never mix and a tiny F keeps its own relative accuracy.
 */
struct RadialTaylor {
  static const int kMax = 400;
  int l;
  double c, k2;
  double r0;
  double lo = 0.0;  // lower end of the interval this series is used on
  int n;
  complex u[kMax];

  void seed(double rSeed, double reach, complex value, complex deriv) {
    r0 = rSeed;
    const double A = l * (l + 1.0) + 2.0 * c * r0 - k2 * r0 * r0;
    const double B = 2.0 * c - 2.0 * k2 * r0;
    const double r02 = r0 * r0;
    u[0] = value;
    u[1] = deriv;
    const double scale = std::max(std::abs(value), std::abs(deriv) * reach);
    double xp = reach;
    n = 2;
    int small = 0;
    for (int m = 0; m + 2 < kMax; m++) {
      complex next = (A - m * (m - 1.0)) * u[m] - 2.0 * r0 * (m + 1.0) * m * u[m + 1];
      if (m >= 1) next += B * u[m - 1];
      if (m >= 2) next -= k2 * u[m - 2];
      next /= r02 * (m + 2.0) * (m + 1.0);
      u[m + 2] = next;
      n = m + 3;
      xp *= reach;
      small = (std::abs(next) * xp < 1.0e-18 * scale) ? small + 1 : 0;
      if (small >= 3) break;
    }
  }
  complex value(double r) const {
    const double x = r - r0;
    complex v = u[n - 1];
    for (int i = n - 2; i >= 0; i--) v = v * x + u[i];
    return v;
  }
  complex derivative(double r) const {
    const double x = r - r0;
    complex v = (double)(n - 1) * u[n - 1];
    for (int i = n - 2; i >= 1; i--) v = v * x + (double)i * u[i];
    return v;
  }
};

/// Panel width at r: <= kPanelPhase radians / e-folds, and <= r/3.
double panelWidth(double p, double kmin, double k2, double c, int l, double r) {
  double q = std::sqrt(std::fabs(localQ2(k2, c, l, r)));
  return std::min(kPanelPhase / (p + std::max(q, kmin)), r / 3.0);
}

/*!
 * Integrates ratio(r) j_l(p r) outwards from r0 to at most rEnd in GL panels,
 * stopping early once
 * restBound(r, |ratio(r)|) < kRelTol |integral|.  Sets rReached.
 */
template <class Ratio, class Bound>
complex realAxis(Ratio &ratio, Bound &restBound, int l, double p, double k2, double c, double r0,
                 double rEnd, double &rReached) {
  const GLRule &gl = glRule();
  complex acc(0.0, 0.0);
  double r = r0;
  const double kmin = std::sqrt(std::max(k2, 0.0));
  for (int panel = 0; panel < kMaxPanels && r < rEnd; panel++) {
    double hi = std::min(r + panelWidth(p, kmin, k2, c, l, r), rEnd);
    double half = 0.5 * (hi - r), mid = 0.5 * (hi + r);
    complex sum(0.0, 0.0), last(0.0, 0.0);
    for (int i = 0; i < kGLOrder; i++) {
      double x = mid + half * gl.x[i];
      complex g = ratio(x);
      sum += gl.w[i] * g * sphJl(l, p * x);
      if (i == gl.last) last = g;
    }
    acc += half * sum;
    r = hi;
    double lastAbs = std::abs(last);
    if (!(lastAbs == lastAbs)) break;  // NaN guard
    if (restBound(r, lastAbs) < kRelTol * std::abs(acc)) break;
  }
  rReached = r;
  return acc;
}

/// Integral of f(t) dt over t in [0, inf) for an integrand decaying like
/// exp(-rate t); panels grow geometrically up to 5/rate.
template <class F>
complex rayIntegral(F &f, double rate, double scale) {
  const GLRule &gl = glRule();
  complex acc(0.0, 0.0);
  double t = 0.0;
  double h = std::min(0.5 / rate, 0.25 * scale);
  for (int panel = 0; panel < 2000; panel++) {
    double half = 0.5 * h, mid = t + half;
    complex sum(0.0, 0.0);
    for (int i = 0; i < kGLOrder; i++) sum += gl.w[i] * f(mid + half * gl.x[i]);
    acc += half * sum;
    t += h;
    double tail = std::abs(f(t)) / rate;
    if (tail < 1.0e-16 * std::abs(acc) || tail == 0.0) break;
    h = std::min(h * 1.6, std::min(5.0 / rate, 0.5 * (scale + t)));
  }
  return acc;
}

/*!
 * WKB barrier integral Int_a^r |q| dr (Langer-corrected) outwards from a,
 * stopped at target or at the turning point.  Returns the radius reached and
 * sets s to the integral there.
 */
double barrierIntegral(double k2, double c, int l, double a, double target, double &s) {
  const double L = (l + 0.5) * (l + 0.5);
  auto kappa = [&](double r) { return std::sqrt(std::max(2.0 * c / r + L / (r * r) - k2, 0.0)); };
  s = 0.0;
  double r = a;
  for (int i = 0; i < 100000 && s < target; i++) {
    double kr = kappa(r);
    if (kr <= 0.0) break;
    double h = std::min(0.05 / kr, 0.05 * r);
    double km = kappa(r + 0.5 * h);
    if (km <= 0.0) break;  // at the turning point
    s += h * km;
    r += h;
  }
  return r;
}

/*!
 * Deep sub-barrier E > 0: Int_a^rc G(r)/G(a) j_l(p r) dr, G integrated inwards
 * from rc (see the header comment).  Values are rescaled on the way in; only
 * the ratio to G(a) is returned.
 */
double inwardBarrier(int l, double c, double k2, double p, double a, double rc) {
  const GLRule &gl = glRule();
  RadialTaylor t;
  t.l = l;
  t.c = c;
  t.k2 = k2;
  double r = rc;
  double y = 1.0;
  double dy = -std::sqrt(std::max(-localQ2(k2, c, l, r), 0.0));  // WKB log-derivative
  double acc = 0.0;
  const double kmin = std::sqrt(std::max(k2, 0.0));
  for (int panel = 0; panel < kMaxPanels && r > a; panel++) {
    double h = panelWidth(p, kmin, k2, c, l, r);
    h = std::min(h, panelWidth(p, kmin, k2, c, l, std::max(a, r - h)));  // steeper inside
    double lo = std::max(r - h, a);
    t.seed(r, r - lo, complex(y, 0.0), complex(dy, 0.0));
    double half = 0.5 * (r - lo), mid = 0.5 * (r + lo);
    double sum = 0.0;
    for (int i = 0; i < kGLOrder; i++) {
      double x = mid + half * gl.x[i];
      sum += gl.w[i] * t.value(x).real() * sphJl(l, p * x);
    }
    acc += half * sum;
    y = t.value(lo).real();
    dy = t.derivative(lo).real();
    r = lo;
    if (std::fabs(y) > 1.0e200) {  // keep y and the integral on the same scale
      y *= 1.0e-200;
      dy *= 1.0e-200;
      acc *= 1.0e-200;
    }
  }
  return acc / y;
}

complex coulombTermUncached(PPair *pair, int l, double E, double p) {
  const int zz = pair->GetZ(1) * pair->GetZ(2);
  if (zz == 0 || !(p > 0.0)) return complex(0.0, 0.0);
  const double mu = pair->GetRedMass() * uconv;  // MeV/c^2
  const double c = zz * fstruc * mu / hbarc;     // eta k, fm^-1 (energy independent)
  const double a = pair->GetChRad();
  complex integral(0.0, 0.0);

  if (E > 0.0) {
    const double k = std::sqrt(2.0 * mu * E) / hbarc;
    const double k2 = k * k;
    const double eta = c / k;
    double s;
    const double rc = barrierIntegral(k2, c, l, a, kBarrierSwitch, s);
    if (s >= kBarrierSwitch) return 2.0 * c * inwardBarrier(l, c, k2, p, a, rc);

    // Matching radius for the tail: the H+ series must converge there.
    double sigma;
    {
      gsl_sf_result lnr, arg;
      gsl_sf_lngamma_complex_e(1.0 + l, eta, &lnr, &arg);
      sigma = arg.val;
    }
    const double rhoTP = eta + std::sqrt(eta * eta + l * (l + 1.0));
    double rhoM = std::max(k * a, std::max(20.0, 1.2 * rhoTP));
    complex hTest;
    for (int it = 0; it < 60 && !hPlusAsymptotic(l, eta, sigma, complex(rhoM, 0.0), hTest); it++) rhoM *= 1.25;
    const double rM = rhoM / k;

    // O_l on [a, rM]: the Coulomb equation integrated INWARDS from the
    // asymptotic series at rM, in Taylor segments (RadialTaylor) of reach
    // X <= r/3 and <= kTaylorPhase/max(|q|, k).  Inwards every solution grows
    // through the barrier at the same rate as G_l, so the error stays at the
    // level of the starting value relative to |O_l| (the part that is lost is
    // the exponentially small F_l, which is irrelevant to O_l = G_l + i F_l).
    // No library Coulomb function is called: the result does not depend on
    // useGSL, and GSL's own G_l/F_l inaccuracies (~2% for l >= 2 at eta ~ 2,
    // rho ~ 1-4) cannot leak in.
    std::vector<RadialTaylor> segs;
    {
      complex o, dO;
      hPlusAsymptotic(l, eta, sigma, complex(rM, 0.0) * k, o, &dO);
      dO *= k;  // d/dr
      double rs = rM;
      while (rs > a && segs.size() < 100000) {
        double q = std::sqrt(std::fabs(localQ2(k2, c, l, rs)));
        double X = std::min(rs / 3.0, kTaylorPhase / std::max(q, k));
        q = std::sqrt(std::fabs(localQ2(k2, c, l, std::max(a, rs - X))));  // steeper inside
        X = std::min(X, kTaylorPhase / std::max(q, k));
        const double lo = std::max(a, rs - X);
        segs.emplace_back();
        RadialTaylor &t = segs.back();
        t.l = l;
        t.c = c;
        t.k2 = k2;
        t.seed(rs, rs - lo, o, dO);
        t.lo = lo;
        o = t.value(lo);
        dO = t.derivative(lo);
        rs = lo;
      }
    }
    complex oa;
    if (segs.empty())
      hPlusAsymptotic(l, eta, sigma, complex(k * a, 0.0), oa);  // a already asymptotic
    else
      oa = segs.back().value(a);
    size_t segIdx = segs.empty() ? 0 : segs.size() - 1;  // segments run from rM down to a
    auto ratio = [&](double r) -> complex {
      while (segIdx > 0 && r > segs[segIdx].r0) segIdx--;
      while (segIdx + 1 < segs.size() && r < segs[segIdx].lo) segIdx++;
      return segs[segIdx].value(r) / oa;
    };
    // |O| falls monotonically outwards and |x j_l(x)| <~ 1: bound on the rest
    // from r, through the barrier to r_tp and the 1/r oscillating tail.
    const double rTP = (c + std::sqrt(c * c + k2 * l * (l + 1.0))) / k2;
    const double dpk = std::max(std::fabs(p - k), 1.0e-12 * p);
    auto bound = [&](double r, double absRatio) {
      double R = std::max(rTP, r);
      return 10.0 * absRatio / p * (std::log(R / r) + 2.0 + 2.0 / (dpk * R));
    };


    double rReached = a;
    if (rM > a) integral = realAxis(ratio, bound, l, p, k2, c, a, rM, rReached);
    if (rReached >= rM) {
      auto h1Part = [&](double t) -> complex {
        complex r(rM, t);
        complex h;
        hPlusAsymptotic(l, eta, sigma, k * r, h);
        return h * sphHankel(l, p * r, +1) * complex(0.0, 1.0);
      };
      const int dir = (p > k) ? -1 : +1;  // h2 half decays like exp(i(k-p)r)
      auto h2Part = [&](double t) -> complex {
        complex r(rM, dir * t);
        complex h;
        hPlusAsymptotic(l, eta, sigma, k * r, h);
        return h * sphHankel(l, p * r, -1) * complex(0.0, (double)dir);
      };
      complex tail = 0.5 * (rayIntegral(h1Part, p + k, rM) + rayIntegral(h2Part, dpk, rM));
      integral += tail / oa;
    }
  } else if (E < 0.0) {
    const double kappa = std::sqrt(-2.0 * mu * E) / hbarc;
    WhitFunc wf(pair);
    double sa;
    const double ma = wf.Scaled(l, a, -E, sa);
    if (ma == 0.0) return complex(0.0, 0.0);
    auto ratio = [&](double r) -> complex {
      double s;
      double m = wf.Scaled(l, r, -E, s);
      double e = s - sa;
      return complex((e < -300.0) ? 0.0 : m / ma * std::pow(10.0, e), 0.0);
    };
    // W decays monotonically, at least as fast as exp(-|q| r) locally.
    auto bound = [&](double r, double absRatio) {
      double q = std::sqrt(-localQ2(-kappa * kappa, c, l, r));
      return 10.0 * absRatio * std::min(1.0, 1.0 / (p * r)) * 2.0 / q;
    };
    double rReached;
    integral = realAxis(ratio, bound, l, p, -kappa * kappa, c, a, 1.0e300, rReached);
  } else {
    // E = 0: the limit of either side, sqrt(r) K_{2l+1}(sqrt(8 c r)).
    const double xa = std::sqrt(8.0 * c * a);
    const double ka = gsl_sf_bessel_Kn_scaled(2 * l + 1, xa);
    auto ratio = [&](double r) -> complex {
      double x = std::sqrt(8.0 * c * r);
      double e = xa - x;
      if (e < -700.0) return complex(0.0, 0.0);
      return complex(std::sqrt(r / a) * gsl_sf_bessel_Kn_scaled(2 * l + 1, x) / ka * std::exp(e), 0.0);
    };
    auto bound = [&](double r, double absRatio) {
      double q = std::sqrt(-localQ2(0.0, c, l, r));
      return 10.0 * absRatio * std::min(1.0, 1.0 / (p * r)) * 2.0 / q;
    };
    double rReached;
    integral = realAxis(ratio, bound, l, p, 0.0, c, a, 1.0e300, rReached);
  }
  return 2.0 * c * integral;
}

}  // namespace

/*
 * EPoint::CalcEDependentValues asks for the term once per (J group, channel),
 * so the same (pair, l, E, p) comes back several times in a row for one point
 * (7x on average in tests/7Li_p_a).  A small per-thread memo of the last
 * results makes the repeats free.  The key holds everything the value depends
 * on, compared exactly.
 */
complex ThmCoulombTerm(PPair *pair, int l, double E, double p, bool useGSL) {
  (void)useGSL;  // O_l is built from its asymptotic expansion, see above
  struct Entry {
    const PPair *pair;
    int zz, l;
    double mu, a, E, p;
    complex value;
  };
  // 64: a spectator-momentum window (ps=...) asks for up to psNodes momenta
  // per (l, E), each a few times.
  static const int kMemo = 64;
  thread_local Entry memo[kMemo];
  thread_local int filled = 0, next = 0;
  const int zz = pair->GetZ(1) * pair->GetZ(2);
  const double mu = pair->GetRedMass(), a = pair->GetChRad();
  for (int i = 0; i < filled; i++) {
    const Entry &e = memo[i];
    if (e.pair == pair && e.zz == zz && e.l == l && e.mu == mu && e.a == a && e.E == E && e.p == p)
      return e.value;
  }
  complex v = coulombTermUncached(pair, l, E, p);
  memo[next] = Entry{pair, zz, l, mu, a, E, p, v};
  next = (next + 1) % kMemo;
  if (filled < kMemo) filled++;
  return v;
}
