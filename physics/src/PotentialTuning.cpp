#include "azure2/physics/PotentialTuning.h"

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Coulomb.h"

namespace azure2 {
namespace physics {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kPi = 3.14159265358979323846;

/// B_2k / (2k) for k = 1 ... 7 (Bernoulli numbers): the coefficients of the asymptotic series of the digamma function,
/// psi(z) ~ ln z - 1/(2z) - sum_k B_2k / (2k z^2k) (Abramowitz and Stegun 6.3.18).
constexpr std::array<double, 7> kBernoulliOver2k = {1.0 / 12.0,          -1.0 / 120.0, 1.0 / 252.0, -1.0 / 240.0,
                                                    1.0 / 132.0,         -691.0 / 32760.0, 1.0 / 12.0};

/// The Sommerfeld parameter and wave number at E > 0.
void Kinematics(int z1, int z2, double mu_MeV, double E_MeV, double &k, double &eta) {
  k = std::sqrt(2.0 * mu_MeV * E_MeV) / kHbarC_MeV_fm;
  eta = z1 * z2 * kFineStructure * mu_MeV / (kHbarC_MeV_fm * k);
}

/// The range of the nuclear term (fm) that sets the expansion's window: the Gaussian width, or R + 2a of a Woods-Saxon.
double NuclearRange(const HybridPotential &p) {
  if (p.central.shape == FormShape::Gaussian) return p.central.width_fm;
  return p.central.radius_fm + 2.0 * p.central.diffuseness_fm;
}

/// Least-squares polynomial of degree `degree` in x through (x_i, y_i), x scaled by `scale`; returns the coefficients of
/// x^0 and x^1 (unscaled).  Normal equations in the scaled variable (well conditioned for x/scale in [0, 1], degree <= 3).
bool FitPolynomial(const std::vector<double> &x, const std::vector<double> &y, int degree, double scale, double &c0,
                   double &c1) {
  const int n = degree + 1;
  std::vector<std::vector<double>> a(n, std::vector<double>(n + 1, 0.0));
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double t = x[i] / scale;
    std::vector<double> pw(2 * n, 1.0);
    for (int m = 1; m < 2 * n; ++m) pw[m] = pw[m - 1] * t;
    for (int r = 0; r < n; ++r) {
      for (int c = 0; c < n; ++c) a[r][c] += pw[r + c];
      a[r][n] += y[i] * pw[r];
    }
  }
  for (int c = 0; c < n; ++c) {  // Gauss-Jordan with partial pivoting
    int piv = c;
    for (int r = c + 1; r < n; ++r)
      if (std::fabs(a[r][c]) > std::fabs(a[piv][c])) piv = r;
    std::swap(a[c], a[piv]);
    if (a[c][c] == 0.0) return false;
    for (int r = 0; r < n; ++r) {
      if (r == c) continue;
      const double f = a[r][c] / a[c][c];
      for (int k = c; k <= n; ++k) a[r][k] -= f * a[c][k];
    }
  }
  c0 = a[0][n] / a[0][0];
  c1 = a[1][n] / a[1][1] / scale;
  return true;
}

}  // namespace

double CoulombH(double eta) {
  if (!(eta > 0.0)) return kNaN;
  if (eta >= 12.0) {
    // h = Re psi(i eta) - ln eta = -sum_k B_2k/(2k) Re (i eta)^-2k = sum_k (-1)^(k+1) B_2k/(2k) eta^-2k: every term > 0
    const double x = 1.0 / (eta * eta);
    double h = 0.0, pw = 1.0;
    for (std::size_t k = 0; k < kBernoulliOver2k.size(); ++k) {
      pw *= x;
      h += (k % 2 == 0 ? 1.0 : -1.0) * kBernoulliOver2k[k] * pw;
    }
    return h;
  }
  // psi(i eta) = psi(i eta + N) - sum_{n=0}^{N-1} 1/(n + i eta); Re 1/(n + i eta) = n/(n^2 + eta^2)
  const int N = 12;
  const std::complex<double> z(N, eta);
  const std::complex<double> z2inv = 1.0 / (z * z);
  std::complex<double> series = 0.0, pw = 1.0;
  for (double c : kBernoulliOver2k) {
    pw *= z2inv;
    series += c * pw;
  }
  const std::complex<double> psi = std::log(z) - 0.5 / z - series;
  double re = psi.real();
  for (int n = 1; n < N; ++n) re -= n / (n * n + eta * eta);
  return re - std::log(eta);
}

double CoulombC0Squared(double eta) {
  if (eta == 0.0) return 1.0;
  const double x = 2.0 * kPi * eta;
  return x / std::expm1(x);
}

double PotentialPhaseShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                           double E_MeV) {
  if (!(E_MeV > 0.0)) throw std::invalid_argument("PotentialPhaseShift: E must be above threshold");
  const double rm = p.matching_radius_fm;
  const InteriorValue v = IntegrateOutward(p, z1, z2, mu_MeV, l, two_l_dot_i, rm, E_MeV);
  double k = 0.0, eta = 0.0;
  Kinematics(z1, z2, mu_MeV, E_MeV, k, eta);
  const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * rm);
  if (w.status != coulomb::Status::ok) return kNaN;
  // u = c (cos d F + sin d G), u' = c k (cos d F' + sin d G') at r_m.  atan of the quotient, not atan2 reduced by pi:
  // a phase of 1e-12 rad (a charged pair far below the barrier) would keep only the ulp of pi (4e-16) of its digits.
  const double num = v.du_dr * w.F - k * v.u * w.Fp, den = k * v.u * w.Gp - v.du_dr * w.G;
  if (den == 0.0) return 0.5 * kPi;
  return std::atan(num / den);
}

BoundStateCount PotentialBoundStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i) {
  BoundStateCount out;
  const double rm = p.matching_radius_fm;
  const InteriorValue v = IntegrateOutward(p, z1, z2, mu_MeV, l, two_l_dot_i, rm, 0.0, &out.nodes);
  out.logDerivative = rm * v.du_dr / v.u;
  const double c = z1 * z2 * kFineStructure * mu_MeV / kHbarC_MeV_fm;  // eta k (fm^-1)
  out.exterior = zero_energy_shift(c, rm, l);
  out.bound = out.nodes + (out.logDerivative < out.exterior ? 1 : 0);
  return out;
}

EffectiveRange PotentialEffectiveRange(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i) {
  EffectiveRange out;
  const bool charged = z1 * z2 != 0;
  if (charged && l > 0) {
    out.why = "the Coulomb-modified effective-range expansion is built for l = 0 only (a charged pair's l = " +
              std::to_string(l) + ")";
    return out;
  }
  const double range = NuclearRange(p);
  if (!(range > 0.0)) {
    out.why = "the potential has no nuclear range";
    return out;
  }
  out.kHigh = 0.5 / range;
  out.kLow = 0.25 * out.kHigh;
  const double c = z1 * z2 * kFineStructure * mu_MeV / kHbarC_MeV_fm;  // eta k (fm^-1)
  if (charged) out.kLow = std::max(out.kLow, c / 30.0);                // eta <= 30 at the lowest energy
  if (!(out.kLow < out.kHigh / 1.5)) {
    out.why = "the window of the expansion is too narrow: eta exceeds 30 above k = 0.5/(1.5 R_N) (a heavy or highly "
              "charged pair)";
    return out;
  }
  const double hc2_2mu = kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * mu_MeV);
  out.eLow = hc2_2mu * out.kLow * out.kLow;
  out.eHigh = hc2_2mu * out.kHigh * out.kHigh;
  out.points = 24;
  std::vector<double> x, y;
  for (int i = 0; i < out.points; ++i) {
    const double k2 = out.kLow * out.kLow + i * (out.kHigh * out.kHigh - out.kLow * out.kLow) / (out.points - 1);
    const double E = hc2_2mu * k2, k = std::sqrt(k2);
    const double d = PotentialPhaseShift(p, z1, z2, mu_MeV, l, two_l_dot_i, E);
    if (!std::isfinite(d) || d == 0.0) {
      out.why = "the phase shift could not be computed at E = " + std::to_string(E) + " MeV";
      return out;
    }
    double K = std::pow(k, 2 * l + 1) / std::tan(d);
    if (charged) {
      const double eta = c / k;
      K = CoulombC0Squared(eta) * k / std::tan(d) + 2.0 * k * eta * CoulombH(eta);
    }
    x.push_back(k2);
    y.push_back(K);
  }
  const double scale = x.back();
  double c0 = 0.0, c1 = 0.0;
  if (!FitPolynomial(x, y, 3, scale, c0, c1)) {
    out.why = "the fit of the effective-range function failed";
    return out;
  }
  out.a = -1.0 / c0;
  out.r = 2.0 * c1;
  // the accuracy of the extrapolation: a quadratic on the same points, a cubic on the lower two thirds
  double q0 = 0.0, q1 = 0.0, l0 = 0.0, l1 = 0.0;
  const std::size_t lower = (2 * x.size()) / 3;
  const std::vector<double> xl(x.begin(), x.begin() + lower), yl(y.begin(), y.begin() + lower);
  if (FitPolynomial(x, y, 2, scale, q0, q1) && FitPolynomial(xl, yl, 3, xl.back(), l0, l1)) {
    out.aSpread = std::max(std::fabs(-1.0 / q0 - out.a), std::fabs(-1.0 / l0 - out.a));
    out.rSpread = std::max(std::fabs(2.0 * q1 - out.r), std::fabs(2.0 * l1 - out.r));
  }
  out.boundStates = PotentialBoundStates(p, z1, z2, mu_MeV, l, two_l_dot_i).bound;
  out.ok = true;
  return out;
}

DepthForLength DepthForScatteringLength(const HybridPotential &p, int z1, int z2, double mu_MeV, int l,
                                        double two_l_dot_i, double a_target) {
  DepthForLength out;
  const double v0 = p.central_depth_MeV;
  if (!(v0 > 0.0)) {
    out.why = "the potential's depth must be positive";
    return out;
  }
  HybridPotential q = p;
  auto at = [&](double v, EffectiveRange &e, int &bound) {
    q.central_depth_MeV = v;
    bound = PotentialBoundStates(q, z1, z2, mu_MeV, l, two_l_dot_i).bound;
    e = PotentialEffectiveRange(q, z1, z2, mu_MeV, l, two_l_dot_i);
    return e.ok;
  };
  EffectiveRange e0;
  int n0 = 0;
  if (!at(v0, e0, n0)) {
    out.why = e0.why;
    return out;
  }
  out.boundStates = n0;
  // a_l decreases with the depth within a branch: deeper when a_l is above the target
  const double dir = e0.a > a_target ? 1.0 : -1.0;
  const double step = 0.005 * v0;
  double vPrev = v0, fPrev = e0.a - a_target;
  for (int i = 1; i <= 600; ++i) {
    const double v = v0 + dir * i * step;
    if (!(v > 0.0)) {
      out.why = "no depth above 0 MeV gives a = " + std::to_string(a_target) + " with " + std::to_string(n0) +
                " bound state(s)";
      return out;
    }
    EffectiveRange e;
    int n = 0;
    if (!at(v, e, n)) {
      out.why = e.why;
      return out;
    }
    if (n != n0) {
      out.why = "the depth that gives a = " + std::to_string(a_target) + " is not in the branch of " + std::to_string(n0) +
                " bound state(s): at " + std::to_string(v) + " MeV the potential binds " + std::to_string(n);
      return out;
    }
    const double f = e.a - a_target;
    if ((f < 0.0) != (fPrev < 0.0)) {
      double lo = vPrev, hi = v, flo = fPrev;
      while (std::fabs(hi - lo) > 1e-10 * v0) {
        const double m = 0.5 * (lo + hi);
        EffectiveRange em;
        int nm = 0;
        if (!at(m, em, nm) || nm != n0) break;
        if (((em.a - a_target) < 0.0) == (flo < 0.0)) {
          lo = m;
          flo = em.a - a_target;
        } else {
          hi = m;
        }
      }
      const double root = 0.5 * (lo + hi);
      EffectiveRange er;
      int nr = 0;
      at(root, er, nr);
      out.found = er.ok && nr == n0;
      out.depth_MeV = root;
      out.a = er.a;
      out.r = er.r;
      q.central_depth_MeV = root;
      out.nodes = PotentialBoundStates(q, z1, z2, mu_MeV, l, two_l_dot_i).nodes;
      if (!out.found) out.why = "the bisection left the branch";
      return out;
    }
    vPrev = v;
    fPrev = f;
  }
  out.why = "no depth within " + std::to_string(dir > 0 ? 4.0 * v0 : 0.0) + " MeV gives a = " + std::to_string(a_target);
  return out;
}

}  // namespace physics
}  // namespace azure2
