#include "azure2/physics/HybridPotential.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <limits>
#include <string>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Whittaker.h"

namespace azure2 {
namespace physics {

namespace {

/// e^2 = alpha hbar c (MeV fm) of the Coulomb term V_C, the constants of the boundary values at r_m.
constexpr double kCoulombE2_MeV_fm = kFineStructure * kHbarC_MeV_fm;

}  // namespace

bool operator==(const FormFactor &a, const FormFactor &b) {
  return a.shape == b.shape && a.radius_fm == b.radius_fm && a.diffuseness_fm == b.diffuseness_fm &&
         a.width_fm == b.width_fm;
}

bool operator==(const SpinOrbitTerm &a, const SpinOrbitTerm &b) {
  return a.depth_MeV == b.depth_MeV && a.form == b.form && a.acts_on == b.acts_on;
}

bool operator==(const HybridPotential &a, const HybridPotential &b) {
  return a.coulomb_radius_fm == b.coulomb_radius_fm && a.central_depth_MeV == b.central_depth_MeV &&
         a.central == b.central && a.matching_radius_fm == b.matching_radius_fm && a.step_fm == b.step_fm &&
         a.spin_orbit == b.spin_orbit;
}

double TwoLDotI(int l, int two_j, int two_i) {
  if (l < 0 || two_i < 0 || two_j < std::abs(2 * l - two_i) || two_j > 2 * l + two_i || (two_j + two_i) % 2 != 0)
    throw std::invalid_argument("TwoLDotI: j = " + std::to_string(two_j) + "/2 is not in |l - i| ... l + i for l = " +
                                std::to_string(l) + ", i = " + std::to_string(two_i) + "/2");
  // j(j+1) - l(l+1) - i(i+1) with j = two_j/2, i = two_i/2: exact in double (multiples of 1/4)
  return (two_j * (two_j + 2) - two_i * (two_i + 2)) / 4.0 - l * (l + 1);
}

SpinOrbitTerm WoodsSaxonSpinOrbit(double depth_MeV, double radius_fm, double diffuseness_fm, SpinOrbitParticle on) {
  if (!(depth_MeV > 0.0) || !(radius_fm > 0.0) || !(diffuseness_fm > 0.0))
    throw std::invalid_argument("spin-orbit Woods-Saxon: V_so, R_so, a_so must be positive");
  SpinOrbitTerm t;
  t.depth_MeV = depth_MeV;
  t.form.shape = FormShape::WoodsSaxon;
  t.form.radius_fm = radius_fm;
  t.form.diffuseness_fm = diffuseness_fm;
  t.acts_on = on;
  return t;
}

HybridPotential WoodsSaxonPotential(double depth_MeV, double radius_fm, double diffuseness_fm) {
  if (depth_MeV <= 0.0 || radius_fm <= 0.0 || diffuseness_fm <= 0.0)
    throw std::invalid_argument("WoodsSaxon: V0, R, a must be positive");
  HybridPotential p;
  p.central_depth_MeV = depth_MeV;
  p.central.shape = FormShape::WoodsSaxon;
  p.central.radius_fm = radius_fm;
  p.central.diffuseness_fm = diffuseness_fm;
  return p;
}

HybridPotential GaussianPotential(double depth_MeV, double width_fm) {
  if (depth_MeV <= 0.0 || width_fm <= 0.0) throw std::invalid_argument("Gaussian: V0, r0 must be positive");
  HybridPotential p;
  p.central_depth_MeV = depth_MeV;
  p.central.shape = FormShape::Gaussian;
  p.central.width_fm = width_fm;
  return p;
}

double NuclearTerm(const HybridPotential &p, double r) {
  if (p.central.shape == FormShape::Gaussian) {
    double x = r / p.central.width_fm;
    return -p.central_depth_MeV * std::exp(-x * x);
  }
  return -p.central_depth_MeV / (1.0 + std::exp((r - p.central.radius_fm) / p.central.diffuseness_fm));
}

double SpinOrbitRadialTerm(const SpinOrbitTerm &t, double r) {
  // f = 1/(1 + e^x), x = (r - R)/a: df/dr = -(1/a) e^x / (1 + e^x)^2, written with e^-|x| (symmetric in x, no overflow)
  const double e = std::exp(-std::fabs((r - t.form.radius_fm) / t.form.diffuseness_fm));
  const double df_dr = -e / (t.form.diffuseness_fm * (1.0 + e) * (1.0 + e));
  return t.depth_MeV * kSpinOrbitLengthSquared_fm2 * df_dr / r;
}

double CoulombSphereTerm(int z1, int z2, double Rc, double r) {
  // Uniformly charged sphere inside R_C, point charge outside.
  if (r < Rc) {
    double x = r / Rc;
    return z1 * z2 * kCoulombE2_MeV_fm / (2.0 * Rc) * (3.0 - x * x);
  }
  return z1 * z2 * kCoulombE2_MeV_fm / r;
}

double RadialPotential(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double r) {
  const double Vc = CoulombSphereTerm(z1, z2, p.coulomb_radius_fm, r);
  double V_nuclear = NuclearTerm(p, r);
  double V_centrifugal = 0.0;
  if (r > 0.0) {
    double hbar2_over_2mu = (kHbarC_MeV_fm * kHbarC_MeV_fm) / (2.0 * mu_MeV);
    V_centrifugal = hbar2_over_2mu * l * (l + 1) / (r * r);
  }
  double V = Vc + V_nuclear + V_centrifugal;
  if (p.spin_orbit && two_l_dot_i != 0.0 && r > 0.0) V += SpinOrbitRadialTerm(*p.spin_orbit, r) * two_l_dot_i;
  return V;
}

std::vector<double> PotentialKey(const std::optional<HybridPotential> &p) {
  if (!p) return {0.0};
  std::vector<double> key{1.0, (double)p->central.shape, p->central_depth_MeV, p->central.radius_fm, p->central.diffuseness_fm,
                          p->central.width_fm, p->coulomb_radius_fm, p->matching_radius_fm, p->step_fm};
  if (p->spin_orbit)
    key.insert(key.end(), {p->spin_orbit->depth_MeV, p->spin_orbit->form.radius_fm, p->spin_orbit->form.diffuseness_fm,
                           p->spin_orbit->acts_on == SpinOrbitParticle::Heavy ? 1.0 : 0.0});
  return key;
}

RadialValue IntegrateInward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm, double E_MeV,
                            const std::function<std::complex<double>(double r_fm)> &outer) {
  // The grid: a step close to step_fm with a on node n_a; the two nodes beyond r_m carry the outer solution.
  const int n_a = std::max(2, static_cast<int>(std::lround(a_fm / p.step_fm)));
  const double dr = a_fm / n_a;
  const int n_top = static_cast<int>(std::ceil(p.matching_radius_fm / dr));
  const int npts = std::max(n_top + 3, n_a + 3);
  // Nodes below n_a - 2 are never read (the difference below reaches n_a - 2), so neither g nor u is made there.
  const int lowest = n_a - 2;

  std::vector<std::complex<double>> u(npts);
  std::vector<double> g(npts);
  for (int i = lowest; i < npts; ++i) {
    const double r = i * dr;
    if (i >= npts - 2) u[i] = outer(r);
    // g(r) = 2 mu / hbar^2 (E - V(r)), fm^-2
    g[i] = 2.0 * mu_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm) * (E_MeV - RadialPotential(p, z1, z2, mu_MeV, l, two_l_dot_i, r));
  }

  // Numerov, inward from r_m.
  RadialValue v;
  const double kRescale = 1e250;
  const double dr2 = dr * dr;
  for (int i = npts - 3; i >= lowest; --i) {
    double denom = 1.0 + (dr2 / 12.0) * g[i];
    if (std::abs(denom) > 1e-20) {
      u[i] = ((2.0 - (5.0 * dr2 / 6.0) * g[i + 1]) * u[i + 1] - (1.0 + (dr2 / 12.0) * g[i + 2]) * u[i + 2]) / denom;
    } else {
      u[i] = u[i + 1];
    }
    if (std::abs(u[i].real()) > kRescale || std::abs(u[i].imag()) > kRescale) {
      for (int j = i; j < npts; ++j) u[j] /= kRescale;
      v.log_scale += std::log(kRescale);
    }
  }

  v.u = u[n_a];
  v.du_dr = (u[n_a - 2] - 8.0 * u[n_a - 1] + 8.0 * u[n_a + 1] - u[n_a + 2]) / (12.0 * dr);
  return v;
}

namespace {

/// ln of the decaying point-Coulomb solution beyond r_m at r (fm) and binding B >= 0 (MeV): W_{-eta_b,l+1/2}(2 kappa r)
/// (unit ANC) for B > 0; sqrt(r) K_{2l+1}(sqrt(8 c r)) (charged) or r^-l (neutral) at B = 0.  NaN if out of range.
double OuterLogDecaying(int l, double c_fm, double mu_MeV, double B_MeV, double r) {
  if (B_MeV > 0.0) {
    const double kappa = std::sqrt(2.0 * mu_MeV * B_MeV) / kHbarC_MeV_fm;  // fm^-1
    const coulomb::WhittakerScaled w = coulomb::whittaker_scaled(l, c_fm / kappa, kappa * r);
    if (w.status != coulomb::Status::ok) return std::numeric_limits<double>::quiet_NaN();
    return std::log(w.mantissa) + (double)w.log10_scale * std::log(10.0);
  }
  if (!(c_fm > 0.0)) return -l * std::log(r);
  const double x = std::sqrt(8.0 * c_fm * r);
  return 0.5 * std::log(r) + std::log(bessel_k_scaled(2 * l + 1, x)) - x;
}

}  // namespace

DecayingValue DecayingSolution(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double r_fm,
                               double B_MeV) {
  const double c = z1 * z2 * kFineStructure * mu_MeV / kHbarC_MeV_fm;  // eta_b kappa (fm^-1)
  DecayingValue v;
  if (r_fm > p.matching_radius_fm) {  // the point-Coulomb solution (B = 0: its log-derivative by a centred difference)
    v.log_W = OuterLogDecaying(l, c, mu_MeV, B_MeV, r_fm);
    if (B_MeV > 0.0) {
      const double kappa = std::sqrt(2.0 * mu_MeV * B_MeV) / kHbarC_MeV_fm;
      v.rho_dlog = coulomb::whittaker_scaled(l, c / kappa, kappa * r_fm).rho_dlog;
    } else {
      const double h = 1e-4 * r_fm;
      v.rho_dlog = r_fm * (OuterLogDecaying(l, c, mu_MeV, B_MeV, r_fm + h) - OuterLogDecaying(l, c, mu_MeV, B_MeV, r_fm - h)) / (2.0 * h);
    }
    return v;
  }
  const double ref = OuterLogDecaying(l, c, mu_MeV, B_MeV, p.matching_radius_fm);
  auto outer = [&](double r) { return std::complex<double>(std::exp(OuterLogDecaying(l, c, mu_MeV, B_MeV, r) - ref), 0.0); };
  const RadialValue u = IntegrateInward(p, z1, z2, mu_MeV, l, two_l_dot_i, r_fm, -B_MeV, outer);
  v.sign = u.u.real() < 0.0 ? -1.0 : 1.0;
  v.log_W = std::log(std::fabs(u.u.real())) + u.log_scale + ref;
  v.rho_dlog = r_fm * u.du_dr.real() / u.u.real();
  return v;
}

double HybridBoundShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                        double B_MeV) {
  // As coulomb_bound_shift: closer to threshold than eta_b = 1000, the linear interpolation between S(-B_min) and S(0).
  const double kEtaMax = 1000.0;
  const double c = z1 * z2 * kFineStructure * mu_MeV / kHbarC_MeV_fm;
  const double bindingMin = (c > 0.0) ? std::pow(kHbarC_MeV_fm * c / kEtaMax, 2.0) / (2.0 * mu_MeV) : 0.0;
  auto S = [&](double B) { return DecayingSolution(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, B).rho_dlog; };
  if (B_MeV >= bindingMin && B_MeV > 0.0) return S(B_MeV);
  const double s0 = S(0.0);
  if (!(B_MeV > 0.0)) return s0;
  return s0 + (S(bindingMin) - s0) * (B_MeV / bindingMin);
}

RegularValue RegularSolution(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                             double E_MeV) {
  const int n_a = std::max(4, static_cast<int>(std::lround(a_fm / p.step_fm)));
  const double dr = a_fm / n_a, dr2 = dr * dr;
  const int npts = n_a + 3;  // nodes 0 .. n_a + 2: the difference at a reads two nodes beyond it
  const double k2 = 2.0 * mu_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  std::vector<double> u(npts, 0.0), g(npts, 0.0);
  for (int i = 1; i < npts; ++i) g[i] = k2 * (E_MeV - RadialPotential(p, z1, z2, mu_MeV, l, two_l_dot_i, i * dr));

  // The first two nodes from the series of the regular solution in a constant potential (the centrifugal term apart).
  const double w = -(g[1] + l * (l + 1) / dr2);  // (2 mu/hbar^2)(V - E) without the centrifugal term, fm^-2
  for (int i = 1; i <= 2; ++i) {
    const double r = i * dr, x = w * r * r;
    double term = 1.0, sum = 1.0;
    for (int k = 1; k <= 6; ++k) {
      term *= x / (2.0 * k * (2 * l + 2 * k + 1));
      sum += term;
    }
    u[i] = std::pow(r, l + 1) * sum;
  }
  const double kRescale = 1e250;
  for (int i = 3; i < npts; ++i) {
    u[i] = ((2.0 - (5.0 * dr2 / 6.0) * g[i - 1]) * u[i - 1] - (1.0 + (dr2 / 12.0) * g[i - 2]) * u[i - 2]) /
           (1.0 + (dr2 / 12.0) * g[i]);
    if (std::fabs(u[i]) > kRescale)
      for (int j = 0; j <= i; ++j) u[j] /= kRescale;
  }

  RegularValue v;
  const double ua = u[n_a];
  const double du = (u[n_a - 2] - 8.0 * u[n_a - 1] + 8.0 * u[n_a + 1] - u[n_a + 2]) / (12.0 * dr);
  v.log_derivative = a_fm * du / ua;
  // int_0^a u^2 dr: Simpson's rule on the nodes; an odd number of intervals takes the first three by the 3/8 rule.
  auto sq = [&](int i) { return u[i] * u[i]; };
  double sum = 0.0;
  int first = 0;
  if (n_a % 2) {
    sum += 3.0 * dr / 8.0 * (sq(0) + 3.0 * sq(1) + 3.0 * sq(2) + sq(3));
    first = 3;
  }
  for (int i = first; i + 2 <= n_a; i += 2) sum += dr / 3.0 * (sq(i) + 4.0 * sq(i + 1) + sq(i + 2));
  v.norm_fm = sum / (ua * ua);
  return v;
}

ChannelTriple InteriorChannelFunctions(const ChannelTriple &outer, double f, double b) {
  if (std::isinf(f)) return outer;
  const std::complex<double> L(outer.S, outer.P);
  const std::complex<double> d = f - L;  // f - L
  const std::complex<double> Lp = b + (L - b) * (f - b) / d;
  ChannelTriple t;
  t.S = Lp.real();
  t.P = Lp.imag();
  const double sign = (f - b) < 0.0 ? -1.0 : 1.0;
  t.exp_minus_i_phi = outer.exp_minus_i_phi * std::conj(d) / std::abs(d) * sign;
  return t;
}

}  // namespace physics
}  // namespace azure2
