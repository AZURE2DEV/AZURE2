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

/// One Numerov step: the solution at the node `next` from the two before it (in the direction of integration), with
/// g = 2 mu (E - V)/hbar^2 at the three nodes and dr2 the squared step.  The one recurrence of both directions.
inline std::complex<double> NumerovNext(double g_next, double g_cur, double g_prev, std::complex<double> u_cur,
                                        std::complex<double> u_prev, double dr2) {
  double denom = 1.0 + (dr2 / 12.0) * g_next;
  if (std::abs(denom) > 1e-20) return ((2.0 - (5.0 * dr2 / 6.0) * g_cur) * u_cur - (1.0 + (dr2 / 12.0) * g_prev) * u_prev) / denom;
  return u_cur;
}

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
         a.spin_orbit == b.spin_orbit && a.waves == b.waves;
}

bool operator==(const InteriorBackground &a, const InteriorBackground &b) {
  return a.boundary == b.boundary && a.subtract_nodes == b.subtract_nodes;
}

bool operator==(const PartialWave &a, const PartialWave &b) {
  return a.l == b.l && a.two_j == b.two_j && a.depth_MeV == b.depth_MeV && a.background == b.background;
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
  // the partial waves (none: the key of a potential without them, as before they existed)
  for (const PartialWave &w : p->waves) {
    key.insert(key.end(), {-1.0, (double)w.l, (double)w.two_j, w.depth_MeV ? 1.0 : 0.0, w.depth_MeV.value_or(0.0),
                           w.background ? 1.0 : 0.0});
    if (w.background) {
      key.insert(key.end(), {w.background->boundary, (double)w.background->subtract_nodes.size()});
      for (int n : w.background->subtract_nodes) key.push_back((double)n);
    }
  }
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
    u[i] = NumerovNext(g[i], g[i + 1], g[i + 2], u[i + 1], u[i + 2], dr2);
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

// ---- Partial waves and the interior background (docs/upgrade/potential_background.md) ----

const PartialWave *FindWave(const HybridPotential &p, int l, int two_j) {
  const PartialWave *any = nullptr;
  for (const PartialWave &w : p.waves) {
    if (w.l != l) continue;
    if (w.two_j == two_j && two_j != kNoSpinOrbit) return &w;
    if (w.two_j == kNoSpinOrbit) any = &w;
  }
  return any;
}

bool DependsOnJ(const HybridPotential &p) {
  if (p.spin_orbit) return true;
  for (const PartialWave &w : p.waves)
    if (w.two_j != kNoSpinOrbit) return true;
  return false;
}

HybridPotential ForPartialWave(const HybridPotential &p, int l, int two_j) {
  HybridPotential q = p;
  const PartialWave *w = FindWave(p, l, two_j);
  if (w && w->depth_MeV) q.central_depth_MeV = *w->depth_MeV;
  return q;
}

const InteriorBackground *BackgroundOf(const HybridPotential &p, int l, int two_j) {
  const PartialWave *w = FindWave(p, l, two_j);
  return (w && w->background) ? &*w->background : nullptr;
}

InteriorValue IntegrateOutward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                               double a_fm, double E_MeV, int *nodes) {
  if (!(a_fm > 0.0) || !(p.step_fm > 0.0)) throw std::invalid_argument("IntegrateOutward: a and the step must be positive");
  // The grid of IntegrateInward: a on node n_a; two nodes beyond it for the five-point difference.
  const int n_a = std::max(2, static_cast<int>(std::lround(a_fm / p.step_fm)));
  const double dr = a_fm / n_a;
  const int npts = n_a + 3;
  const double dr2 = dr * dr;
  const double k2 = 2.0 * mu_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  std::vector<double> g(npts, 0.0);
  for (int i = 1; i < npts; ++i) g[i] = k2 * (E_MeV - RadialPotential(p, z1, z2, mu_MeV, l, two_l_dot_i, i * dr));
  std::vector<std::complex<double>> u(npts);
  u[0] = 0.0;
  u[1] = std::pow(dr, l + 1);
  // (1 + dr2/12 g) u at the origin, the limit of g u there (u = c r^(l+1), c = u(h)/h^(l+1)): l = 1, g u -> -l(l+1) c
  // = -2 u(h)/h^2; l = 0 with a point charge (R_C = 0), g u -> -k2 Z1 Z2 e^2 c = -k2 Z1 Z2 e^2 u(h)/h; otherwise 0 (V
  // finite at the origin, or u vanishing faster than the centrifugal term grows).
  double origin = 0.0;
  if (l == 1) origin = -(dr2 / 12.0) * 2.0 * u[1].real() / dr2;
  if (l == 0 && !(p.coulomb_radius_fm > 0.0)) origin = -(dr2 / 12.0) * k2 * z1 * z2 * kCoulombE2_MeV_fm * u[1].real() / dr;
  InteriorValue v;
  const double kRescale = 1e250;
  for (int i = 1; i + 1 < npts; ++i) {
    if (i == 1)
      u[2] = NumerovNext(g[2], g[1], 0.0, u[1], std::complex<double>(origin, 0.0), dr2);
    else
      u[i + 1] = NumerovNext(g[i + 1], g[i], g[i - 1], u[i], u[i - 1], dr2);
    if (std::abs(u[i + 1].real()) > kRescale) {
      for (int j = 0; j <= i + 1; ++j) u[j] /= kRescale;
      v.log_scale += std::log(kRescale);
    }
  }
  v.u = u[n_a].real();
  v.du_dr = ((u[n_a - 2] - 8.0 * u[n_a - 1] + 8.0 * u[n_a + 1] - u[n_a + 2]) / (12.0 * dr)).real();
  if (nodes) {
    int n = 0;
    for (int i = 2; i < n_a; ++i)
      if ((u[i].real() < 0.0) != (u[i - 1].real() < 0.0) && u[i].real() != 0.0) ++n;
    *nodes = n;
  }
  return v;
}

double InteriorLogDerivative(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                             double a_fm, double E_MeV) {
  const InteriorValue v = IntegrateOutward(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  return a_fm * v.du_dr / v.u;
}

double InteriorRFunction(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                         double E_MeV, double boundary) {
  const InteriorValue v = IntegrateOutward(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  return v.u / (a_fm * v.du_dr - boundary * v.u);
}

namespace {

/// f(E) = a u'(a) - B u(a) of the regular solution, divided by exp(log_scale) > 0: its sign and zeros are f's.
double BoxCondition(const HybridPotential &p, int z1, int z2, double mu, int l, double ldoti, double a, double B,
                    double E, int *nodes = nullptr) {
  const InteriorValue v = IntegrateOutward(p, z1, z2, mu, l, ldoti, a, E, nodes);
  return a * v.du_dr - B * v.u;
}

}  // namespace

std::vector<BoxState> InteriorBoxStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                        double a_fm, double boundary, double E_lo_MeV, double E_hi_MeV, double dE_MeV) {
  if (!(dE_MeV > 0.0)) throw std::invalid_argument("InteriorBoxStates: the scan step must be positive");
  // No eigenvalue lies below the minimum of V on (0, a].
  const int n_a = std::max(2, static_cast<int>(std::lround(a_fm / p.step_fm)));
  double vmin = std::numeric_limits<double>::infinity();
  for (int i = 1; i <= n_a; ++i) vmin = std::min(vmin, RadialPotential(p, z1, z2, mu_MeV, l, two_l_dot_i, i * a_fm / n_a));
  double E = std::max(E_lo_MeV, vmin);
  std::vector<BoxState> states;
  double f0 = BoxCondition(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, boundary, E);
  while (E < E_hi_MeV) {
    const double E1 = std::min(E + dE_MeV, E_hi_MeV);
    const double f1 = BoxCondition(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, boundary, E1);
    if ((f0 < 0.0) != (f1 < 0.0)) {
      double lo = E, hi = E1, flo = f0;
      while (hi - lo > 1e-11) {
        const double mid = 0.5 * (lo + hi);
        const double fm = BoxCondition(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, boundary, mid);
        if ((fm < 0.0) == (flo < 0.0)) {
          lo = mid;
          flo = fm;
        } else {
          hi = mid;
        }
      }
      BoxState s;
      s.E_MeV = 0.5 * (lo + hi);
      BoxCondition(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, boundary, s.E_MeV, &s.nodes);
      // gamma^2 = -1/(dl/dE): Richardson on the central difference of l(E) = a u'/u
      auto ell = [&](double e) { return InteriorLogDerivative(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, e); };
      const double h = 1e-3;
      const double d1 = (ell(s.E_MeV + h) - ell(s.E_MeV - h)) / (2.0 * h);
      const double d2 = (ell(s.E_MeV + 0.5 * h) - ell(s.E_MeV - 0.5 * h)) / h;
      s.gamma2_MeV = -1.0 / ((4.0 * d2 - d1) / 3.0);
      states.push_back(s);
    }
    E = E1;
    f0 = f1;
  }
  return states;
}

double BackgroundR(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                   double E_MeV, const InteriorBackground &bg) {
  double R = InteriorRFunction(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV, bg.boundary);
  if (bg.subtract_nodes.empty()) return R;
  const int highest = *std::max_element(bg.subtract_nodes.begin(), bg.subtract_nodes.end());
  // The box states up to the highest named: scan upward in windows until it is found.
  std::vector<BoxState> states;
  for (double hi = 0.0; hi <= 200.0; hi += 25.0) {
    states = InteriorBoxStates(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, bg.boundary, -1e300, hi);
    if (!states.empty() && states.back().nodes >= highest) break;
  }
  for (int n : bg.subtract_nodes) {
    const BoxState *s = nullptr;
    for (const BoxState &b : states)
      if (b.nodes == n) s = &b;
    if (!s) throw std::runtime_error("BackgroundR: no box state with " + std::to_string(n) + " nodes below 200 MeV");
    R -= s->gamma2_MeV / (s->E_MeV - E_MeV);
  }
  return R;
}

WavesAtRadius DressWaves(const WavesAtRadius &w, double rho, double R0, double B) {
  WavesAtRadius d;
  const double jF = rho * w.Fp - B * w.F, jG = rho * w.Gp - B * w.G;  // rho u' - B u
  d.F = w.F - R0 * jF;
  d.G = w.G - R0 * jG;
  d.Fp = (rho * w.Fp - B * R0 * jF) / rho;
  d.Gp = (rho * w.Gp - B * R0 * jG) / rho;
  return d;
}

}  // namespace physics
}  // namespace azure2
