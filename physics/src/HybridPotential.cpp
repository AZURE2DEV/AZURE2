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
  return a.boundary == b.boundary && a.subtract_nodes == b.subtract_nodes && a.capture == b.capture;
}

bool operator==(const PartialWave &a, const PartialWave &b) {
  return a.l == b.l && a.two_j == b.two_j && a.depth_MeV == b.depth_MeV && a.width_fm == b.width_fm &&
         a.background == b.background;
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
                           w.width_fm ? 1.0 : 0.0, w.width_fm.value_or(0.0), w.background ? 1.0 : 0.0});
    if (w.background) {
      key.insert(key.end(), {w.background->boundary, (double)w.background->subtract_nodes.size()});
      for (int n : w.background->subtract_nodes) key.push_back((double)n);
      // interior capture: a mark only when asked, so that the key of a background without it is the one before it existed
      if (w.background->capture) key.insert(key.end(), {-2.0, 1.0});
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
  if (w && w->width_fm) q.central.width_fm = *w->width_fm;
  return q;
}

HybridPotential BelowThresholdPotential(const HybridPotential &p, int l, int two_j) {
  if (BackgroundOf(p, l, two_j))
    throw std::runtime_error("the partial wave l = " + std::to_string(l) +
                             (two_j >= 0 ? ", 2J = " + std::to_string(two_j) : std::string()) +
                             " of a hybrid pair has an interior background and is evaluated at or below threshold (a level "
                             "below the pair's threshold, a closed channel or a bound final state of external capture): "
                             "refused, the background's decaying solution (ANC and capture final state) is not built "
                             "(docs/upgrade/potential_background.md section 6)");
  return ForPartialWave(p, l, two_j);
}

const InteriorBackground *BackgroundOf(const HybridPotential &p, int l, int two_j) {
  const PartialWave *w = FindWave(p, l, two_j);
  return (w && w->background) ? &*w->background : nullptr;
}

namespace {

/// The regular solution on the grid of IntegrateInward (a on node n_a, two nodes beyond it), divided by exp(log_scale).
struct OutwardGrid {
  int n_a = 0;
  double dr = 0.0;
  std::vector<std::complex<double>> u;
  double log_scale = 0.0;
};

OutwardGrid IntegrateOutwardGrid(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                 double a_fm, double E_MeV) {
  if (!(a_fm > 0.0) || !(p.step_fm > 0.0)) throw std::invalid_argument("IntegrateOutward: a and the step must be positive");
  OutwardGrid out;
  // The grid of IntegrateInward: a on node n_a; two nodes beyond it for the five-point difference.
  const int n_a = std::max(2, static_cast<int>(std::lround(a_fm / p.step_fm)));
  const double dr = a_fm / n_a;
  const int npts = n_a + 3;
  const double dr2 = dr * dr;
  const double k2 = 2.0 * mu_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  std::vector<double> g(npts, 0.0);
  for (int i = 1; i < npts; ++i) g[i] = k2 * (E_MeV - RadialPotential(p, z1, z2, mu_MeV, l, two_l_dot_i, i * dr));
  std::vector<std::complex<double>> &u = out.u;
  u.resize(npts);
  u[0] = 0.0;
  u[1] = std::pow(dr, l + 1);
  // (1 + dr2/12 g) u at the origin, the limit of g u there (u = c r^(l+1), c = u(h)/h^(l+1)): l = 1, g u -> -l(l+1) c
  // = -2 u(h)/h^2; l = 0 with a point charge (R_C = 0), g u -> -k2 Z1 Z2 e^2 c = -k2 Z1 Z2 e^2 u(h)/h; otherwise 0 (V
  // finite at the origin, or u vanishing faster than the centrifugal term grows).
  double origin = 0.0;
  if (l == 1) origin = -(dr2 / 12.0) * 2.0 * u[1].real() / dr2;
  if (l == 0 && !(p.coulomb_radius_fm > 0.0)) origin = -(dr2 / 12.0) * k2 * z1 * z2 * kCoulombE2_MeV_fm * u[1].real() / dr;
  const double kRescale = 1e250;
  for (int i = 1; i + 1 < npts; ++i) {
    if (i == 1)
      u[2] = NumerovNext(g[2], g[1], 0.0, u[1], std::complex<double>(origin, 0.0), dr2);
    else
      u[i + 1] = NumerovNext(g[i + 1], g[i], g[i - 1], u[i], u[i - 1], dr2);
    if (std::abs(u[i + 1].real()) > kRescale) {
      for (int j = 0; j <= i + 1; ++j) u[j] /= kRescale;
      out.log_scale += std::log(kRescale);
    }
  }
  out.n_a = n_a;
  out.dr = dr;
  return out;
}

/// Simpson's rule over the nodes 0 .. n of step h (the 3/8 rule on the last three intervals when n is odd).
double SimpsonNodes(const std::vector<double> &f, int n, double h) {
  if (n < 2) throw std::invalid_argument("SimpsonNodes: fewer than two intervals");
  double s = 0.0;
  int m = n;
  if (n % 2) {
    if (n < 3) throw std::invalid_argument("SimpsonNodes: one interval");
    m = n - 3;
    s += 3.0 * h / 8.0 * (f[m] + 3.0 * f[m + 1] + 3.0 * f[m + 2] + f[m + 3]);
  }
  if (m > 0) {
    double t = f[0] + f[m];
    for (int i = 1; i < m; ++i) t += (i % 2 ? 4.0 : 2.0) * f[i];
    s += h / 3.0 * t;
  }
  return s;
}

/// The regular solution's values at the nodes 0 .. n_a divided by its value at a.
std::vector<double> NormalisedAtA(const OutwardGrid &g) {
  std::vector<double> phi(g.n_a + 1);
  const double ua = g.u[g.n_a].real();
  for (int i = 0; i <= g.n_a; ++i) phi[i] = g.u[i].real() / ua;
  return phi;
}

}  // namespace

InteriorValue IntegrateOutward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                               double a_fm, double E_MeV, int *nodes) {
  const OutwardGrid grid = IntegrateOutwardGrid(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  const std::vector<std::complex<double>> &u = grid.u;
  const int n_a = grid.n_a;
  const double dr = grid.dr;
  InteriorValue v;
  v.log_scale = grid.log_scale;
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

std::vector<BoxState> SubtractedBoxStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                          double a_fm, const InteriorBackground &bg) {
  std::vector<BoxState> named;
  if (bg.subtract_nodes.empty()) return named;
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
    named.push_back(*s);
  }
  return named;
}

double BackgroundR(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                   double E_MeV, const InteriorBackground &bg) {
  double R = InteriorRFunction(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV, bg.boundary);
  if (bg.subtract_nodes.empty()) return R;
  for (const BoxState &s : SubtractedBoxStates(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, bg)) R -= s.gamma2_MeV / (s.E_MeV - E_MeV);
  return R;
}

bool HasInteriorCapture(const HybridPotential &p) {
  for (const PartialWave &w : p.waves)
    if (w.background && w.background->capture) return true;
  return false;
}

const InteriorBackground *InteriorCaptureOf(const HybridPotential &p, int l, int two_j) {
  const InteriorBackground *bg = BackgroundOf(p, l, two_j);
  return (bg && bg->capture) ? bg : nullptr;
}

double InteriorCaptureIntegral(const HybridPotential &pi, int li, double two_l_dot_i_i, double E_MeV,
                               const InteriorBackground &bg, const HybridPotential &pf, int lf, double two_l_dot_i_f,
                               double Ef_MeV, int z1, int z2, double mu_MeV, double a_fm,
                               const std::function<double(double r_fm)> &weight) {
  // The entrance and final waves on one grid: the step of IntegrateOutward is set by a and step_fm alone.
  if (pi.step_fm != pf.step_fm) throw std::invalid_argument("InteriorCaptureIntegral: the two potentials have different steps");
  const OutwardGrid gi = IntegrateOutwardGrid(pi, z1, z2, mu_MeV, li, two_l_dot_i_i, a_fm, E_MeV);
  const OutwardGrid gf = IntegrateOutwardGrid(pf, z1, z2, mu_MeV, lf, two_l_dot_i_f, a_fm, Ef_MeV);
  const int n = gi.n_a;
  const double h = gi.dr;
  // R_pot(E) as InteriorRFunction forms it (the same grid and difference): g_E(a) = R0, the background of the dressed waves
  const InteriorValue vi = IntegrateOutward(pi, z1, z2, mu_MeV, li, two_l_dot_i_i, a_fm, E_MeV);
  const double Rpot = vi.u / (a_fm * vi.du_dr - bg.boundary * vi.u);
  std::vector<double> g = NormalisedAtA(gi);
  for (double &x : g) x *= Rpot;
  for (const BoxState &s : SubtractedBoxStates(pi, z1, z2, mu_MeV, li, two_l_dot_i_i, a_fm, bg)) {
    const std::vector<double> phi = NormalisedAtA(IntegrateOutwardGrid(pi, z1, z2, mu_MeV, li, two_l_dot_i_i, a_fm, s.E_MeV));
    const double c = s.gamma2_MeV / (s.E_MeV - E_MeV);
    for (int i = 0; i <= n; ++i) g[i] -= c * phi[i];
  }
  const std::vector<double> phif = NormalisedAtA(gf);
  std::vector<double> f(n + 1);
  for (int i = 0; i <= n; ++i) f[i] = i == 0 ? 0.0 : phif[i] * weight(i * h) * g[i];
  return SimpsonNodes(f, n, h);
}

double BoundMatching(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                     double E_MeV) {
  const InteriorValue u = IntegrateOutward(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  const DecayingValue w = DecayingSolution(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, -E_MeV);
  return w.sign * (w.rho_dlog * u.u - a_fm * u.du_dr);
}

PotentialEigenvalue NearestEigenvalue(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                      double a_fm, double E_MeV, double window_MeV) {
  PotentialEigenvalue best;
  if (!(E_MeV < 0.0) || !(window_MeV > 0.0)) throw std::invalid_argument("NearestEigenvalue: E < 0 and a positive window");
  const double lo = E_MeV - window_MeV, hi = std::min(E_MeV + window_MeV, -1e-6);
  const int n = 100;
  const double step = (hi - lo) / n;
  auto f = [&](double e) { return BoundMatching(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, e); };
  double e0 = lo, f0 = f(e0);
  for (int k = 1; k <= n; ++k) {
    const double e1 = lo + k * step, f1 = f(e1);
    if ((f0 < 0.0) != (f1 < 0.0)) {
      double a = e0, b = e1, fa = f0;
      while (b - a > 1e-10) {
        const double m = 0.5 * (a + b), fm = f(m);
        if ((fm < 0.0) == (fa < 0.0)) {
          a = m;
          fa = fm;
        } else {
          b = m;
        }
      }
      const double root = 0.5 * (a + b);
      if (!best.found || std::fabs(root - E_MeV) < std::fabs(best.E_MeV - E_MeV)) best = {true, root};
    }
    e0 = e1;
    f0 = f1;
  }
  return best;
}

DepthAtEnergy DepthBindingAt(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                             double E_MeV) {
  DepthAtEnergy best;
  const double v0 = p.central_depth_MeV;
  auto f = [&](double v) {
    HybridPotential q = p;
    q.central_depth_MeV = v;
    return BoundMatching(q, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  };
  const int n = 120;  // 0.25 % steps over +-30 %
  double x0 = 0.7 * v0, f0 = f(x0);
  for (int k = 1; k <= 2 * n; ++k) {
    const double x1 = 0.7 * v0 + k * (0.6 * v0 / (2 * n)), f1 = f(x1);
    if ((f0 < 0.0) != (f1 < 0.0)) {
      double a = x0, b = x1, fa = f0;
      while (b - a > 1e-9 * v0) {
        const double m = 0.5 * (a + b), fm = f(m);
        if ((fm < 0.0) == (fa < 0.0)) {
          a = m;
          fa = fm;
        } else {
          b = m;
        }
      }
      const double root = 0.5 * (a + b);
      if (!best.found || std::fabs(root - v0) < std::fabs(best.depth_MeV - v0)) best = {true, root};
    }
    x0 = x1;
    f0 = f1;
  }
  return best;
}

double InteriorNormIntegral(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                            double E_MeV) {
  const OutwardGrid grid = IntegrateOutwardGrid(p, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV);
  std::vector<double> phi = NormalisedAtA(grid);
  for (double &x : phi) x *= x;
  return SimpsonNodes(phi, grid.n_a, grid.dr);
}

double DressShift(double S, double R0, double B) { return B + (S - B) / (1.0 - R0 * (S - B)); }

double ChannelBoundShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, int two_j, double two_l_dot_i,
                         double a_fm, double B_MeV) {
  const HybridPotential q = ForPartialWave(p, l, two_j);
  const double S = HybridBoundShift(q, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, B_MeV);
  const InteriorBackground *bg = BackgroundOf(p, l, two_j);
  if (!bg) return S;
  return DressShift(S, BackgroundR(q, z1, z2, mu_MeV, l, two_l_dot_i, a_fm, -B_MeV, *bg), bg->boundary);
}

double BackgroundWFactor(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, int two_j, double two_l_dot_i,
                         double a_fm, double E_MeV, double S_dressed) {
  const InteriorBackground *bg = BackgroundOf(p, l, two_j);
  if (!bg) return 1.0;
  const double R0 = BackgroundR(ForPartialWave(p, l, two_j), z1, z2, mu_MeV, l, two_l_dot_i, a_fm, E_MeV, *bg);
  return 1.0 / (1.0 + R0 * (S_dressed - bg->boundary));
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
