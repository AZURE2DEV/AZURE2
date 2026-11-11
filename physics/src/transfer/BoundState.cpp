#include "azure2/physics/transfer/BoundState.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/coulomb/Whittaker.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;

/// 2 mu / hbar^2 (MeV^-1 fm^-2) of a pair.
double TwoMuOverHbar2(const ChannelKinematics &pair) {
  return 2.0 * pair.mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
}

/// ln W_{-eta, l+1/2}(2 kappa r) (unit ANC).
double LogWhittaker(int l, double eta, double kappa, double r) {
  const coulomb::WhittakerScaled w = coulomb::whittaker_scaled(l, eta, kappa * r);
  if (w.status != coulomb::Status::ok)
    throw std::runtime_error("bound state: Whittaker function failed at r = " + std::to_string(r) + " fm");
  return std::log(w.mantissa) + static_cast<double>(w.log10_scale) * std::log(10.0);
}

void CheckRequest(const BoundStateRequest &q) {
  if (!(q.binding_MeV > 0.0)) throw std::invalid_argument("bound state: B must be positive");
  if (!(q.pair.mu_u > 0.0) || q.pair.relativistic) throw std::invalid_argument("bound state: mu > 0, non-relativistic");
  if (q.l < 0 || q.nodes < 0) throw std::invalid_argument("bound state: l, nodes >= 0");
  if (!(q.step_fm > 0.0) || !(q.r_max_fm > 20.0 * q.step_fm)) throw std::invalid_argument("bound state: bad grid");
  if (q.two_s == 0 && q.two_j != 2 * q.l) throw std::invalid_argument("bound state: j = l without spin");
}

/// Everything of Eq. (3.1) that does not depend on V0, on the grid.
struct Problem {
  int n = 0, ic = 0;  // nodes; matching node
  double h = 0.0, k2mu = 0.0, energy = 0.0;
  std::vector<double> form, fixed, centrifugal;
};

/// q_i = (2 mu/hbar^2)(E - V) - l(l+1)/r^2 at depth V0, and q_0 of the series start.
void Coefficients(const Problem &p, double v0, std::vector<cplx> &q, cplx &q0) {
  q.assign(p.n, 0.0);
  for (int i = 1; i < p.n; ++i) q[i] = p.k2mu * (p.energy + v0 * p.form[i] - p.fixed[i]) - p.centrifugal[i];
  q0 = p.k2mu * (p.energy + v0 * p.form[1] - p.fixed[1]);
}

int SignChanges(const std::vector<cplx> &u, int last) {
  int changes = 0;
  for (int i = 2; i <= last; ++i)
    if ((u[i].real() < 0.0) != (u[i - 1].real() < 0.0) && u[i].real() != 0.0) ++changes;
  return changes;
}

}  // namespace

BoundState SolveBoundStateOfForm(const std::function<double(double)> &form, const std::function<double(double)> &fixed,
                                 const BoundStateRequest &q) {
  CheckRequest(q);
  if (!(q.matching_radius_fm > 0.0) || q.matching_radius_fm > q.r_max_fm - 5.0 * q.step_fm)
    throw std::invalid_argument("bound state: the matching radius must be inside the grid");
  Problem p;
  p.h = q.step_fm;
  p.n = static_cast<int>(std::lround(q.r_max_fm / p.h)) + 1;
  p.ic = static_cast<int>(std::lround(q.matching_radius_fm / p.h));
  p.k2mu = TwoMuOverHbar2(q.pair);
  p.energy = -q.binding_MeV;
  p.form.assign(p.n, 0.0);
  p.fixed.assign(p.n, 0.0);
  p.centrifugal.assign(p.n, 0.0);
  for (int i = 1; i < p.n; ++i) {
    const double r = i * p.h;
    p.form[i] = form(r);
    p.fixed[i] = fixed(r);
    p.centrifugal[i] = q.l * (q.l + 1) / (r * r);
  }

  BoundState s;
  s.l = q.l;
  s.two_j = q.two_j;
  s.two_s = q.two_s;
  s.binding_MeV = q.binding_MeV;
  s.kappa_fm = rho_bound(q.pair, 1.0, q.binding_MeV);
  s.eta = eta_bound(q.pair, q.binding_MeV);
  s.step_fm = p.h;

  // The inward start: the Whittaker function at the last two nodes, scaled to 1 at the last.
  const double logw_last = LogWhittaker(q.l, s.eta, s.kappa_fm, (p.n - 1) * p.h);
  const double logw_before = LogWhittaker(q.l, s.eta, s.kappa_fm, (p.n - 2) * p.h);
  const cplx in_before = std::exp(logw_before - logw_last);

  std::vector<cplx> qv;
  cplx q0;
  auto count = [&](double v0) {
    Coefficients(p, v0, qv, q0);
    return SignChanges(NumerovOutward(qv, p.h, q.l, q0), p.n - 1);
  };
  // The mismatch: the Casoratian of z = (1 + h^2 q/12) u of the two solutions at nodes ic, ic + 1, over their norms
  // (continuous in V0, zero exactly where -B is an eigenvalue).
  std::vector<cplx> out, in;
  auto mismatch = [&](double v0) {
    Coefficients(p, v0, qv, q0);
    out = NumerovOutward(qv, p.h, q.l, q0);
    in = NumerovInward(qv, p.h, in_before, 1.0, p.ic);
    const double h2 = p.h * p.h / 12.0;
    const double zo0 = ((1.0 + h2 * qv[p.ic]) * out[p.ic]).real(), zo1 = ((1.0 + h2 * qv[p.ic + 1]) * out[p.ic + 1]).real();
    const double zi0 = ((1.0 + h2 * qv[p.ic]) * in[p.ic]).real(), zi1 = ((1.0 + h2 * qv[p.ic + 1]) * in[p.ic + 1]).real();
    return (zo0 * zi1 - zo1 * zi0) / (std::hypot(zo0, zo1) * std::hypot(zi0, zi1));
  };

  // Bracket the depth by the Sturm count, then make sure the mismatch changes sign in the bracket (the zero of a
  // weakly bound solution may lie beyond r_max just above the eigen-depth).
  double hi = 10.0;
  while (count(hi) <= q.nodes) {
    hi *= 2.0;
    if (hi > 1e5) throw std::runtime_error("bound state: no depth below 1e5 MeV binds the state");
  }
  double lo = 0.0;
  if (count(lo) > q.nodes) throw std::runtime_error("bound state: bound without a well (check the inputs)");
  for (int it = 0; it < 200 && hi - lo > 1e-6 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    (count(mid) <= q.nodes ? lo : hi) = mid;
  }
  const double f_hi = mismatch(hi);
  double f_lo = mismatch(lo);
  for (int k = 1; f_lo * f_hi > 0.0; ++k) {
    lo = hi - std::pow(2.0, k) * (hi - lo);
    if (lo < 0.0 || k > 60 || count(lo) < q.nodes) throw std::runtime_error("bound state: depth not bracketed");
    f_lo = mismatch(lo);
  }
  // Illinois (modified regula falsi) on the mismatch.
  double a = lo, b = hi, fa = f_lo, fb = f_hi;
  int side = 0;
  double v0 = 0.5 * (a + b);
  for (int it = 0; it < 200; ++it) {
    v0 = (a * fb - b * fa) / (fb - fa);
    const double fv = mismatch(v0);
    if (fv == 0.0 || std::fabs(b - a) < 1e-14 * std::fabs(v0)) break;
    if (fv * fb > 0.0) {
      b = v0, fb = fv;
      if (side == -1) fa *= 0.5;
      side = -1;
    } else {
      a = v0, fa = fv;
      if (side == +1) fb *= 0.5;
      side = +1;
    }
  }
  mismatch(v0);  // the solutions at the depth found
  s.depth_MeV = v0;

  // Join: outward below the matching node, the inward solution scaled to it above (least squares on two nodes).
  const double scale = (out[p.ic].real() * in[p.ic].real() + out[p.ic + 1].real() * in[p.ic + 1].real()) /
                       (in[p.ic].real() * in[p.ic].real() + in[p.ic + 1].real() * in[p.ic + 1].real());
  s.u.assign(p.n, 0.0);
  for (int i = 0; i < p.n; ++i) s.u[i] = i <= p.ic ? out[i].real() : scale * in[i].real();
  std::vector<double> u2(p.n), r2u2(p.n);
  for (int i = 0; i < p.n; ++i) u2[i] = s.u[i] * s.u[i];
  double norm = std::sqrt(Simpson(u2, p.h));
  if (s.u[1] < 0.0) norm = -norm;  // positive next to the origin
  for (double &x : s.u) x /= norm;
  for (int i = 0; i < p.n; ++i) r2u2[i] = (i * p.h) * (i * p.h) * s.u[i] * s.u[i];
  s.rms_radius_fm = std::sqrt(Simpson(r2u2, p.h));
  // u = scale in / norm and in(r_{n-1}) = 1 = W(r_{n-1}) exp(-logw_last): C = scale exp(-logw_last) / norm.
  s.anc_fm_half = scale / norm * std::exp(-logw_last);

  s.potential_MeV.assign(p.n, 0.0);
  for (int i = 1; i < p.n; ++i) s.potential_MeV[i] = -v0 * p.form[i] + p.fixed[i];
  s.potential_MeV[0] = s.potential_MeV[1];  // multiplies u(0) = 0; the spin-orbit form is singular at r = 0
  s.coulomb_MeV.assign(p.n, 0.0);
  return s;
}

BoundState SolveBoundState(const OpticalPotential &well, const BoundStateRequest &request) {
  if (well.imaginary_volume.depth_MeV != 0.0 || well.imaginary_surface.depth_MeV != 0.0 ||
      well.imaginary_spin_orbit.depth_MeV != 0.0)
    throw std::invalid_argument("bound state: the well must be real");
  if (well.real_volume.depth_MeV == 0.0) throw std::invalid_argument("bound state: the well needs a real volume term");
  OpticalPotential unit = NoPotential();
  unit.convention = well.convention;
  unit.target_mass_number = well.target_mass_number;
  unit.projectile_mass_number = well.projectile_mass_number;
  unit.real_volume = well.real_volume;
  unit.real_volume.depth_MeV = -1.0;  // volume_MeV = +f
  OpticalPotential rest = well;
  rest.real_volume.depth_MeV = 0.0;
  const int z1z2 = request.pair.z1 * request.pair.z2;
  const double two_l_dot_s = request.two_s > 0 ? TwoLDotI(request.l, request.two_j, request.two_s) : 0.0;
  BoundStateRequest q = request;
  if (q.matching_radius_fm == 0.0) q.matching_radius_fm = RadiusOf(well, well.real_volume.r0_fm);
  auto form = [&](double r) { return TermsAt(unit, 0, r).volume_MeV.real(); };
  auto fixed = [&](double r) { return PotentialAt(rest, z1z2, two_l_dot_s, r).real(); };
  BoundState s = SolveBoundStateOfForm(form, fixed, q);
  for (std::size_t i = 0; i < s.coulomb_MeV.size(); ++i)
    s.coulomb_MeV[i] = TermsAt(well, z1z2, std::max(i, std::size_t(1)) * s.step_fm).coulomb_MeV;
  return s;
}

double ZeroRangeConstant(const BoundState &s) {
  if (s.l != 0) throw std::invalid_argument("ZeroRangeConstant: an l = 0 state");
  std::vector<double> f(s.u.size());
  for (std::size_t i = 0; i < f.size(); ++i) f[i] = (i * s.step_fm) * s.potential_MeV[i] * s.u[i];
  return std::sqrt(4.0 * std::acos(-1.0)) * Simpson(f, s.step_fm);
}

BoundState AnalyticTail(TailShape shape, const BoundStateRequest &q, double anc, double r_min) {
  CheckRequest(q);
  BoundState s;
  s.l = q.l;
  s.two_j = q.two_j;
  s.two_s = q.two_s;
  s.binding_MeV = q.binding_MeV;
  s.kappa_fm = rho_bound(q.pair, 1.0, q.binding_MeV);
  s.eta = shape == TailShape::Yukawa ? 0.0 : eta_bound(q.pair, q.binding_MeV);
  s.step_fm = q.step_fm;
  s.anc_fm_half = anc;
  const int n = static_cast<int>(std::lround(q.r_max_fm / q.step_fm)) + 1;
  s.u.assign(n, 0.0);
  s.potential_MeV.assign(n, 0.0);
  s.coulomb_MeV.assign(n, 0.0);
  for (int i = 1; i < n; ++i) {
    const double r = i * q.step_fm;
    if (r < r_min) continue;
    s.u[i] = shape == TailShape::Yukawa ? anc * std::exp(-s.kappa_fm * r)
                                        : anc * std::exp(LogWhittaker(q.l, s.eta, s.kappa_fm, r));
  }
  return s;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
