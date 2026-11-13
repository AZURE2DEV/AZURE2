#include "azure2/physics/transfer/DistortedWaves.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

/// H^+-(eta, rho) = G +- iF.
void Hankel(int l, double eta, double rho, cplx &h_plus, cplx &h_minus) {
  const coulomb::CoulombWaves c = CoulombForMatching(l, eta, rho);
  h_plus = {c.G, c.F};
  h_minus = {c.G, -c.F};
}

bool HasSpinOrbit(const OpticalPotential &p) {
  return p.real_spin_orbit.depth_MeV != 0.0 || p.imaginary_spin_orbit.depth_MeV != 0.0;
}

/// f_m'm(theta) of Eq. (4.3) (fm), m = two_m/2, m' = two_mp/2.
cplx Amplitude(const DistortedWaves &w, const std::vector<double> &ylm, int m_stride, double theta, int two_mp,
               int two_m) {
  cplx f = 0.0;
  if (two_mp == two_m && w.eta != 0.0) {
    const double s2 = std::pow(std::sin(0.5 * theta), 2);
    f = -w.eta / (2.0 * w.k_fm * s2) * std::exp(cplx(0.0, -w.eta * std::log(s2) + 2.0 * w.sigma[0]));
  }
  const int two_s = w.two_s;
  const int mu = (two_m - two_mp) / 2;  // m - m'
  for (const PartialWave &pw : w.waves) {
    const int l = pw.l;
    if (std::abs(mu) > l) continue;
    // Without spin-orbit one wave stands for every j: the sum over j of the two Clebsch-Gordan coefficients.
    const int j_lo = w.spin_orbit ? pw.two_j : std::abs(2 * l - two_s);
    const int j_hi = w.spin_orbit ? pw.two_j : 2 * l + two_s;
    double coupling = 0.0;
    for (int two_j = j_lo; two_j <= j_hi; two_j += 2)
      coupling += clebsch_gordan(2 * l, 0, two_s, two_m, two_j, two_m) *
                  clebsch_gordan(2 * l, 2 * mu, two_s, two_mp, two_j, two_m);
    if (coupling == 0.0) continue;
    const double y = ylm[l * m_stride + std::abs(mu)] * ((mu < 0 && (mu % 2 != 0)) ? -1.0 : 1.0);
    f += std::sqrt(4.0 * kPi * (2 * l + 1)) / (cplx(0.0, 2.0) * w.k_fm) * std::exp(cplx(0.0, 2.0 * w.sigma[l])) *
         (pw.S - 1.0) * coupling * y;
  }
  (void)theta;
  return f;
}

double SumSquares(const DistortedWaves &w, double theta_deg, double &rutherford) {
  if (!(theta_deg > 0.0) || theta_deg > 180.0) throw std::invalid_argument("elastic: 0 < theta <= 180 deg");
  const double theta = theta_deg * kPi / 180.0;
  int l_max = 0;
  for (const PartialWave &pw : w.waves) l_max = std::max(l_max, pw.l);
  const int m_max = w.two_s;  // |m - m'| <= 2s
  const std::vector<double> ylm = SphericalHarmonicsAtPhiZero(l_max, m_max, std::cos(theta));
  double sum = 0.0;
  for (int two_m = -w.two_s; two_m <= w.two_s; two_m += 2)
    for (int two_mp = -w.two_s; two_mp <= w.two_s; two_mp += 2)
      sum += std::norm(Amplitude(w, ylm, m_max + 1, theta, two_mp, two_m));
  sum /= (w.two_s + 1);
  rutherford = std::pow(w.eta / (2.0 * w.k_fm * std::pow(std::sin(0.5 * theta), 2)), 2);
  return sum;
}

}  // namespace

namespace {

void CheckRequest(const ScatteringRequest &q) {
  if (!(q.energy_MeV > 0.0)) throw std::invalid_argument("distorted waves: E must be positive");
  if (q.two_s < 0 || q.two_s > 2) throw std::invalid_argument("distorted waves: projectile spin 0, 1/2 or 1");
  if (!(q.pair.mu_u > 0.0) || q.pair.relativistic) throw std::invalid_argument("distorted waves: mu > 0, non-relativistic");
  if (!(q.step_fm > 0.0) || q.l_max < 0) throw std::invalid_argument("distorted waves: bad grid or l_max");
}

/// The grid and the kinematics of a channel reaching the node nearest radius_fm.
RadialChannel Grid(const ScatteringRequest &q, double radius_fm) {
  CheckRequest(q);
  RadialChannel c;
  c.k_fm = rho(q.pair, 1.0, q.energy_MeV);
  c.eta = eta(q.pair, q.energy_MeV);
  c.step_fm = q.step_fm;
  c.two_s = q.two_s;
  c.k2mu_per_MeV_fm2 = 2.0 * q.pair.mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  const int n = static_cast<int>(std::lround(radius_fm / q.step_fm)) + 1;
  if (n < 12) throw std::invalid_argument("distorted waves: fewer than 12 nodes");
  c.central_MeV.assign(n, 0.0);
  return c;
}

/// The table of central and spin-orbit terms from two functions of r.
void Fill(RadialChannel &c, const std::function<cplx(double)> &central_of, const std::function<cplx(double)> &so_of) {
  const int n = static_cast<int>(c.central_MeV.size());
  if (c.spin_orbit) c.spin_orbit_MeV.assign(n, 0.0);
  for (int i = 1; i < n; ++i) {
    c.central_MeV[i] = central_of(i * c.step_fm);
    if (c.spin_orbit) c.spin_orbit_MeV[i] = so_of(i * c.step_fm);
  }
  if (c.spin_orbit) c.spin_orbit_at_origin_MeV = c.spin_orbit_MeV[1];
}

/// u_l = x j_l(x), l = 0 ... l_max, at x >= 0 into u (l_max + 1 values): upward from j_0, j_1 where l_max < x,
/// else Miller's downward recurrence from far above both, normalised with sum (2l+1) j_l^2 = 1.
void RiccatiBessel(int l_max, double x, double *u) {
  if (!(x > 0.0)) {
    for (int l = 0; l <= l_max; ++l) u[l] = 0.0;
    return;
  }
  if (l_max < x) {
    double jm = std::sin(x) / x, j = std::sin(x) / (x * x) - std::cos(x) / x;
    u[0] = x * jm;
    for (int l = 1; l <= l_max; ++l) {
      u[l] = x * j;
      const double next = (2 * l + 1) / x * j - jm;
      jm = j;
      j = next;
    }
    return;
  }
  const double top = std::max<double>(l_max, x);
  const int start = static_cast<int>(top + 30.0 + std::sqrt(40.0 * top));
  std::vector<double> j(start + 2, 0.0);
  j[start] = 1.0;
  for (int l = start; l >= 1; --l) {
    j[l - 1] = (2 * l + 1) / x * j[l] - j[l + 1];
    if (std::fabs(j[l - 1]) > 1e100)  // kept below 1e100, so that the sum of squares cannot overflow
      for (int m = l - 1; m <= start; ++m) j[m] *= 1e-100;
  }
  double sum = 0.0;
  for (int l = start; l >= 0; --l) sum += (2 * l + 1) * j[l] * j[l];
  // The sign from the closed form of j_0 or, where it is small, j_1.
  const double j0 = std::sin(x) / x, j1 = std::sin(x) / (x * x) - std::cos(x) / x;
  const double sign = std::fabs(j0) > std::fabs(j1) ? (j0 * j[0] < 0.0 ? -1.0 : 1.0) : (j1 * j[1] < 0.0 ? -1.0 : 1.0);
  const double scale = sign / std::sqrt(sum);
  for (int l = 0; l <= l_max; ++l) u[l] = x * j[l] * scale;
}

/// u_l(r_i) = F_l(0, k r_i), l = 0 ... l_max, i < n: [l][i].
std::vector<std::vector<double>> RiccatiTable(int l_max, double k, double h, int n) {
  std::vector<std::vector<double>> u(l_max + 1, std::vector<double>(n, 0.0));
  std::vector<double> ric(l_max + 1);
  for (int i = 1; i < n; ++i) {
    RiccatiBessel(l_max, k * i * h, ric.data());
    for (int l = 0; l <= l_max; ++l) u[l][i] = ric[l];
  }
  return u;
}

/// R + 12 a of the largest term of a potential and its Coulomb radius (fm).
double PotentialRange(const OpticalPotential &p) {
  double range = RadiusOf(p, p.coulomb_r0_fm);
  for (const WoodsSaxonTerm *t : {&p.real_volume, &p.imaginary_volume, &p.imaginary_surface, &p.real_spin_orbit,
                                  &p.imaginary_spin_orbit})
    if (t->depth_MeV != 0.0) range = std::max(range, RadiusOf(p, t->r0_fm) + 12.0 * t->diffuseness_fm);
  return range;
}

/// The outer matching radius (fm) of wave l and the inner node's distance in nodes.
double OuterRadius(const RadialChannel &c, int l, int stored_nodes, const MatchingBeyondTurningPoint &rule,
                   int &inner_offset) {
  const double turning = (c.eta + std::sqrt(c.eta * c.eta + l * (l + 1.0))) / c.k_fm;
  const double quarter = 0.5 * kPi / c.k_fm;
  inner_offset = std::max(2, static_cast<int>(std::lround(quarter / c.step_fm)));
  return std::max(std::max((stored_nodes - 1) * c.step_fm, c.potential_range_fm), turning) +
         std::max(rule.quarter_wavelengths * quarter, rule.min_fm);
}

}  // namespace

RadialChannel TabulateChannel(const OpticalPotential &potential, const ScatteringRequest &request, double radius_fm) {
  RadialChannel c = Grid(request, radius_fm);
  const int z1z2 = request.pair.z1 * request.pair.z2;
  c.spin_orbit = request.two_s > 0 && HasSpinOrbit(potential);
  const double rc = RadiusOf(potential, potential.coulomb_r0_fm);
  c.point_charge_MeV_fm = rc == 0.0 ? z1z2 * kFineStructure * kHbarC_MeV_fm : 0.0;
  c.potential_range_fm = PotentialRange(potential);
  Fill(
      c,
      [&](double r) {
        const PotentialTerms t = TermsAt(potential, z1z2, r);
        return t.coulomb_MeV + t.volume_MeV + t.surface_MeV;
      },
      [&](double r) { return TermsAt(potential, z1z2, r).spin_orbit_form_MeV; });
  // At the origin: the terms at r = 0, with the sphere's 3/2 Z1 Z2 e^2/R_C (a point charge is C/r, not in it).
  const PotentialTerms t0 = TermsAt(potential, rc > 0.0 ? z1z2 : 0, 0.0);
  c.central_at_origin_MeV = t0.coulomb_MeV + t0.volume_MeV + t0.surface_MeV;
  return c;
}

RadialChannel PlaneWaveChannel(const ScatteringRequest &request, double radius_fm) {
  ScatteringRequest q = request;
  q.pair.z1 = q.pair.z2 = 0;
  RadialChannel c = Grid(q, radius_fm);
  c.plane_wave = true;
  return c;
}

WaveMatching MatchingAtLastNode(const RadialChannel &c) {
  const int n = static_cast<int>(c.central_MeV.size());
  return {n - 1, n - 1 - std::max(5, static_cast<int>(std::lround(1.0 / c.step_fm)))};
}

WaveMatching MatchingFor(const RadialChannel &c, int l, int stored_nodes, const MatchingBeyondTurningPoint &rule) {
  int d = 0;
  int outer = static_cast<int>(std::ceil(OuterRadius(c, l, stored_nodes, rule, d) / c.step_fm));
  if (outer - d < stored_nodes) outer = stored_nodes + d;
  return {outer, outer - d};
}

double MatchingRadius(const RadialChannel &c, int l_max, int stored_nodes, const MatchingBeyondTurningPoint &rule) {
  int top = 0;
  for (int l = 0; l <= l_max; ++l) top = std::max(top, MatchingFor(c, l, stored_nodes, rule).outer_node);
  return (top + 1) * c.step_fm;
}

PartialWave SolvePartialWave(const RadialChannel &c, int l, int two_j, const WaveMatching &m, int stored_nodes) {
  const int n_table = static_cast<int>(c.central_MeV.size());
  PartialWave pw;
  pw.l = l;
  pw.two_j = c.spin_orbit ? two_j : 2 * l;
  if (c.plane_wave) {
    if (stored_nodes > n_table) throw std::invalid_argument("plane wave: more stored nodes than the channel has");
    pw.S = 1.0;
    const std::vector<double> u = RiccatiTable(l, c.k_fm, c.step_fm, stored_nodes)[l];
    pw.u.assign(u.begin(), u.end());
    return pw;
  }
  if (!(m.outer_node < n_table) || !(m.inner_node > 2) || !(m.inner_node < m.outer_node) ||
      stored_nodes > m.outer_node + 1)
    throw std::invalid_argument("distorted waves: the matching nodes lie outside the channel's grid");
  const double h = c.step_fm, k2 = c.k_fm * c.k_fm, k2mu = c.k2mu_per_MeV_fm2;
  const double two_ls = c.spin_orbit ? TwoLDotI(l, two_j, c.two_s) : 0.0;
  const int n = m.outer_node + 1;
  std::vector<cplx> qv(n);
  for (int i = 1; i < n; ++i) {
    const cplx so = c.spin_orbit ? c.spin_orbit_MeV[i] * two_ls : cplx(0.0);
    qv[i] = k2 - k2mu * (c.central_MeV[i] + so) - l * (l + 1) / ((i * h) * (i * h));
  }
  const cplx q0 = k2 - k2mu * (c.central_at_origin_MeV + c.spin_orbit_at_origin_MeV * two_ls);
  std::vector<cplx> u = NumerovOutward(qv, h, l, q0, -k2mu * c.point_charge_MeV_fm);
  cplx hp1, hm1, hp2, hm2;
  Hankel(l, c.eta, c.k_fm * m.outer_node * h, hp1, hm1);
  Hankel(l, c.eta, c.k_fm * m.inner_node * h, hp2, hm2);
  const cplx u1 = u[m.outer_node], u2 = u[m.inner_node];
  pw.S = (hm1 * u2 - hm2 * u1) / (hp1 * u2 - hp2 * u1);
  const cplx scale = cplx(0.0, 0.5) * (hm1 - pw.S * hp1) / u1;
  if (!std::isfinite(std::abs(pw.S)) || !std::isfinite(std::abs(scale)) || !(std::abs(scale) > 0.0))
    throw std::runtime_error("distorted waves: the wave l = " + std::to_string(l) + " cannot be normalised");
  u.resize(stored_nodes);
  for (cplx &x : u) x *= scale;
  pw.u = std::move(u);
  return pw;
}

std::vector<double> CoulombPhases(const RadialChannel &c, int l_max) {
  std::vector<double> sigma;
  for (int l = 0; l <= l_max; ++l) sigma.push_back(c.plane_wave ? 0.0 : coulomb::coulomb_phase(l, c.eta));
  return sigma;
}

namespace {

/// Every wave l <= l_max (and every j) of a channel, matched at its last node.
DistortedWaves SolveAll(const RadialChannel &c, const ScatteringRequest &q) {
  DistortedWaves w;
  w.k_fm = c.k_fm;
  w.eta = c.eta;
  w.step_fm = c.step_fm;
  w.two_s = c.two_s;
  w.spin_orbit = c.spin_orbit;
  w.sigma = CoulombPhases(c, q.l_max);
  const int n = static_cast<int>(c.central_MeV.size());
  if (c.plane_wave) {
    const std::vector<std::vector<double>> u = RiccatiTable(q.l_max, c.k_fm, c.step_fm, n);
    for (int l = 0; l <= q.l_max; ++l) {
      PartialWave pw;
      pw.l = l;
      pw.two_j = 2 * l;
      pw.S = 1.0;
      if (q.keep_waves) pw.u.assign(u[l].begin(), u[l].end());
      w.waves.push_back(std::move(pw));
    }
    return w;
  }
  const WaveMatching m = MatchingAtLastNode(c);
  for (int l = 0; l <= q.l_max; ++l) {
    const int j_lo = w.spin_orbit ? std::abs(2 * l - q.two_s) : 2 * l;
    const int j_hi = w.spin_orbit ? 2 * l + q.two_s : 2 * l;
    for (int two_j = j_lo; two_j <= j_hi; two_j += 2) {
      PartialWave pw = SolvePartialWave(c, l, two_j, m, n);
      if (!q.keep_waves) pw.u.clear();
      w.waves.push_back(std::move(pw));
    }
  }
  return w;
}

void CheckMatchingRadius(const ScatteringRequest &q) {
  if (!(q.matching_radius_fm > 2.0)) throw std::invalid_argument("distorted waves: bad grid or l_max");
}

}  // namespace

DistortedWaves SolveDistortedWavesInPotential(const std::function<cplx(double)> &central_of,
                                              const std::function<cplx(double)> &so_of, bool has_spin_orbit,
                                              double point_charge, const ScatteringRequest &q) {
  CheckMatchingRadius(q);
  RadialChannel c = Grid(q, q.matching_radius_fm);
  c.spin_orbit = q.two_s > 0 && has_spin_orbit;
  c.point_charge_MeV_fm = point_charge;
  c.potential_range_fm = q.matching_radius_fm;
  Fill(c, central_of, so_of);
  // The finite part at the origin, taken at r = h (the point charge's C/r removed).
  c.central_at_origin_MeV = c.central_MeV[1] - point_charge / q.step_fm;
  return SolveAll(c, q);
}

DistortedWaves SolveDistortedWaves(const OpticalPotential &potential, const ScatteringRequest &q) {
  CheckMatchingRadius(q);
  return SolveAll(TabulateChannel(potential, q, q.matching_radius_fm), q);
}

DistortedWaves PlaneWaves(const ScatteringRequest &q) {
  CheckMatchingRadius(q);
  return SolveAll(PlaneWaveChannel(q, q.matching_radius_fm), q);
}

RadialChannel ChannelFor(const WaveChannel &ch, double step_fm, int stored_nodes, int l_max,
                         const MatchingBeyondTurningPoint &rule) {
  ScatteringRequest q;
  q.pair = ch.pair;
  q.energy_MeV = ch.energy_MeV;
  q.step_fm = step_fm;
  const double stored = (stored_nodes - 1) * step_fm;
  if (ch.plane_wave) return PlaneWaveChannel(q, std::max(stored, 12.0 * step_fm));
  // The kinematics and the potential's range from a short table, then the table as far as the waves need.
  const RadialChannel probe = TabulateChannel(ch.potential, q, 12.0 * step_fm);
  return TabulateChannel(ch.potential, q, MatchingRadius(probe, l_max, stored_nodes, rule));
}

DistortedWaves SolveStoredWaves(const WaveChannel &ch, double step_fm, int stored_nodes, int l_max,
                                const MatchingBeyondTurningPoint &rule) {
  const RadialChannel c = ChannelFor(ch, step_fm, stored_nodes, l_max, rule);
  DistortedWaves w;
  w.k_fm = c.k_fm;
  w.eta = c.eta;
  w.step_fm = c.step_fm;
  w.sigma = CoulombPhases(c, l_max);
  if (c.plane_wave) {
    const std::vector<std::vector<double>> u = RiccatiTable(l_max, c.k_fm, c.step_fm, stored_nodes);
    for (int l = 0; l <= l_max; ++l) {
      PartialWave pw;
      pw.l = l;
      pw.two_j = 2 * l;
      pw.S = 1.0;
      pw.u.assign(u[l].begin(), u[l].end());
      w.waves.push_back(std::move(pw));
    }
    return w;
  }
  for (int l = 0; l <= l_max; ++l)
    w.waves.push_back(SolvePartialWave(c, l, 2 * l, MatchingFor(c, l, stored_nodes, rule), stored_nodes));
  return w;
}

const PartialWave &WaveOf(const DistortedWaves &w, int l, int two_j) {
  for (const PartialWave &pw : w.waves)
    if (pw.l == l && (!w.spin_orbit || pw.two_j == two_j)) return pw;
  throw std::out_of_range("distorted waves: l = " + std::to_string(l) + " not solved");
}

double ElasticOverRutherford(const DistortedWaves &w, double theta_deg) {
  double rutherford = 0.0;
  const double sum = SumSquares(w, theta_deg, rutherford);
  return w.eta != 0.0 ? sum / rutherford : 10.0 * sum;
}

double ElasticCrossSection(const DistortedWaves &w, double theta_deg) {
  double rutherford = 0.0;
  return 10.0 * SumSquares(w, theta_deg, rutherford);
}

double ReactionCrossSection(const DistortedWaves &w) {
  double sum = 0.0;
  for (const PartialWave &pw : w.waves) {
    const int l = pw.l;
    const double weight = w.spin_orbit ? (pw.two_j + 1.0) : (2 * l + 1.0) * (w.two_s + 1.0);
    sum += weight / (w.two_s + 1.0) * (1.0 - std::norm(pw.S));
  }
  return 10.0 * kPi / (w.k_fm * w.k_fm) * sum;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
