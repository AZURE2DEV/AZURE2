#include "azure2/physics/transfer/DistortedWaves.h"

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

DistortedWaves SolveDistortedWavesInPotential(const std::function<cplx(double)> &central_of,
                                              const std::function<cplx(double)> &so_of, bool has_spin_orbit,
                                              double point_charge, const ScatteringRequest &q) {
  if (!(q.energy_MeV > 0.0)) throw std::invalid_argument("distorted waves: E must be positive");
  if (q.two_s < 0 || q.two_s > 2) throw std::invalid_argument("distorted waves: projectile spin 0, 1/2 or 1");
  if (!(q.pair.mu_u > 0.0) || q.pair.relativistic) throw std::invalid_argument("distorted waves: mu > 0, non-relativistic");
  if (!(q.step_fm > 0.0) || !(q.matching_radius_fm > 2.0) || q.l_max < 0)
    throw std::invalid_argument("distorted waves: bad grid or l_max");
  DistortedWaves w;
  w.k_fm = rho(q.pair, 1.0, q.energy_MeV);
  w.eta = eta(q.pair, q.energy_MeV);
  w.step_fm = q.step_fm;
  w.two_s = q.two_s;
  w.spin_orbit = q.two_s > 0 && has_spin_orbit;
  const double h = q.step_fm;
  const int n = static_cast<int>(std::lround(q.matching_radius_fm / h)) + 1;
  if (n < 12) throw std::invalid_argument("distorted waves: fewer than 12 nodes");
  const int i1 = n - 1, i2 = n - 1 - std::max(5, static_cast<int>(std::lround(1.0 / h)));
  const double k2mu = 2.0 * q.pair.mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  const double k2 = w.k_fm * w.k_fm;
  const double q_m1 = -k2mu * point_charge;  // the 1/r coefficient of q (fm^-1)

  std::vector<cplx> central(n, 0.0), so(n, 0.0);
  for (int i = 1; i < n; ++i) {
    central[i] = central_of(i * h);
    if (w.spin_orbit) so[i] = so_of(i * h);
  }
  for (int l = 0; l <= q.l_max; ++l) w.sigma.push_back(coulomb::coulomb_phase(l, w.eta));

  std::vector<cplx> qv(n);
  for (int l = 0; l <= q.l_max; ++l) {
    const int j_lo = w.spin_orbit ? std::abs(2 * l - q.two_s) : 2 * l;
    const int j_hi = w.spin_orbit ? 2 * l + q.two_s : 2 * l;
    cplx hp1, hm1, hp2, hm2;
    Hankel(l, w.eta, w.k_fm * i1 * h, hp1, hm1);
    Hankel(l, w.eta, w.k_fm * i2 * h, hp2, hm2);
    for (int two_j = j_lo; two_j <= j_hi; two_j += 2) {
      const double two_ls = w.spin_orbit ? TwoLDotI(l, two_j, q.two_s) : 0.0;
      for (int i = 1; i < n; ++i)
        qv[i] = k2 - k2mu * (central[i] + so[i] * two_ls) - l * (l + 1) / ((i * h) * (i * h));
      // q_0: the finite part at the origin, taken at r = h (the point charge's C/r removed)
      const cplx q0 = k2 - k2mu * (central[1] - point_charge / h + so[1] * two_ls);
      std::vector<cplx> u = NumerovOutward(qv, h, l, q0, q_m1);
      PartialWave pw;
      pw.l = l;
      pw.two_j = two_j;
      pw.S = (hm1 * u[i2] - hm2 * u[i1]) / (hp1 * u[i2] - hp2 * u[i1]);
      if (q.keep_waves) {
        const cplx c = cplx(0.0, 0.5) * (hm1 - pw.S * hp1) / u[i1];
        for (cplx &x : u) x *= c;
        pw.u = std::move(u);
      }
      w.waves.push_back(std::move(pw));
    }
  }
  return w;
}

DistortedWaves SolveDistortedWaves(const OpticalPotential &potential, const ScatteringRequest &q) {
  const int z1z2 = q.pair.z1 * q.pair.z2;
  auto central = [&](double r) {
    const PotentialTerms t = TermsAt(potential, z1z2, r);
    return t.coulomb_MeV + t.volume_MeV + t.surface_MeV;
  };
  auto so = [&](double r) { return TermsAt(potential, z1z2, r).spin_orbit_form_MeV; };
  const double point = RadiusOf(potential, potential.coulomb_r0_fm) == 0.0 ? z1z2 * kFineStructure * kHbarC_MeV_fm : 0.0;
  return SolveDistortedWavesInPotential(central, so, HasSpinOrbit(potential), point, q);
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
