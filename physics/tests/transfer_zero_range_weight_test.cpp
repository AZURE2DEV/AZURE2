// transfer_zero_range_weight (docs/upgrade/transfer.md section 24): the pieces the engine's THM distortion factor
// calls, against answers that do not come from them.
//
//  1. PlaneWaves: u_l = F_l(0, kr) = kr j_l(kr) for l <= 60 at kr from 1e-3 to 300: beyond kr = l/2 against the
//     Coulomb functions of physics/coulomb at eta = 0 (COULFG/COULCC through CoulombForMatching), 1e-12 of
//     max(1, |F_l|); below, against the power series x^{l+1}/(2l+1)!! sum_m (-x^2/2)^m/(m! (2l+3)...(2l+2m+1)), 1e-12 relative.
//  2. SolvePartialWave matched beyond the turning point (MatchingFor) against SolveDistortedWaves matched at the last
//     node of a longer grid, for a complex Woods-Saxon well with a Coulomb sphere, l <= 25: S and u agree to 1e-7
//     (both carry Numerov's error, which grows with the matching radius: h = 0.005 fm).
//  3. Plane waves in both channels with a Yukawa tail e^{-kappa r}: the weight is the Fourier transform of
//     phi = e^{-kappa r}/r, M = sqrt(4 pi) O_0 = 4 pi/(kappa^2 + q^2), q = |k_f - c_i k_i| (c_f = 1), at three
//     directions; ZeroRangeMomentumAmplitude the same; and with r_min = 1.5 fm both against the closed form
//     (4 pi/q) Im[e^{-(kappa - iq) r_min}/(kappa - iq)]; 1e-6 relative.
//  4. The sum's bound: twenty more partial waves after the stop change O_0 by less than 1e-12 relative.
//  5. The two implementations of Eq. (6.7) for an l = 0 overlap (ZeroRangeOrbitalIntegral, ZeroRangeWeightSum) with the
//     same waves, point Coulomb and optical, three directions: 1e-10 relative.
//  6. FRESCO, through the zero-range 12C(6Li,d)16O benchmark (tests/data/fresco/li6d_zero_range): the sum with
//     c_i = 1, c_f = alpha = m_A/m_F and the alpha-12C bound state as the l = 0 overlap gives
//     dsigma/dOmega = 10 mu_i mu_f/(2 pi hbar^2)^2 (k_f/k_i) |D0 O_0|^2 (all spins couple trivially for
//     l_a = l_F = 0), against FRESCO's reference to 1 % (the test of ZeroRangeCrossSection asks the same).
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "transfer_fresco_common.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/RadialNumerics.h"
#include "azure2/physics/transfer/Transfer.h"
#include "azure2/physics/transfer/ZeroRangeWeight.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

using cplx = std::complex<double>;
const double kPi = std::acos(-1.0);

void TestPlaneWaves() {
  ScatteringRequest q;
  q.pair.mu_u = 1.0;
  q.energy_MeV = 2.0;
  q.step_fm = 0.05;
  q.l_max = 60;
  const double k = rho(q.pair, 1.0, q.energy_MeV);
  q.matching_radius_fm = 300.0 / k;
  const DistortedWaves w = PlaneWaves(q);
  double worst = 0.0;
  for (int l = 0; l <= 60; l += 3)
    for (double x : {1e-3, 0.05, 0.7, 3.0, 11.0, 29.9, 45.0, 60.5, 140.0, 299.9}) {
      const int i = static_cast<int>(std::lround(x / (k * q.step_fm)));
      if (i == 0) continue;
      const double xi = k * i * q.step_fm;
      // Below the turning point F_l is small: relative to itself (its series where COULFG/COULCC underflow).
      double F = 0.0;
      if (xi < 0.5 * l) {
        double term = 1.0;
        for (int m = 1; m <= l; ++m) term *= xi / (2 * m + 1);
        term *= xi;  // x^{l+1}/(2l+1)!!
        for (int m = 0; m < 200 && std::fabs(term) > 1e-18 * std::fabs(F); ++m) {
          F += term;
          term *= -xi * xi / (2.0 * (m + 1) * (2 * l + 2 * m + 3));
        }
      } else {
        F = CoulombForMatching(l, 0.0, xi).F;
      }
      const double scale = xi < 0.5 * l ? std::fabs(F) : std::max(1.0, std::fabs(F));
      worst = std::max(worst, std::abs(WaveOf(w, l, 2 * l).u[i] - F) / scale);
    }
  char buf[120];
  std::snprintf(buf, sizeof buf, "plane waves against F_l(0, kr): %.2e", worst);
  Check(worst < 1e-12, buf);
  std::printf("%s\n", buf);
}

void TestMatching() {
  OpticalPotential p;
  // Diffuseness 0.3 fm: beyond R + 12 a and the margin the nuclear terms are below e^-17 of their depth, so that
  // what the two matchings differ by is the matching itself (the tail beyond is not in the closer one).
  p.real_volume = {50.0, 5.5, 0.3};
  p.imaginary_volume = {20.0, 5.6, 0.3};
  p.coulomb_r0_fm = 5.5;
  ScatteringRequest q;
  q.pair.z1 = 7, q.pair.z2 = 6, q.pair.mu_u = 14.0 * 12.0 / 26.0;
  q.energy_MeV = 13.8;
  q.step_fm = 0.005;
  q.matching_radius_fm = 60.0;
  q.l_max = 25;
  const DistortedWaves all = SolveDistortedWaves(p, q);
  WaveChannel ch;
  ch.potential = p;
  ch.pair = q.pair;
  ch.energy_MeV = q.energy_MeV;
  const int stored = 2001;  // to 10 fm
  const DistortedWaves w = SolveStoredWaves(ch, q.step_fm, stored, q.l_max, MatchingBeyondTurningPoint());
  double worst_s = 0.0, worst_u = 0.0;
  for (int l = 0; l <= q.l_max; ++l) {
    worst_s = std::max(worst_s, std::abs(w.waves[l].S - all.waves[l].S));
    double top = 0.0;
    for (int i = 0; i < stored; ++i) top = std::max(top, std::abs(all.waves[l].u[i]));
    for (int i = 0; i < stored; ++i)
      worst_u = std::max(worst_u, std::abs(w.waves[l].u[i] - all.waves[l].u[i]) / std::max(top, 1e-300));
  }
  char buf[120];
  std::snprintf(buf, sizeof buf, "matching beyond the turning point = at 60 fm: S %.2e, u %.2e", worst_s, worst_u);
  Check(worst_s < 1e-7 && worst_u < 1e-7, buf);
  std::printf("%s\n", buf);
}

/// The Yukawa tail of the s-x pair of 12C(14N,d) (unit ANC, r_min) on h = 0.02 fm to r_min + 50/kappa.
RadialOverlap Yukawa(double r_min, double &kappa) {
  BoundStateRequest b;
  b.pair.z1 = 1, b.pair.z2 = 6, b.pair.mu_u = 2.0 * 12.0 / 14.0;
  b.binding_MeV = 10.27;
  b.step_fm = 0.02;
  kappa = rho_bound(b.pair, 1.0, b.binding_MeV);
  b.r_max_fm = r_min + 2.0 * std::ceil(25.0 / kappa / b.step_fm) * b.step_fm;
  return RadialOverlapOf(AnalyticTail(TailShape::Yukawa, b, 1.0, r_min));
}

void TestPlaneWaveWeight() {
  WaveChannel in, out;
  in.plane_wave = out.plane_wave = true;
  in.pair.mu_u = 14.0 * 12.0 / 26.0, in.energy_MeV = 13.85;
  out.pair.mu_u = 2.0 * 24.0 / 26.0, out.energy_MeV = 2.6;
  ZeroRangeWeightRequest q;
  q.c_i = 2.0 / 14.0;
  double worst = 0.0;
  for (double r_min : {0.0, 1.5}) {
    double kappa = 0.0;
    const RadialOverlap o = Yukawa(r_min, kappa);
    const ZeroRangeEntrance e = PrepareZeroRangeEntrance(o, in, q);
    const std::vector<double> x = {1.0, 0.3, -1.0};
    const ZeroRangeWeight w = ZeroRangeWeightSum(o, e, out, x, q);
    const double ki = rho(in.pair, 1.0, in.energy_MeV) * q.c_i, kf = rho(out.pair, 1.0, out.energy_MeV);
    for (std::size_t d = 0; d < x.size(); ++d) {
      const double qq = std::sqrt(kf * kf + ki * ki - 2.0 * kf * ki * x[d]);
      const cplx z(kappa, -qq);
      const double exact = r_min == 0.0 ? 4.0 * kPi / (kappa * kappa + qq * qq)
                                        : 4.0 * kPi / qq * std::imag(std::exp(-z * r_min) / z);
      const double m = std::abs(std::sqrt(4.0 * kPi) * w.O[d]);
      const double pw = ZeroRangeMomentumAmplitude(o, qq);
      worst = std::max(worst, std::max(std::fabs(m - exact), std::fabs(pw - exact)) / exact);
      Check(std::fabs(m - exact) < 1e-6 * exact && std::fabs(pw - exact) < 1e-6 * exact,
            "plane waves, r_min " + std::to_string(r_min) + ", x " + std::to_string(x[d]) + ": M " + std::to_string(m) +
                ", M_PW " + std::to_string(pw) + ", closed form " + std::to_string(exact));
    }
  }
  std::printf("plane waves, Yukawa tail, r_min 0 and 1.5 fm: M and M_PW against the closed form, worst %.2e\n", worst);
}

void TestBound() {
  WaveChannel in, out;
  in.potential.real_volume = {50.0, 5.5, 0.6};
  in.potential.imaginary_volume = {20.0, 5.6, 0.6};
  in.potential.coulomb_r0_fm = 5.5;
  in.pair.z1 = 7, in.pair.z2 = 6, in.pair.mu_u = 14.0 * 12.0 / 26.0, in.energy_MeV = 13.85;
  out.potential.real_volume = {90.0, 3.4, 0.75};
  out.potential.imaginary_surface = {10.0, 3.8, 0.65};
  out.potential.coulomb_r0_fm = 3.8;
  out.pair.z1 = 1, out.pair.z2 = 12, out.pair.mu_u = 2.0 * 24.0 / 26.0, out.energy_MeV = 2.6;
  double kappa = 0.0;
  const RadialOverlap o = Yukawa(0.0, kappa);
  ZeroRangeWeightRequest q;
  q.c_i = 2.0 / 14.0;
  const ZeroRangeEntrance e = PrepareZeroRangeEntrance(o, in, q);
  const ZeroRangeWeight w = ZeroRangeWeightSum(o, e, out, {1.0}, q);
  q.extra_partial_waves = 20;
  const ZeroRangeWeight more = ZeroRangeWeightSum(o, e, out, {1.0}, q);
  const double change = std::abs(more.O[0] - w.O[0]) / std::abs(w.O[0]);
  char buf[120];
  std::snprintf(buf, sizeof buf, "twenty more partial waves after l = %d (entrance l <= %d): %.2e", w.l_max,
                static_cast<int>(e.u.size()) - 1, change);
  Check(more.l_max == std::min(w.l_max + 20, static_cast<int>(e.u.size()) - 1) && change < 1e-12, buf);
  std::printf("%s\n", buf);
}

/// The two implementations of Eq. (6.7) for an l = 0 overlap: ZeroRangeOrbitalIntegral (Clebsch-Gordan coefficients,
/// l <= 40, waves read by interpolation) and ZeroRangeWeightSum (Legendre recurrence, the sum stopped by its bound), at
/// c_i = 2/14, c_f = 1 with the same waves (SolveStoredWaves on the overlap's grid, c_i h for the entrance), point
/// Coulomb only and optical, three directions.
void TestTwoImplementations() {
  double worst = 0.0;
  for (bool optical : {false, true}) {
    WaveChannel in, out;
    in.pair.z1 = 7, in.pair.z2 = 6, in.pair.mu_u = 14.0 * 12.0 / 26.0, in.energy_MeV = 2.0;
    out.pair.z1 = 1, out.pair.z2 = 12, out.pair.mu_u = 2.0 * 24.0 / 26.0, out.energy_MeV = 1.0;
    if (optical) {
      in.potential.real_volume = {50.0, 5.5, 0.6};
      in.potential.imaginary_volume = {20.0, 5.6, 0.6};
      in.potential.coulomb_r0_fm = 5.5;
      out.potential.real_volume = {90.0, 3.4, 0.75};
      out.potential.imaginary_surface = {10.0, 3.8, 0.65};
      out.potential.coulomb_r0_fm = 3.8;
    }
    double kappa = 0.0;
    const RadialOverlap o = Yukawa(0.0, kappa);
    ZeroRangeWeightRequest q;
    q.c_i = 2.0 / 14.0;
    const int n = static_cast<int>(o.u.size());
    const DistortedWaves wi = SolveStoredWaves(in, q.c_i * o.step_fm, n, 40, q.matching);
    const DistortedWaves wf = SolveStoredWaves(out, o.step_fm, n, 40, q.matching);
    const ZeroRangeEntrance e = PrepareZeroRangeEntrance(o, in, q);
    const std::vector<double> theta = {0.0, 70.0, 180.0};
    std::vector<double> x;
    for (double th : theta) x.push_back(std::cos(th * kPi / 180.0));
    const ZeroRangeWeight w = ZeroRangeWeightSum(o, e, out, x, q);
    for (std::size_t d = 0; d < x.size(); ++d) {
      const cplx O = ZeroRangeOrbitalIntegral(o, wi, wf, q.c_i, 1.0, theta[d])[0];
      worst = std::max(worst, std::abs(O - w.O[d]) / std::abs(w.O[d]));
    }
  }
  char buf[140];
  std::snprintf(buf, sizeof buf,
                "ZeroRangeOrbitalIntegral = ZeroRangeWeightSum (Coulomb only and optical, 3 directions): %.2e", worst);
  Check(worst < 1e-10, buf);
  std::printf("%s\n", buf);
}

void TestFresco() {
  const auto runs = ReadCase("li6d_zero_range");
  if (runs.empty()) return;
  const ScatteringRequest in = Scattering(3, 6, 6.0, 12.0, 0, 20.0 * 12.0 / 18.0, 0.0, 38);
  const ScatteringRequest out = Scattering(1, 8, 2.0, 16.0, 0, 20.0 * 12.0 / 18.0, 5.6883, 38);
  const RadialOverlap residual =
      RadialOverlapOf(SolveBoundState(AlphaC12(), Bound(2, 6, 12.0, 4.0, 0, 0, 0, 7.162, 2)));
  WaveChannel wi, wf;
  wi.potential = Li6C12(), wi.pair = in.pair, wi.energy_MeV = in.energy_MeV;
  wf.potential = D16O(), wf.pair = out.pair, wf.energy_MeV = out.energy_MeV;
  ZeroRangeWeightRequest q;
  q.c_i = 1.0;
  q.c_f = 12.0 / 16.0;
  const ZeroRangeEntrance e = PrepareZeroRangeEntrance(residual, wi, q);
  std::vector<double> x;
  for (double th : Angles()) x.push_back(std::cos(th * kPi / 180.0));
  const ZeroRangeWeight w = ZeroRangeWeightSum(residual, e, wf, x, q);
  const double mu_i = in.pair.mu_u * kAtomicMassUnit_MeV, mu_f = out.pair.mu_u * kAtomicMassUnit_MeV;
  const double ki = rho(in.pair, 1.0, in.energy_MeV), kf = rho(out.pair, 1.0, out.energy_MeV);
  const double hbar2 = kHbarC_MeV_fm * kHbarC_MeV_fm, d0 = 68.3056;
  const double factor = 10.0 * mu_i * mu_f / std::pow(2.0 * kPi * hbar2, 2) * (kf / ki) * d0 * d0;
  std::vector<double> ours;
  for (const cplx &o : w.O) ours.push_back(factor * std::norm(o));
  CheckTransfer("li6d_zero_range", ours, 0.01);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  TestPlaneWaves();
  TestMatching();
  TestPlaneWaveWeight();
  TestBound();
  TestTwoImplementations();
  TestFresco();
  if (failures) {
    std::printf("transfer_zero_range_weight: %d check(s) failed\n", failures);
    return 1;
  }
  std::printf("transfer_zero_range_weight: all checks passed\n");
  return 0;
}
