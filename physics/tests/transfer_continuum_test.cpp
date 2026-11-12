// transfer_continuum (docs/upgrade/transfer.md section 14): the scattering states and bins of x + A against answers
// that do not come from the library's integration.
//
//  1. No nuclear potential: delta = 0 and u_E = F_l(eta, kr) to 1e-10 relative to max|F| (COULFG beyond the turning
//     point, COULCC inside it), neutral and point charge (eta = 0.49), l = 0 ... 6, every node (h = 0.004 fm).
//  2. A neutral square well (depth V0, radius R), l = 0: delta from K cot(KR) = k cot(kR + delta) and u from
//     sin(kr + delta) outside, A sin(Kr) inside, solved here; 1e-6 at h = 2e-4 fm (the step sits on a node).
//  3. The delta(E - E') normalisation of Eq. (14.2): a bin of width dE -> 0 is sqrt(dE) u_E^(E) (Eq. (14.3) with the
//     amplitude-1 states), checked at dE = 1e-4 MeV to 1e-6 inside the nucleus.
//  4. Orthonormality in one neutral Woods-Saxon well with a bound state: the bound state and a bin are orthogonal
//     (1e-6); a bin's norm on [0, R] is 1 - 2/(pi dk R) and two adjacent bins overlap by 1/(pi sqrt(dk dk') R) on
//     [0, R] (the 1/R tails of Eq. (14.3) beyond R, derived in transfer.md section 14), and two bins with a gap are
//     orthogonal: to 1e-5 at R = 1000 fm (1201 nodes in k) (the next terms are O(1/(dk^2 R^2))).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/transfer/BoundState.h"
#include "azure2/physics/transfer/Continuum.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;

namespace {

int failures = 0;
const double kPi = std::acos(-1.0);

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
  std::printf("%s\n", what.c_str());
}

double TwoMuOverHbar2(double mu_u) { return 2.0 * mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm); }

// ---------------------------------------------------------------- 1. free and point-Coulomb waves
void TestFree() {
  for (int charged = 0; charged <= 1; ++charged) {
    ContinuumRequest q;
    q.pair.z1 = charged ? 1 : 0;
    q.pair.z2 = 7;
    q.pair.mu_u = 15.0 / 16.0;
    q.step_fm = 0.004;
    q.r_max_fm = 20.0;
    double worst_u = 0.0, worst_d = 0.0;
    for (int l = 0; l <= 6; ++l) {
      q.l = l;
      q.two_j = 2 * l;
      const ScatteringState s = SolveScatteringState(NoPotential(), q, 4.6875);
      worst_d = std::max(worst_d, std::fabs(s.phase_shift_rad));
      double fmax = 0.0, dmax = 0.0;
      for (std::size_t i = 1; i < s.u.size(); ++i) {
        const double rho = s.k_fm * i * q.step_fm, tp = s.eta + std::sqrt(s.eta * s.eta + l * (l + 1.0));
        const double F =
            coulomb::coulomb(l, s.eta, rho, rho >= tp ? coulomb::Method::coulfg : coulomb::Method::coulcc).F;
        fmax = std::max(fmax, std::fabs(F));
        dmax = std::max(dmax, std::fabs(s.u[i] - F));
      }
      worst_u = std::max(worst_u, dmax / fmax);
    }
    char buf[200];
    std::snprintf(buf, sizeof buf, "no potential (%s): |delta| %.1e, u - F_l %.1e of max|F| (l <= 6)",
                  charged ? "eta = 0.49" : "neutral", worst_d, worst_u);
    Check(worst_d < 1e-10 && worst_u < 1e-10, buf);
  }
}

// ---------------------------------------------------------------- 2. square well
void TestSquareWell() {
  const double V0 = 30.0, R = 3.0, mu = 1.0, E = 5.0;
  ContinuumRequest q;
  q.pair.mu_u = mu;
  q.step_fm = 2e-4;
  q.r_max_fm = 12.0;
  RealPotential v;
  v.total_MeV = [&](double r) { return r < R - 1e-9 ? -V0 : (r > R + 1e-9 ? 0.0 : -0.5 * V0); };
  v.coulomb_MeV = [](double) { return 0.0; };
  const ScatteringState s = SolveScatteringStateIn(v, q, E);
  const double c = TwoMuOverHbar2(mu), k = std::sqrt(c * E), K = std::sqrt(c * (E + V0));
  const double delta = std::atan(k * std::tan(K * R) / K) - k * R;  // mod pi
  const double dd = std::fabs(std::sin(s.phase_shift_rad - delta));
  // u with the library's sign choice: delta in (-pi/2, pi/2]
  const double A = std::sin(k * R + s.phase_shift_rad) / std::sin(K * R);
  double worst = 0.0;
  for (std::size_t i = 1; i < s.u.size(); ++i) {
    const double r = i * q.step_fm;
    const double u = r < R ? A * std::sin(K * r) : std::sin(k * r + s.phase_shift_rad);
    worst = std::max(worst, std::fabs(s.u[i] - u));
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "square well: |sin(delta - delta_exact)| %.1e, u %.1e", dd, worst);
  Check(dd < 1e-6 && worst < 1e-6, buf);
}

// ---------------------------------------------------------------- 3. delta(E - E') normalisation
void TestNormalisation() {
  OpticalPotential p;
  p.convention = RadiusConvention::TargetMass;
  p.target_mass_number = 15;
  p.coulomb_r0_fm = 1.25;
  p.real_volume = {50.0, 1.25, 0.65};
  ContinuumRequest q;
  q.pair.z1 = 1;
  q.pair.z2 = 7;
  q.pair.mu_u = 15.0 / 16.0;
  q.l = 1;
  q.two_j = 1;
  q.two_s = 1;
  q.step_fm = 0.02;
  q.r_max_fm = 20.0;
  const double E = 2.0, dE = 1e-4;
  const ContinuumBin b = MakeBin(p, q, E, dE, 5, BinWeight::Flat, BinRule::Simpson);
  const ScatteringState s = SolveScatteringState(p, q, E);
  const double f = NormalisationFactor(s, q.pair.mu_u, ContinuumNormalisation::DeltaE);
  double worst = 0.0;
  for (int i = 50; i <= 250; i += 50)  // r = 1 ... 5 fm
    worst = std::max(worst, std::fabs(b.u[i] / (std::sqrt(dE) * f * s.u[i]) - 1.0));
  char buf[200];
  std::snprintf(buf, sizeof buf, "narrow bin = sqrt(dE) u^(E) at r = 1-5 fm: %.1e", worst);
  Check(worst < 1e-6, buf);
}

// ---------------------------------------------------------------- 4. orthonormality
double Overlap(const std::vector<double> &a, const std::vector<double> &b, double h) {
  const std::size_t n = std::min(a.size(), b.size());
  double s = 0.0;
  for (std::size_t i = 1; i < n; ++i) s += a[i] * b[i] * (i + 1 == n ? 0.5 : 1.0);
  return s * h;
}

void TestOrthonormality() {
  OpticalPotential p;
  p.convention = RadiusConvention::Absolute;
  p.real_volume = {40.0, 3.0, 0.6};
  BoundStateRequest bq;
  bq.pair.mu_u = 16.0 / 17.0;
  bq.binding_MeV = 3.0;
  bq.step_fm = 0.05;
  bq.r_max_fm = 40.0;
  const BoundState B = SolveBoundState(p, bq);
  p.real_volume.depth_MeV = B.depth_MeV;  // the continuum of the same well
  ContinuumRequest q;
  q.pair.mu_u = bq.pair.mu_u;
  q.step_fm = 0.05;
  q.r_max_fm = 1000.0;
  const double c = TwoMuOverHbar2(q.pair.mu_u);
  auto e_of = [&](double k) { return k * k / c; };
  // three bins in k: [0.8, 1.2], [1.2, 1.6] (adjacent), [1.8, 2.2] (a gap)
  const double k[4] = {0.8, 1.2, 1.6, 1.8};
  std::vector<ContinuumBin> bins;
  for (int i = 0; i < 3; ++i) {
    const double k1 = i < 2 ? k[i] : k[3], k2 = k1 + 0.4;
    bins.push_back(MakeBin(p, q, 0.5 * (e_of(k1) + e_of(k2)), e_of(k2) - e_of(k1), 1201, BinWeight::Flat,
                           BinRule::Simpson));
  }
  const double h = q.step_fm, R = q.r_max_fm, dk = 0.4;
  const double bound = std::fabs(Overlap(B.u, bins[0].u, h));
  const double norm = std::fabs(Overlap(bins[0].u, bins[0].u, h) + 2.0 / (kPi * dk * R) - 1.0);
  const double adjacent = std::fabs(Overlap(bins[0].u, bins[1].u, h) - 1.0 / (kPi * dk * R));
  const double gap = std::fabs(Overlap(bins[0].u, bins[2].u, h));
  char buf[300];
  std::snprintf(buf, sizeof buf,
                "orthonormality: <bound|bin> %.1e; bin norm on [0, R] - (1 - 2/(pi dk R)) %.1e; adjacent bins - "
                "1/(pi dk R) %.1e; bins with a gap %.1e (R = %.0f fm)",
                bound, norm, adjacent, gap, R);
  Check(bound < 1e-6 && norm < 1e-5 && adjacent < 1e-5 && gap < 1e-5, buf);
}

}  // namespace

int main() {
  TestFree();
  TestSquareWell();
  TestNormalisation();
  TestOrthonormality();
  if (failures) {
    std::printf("transfer_continuum: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_continuum: all checks passed\n");
  return 0;
}
