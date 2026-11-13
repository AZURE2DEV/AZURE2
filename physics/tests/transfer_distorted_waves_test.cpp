// transfer_distorted_waves (docs/upgrade/transfer.md sections 4 and 9): the distorted waves against answers that do
// not come from the library's integration.
//
//  1. No nuclear potential: S = 1 and u_l = F_l(eta, kr) to 1e-10 relative to max|F| (F of COULFG beyond the turning
//     point, of COULCC inside it), for a neutral pair and for a point charge (eta = 0.49), l = 0 ... 8, every node
//     (h = 0.004 fm: at smaller steps the round-off of 5000+ Numerov steps passes 1e-10).
//  2. A neutral square well (complex depth -V - iW, R = 3 fm), l = 0: S against e^{2i delta} from
//     K cot(KR) = k cot(kR + delta), K^2 = k^2 + 2 mu (V + iW)/hbar^2, solved here; 1e-6 at h = 2e-4 fm (O(h^2): the
//     step sits on a node, with the mean of the two sides).
//  3. A real potential (Woods-Saxon, spin-orbit, Coulomb sphere; spin 1/2 and 1): |S_lj| = 1 to 1e-10.
//  4. The reaction cross section of a pure absorber equals (pi/k^2) sum (2j+1)/(2s+1) (1 - |S|^2) written out here,
//     and sigma/sigma_R -> 1 for a potential that does not reach the matching radius at forward angles (5 deg).
//  5. The Coulomb sphere at the origin: TermsAt's V_C is 3 Z1 Z2 e^2/(2 R_C) at r = 0 and at r = 1e-300 fm (finite),
//     continuous at r = 1e-8 fm, and Z1 Z2 e^2/r outside (1e-14 relative).
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/transfer/DistortedWaves.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using cplx = std::complex<double>;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

const double kPi = std::acos(-1.0);

// ---------------------------------------------------------------- 1. free and point-Coulomb waves
void TestFreeWaves() {
  for (int charged = 0; charged <= 1; ++charged) {
    ScatteringRequest q;
    q.pair.z1 = charged ? 1 : 0;
    q.pair.z2 = 7;
    q.pair.mu_u = 15.0 / 16.0;
    q.energy_MeV = 4.6875;
    q.step_fm = 0.004;
    q.matching_radius_fm = 20.0;
    q.l_max = 8;
    const DistortedWaves w = SolveDistortedWaves(NoPotential(), q);
    double worst_s = 0.0, worst_u = 0.0;
    for (const PartialWave &pw : w.waves) {
      worst_s = std::max(worst_s, std::abs(pw.S - 1.0));
      double fmax = 0.0;
      std::vector<double> F(pw.u.size(), 0.0);
      for (std::size_t i = 1; i < pw.u.size(); ++i) {
        // the reference: COULFG beyond the turning point (6e-14), COULCC inside it (1.5e-12; 4.5e-13 at eta = 0)
        const double rho = w.k_fm * i * q.step_fm, tp = w.eta + std::sqrt(w.eta * w.eta + pw.l * (pw.l + 1.0));
        F[i] = coulomb::coulomb(pw.l, w.eta, rho, rho >= tp ? coulomb::Method::coulfg : coulomb::Method::coulcc).F;
        fmax = std::max(fmax, std::fabs(F[i]));
      }
      for (std::size_t i = 1; i < pw.u.size(); ++i) worst_u = std::max(worst_u, std::abs(pw.u[i] - F[i]) / fmax);
    }
    char buf[200];
    std::snprintf(buf, sizeof buf, "%s, no nuclear potential: |S - 1| %.2e, |u - F|/max|F| %.2e (eta %.3f)",
                  charged ? "point charge" : "neutral", worst_s, worst_u, w.eta);
    Check(worst_s < 1e-10 && worst_u < 1e-10, buf);
  }
}

// ---------------------------------------------------------------- 2. neutral square well
void TestSquareWell() {
  const double R = 3.0, V = 40.0, W = 5.0, mu_u = 15.0 / 16.0, E = 5.0;
  const double k2mu = 2.0 * mu_u * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  ScatteringRequest q;
  q.pair.z1 = 0;
  q.pair.z2 = 7;
  q.pair.mu_u = mu_u;
  q.energy_MeV = E;
  q.step_fm = 0.0002;
  q.matching_radius_fm = 8.0;
  q.l_max = 0;
  auto well = [&](double r) { return std::fabs(r - R) < 1e-9 ? cplx(-V, -W) * 0.5 : (r < R ? cplx(-V, -W) : cplx(0.0)); };
  const DistortedWaves w = SolveDistortedWavesInPotential(well, [](double) { return cplx(0.0); }, false, 0.0, q);
  const double k = std::sqrt(k2mu * E);
  const cplx K = std::sqrt(cplx(k * k + k2mu * V, k2mu * W));
  // u = A sin(Kr) inside, sin(kr + delta) outside: delta = atan(k tan(KR)/K) - kR (complex), S = e^{2 i delta}
  const cplx delta = std::atan(k * std::tan(K * R) / K) - k * R;
  const cplx S = std::exp(cplx(0.0, 2.0) * delta);
  char buf[200];
  std::snprintf(buf, sizeof buf, "square well l = 0: S = %.10f%+.10fi against %.10f%+.10fi", w.waves[0].S.real(),
                w.waves[0].S.imag(), S.real(), S.imag());
  Check(std::abs(w.waves[0].S - S) < 1e-6, buf);
  Check(std::fabs(w.k_fm - k) < 1e-14 * k, "k = sqrt(2 mu E)/hbar");
}

// ---------------------------------------------------------------- 3. real potential: |S| = 1
void TestUnitarity() {
  for (int two_s = 1; two_s <= 2; ++two_s) {
    OpticalPotential p;
    p.convention = RadiusConvention::TargetMass;
    p.target_mass_number = 15;
    p.coulomb_r0_fm = 1.25;
    p.real_volume = {50.0, 1.2, 0.65};
    p.real_spin_orbit = {6.0, 1.1, 0.6};
    ScatteringRequest q;
    q.pair.z1 = 1;
    q.pair.z2 = 7;
    q.pair.mu_u = two_s == 1 ? 15.0 / 16.0 : 30.0 / 17.0;
    q.two_s = two_s;
    q.energy_MeV = 6.0;
    q.step_fm = 0.05;
    q.matching_radius_fm = 20.0;
    q.l_max = 15;
    const DistortedWaves w = SolveDistortedWaves(p, q);
    double worst = 0.0;
    for (const PartialWave &pw : w.waves) worst = std::max(worst, std::fabs(std::abs(pw.S) - 1.0));
    char buf[160];
    std::snprintf(buf, sizeof buf, "real potential, 2s = %d: max ||S| - 1| = %.2e over %zu waves", two_s, worst,
                  w.waves.size());
    Check(worst < 1e-10 && w.waves.size() == static_cast<std::size_t>(16 * (two_s + 1) - two_s), buf);
  }
}

// ---------------------------------------------------------------- 4. reaction cross section, forward limit
void TestObservables() {
  OpticalPotential p;
  p.convention = RadiusConvention::TargetMass;
  p.target_mass_number = 15;
  p.coulomb_r0_fm = 1.3;
  p.real_volume = {90.0, 1.15, 0.75};
  p.imaginary_surface = {10.0, 1.35, 0.6};
  p.real_spin_orbit = {6.0, 1.07, 0.66};
  ScatteringRequest q;
  q.pair.z1 = 1;
  q.pair.z2 = 7;
  q.pair.mu_u = 30.0 / 17.0;
  q.two_s = 2;
  q.energy_MeV = 10.0 * 15.0 / 17.0;
  q.step_fm = 0.05;
  q.matching_radius_fm = 25.0;
  q.l_max = 25;
  const DistortedWaves w = SolveDistortedWaves(p, q);
  double sum = 0.0;
  for (const PartialWave &pw : w.waves) sum += (pw.two_j + 1.0) / 3.0 * (1.0 - std::norm(pw.S));
  const double expect = 10.0 * kPi / (w.k_fm * w.k_fm) * sum;
  Check(std::fabs(ReactionCrossSection(w) - expect) < 1e-12 * expect, "reaction cross section, spin 1");
  const double ratio = ElasticOverRutherford(w, 1.0);
  char buf[120];
  std::snprintf(buf, sizeof buf, "sigma/sigma_R at 1 deg = %.6f", ratio);
  Check(std::fabs(ratio - 1.0) < 1e-3, buf);
}

}  // namespace

void TestCoulombSphereAtOrigin() {
  OpticalPotential p;
  p.coulomb_r0_fm = 3.0;
  const int z1z2 = 6;
  const double e2 = kFineStructure * kHbarC_MeV_fm, centre = 1.5 * z1z2 * e2 / 3.0;
  double worst = 0.0;
  for (double r : {0.0, 1e-300, 1e-8}) {
    const double v = TermsAt(p, z1z2, r).coulomb_MeV;
    worst = std::max(worst, std::isfinite(v) ? std::fabs(v / centre - 1.0) : 1.0);
  }
  worst = std::max(worst, std::fabs(TermsAt(p, z1z2, 4.5).coulomb_MeV / (z1z2 * e2 / 4.5) - 1.0));
  Check(worst < 1e-14, "Coulomb sphere at r = 0, 1e-300, 1e-8 fm and outside: " + std::to_string(worst));
  std::printf("Coulomb sphere at the origin: largest relative difference %.2e\n", worst);
}

int main() {
  TestCoulombSphereAtOrigin();
  TestFreeWaves();
  TestSquareWell();
  TestUnitarity();
  TestObservables();
  if (failures) {
    std::printf("transfer_distorted_waves: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_distorted_waves: all checks passed\n");
  return 0;
}
